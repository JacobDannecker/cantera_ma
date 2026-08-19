from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import special
import cantera as ct
import time
from colorama import Fore, Back, Style

def add_attributes(f, file_path, name, wall_params, z_stoich, runtime, enthalpy_refinement, enthalpy_curve):
    z_array = f.mixture_fraction(wall_params["mix_frac"])
    h_mass_array = f.enthalpy_mass
    h_mole_array = f.enthalpy_mole
    hdf5_file = h5py.File(file_path, "a")
    group_z = name + "/flame/z"
    group_h_mass = name + "/flame/h_mass"
    group_h_mole = name + "/flame/h_mole"
    group_flame = name + "/flame"
    group_refine = name + "/flame/refine-criteria"
    # Data
    hdf5_file.create_dataset(name=group_z, data=z_array)
    hdf5_file.create_dataset(name=group_h_mass, data=h_mass_array)
    hdf5_file.create_dataset(name=group_h_mole, data=h_mole_array)
    # Attributes
    hdf5_file[group_flame].add_attr = "runtime"
    hdf5_file[group_z].add_attr = "mix_frac"
    hdf5_file[group_z].add_attr = "fuel"
    hdf5_file[group_z].add_attr = "oxidizer"
    hdf5_file[group_z].add_attr = "basis"
    hdf5_file[group_z].add_attr = "chi_st"
    hdf5_file[group_z].add_attr = "z_wall"
    hdf5_file[group_z].add_attr = "z_wall_width"
    hdf5_file[group_refine].add_attr = "enthalpy_refinement"
    hdf5_file[group_refine].add_attr = "enthalpy_curve"

    hdf5_file[group_z].add_attr = "factor"

    hdf5_file[group_z].attrs["mix_frac"] = wall_params["mix_frac"]
    hdf5_file[group_z].attrs["fuel"] = wall_params["fuel"]
    hdf5_file[group_z].attrs["oxidizer"] = wall_params["oxidizer"]
    hdf5_file[group_z].attrs["basis"] = wall_params["basis"]
    hdf5_file[group_z].attrs["chi_st"] = chi_stoich(f, z_stoich)
    hdf5_file[group_z].attrs["z_wall"] = wall_params["Z_wall"]
    hdf5_file[group_z].attrs["factor"] = wall_params["factor"]

    hdf5_file[group_z].attrs["z_wall_width"] = wall_params["Z_wall_width"]
    hdf5_file[group_refine].attrs["enthalpy_refinement"] = enthalpy_refinement
    hdf5_file[group_refine].attrs["enthalpy_curve"] = enthalpy_curve
    
    if runtime:
        hdf5_file[group_flame].attrs["runtime"]= runtime
    else:
        hdf5_file[group_flame].attrs["runtime"]= "None"


    hdf5_file.close()

def save_with_attributes(f, file_path, name, wall_params, z_stoich, enthalpy_refinement, enthalpy_curve, info=True, runtime=False):
    # Save state of flame to hdf5 file. Add relevant data.
    z_wall = wall_params["Z_wall"]
    if info:
        print("##############################################################")
        print_g(f"Solved at z_wall: {z_wall}, saving flame")
        print("##############################################################\n")
    f.save(file_path, name=name, overwrite=True)
    add_attributes(f, file_path, name, wall_params, z_stoich, runtime, enthalpy_refinement, enthalpy_curve)

def chi_stoich(f, z_stoich):
    a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
#    a = np.mean(np.abs(np.gradient(f.velocity / np.abs(np.gradient(f.mixture_fraction("H")))+0.0001)))
    chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
    return chi_stoich

def get_delta_T(f, wall_params):
    # Worst-case (largest) remaining temperature deviation from T_wall among
    # grid points actually inside the wall's influence (Z >= Z_wall). Using
    # only the single nearest-to-Z_wall point instead can pick a point just
    # outside a narrow blending window (e.g. Z_wall close to 1), where the
    # sink has no effect and increasing "factor" can never converge it.
    Z = f.mixture_fraction(wall_params["mix_frac"])
    T = f.T
    mask = Z >= wall_params["Z_wall"]
    if not np.any(mask):
        # No point currently classified as wall-affected; fall back to the
        # nearest point as the best available proxy.
        idx_wall = np.abs(Z - wall_params["Z_wall"]).argmin()
        return f.T[idx_wall] - wall_params["T_wall"]
    return np.max(T[mask] - wall_params["T_wall"])

def get_z_stoich(gas, wall_params, reaction_yaml):
    fuel = wall_params["fuel"]
    oxidizer = wall_params["oxidizer"]
    mix_frac = wall_params["mix_frac"]
    gas.set_equivalence_ratio(1.0, fuel, oxidizer)
    z_st = gas.mixture_fraction(fuel, oxidizer, element=mix_frac)
    return z_st

class SolveFailure(Exception):
    """Carries a failure_type classifying why solve_with_wall failed."""
    def __init__(self, failure_type, message):
        self.failure_type = failure_type
        super().__init__(message)

def classify_failure(msg):
    msg_lower = msg.lower()
    if "max number of grid points" in msg_lower:
        return "grid_limit"
    elif "bandmatrix" in msg_lower or ("matrix" in msg_lower and "singular" in msg_lower):
        return "singular"
    elif "newton steady-state solve failed" in msg_lower:
        return "convergence"
    elif "timestep" in msg_lower or "time step" in msg_lower:
        return "timestep"
    elif "enthalpy" in msg_lower:
        return "enthalpy_refinement"
    else:
        return "unknown"

#: Cold-start factor seed used when no factor_last_working is supplied.
DEFAULT_FACTOR_SEED = 1e6
_MAX_LOG_STEP = np.log(10.0)   # cap a single secant step to 10x in factor
_MIN_LOG_STEP = np.log(1.05)   # minimum forward progress per secant step

def get_next_factor(history, current_factor, factor_increase, delta_T_max, max_factor):
    """Propose the next wall 'factor' via a bounded secant step on
    ln(factor) targeting delta_T_wall == delta_T_max, given the
    (ln(factor), delta_T_wall) history tried so far. Falls back to geometric
    doubling when there's not enough history or the secant estimate isn't
    trustworthy (non-finite, or pointing the wrong way)."""
    if len(history) < 2:
        next_factor = current_factor * factor_increase
    else:
        (x0, y0), (x1, y1) = history[-2], history[-1]
        y0 -= delta_T_max
        y1 -= delta_T_max
        if y1 == y0 or not np.isfinite(y0) or not np.isfinite(y1):
            next_factor = current_factor * factor_increase
        else:
            step = x1 - y1 * (x1 - x0) / (y1 - y0) - x1
            if not np.isfinite(step) or step <= 0:
                next_factor = current_factor * factor_increase
            else:
                step = min(max(step, _MIN_LOG_STEP), _MAX_LOG_STEP)
                next_factor = np.exp(x1 + step)
    return min(next_factor, max_factor)

def solve_with_wall(f, wall_params, delta_T_max=1., factor_last_working=False,
                    factor_seed=DEFAULT_FACTOR_SEED, factor_increase=2,
                    factor_decrease=0.9, max_plain_errors=3, max_auto_errors=2,
                    loglevel=0, refine_grid=True, auto=False):
    """Solve with the permeable wall active, searching for the smallest
    'factor' (sink stiffness) that pulls delta_T_wall below delta_T_max: a
    secant search on ln(factor), escalating to auto=True after repeated
    plain-solve failures.

    :return: (runtime, factor_last_working)
    :raises SolveFailure: if no converged solution is found.
    """
    start_time = time.time()
    f.max_grid_points = 10000
    max_factor = 1e19

    wall_params["factor"] = factor_last_working if factor_last_working else factor_seed

    history = []  # (ln(factor), delta_T_wall) pairs tried at this step
    plain_errors = 0
    auto_errors = 0
    use_auto = auto
    last_error_msg = ""

    while True:
        try:
            print(f"Before solve_with_wall mdot f, o : {f.fuel_inlet.mdot}, {f.oxidizer_inlet.mdot}")
            print(f"Factor before solve: {wall_params['factor']} (auto={use_auto})")
            print(f"Grid refinement status: {refine_grid}")
            print(f"Gridpoints: {f.grid.shape}")

            f.flame.set_non_adiabatic_wall(wall_params)
            f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=use_auto)
            delta_T_wall = get_delta_T(f, wall_params)
            print_c(f"Delta T wall: {delta_T_wall}")
            history.append((np.log(wall_params["factor"]), delta_T_wall))
            if delta_T_wall < delta_T_max:
                strain_max = f.strain_rate("max")
                print_c(f"Solved with \n m_f: {f.fuel_inlet.mdot}, \n m_o: {f.oxidizer_inlet.mdot}, "
                        f"\n delta_T_wall: {delta_T_wall}, \n n_grid: {f.grid.shape}")
                print_m(f"\n strain_max: {strain_max}, \n T_max = {np.max(f.T)}")
                runtime = time.time() - start_time
                return runtime, wall_params["factor"]
            correct_enthalpy_flame(f, "H2O")                                         

            # A successful (if insufficiently strong) solve resets the
            # plain-solve error streak and drops back to the caller's auto
            # setting; if factor is already pinned at its cap, further
            # increases can't help, so stop instead of looping forever.
            if wall_params["factor"] >= max_factor:
                raise SolveFailure("max factor",
                    f"Reached max factor ({max_factor:.3g}) without reaching "
                    f"delta_T_max={delta_T_max} (stuck at delta_T_wall="
                    f"{delta_T_wall:.4g}).")
            plain_errors = 0
            use_auto = auto
            wall_params["factor"] = get_next_factor(
                history, wall_params["factor"], factor_increase, delta_T_max, max_factor)

        except ct.CanteraError as err:
            last_error_msg = str(err)
            print(err)

            if not use_auto and plain_errors < max_plain_errors:
                plain_errors += 1
                wall_params["factor"] /= factor_increase
                factor_increase = max(factor_increase * factor_decrease, 1.01)
                print_y(f"Plain solve failed ({plain_errors}/{max_plain_errors}). "
                        f"New factor {wall_params['factor']}")
            elif auto_errors < max_auto_errors:
                auto_errors += 1
                use_auto = True
                print_y(f"Falling back to auto continuation "
                        f"({auto_errors}/{max_auto_errors}) at factor "
                        f"{wall_params['factor']}")
            else:
                failure_type = classify_failure(last_error_msg)
                raise SolveFailure(failure_type,
                    f"Failed solve_with_wall() at factor={wall_params['factor']:.4g}: "
                    f"{last_error_msg}") from err

            if wall_params["factor"] >= max_factor:
                raise SolveFailure("max factor", "Max factor reached.") from err

def runtime(func):
    def wrapper():
        start = time.time()
        val = func()
        end = time.time()
        total = end - start
        print(f"Runtime: {total}")
        return val
    return wrapper

def print_r(a):
    print(Fore.RED + str(a) + Style.RESET_ALL)

def print_g(a):
    print(Fore.GREEN + str(a) + Style.RESET_ALL)

def print_c(a):
    print(Fore.CYAN + str(a) + Style.RESET_ALL)

def print_y(a):
    print(Fore.YELLOW + str(a) + Style.RESET_ALL)

def print_m(a):
    print(Fore.MAGENTA + str(a) + Style.RESET_ALL)

def load_data(file_path, name, C):                                              
    reaction_mechanism = "h2o2.yaml"                                            
    gas = ct.Solution(reaction_mechanism)                                       
    width = 18e-3                                                               
    grid = np.linspace(0, width, 250)                                           
    f = ct.CounterflowDiffusionFlame(gas, grid=grid)                            
    f.P = 1.e5                                                                  
    f.fuel_inlet.mdot = 0.5                                                     
    f.oxidizer_inlet.mdot = 3.0                                                 
    f.fuel_inlet.X = "H2:1"                                                     
    f.oxidizer_inlet.X = "O2:1"                                                 
    f.fuel_inlet.T = 300                                                        
    f.oxidizer_inlet.T = 300                                                    
    f.transport_model = "unity-Lewis-number"                                    
    f.restore(file_path, name)                                                  
    return f.grid, f.T, f.Y, f.strain_rate("max"), f.gas.species_index(C), reaction_mechanism
                          
def get_species_index(species):
    reaction_mechanism = "h2o2.yaml"                                            
    gas = ct.Solution(reaction_mechanism)                                       
    return gas.species_index(species)

def compute_enthalpy_and_Z(gas, T, Y, fuel_idx=0, oxidizer_idx=3):              
    n = len(T)                                                                  
    h = np.empty(n)                                                             
    Z = np.empty(n)                                                             
    fuel = np.zeros(gas.n_species)                                              
    oxidizer = np.zeros(gas.n_species)                                          
    fuel[fuel_idx] = 1.0                                                        
    oxidizer[oxidizer_idx] = 1.0                                                
    for j in range(n):                                                          
        gas.TPY = T[j], gas.P, Y.T[j]                                           
        h[j] = gas.enthalpy_mass                                                
        Z[j] = gas.mixture_fraction(fuel, oxidizer, basis="mass")               
    return h, Z        

def compute_cold_enthalpy_and_Z(gas, Y, T_cold=300.0, fuel_idx=0, oxidizer_idx=3):
    n = Y.shape[1]
    h = np.empty(n)
    Z = np.empty(n)
    fuel = np.zeros(gas.n_species);  fuel[fuel_idx] = 1.0
    oxidizer = np.zeros(gas.n_species);  oxidizer[oxidizer_idx] = 1.0
    for j in range(n):
        gas.TPY = T_cold, gas.P, Y.T[j]
        h[j] = gas.enthalpy_mass
        Z[j] = gas.mixture_fraction(fuel, oxidizer, basis="mass")
    return h, Z

def perfect_v_shape(Z, h, zero_ends=False):
    i_tip = np.argmin(h)
    z_min = Z[i_tip]
    if zero_ends:
        h_left, h_right = 0.0, 0.0
    else:
        h_left, h_right = h[0], h[-1]

    z_range = Z[-1] - Z[0]
    tol = 1e-10 * max(abs(z_range), 1e-12)
    left_span = Z[i_tip] - Z[0]
    right_span = Z[-1] - Z[i_tip]
    left_slope = (h[i_tip] - h_left) / left_span if abs(left_span) > tol else 0.0
    right_slope = (h_right - h[i_tip]) / right_span if abs(right_span) > tol else 0.0
    h_v = np.where(
        np.arange(len(Z)) <= i_tip,
        h_left + left_slope * (Z - Z[0]),
        h[i_tip] + right_slope * (Z - Z[i_tip])
    )
    return h_v, i_tip
                                                          
def perfect_line(Z, h, zero_ends=False):
    if zero_ends:
        z_line = np.array([Z[0], Z[-1]])
        h_line = np.array([0.0, 0.0])
    else:
        z_line = np.array([Z[0], Z[-1]])
        h_line = np.array([h[0], h[-1]])
    slope = (h_line[1] - h_line[0]) / (z_line[1] - z_line[0])
    intercept = h_line[0] - slope * z_line[0]
    h_line = slope * Z + intercept
    return h_line

def temperature_from_HPY(gas, h, Y, P=None):
   if P is None:                                                                
       P = gas.P                                                                
   n = len(h)                                                                   
   T = np.empty(n)                                                              
   for j in range(n):                                                           
       gas.HPY = h[j], P, Y.T[j]                                                
       T[j] = gas.T                                                             
   return T                                                                     
                                                                                
def correct_enthalpy(file_path, name, C, style="vshape", active=True, zero_ends=False):                                                
   grid, T_orig, Y, a_max, idx_C, mech = load_data(file_path, name, C)          
   gas_i = ct.Solution(mech)                                                    
   P = gas_i.P                                                                  
   h_orig, Z = compute_enthalpy_and_Z(gas_i, T_orig, Y)                         
   if not active:
       return T_orig, Y, h_orig, Z, a_max, idx_C, 0.0, 0.0
   if style == "vshape":
       h_v, i_tip = perfect_v_shape(Z, h_orig, zero_ends=zero_ends)                     
   elif style == "line":
       h_v = perfect_line(Z, h_orig, zero_ends=zero_ends)                     
   T_v = temperature_from_HPY(gas_i, h_v, Y)                                    
   max_dh = np.max(np.abs(h_orig - h_v))                                        
   max_dT = np.max(np.abs(T_orig - T_v))                                        
   print(f"Max dh: {max_dh},Max dt: {max_dT}")                                  
   return T_v, Y, h_v, Z, a_max, idx_C, max_dh, max_dT                                          

def correct_enthalpy_flame(f, C, style="vshape"):                                                
    grid = f.grid
    T_orig = f.T
    Y = f.Y
    gas_i = ct.Solution("h2o2.yaml")                                                    
    P = gas_i.P                                                                  
    print(grid.shape, T_orig.shape, Y.shape)                                     
    h_orig, Z = compute_enthalpy_and_Z(gas_i, T_orig, Y)                         
    if style == "vshape":
        h_v, i_tip = perfect_v_shape(Z, h_orig, zero_ends=False)                     
    elif style == "line":
        h_v = perfect_line(Z, h_orig, zero_ends=False)                     
    T_v = temperature_from_HPY(gas_i, h_v, Y)                                    
    max_dh = np.max(np.abs(h_orig - h_v))                                        
    max_dT = np.max(np.abs(T_orig - T_v))                                        
    f.flame.set_values("T", T_v)                                              
                                                
def rms(data):
    return np.sqrt(np.mean(data ** 2))
