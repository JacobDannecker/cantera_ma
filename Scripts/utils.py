from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import special
import cantera as ct
import time
from colorama import Fore, Back, Style

def add_attributes(f, file_path, name, wall_params, z_stoich, runtime):
    z_array = f.mixture_fraction(wall_params["mix_frac"])
    h_mass_array = f.enthalpy_mass
    h_mole_array = f.enthalpy_mole
    hdf5_file = h5py.File(file_path, "a")
    group_z = name + "/flame/z"
    group_h_mass = name + "/flame/h_mass"
    group_h_mole = name + "/flame/h_mole"
    group_flame = name + "/flame"
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
    hdf5_file[group_z].add_attr = "factor"

    hdf5_file[group_z].attrs["mix_frac"] = wall_params["mix_frac"]
    hdf5_file[group_z].attrs["fuel"] = wall_params["fuel"]
    hdf5_file[group_z].attrs["oxidizer"] = wall_params["oxidizer"]
    hdf5_file[group_z].attrs["basis"] = wall_params["basis"]
    hdf5_file[group_z].attrs["chi_st"] = chi_stoich(f, z_stoich)
    hdf5_file[group_z].attrs["z_wall"] = wall_params["Z_wall"]
    hdf5_file[group_z].attrs["factor"] = wall_params["factor"]
       
    if runtime:
        hdf5_file[group_flame].attrs["runtime"]= runtime
    else:
        hdf5_file[group_flame].attrs["runtime"]= "None"


    hdf5_file.close()

def save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True, runtime=False):
    # Save state of flame to hdf5 file. Add relevant data.
    z_wall = wall_params["Z_wall"]
    if info:
        print("##############################################################")
        print_g(f"Solved at z_wall: {z_wall}, saving flame")
        print("##############################################################\n")
    f.save(file_path, name=name, overwrite=True)
    add_attributes(f, file_path, name, wall_params, z_stoich, runtime)

def chi_stoich(f, z_stoich):
    a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
#    a = np.mean(np.abs(np.gradient(f.velocity / np.abs(np.gradient(f.mixture_fraction("H")))+0.0001)))
    chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
    return chi_stoich

def get_delta_T(f, wall_params):
    # Returns True if delta_T smaller than delta_T_max else reurns False
    idx_wall = np.abs(f.mixture_fraction(wall_params["mix_frac"]) - wall_params["Z_wall"]).argmin()
    delta_T_wall =  f.T[idx_wall] - wall_params["T_wall"]
    return delta_T_wall

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

def solve_with_wall(f, wall_params, delta_T_max=1.,
                    factor_last_working=False, factor_increase=2, factor_decrease=0.9, loglevel=0, refine_grid=True, auto=True):
    start_time = time.time()
    f.max_grid_points = 10000                                                    
    error_counter = 0
    z_wall = wall_params["Z_wall"]
    delta_T_ok = False
    failed = False
    wall_params["factor"] = 100.0
    max_factor = 1e19
    set_factor = False
    last_error_msg = ""
    while not delta_T_ok and not failed:
        try:
            if factor_last_working and not set_factor:
                wall_params["factor"] = factor_last_working
                set_factor = True
            else:
                wall_params["factor"] = min(
                    wall_params["factor"] * factor_increase, max_factor
                )
                if wall_params["factor"] >= max_factor:
                    raise SolveFailure("max factor", "Max factor reached.")
            
            print(f"Before ut.solve mdot f, o : {f.fuel_inlet.mdot}, {f.oxidizer_inlet.mdot}")
            print(f"Factor before ut.solve: {wall_params['factor']}")
            print(f"Grid refinement status: {refine_grid}")
            print(f"Gridpoints: {f.grid.shape}")
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=auto)                           
            delta_T_wall = get_delta_T(f, wall_params)
            print_c(f"Delta T wall: {delta_T_wall}")
            if delta_T_wall < delta_T_max:
                strain_max = f.strain_rate("max")
                delta_T_ok = True
                print_c(f"Solved with \n m_f: {f.fuel_inlet.mdot}, \n m_o: {f.oxidizer_inlet.mdot}, \n delta_T_wall: {delta_T_wall}, \n n_grid: {f.grid.shape}")
                print_m(f"\n strain_max: {strain_max}, \n T_max = {np.max(f.T)}")
                factor_last_working = wall_params["factor"]
                end_time = time.time()
                runtime = end_time - start_time
                return runtime, factor_last_working

        except ct.CanteraError as err:
            last_error_msg = str(err)
            print(err)
            error_counter += 1
            print_r(f"Had an exception in solve_with_wall errors: {error_counter}")
            if error_counter <= 3:
                if set_factor == True: 
                    pass
                else:
                    wall_params["factor"] /= factor_increase
                print_y(f"New factor {wall_params['factor']}")

                if factor_increase > 1.2:
                    factor_increase *= factor_decrease
            else:
                print_r("No solution found. Leaving solve_with_wall()")
                failed = True

    if failed:
        failure_type = classify_failure(last_error_msg)
        raise SolveFailure(failure_type,
            f"Failed solve_with_wall() at factor={wall_params['factor']:.1f}: "
            f"{last_error_msg}")


def runtime(func):
    def wrapper():
        start = time.time()
        end = time.time()
        val = func()
        total = end - start
        print(f"Runtime: {total}")
        return val
    return wrapper



def print_r(a):
    print(Fore.RED + a + Style.RESET_ALL)

def print_g(a):
    print(Fore.GREEN + a + Style.RESET_ALL)

def print_c(a):
    print(Fore.CYAN + a + Style.RESET_ALL)

def print_y(a):
    print(Fore.YELLOW + a + Style.RESET_ALL)

def print_m(a):
    print(Fore.MAGENTA + a + Style.RESET_ALL)


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

def perfect_v_shape(Z, h, zero_ends=False):                                     
    i_tip = np.argmin(h)                                                        
    if zero_ends:                                                               
        left_z = np.array([Z[0], Z[i_tip]])                                     
        left_h = np.array([0.0, h[i_tip]])                                      
        right_z = np.array([Z[i_tip], Z[-1]])                                   
        right_h = np.array([h[i_tip], 0.0])                                     
    else:                                                                       
        left_z = Z[: i_tip + 1]                                                 
        left_h = h[: i_tip + 1]                                                 
        right_z = Z[i_tip:]                                                     
        right_h = h[i_tip:]                                                     
    left = np.polyfit(left_z, left_h, 1)                                        
    right = np.polyfit(right_z, right_h, 1)                                     
    h_v = np.where(np.arange(len(Z)) <= i_tip, np.polyval(left, Z),             
                   np.polyval(right, Z))                                        
    return h_v, i_tip                                                          


def perfect_line(Z, h, zero_ends=False):                                        
    if zero_ends:                                                               
        z_line = np.array([Z[0], Z[-1]])                                        
        h_line = np.array([0.0, 0.0])                                           
    else:                                                                       
        z_line = Z                                                               
        h_line = h                                                               
    coeffs = np.polyfit(z_line, h_line, 1)                                       
    h_line = np.polyval(coeffs, Z)                                               
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
                                                                                
def correct_enthalpy(file_path, name, C, style="vshape"):                                                
   grid, T_orig, Y, a_max, idx_C, mech = load_data(file_path, name, C)          
   gas_i = ct.Solution(mech)                                                    
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
   print(f"Max dh: {max_dh},Max dt: {max_dT}")                                  
   return T_v, Y, h_v, Z, a_max, idx_C                                          
                                                  
