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
        print_g(f"Solved at z_wall: {z_wall}")
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
    f.max_grid_points = 4000                                                    
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
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=loglevel, refine_grid=refine_grid, auto=auto)                           
            delta_T_wall = get_delta_T(f, wall_params)
            if delta_T_wall < delta_T_max:
                strain_max = f.strain_rate("max")
                delta_T_ok = True
                print_c(f"Solved with \n m_f: {f.fuel_inlet.mdot}, \n m_o: {f.oxidizer_inlet.mdot}, \n delta_T_wall: {delta_T_wall}, \n n_grid: {f.grid.shape}")
                print_m(f"\n strain_max: {strain_max}, \n T_max = {np.max(f.T)}")
                print("##############################################################")
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














