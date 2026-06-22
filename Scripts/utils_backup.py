from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import special
import cantera as ct

def add_attributes(f, file_path, name, wall_params, z_stoich):
    z_array = f.mixture_fraction(wall_params["mix_frac"])
    h_mass_array = f.enthalpy_mass
    h_mole_array = f.enthalpy_mole
    hdf5_file = h5py.File(file_path, "a")
    group_z = name + "/flame/z"
    group_h_mass = name + "/flame/h_mass"
    group_h_mole = name + "/flame/h_mole"
    # Data
    hdf5_file.create_dataset(name=group_z, data=z_array)
    hdf5_file.create_dataset(name=group_h_mass, data=h_mass_array)
    hdf5_file.create_dataset(name=group_h_mole, data=h_mole_array)
    # Attributes
    hdf5_file[group_z].add_attr = "mix_frac"
    hdf5_file[group_z].add_attr = "fuel"
    hdf5_file[group_z].add_attr = "oxidizer"
    hdf5_file[group_z].add_attr = "basis"
    hdf5_file[group_z].add_attr = "chi_st"
    hdf5_file[group_z].add_attr = "z_wall"
    hdf5_file[group_z].attrs["mix_frac"] = wall_params["mix_frac"]
    hdf5_file[group_z].attrs["fuel"] = wall_params["fuel"]
    hdf5_file[group_z].attrs["oxidizer"] = wall_params["oxidizer"]
    hdf5_file[group_z].attrs["basis"] = wall_params["basis"]
    hdf5_file[group_z].attrs["chi_st"] = chi_stoich(f, z_stoich)
    hdf5_file[group_z].attrs["z_wall"] = wall_params["z_wall"]

    hdf5_file.close()

def save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True):
    # Save state of flame to hdf5 file. Add relevant data.
    z_wall = wall_params["z_wall"]
    if info:
        print(f"##### Saved solution at z_wall: {z_wall} ##### \n")
    f.save(file_path, name=name, overwrite=True)
    add_attributes(f, file_path, name, wall_params, z_stoich)

def chi_stoich(f, z_stoich):
    a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
    chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
    return chi_stoich

def get_delta_T(f, wall_params):
    # Returns True if delta_T smaller than delta_T_max else reurns False
    idx_wall = np.abs(f.mixture_fraction(wall_params["mix_frac"]) - wall_params["z_wall"]).argmin()
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

def solve_with_wall(f, wall_params, delta_T_max=1., factor_increase=2, loglevel=0):
    error_counter = 0
    z_wall = wall_params["z_wall"]
    failed = False
    max_factor = 1e25
    last_error_msg = ""
    print(f"Start solving with wall at z: {z_wall}")
    while True:
        try:
            wall_params["factor"] = min(wall_params["factor"] * factor_increase, max_factor)
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=loglevel, refine_grid=True, auto=True)                           
            delta_T_wall = get_delta_T(f, wall_params)
            if delta_T_wall < delta_T_max:
                print(f"Solution found at delta_T_wall: {delta_T_wall}")
                print(f"mdot f, o : {f.fuel_inlet.mdot}, {f.oxidizer_inlet.mdot}")
                print(f"Strain max: {f.strain_rate("max")}")
                return True
        except ct.CanteraError as err:
            last_error_msg = str(err)
            print(err)
            error_counter += 1
            print(f"+++++ Had an exception in solve_with_wall. Error count {error_counter}")
            if error_counter <= 3:
                # Reset factor and reduce factor_increase
                wall_params["factor"] /= factor_increase
                if factor_increase > 1.2:
                    factor_increase *= 0.9
            else:
                print("+++++ No solution found :/ Leaving solve_with_wall()")
                raise SolveFailure(failure_type,                                        
                f"Failed solve_with_wall() at factor={wall_params['factor']:.1f}: " 
                f"{last_error_msg}")  
                return False 






















