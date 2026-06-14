# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
Diffusion flame extinction strain rate
======================================

This example computes the extinction point of a counterflow diffusion flame.
A hydrogen-oxygen diffusion flame at 1 bar is studied.

The tutorial makes use of the scaling rules derived by Fiala and Sattelmayer
(doi:10.1155/2014/484372). Please refer to this publication for a detailed
explanation. Also, please don't forget to cite it if you make use of it.

Requires: cantera >= 3.2, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, diffusion flame, strained flame, extinction,
          saving output, plotting
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import special
import cantera as ct

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
    hdf5_file[group_z].attrs["z_wall"] = wall_params["Z_wall"]

    hdf5_file.close()

def chi_stoich(f, z_stoich):
    #a = f.strain_rate("mean")
    a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
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

def save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True):
    # Save state of flame to hdf5 file. Add relevant data.
    z_wall = wall_params["Z_wall"]
    if info:
        print("\n##############################################################")
        print(f"Solved at z_wall: {z_wall}")
        print("##############################################################")
    f.save(file_path, name=name, overwrite=True)
    add_attributes(f, file_path, name, wall_params, z_stoich)


def solve_with_wall(f, wall_params, name_fallback="initial", delta_T_max=1.,
                    factor_last_working=False, factor_increase=2, loglevel=0):
    error_counter = 0
    z_wall = wall_params["Z_wall"]
    delta_T_ok = False
    failed = False
    wall_params["factor"] = 100.0
    max_factor = 1e25
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
            print(wall_params["factor"])
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=loglevel, refine_grid=True, auto=True)                           
            delta_T_wall = get_delta_T(f, wall_params)
            if delta_T_wall < delta_T_max:
                print(f"mdot f, o : {f.fuel_inlet.mdot}, {f.oxidizer_inlet.mdot}")
                strain_max = f.strain_rate("max")
                print(f"Strain max: {strain_max}")
                delta_T_ok = True
                print(f"Solution found at delta_T_wall: {delta_T_wall}")
                return factor_last_working

        except ct.CanteraError as err:
            last_error_msg = str(err)
            print(err)
            error_counter += 1
            print(f"Error count {error_counter}")
            print("======================Had an exception in solve_with_wall")
            if error_counter <= 3:
                # Reset factor and reduce factor_increase
                wall_params["factor"] /= factor_increase
                if factor_increase > 1.2:
                    factor_increase *= 0.9
            else:
                print("No solution found. Leaving solve_with_wall()")
                failed = True

    if failed:
        failure_type = classify_failure(last_error_msg)
        raise SolveFailure(failure_type,
            f"Failed solve_with_wall() at factor={wall_params['factor']:.1f}: "
            f"{last_error_msg}")

# Flame settings                                                            
reaction_mechanism = "h2o2.yaml"                                            
gas = ct.Solution(reaction_mechanism)                                       
width = 18e-3                                                               
grid = np.linspace(0, width, 250)                                           
f = ct.CounterflowDiffusionFlame(gas, grid=grid)                            
f.P = 1.e5                                                                  
f.fuel_inlet.mdot = 0.1 
f.oxidizer_inlet.mdot = 0.5 
f.fuel_inlet.X = "H2:1"                                                     
f.oxidizer_inlet.X = "O2:1"                                                 
f.fuel_inlet.T = 300                                                        
f.oxidizer_inlet.T = 300                                                    
f.transport_model = "unity-Lewis-number"                            
f.set_refine_criteria(ratio=3, slope=0.5, curve=0.05, prune=0.04,                 
                    enthalpy=False, enthalpy_curve=0.05) 
# Wall                                                                      
wall_params = {                                                             
'Z_wall': 0.3,                                                                
'T_wall': 300.0,                                                            
'factor': 1000,                                                                
'mix_frac': 'H',                                                       
'fuel': 'H2',                                                               
'oxidizer': 'O2',                                                           
'basis': 'mass'                                                             
}                                                                           
 
z_stoich = 0.111 
output_path = Path("Scripts/Data")
file_path = str(output_path / "stable_test03no2.h5")
csv_path = str(output_path / "stable_test03no2.csv")
fig_path = str(output_path / "stable_test03no2.png")
# Names                                                                     
name = "initial"                                                            
names = [name,]  



# Define a limit for the maximum temperature below which the flame is
# considered as extinguished and the computation is aborted
temperature_limit_extinction = max(f.oxidizer_inlet.T, f.fuel_inlet.T)

# Initialize and solve
#f.solve(loglevel=0, refine_grid=True, auto=True)
f.set_initial_guess(data="Scripts/Data/enthalpy.h5", group="initial")

#factor_last_working = 0 
#factor_last_working = solve_with_wall(f, wall_params, name_fallback=names[-1],factor_last_working=factor_last_working, delta_T_max=1.0, loglevel=0)

save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)

# %%
# Compute Extinction Strain Rate
# ------------------------------
# Exponents for the initial solution variation with changes in strain rate
# Taken from Fiala and Sattelmayer (2014)
exp_d_a = - 1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

# Set normalized initial strain rate
alpha = [1.]
# Initial relative strain rate increase
delta_alpha = 5.
# Factor of refinement of the strain rate increase
delta_alpha_factor = 50.
# Limit of the refinement: Minimum normalized strain rate increase
delta_alpha_min = .001
# Limit of the Temperature decrease
delta_T_min = 1  # K

# Iteration indicator
n = 0
# Indicator of the latest flame still burning
n_last_burning = 0
# List of peak temperatures
T_max = [np.max(f.T)]
# List of maximum axial velocity gradients
a_max = [np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))]

# %%
# Simulate counterflow flames at increasing strain rates until the flame is
# extinguished. To achieve a fast simulation, an initial coarse strain rate
# increase is set. This increase is reduced after an extinction event and
# the simulation is again started based on the last burning solution.
# The extinction point is considered to be reached if the abortion criteria
# on strain rate increase and peak temperature decrease are fulfilled.
while True:
    n += 1
    # Update relative strain rates
    alpha.append(alpha[n_last_burning] + delta_alpha)
    strain_factor = alpha[-1] / alpha[n_last_burning]
    # Create an initial guess based on the previous solution
    # Update grid
    # Note that grid scaling changes the diffusion flame width
    f.flame.grid *= strain_factor ** exp_d_a
    # Update mass fluxes
    f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
    f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
    # Update velocities
    f.flame.set_values("velocity", f.flame.velocity * strain_factor ** exp_u_a)
    f.flame.set_values("spreadRate", f.flame.spread_rate * strain_factor ** exp_V_a)
    # Update pressure curvature
    f.flame.set_values(
        "Lambda", f.flame.radial_pressure_gradient * strain_factor ** exp_lam_a)

    solved = False
    failure_type = None
    try:
        factor_last_working = 1000
        solve_with_wall(f, wall_params, name_fallback=names[-1],
                        factor_last_working=factor_last_working,
                        delta_T_max=1.0, loglevel=0)
        solved = True
    except (SolveFailure, ct.CanteraError) as e:
        failure_type = getattr(e, 'failure_type', classify_failure(str(e)))
        print(f"Error: Did not converge at n = {n}, type={failure_type}")

    if solved:
        t_max_val = float(np.max(f.T))
        a_max_val = float(np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid))))
        T_max.append(t_max_val)
        a_max.append(a_max_val)
        print(f"MAX Temp {t_max_val}")

        if not np.isclose(t_max_val, temperature_limit_extinction):
            # Flame is still burning, so proceed to next strain rate
            n_last_burning = n
            name = f"extinction/{n:04d}"
            names.append(name)
            save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)
            print(f'Flame burning at alpha = {alpha[-1]:8.4f}. Proceeding, '
                  f'delta_alpha = {delta_alpha}')
            continue
        else:
            print(f"Flame extinguished (solved, T~ambient) at alpha = {alpha[-1]:8.4f}")
    else:
        # Solver failed — f is in undefined state; do NOT use it for T/a_max
        print(f"Flame extinguished (solver failed: {failure_type}) at alpha = {alpha[-1]:8.4f}")

    # --- Extinction handling (arrive here from both branches) ---
    if delta_alpha < delta_alpha_min:
        # Converged: step is refined to minimum, flame is out
        name = f"extinction/{n:04d}"
        names.append(name)
        if solved:
            save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)
        else:
            print("  (not saving — solver failed, no valid solution)")
        print(f'Flame extinguished at alpha = {alpha[-1]:8.4f}. '
              'Abortion criterion satisfied.')
        break
    else:
        # Refine step and retry from last burning solution
        delta_alpha = delta_alpha / delta_alpha_factor
        alpha.pop()  # discard the unphysical alpha entry
        n = n-1
        print(f'Flame extinguished at alpha = {alpha[-1]:8.4f} (discarded). '
              f'Restoring alpha = {alpha[n_last_burning]:8.4f} and '
              f'trying delta_alpha = {delta_alpha}')
        print(f"Setting initial guess {names[-1]}")
        f.set_initial_guess(data=file_path, group=names[-1])


# %%
# Results
# -------
# Print some parameters at the extinction point, after restoring the last burning
# solution
name = f"extinction/{n_last_burning:04d}"
f.restore(file_path, name)

print('----------------------------------------------------------------------')
print('Parameters at the extinction point:')
print('Pressure p={0} bar'.format(f.P / 1e5))
print('Peak temperature T={0:4.0f} K'.format(np.max(f.T)))
print('Mean axial strain rate a_mean={0:.2e} 1/s'.format(f.strain_rate('mean')))
print('Maximum axial strain rate a_max={0:.2e} 1/s'.format(f.strain_rate('max')))
print('Fuel inlet potential flow axial strain rate a_fuel={0:.2e} 1/s'.format(
      f.strain_rate('potential_flow_fuel')))
print('Oxidizer inlet potential flow axial strain rate a_ox={0:.2e} 1/s'.format(
      f.strain_rate('potential_flow_oxidizer')))
print('Axial strain rate at stoichiometric surface a_stoich={0:.2e} 1/s'.format(
      f.strain_rate('stoichiometric', fuel='H2')))

# %%
# Plot the maximum temperature over the maximum axial velocity gradient
plt.figure()
plt.semilogx(a_max, T_max)
plt.xlabel(r'$a_{max}$ [1/s]')
plt.ylabel(r'$T_{max}$ [K]')
plt.savefig(output_path / "figure_T_max_a_max.png")
plt.show()

