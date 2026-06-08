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
        print("\n############################################################")
        print(f"Solved at z_wall: {z_wall}")
        print("##############################################################")
    f.save(file_path, name=name, overwrite=True)
    add_attributes(f, file_path, name, wall_params, z_stoich)


def solve_with_wall(f, wall_params, name_fallback="initial", delta_T_max=1., factor_last_working=False, factor_increase=2, loglevel=0):
    error_counter = 0
    z_wall = wall_params["Z_wall"]
    delta_T_ok = False
    failed = False
    wall_params["factor"] = 100.0
    max_factor = 1e25
    set_factor = False
    while not delta_T_ok and not failed:
        try:
            if factor_last_working and not set_factor:
                wall_params["factor"] = factor_last_working
                set_factor = True
            else:
                wall_params["factor"] = min(
                    wall_params["factor"] * factor_increase, max_factor
                )
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=loglevel, refine_grid=True)                           
            delta_T_wall = get_delta_T(f, wall_params)
            if delta_T_wall < delta_T_max:
                print(f"mdot f, o : {f.fuel_inlet.mdot}, {f.oxidizer_inlet.mdot}")
                print(f"Strain max: {f.strain_rate("max")}")
                delta_T_ok = True
                print(f"Solution found at delta_T_wall: {delta_T_wall}")
                return factor_last_working


        except ct.CanteraError as err:      
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
                raise ct.CanteraError("Failed solve_with_wall()")

           # if  error_counter == 4:
           #     # Start from initial solution wit wall
           #     # reset factor and factor_increase
           #     print("Try with initial solution as ininital guess.")       
           #     f.set_initial_guess(data=file_path, group=name_fallback)         
           #     #f.from_array(name_fallback)
           #     wall_params["factor"] = 10                                  
           #     factor_increase = 2                                         

           # if error_counter > 4:
           #     print("Reached max Errors!!!!")
           #     # Reset factor and reduce factor_increase
           #     wall_params["factor"] /= factor_increase                       
           #     if factor_increase > 1.2:                                       
           #         factor_increase *= 0.9                                      

        #if error_counter > 5:
        #    print("No solution found. error_counter: {error_counter}")
        #    failed = True 
        #    raise ct.CanteraError("No solution with wall found!!!!!!!!!!!!")
         

# Flame settings                                                            
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
f.set_refine_criteria(ratio=3, slope=0.5, curve=0.5, prune=0.0,                 
                    enthalpy=True, enthalpy_curve=0.1) 
# Wall                                                                      
wall_params = {                                                             
'Z_wall': 0.8,                                                                
'T_wall': 300.0,                                                            
'factor': 10000,                                                                
'mix_frac': 'Bilger',                                                       
'fuel': 'H2',                                                               
'oxidizer': 'O2',                                                           
'basis': 'mass'                                                             
}                                                                           
 
z_stoich = 0.111 
file_path = f"Scripts/Data/stable_080.h5"                     
csv_path = f"Scripts/Data/stable_080.csv"
fig_path = f"Scripts/Data/stable_080.png"
# Names                                                                     
name = "initial"                                                            
names = [name,]  



# Define a limit for the maximum temperature below which the flame is
# considered as extinguished and the computation is aborted
temperature_limit_extinction = max(f.oxidizer_inlet.T, f.fuel_inlet.T)

# Initialize and solve
print('Creating the initial solution')

#f.solve(loglevel=0, auto=True)
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
delta_alpha = 1.
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

    try:                                                                        
        factor_last_working = 1000
        solve_with_wall(f, wall_params, name_fallback="initial",factor_last_working=factor_last_working, delta_T_max=1.0, loglevel=0)
    except ct.CanteraError as e:
        print('Error: Did not converge at n =', n, e)

    T_max.append(np.max(f.T))
    a_max.append(np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid))))
    print(f"MAX Temp {np.max(f.T)}")
    if not np.isclose(np.max(f.T), temperature_limit_extinction):
        # Flame is still burning, so proceed to next strain rate
        n_last_burning = n
        name = f"extinction/{n:04d}" 
        names.append(name)
        save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)
        print('Flame burning at alpha = {:8.4F}. Proceeding to the next iteration, '
              'with delta_alpha = {}'.format(alpha[-1], delta_alpha))
    elif ((T_max[-2] - T_max[-1] < delta_T_min) and (delta_alpha < delta_alpha_min)):
        # If the temperature difference is too small and the minimum relative
        # strain rate increase is reached, save the last, non-burning, solution
        # to the output file and break the loop
        name = f"extinction/{n:04d}" 
        names.append(name)
        save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)

        print('Flame extinguished at alpha = {0:8.4F}.'.format(alpha[-1]),
              'Abortion criterion satisfied.')
        break
    else:
        # Procedure if flame extinguished but abortion criterion is not satisfied
        # Reduce relative strain rate increase
        delta_alpha = delta_alpha / delta_alpha_factor

        print('Flame extinguished at alpha = {0:8.4F}. Restoring alpha = {1:8.4F} and '
              'trying delta_alpha = {2}'.format(
                  alpha[-1], alpha[n_last_burning], delta_alpha))

        # Restore last burning solution
        f.set_initial_guess(data="Scripts/Data/enthalpy.h5", group="initial")

# %%
# Results
# -------
# Print some parameters at the extinction point, after restoring the last burning
# solution
name = names(f"extinction/{n_last_burning:04d}")
f.restore(file_name, entry)

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

