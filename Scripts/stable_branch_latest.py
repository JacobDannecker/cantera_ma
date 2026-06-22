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
import utils as ut



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

#f.max_time_step_count = 500

# Wall                                                                      
wall_params = {                                                             
'z_wall': 0.3,                                                                
'T_wall': 300.0,                                                            
'factor': 1000,                                                                
'mix_frac': 'H',                                                       
'fuel': 'H2',                                                               
'oxidizer': 'O2',                                                           
'basis': 'mass'                                                             
}                                                                           
 


z_stoich = ut.get_z_stoich(gas, wall_params, reaction_mechanism)
output_path = Path("Scripts/Data")
file_path = str(output_path / "standaloneNo4.h5")
csv_path = str(output_path / "standaloneNo4.csv")
fig_path = str(output_path / "standaloneNo4.png")
# Names                                                                     
name = "initial"
names = [name,]  

# Define a limit for the maximum temperature below which the flame is
# considered as extinguished and the computation is aborted
temperature_limit_extinction = max(f.oxidizer_inlet.T, f.fuel_inlet.T)

f.set_initial_guess(data="Scripts/Data/enthalpy.h5", group="initial")
ut.save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)      
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
delta_alpha = 5 
# Factor of refinement of the strain rate increase
delta_alpha_factor = 50.
# Limit of the refinement: Minimum normalized strain rate increase
delta_alpha_min = .001
# Limit of the Temperature decrease

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
        wall_params["factor"] = 1000
        solved = ut.solve_with_wall(f, wall_params, delta_T_max=1.0, loglevel=0)
    except (ut.SolveFailure, ct.CanteraError) as e:
        print(e)
        failure_type = getattr(e, 'failure_type', ut.classify_failure(str(e)))
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
            name = f"{n:04d}"
            names.append(name)
            ut.save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)
            print(f'Flame burning at alpha = {alpha[n_last_burning]:8.4f}. Proceeding, '
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
        name = f"{n:04d}"
        names.append(name)
        if solved:
            ut.save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)
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
name = f"{n_last_burning:04d}"
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

