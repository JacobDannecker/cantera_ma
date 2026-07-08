import numpy as np                                                              
import cantera as ct                                                            
from scipy import special                                                       
from matplotlib import pyplot as plt                                            
import time                                                                     
import h5py
import logging
import sys
import pandas as pd
import utils as ut

# Flame settings                                                            
reaction_mechanism = "h2o2.yaml"                                            
gas = ct.Solution(reaction_mechanism)                                       
width = 18e-3                                                               
grid = np.linspace(0, width, 300)                                           
f = ct.CounterflowDiffusionFlame(gas, grid=grid)                            
f.P = 1.e5                                                                  
f.fuel_inlet.X = "H2:1"
f.oxidizer_inlet.X = "O2:1"
f.fuel_inlet.T = 300
f.oxidizer_inlet.T = 300

# Save data in:
file_path = "Data/unstable_temp_change_3_lower_step_Z0.5000.h5"


# Initialize
start_file = "Data/stable_Z0.5000.h5"
h5_file = h5py.File(start_file, "r")
keys = ["extinction/" + key for key in h5_file["extinction"].keys()]

last_run = keys[-1]
first_run = keys[0]

# Restore first solution
f.restore(start_file, first_run)

a_max = strain_rate = f.strain_rate("max")

wall_params = {
    'Z_wall': h5_file[last_run]["flame"]["z"].attrs["z_wall"],
    'T_wall': 300.0,
    'factor': h5_file[last_run]["flame"]["z"].attrs["factor"],
    'mix_frac': 'H',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}

z_stoich = ut.get_z_stoich(gas, wall_params, reaction_mechanism)

# Restore last solution
f.restore(start_file, last_run)

f.set_refine_criteria(ratio=3.0, slope=0.5, curve=0.05, prune=0.01, enthalpy=False, enthalpy_curve=0.5)
f.transport_model = "unity-Lewis-number"

# Save starting solution in new file
# ut.save_with_attributes(f, file_path, "initial", wall_params, z_stoich, info=True)


##############################################################################
# Flame Continuation
trapezoid = getattr(np, "trapezoid", None) or np.trapz
# Maximum number of steps to take
n_max = 500

# Relative temperature defining control point locations, with 1 being the peak
# temperature and 0 being the inlet temperature. Lower values tend to avoid solver
# failures early on, while using higher values on the unstable branch tend to help
# with finding solutions where the peak temperature is very low.
initial_spacing = 0.6
unstable_spacing =  0.95

# Amount to adjust temperature at the control point each step [K]
temperature_increment = 50
max_increment = 100

# Try to keep T_max from changing more than this much each step [K]
target_delta_T_max = 10.

# Stop after this many successive errors
max_error_count = 5
error_count = 0

# Stop after any failure if the strain rate has dropped to this fraction of the maximum
strain_rate_tol = 0.1

f.two_point_control_enabled = True
#f.set_max_jac_age(1, 1)
# Prevent two point control from finding solutions with negative inlet velocities
f.flame.set_bounds(spread_rate=(-1e-5, 1e20))
f.max_time_step_count = 1000 
T_max = max(f.T)
#a_max = strain_rate = f.strain_rate('max')
print("tmax: ", T_max)
data = []  # integral output quantities for each step
solved = False

factor_last_working = wall_params["factor"]
refine_grid = False
auto = False

for i in range(n_max):
    print(f"i = {i}")
    if strain_rate > 0.98 * a_max:
        spacing = initial_spacing
    else:
        spacing = unstable_spacing

    backup_state = f.to_array()
    
    control_temperature = np.min(f.T) + spacing*(np.max(f.T) - np.min(f.T))

    # Store the flame state in case the iteration fails and we need to roll back
    backup_state = f.to_array()

    print("In starin loop:")
    print(f"T_max: {np.max(f.T)}")
    print(f'Iteration {i}: Control temperature = {control_temperature:.2f} K')

    f.set_left_control_point(control_temperature)
    f.set_right_control_point(control_temperature)

    print(f"Temperature increment = {temperature_increment}")

    ## This decrement is what drives the two-point control. If failure
    ## occurs, try decreasing the decrement.
    if i == 1:
        f.right_control_point_temperature -= 1.
        f.left_control_point_temperature -= 1.
    else:
        f.right_control_point_temperature -= temperature_increment
        f.left_control_point_temperature -= temperature_increment
    f.clear_stats()

    if (f.left_control_point_temperature < f.fuel_inlet.T + 100
        or f.right_control_point_temperature < f.oxidizer_inlet.T + 100
    ):
        print("SUCCESS! Stopping because control point temperature is "
                    "sufficiently close to inlet temperature.")
        break

    try:
        runtime, factor_last_working = ut.solve_with_wall(f, wall_params, factor_last_working=factor_last_working, 
                                                 factor_increase=1.1, factor_decrease=0.9, delta_T_max=1.0, loglevel=1, refine_grid=refine_grid, auto=auto)

        print(f"Factor last working: {factor_last_working}")
        if abs(max(f.T) - T_max) < 0.8 * target_delta_T_max:
            # Max temperature is changing slowly. Try a larger increment next step
            temperature_increment = min(temperature_increment + 3, max_increment)
        elif abs(max(f.T) - T_max) > target_delta_T_max:
            # Max temperature is changing quickly. Scale down increment for next step
            temperature_increment *= 0.9 * target_delta_T_max / (abs(max(f.T) - T_max))
        solved = True
        error_count = 0
    except (ut.SolveFailure, ct.CanteraError) as err:
        print(err)
        if strain_rate / a_max < strain_rate_tol:
            print('SUCCESS! Traversed unstable branch down to '
                        f'{100 * strain_rate / a_max:.2f}% of the maximum strain rate.')
            break
        # Restore the previous solution and try a smaller temperature increment for the
        # next iteration
        f.from_array(backup_state)
        temperature_increment = 0.7 * temperature_increment
        error_count += 1
        print(f"Did not converge on iteration {i}. New dT = {temperature_increment:.2f}")


    print(f"Saving on iteration{i}")
    name = f"iteration{i}"
    ut.save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)

    # Collect output stats
    T_max = max(f.T)
    T_mid = 0.5 * (min(f.T) + max(f.T))
    s = np.where(f.T > T_mid)[0][[0,-1]]
    width = f.grid[s[1]] - f.grid[s[0]]
    strain_rate = f.strain_rate('max')
    a_max = max(strain_rate, a_max)

    data.append({
        'T_max': max(f.T),
        'strain_rate': strain_rate,
        'heat_release_rate': trapezoid(f.heat_release_rate, f.grid),
        'n_points': len(f.grid),
        'flame_width': width,
        'Tc_increment': temperature_increment,
        'time_steps': sum(f.time_step_stats),
        'eval_count': sum(f.eval_count_stats),
        'cpu_time': sum(f.jacobian_time_stats + f.eval_time_stats),
        'errors': error_count
    })
    
    df = pd.DataFrame.from_records(data)

    if error_count >= max_error_count:
        print(f'FAILURE! Stopping after {error_count} successive solver '
                       'errors.')
        break

print(f'Stopped after {i} iterations')

plt.figure()
plt.plot(df.strain_rate, df.T_max)
plt.xlabel('a_max')
plt.ylabel('T_max')
plt.show()

                                
