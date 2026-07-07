import numpy as np                                                              
import cantera as ct                                                            
from scipy import special                                                       
from matplotlib import pyplot as plt                                            
import time                                                                     
import h5py
import logging
import sys
import pandas as pd
import utils_unstable as ut


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
f.transport_model = "unity-Lewis-number"

# Save data in:
file_path = "Scripts/Data/unstable_Run31_Z0.5000.h5"
csv_path = "Scripts/Data/unstable_Run31_Z0.5000.csv"
fig_path = "Scripts/Data/unstable_Run31_Z0.5000.png"
# Initialize flame from enthalpy guess and solve with wall


start_file = "Scripts/Data/Run30_final_Z0.5000.h5"
h5_file = h5py.File(start_file, "r")
keys = ["extinction/" + key for key in h5_file["extinction"].keys()]
second_last_run = keys[-1]
first_run = keys[0]

f.restore(start_file, first_run)
a_max = strain_rate = f.strain_rate("max")

f.fuel_inlet.mdot = h5_file[second_last_run]["fuel_inlet"].attrs["mass-flux"]
f.oxidizer_inlet.mdot = h5_file[second_last_run]["oxidizer_inlet"].attrs["mass-flux"]

z_wall = h5_file[second_last_run]["flame"]["z"].attrs["z_wall"]

wall_params = {
    'Z_wall': z_wall,
    'T_wall': 300.0,
    'factor': 1000,
    'mix_frac': 'H',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}

z_stoich = ut.get_z_stoich(gas, wall_params, reaction_mechanism)

f.restore(start_file, second_last_run)
f.set_initial_guess(start_file, second_last_run)
f.set_refine_criteria(ratio=3.0, slope=0.5, curve=0.05, prune=0.04, enthalpy=False, enthalpy_curve=0.5)

ut.save_with_attributes(f, file_path, "initial", wall_params, z_stoich, info=True)


names = ["initial"]


##############################################################################
# Flame Continuation
trapezoid = getattr(np, "trapezoid", None) or np.trapz
# Maximum number of steps to take
n_max = 500

# Relative temperature defining control point locations, with 1 being the peak
# temperature and 0 being the inlet temperature. Lower values tend to avoid solver
# failures early on, while using higher values on the unstable branch tend to help
# with finding solutions where the peak temperature is very low.
initial_spacing = 0.9
unstable_spacing =  0.5

# Amount to adjust temperature at the control point each step [K]
temperature_increment = 20
max_increment = 100

# Try to keep T_max from changing more than this much each step [K]
target_delta_T_max = 100

# Stop after this many successive errors
max_error_count = 5
error_count = 0

# Stop after any failure if the strain rate has dropped to this fraction of the maximum
strain_rate_tol = 0.001

f.two_point_control_enabled = True

# Prevent two point control from finding solutions with negative inlet velocities
f.flame.set_bounds(spread_rate=(-1e-5, 1e20))
f.max_time_step_count = 1000 
T_max = max(f.T)
#a_max = strain_rate = f.strain_rate('max')
data = []  # integral output quantities for each step
solved = False
factor_last_working = 7730941132800.0
refine_grid = True

for i in range(n_max):
    print(f"i = {i}")
    #if strain_rate > 0.98 * a_max:
    #    spacing = initial_spacing
    #else:
    spacing = unstable_spacing
    control_temperature = np.min(f.T) + spacing*(np.max(f.T) - np.min(f.T))

    # Store the flame state in case the iteration fails and we need to roll back
    backup_state = f.to_array()
    print("In starin loop:")
    print(f"T_max: {np.max(f.T)}")
    
    print(f'Iteration {i}: Control temperature = {control_temperature:.2f} K')
    f.set_left_control_point(control_temperature)
    f.set_right_control_point(control_temperature)
    print(f"Temperature increment = {temperature_increment}")
    #print(control_temperature, f.left_control_point_temperature, f.right_control_point_temperature)
    ## This decrement is what drives the two-point control. If failure
    ## occurs, try decreasing the decrement.
    #delta_left =  f.left_control_point_temperature - control_temperature 
    #delta_right = f.right_control_point_temperature - control_temperature 
    #if delta_left > 0:
    #    f.left_control_point_temperature -= (temperature_increment*0.5 + delta_left)
    #else: 
    #    f.left_control_point_temperature -= temperature_increment

    #if delta_right > 0:
    #    f.right_control_point_temperature -= (temperature_increment*0.5 + delta_right)
    #else: 
    #    f.right_control_point_temperature -= temperature_increment
    #
    #print(control_temperature, f.left_control_point_temperature, f.right_control_point_temperature)
   
    #idx = np.argmin(f.mixture_fraction("H") - 0.5)
    #T_at_wall = f.T[idx]
    #if control_temperature < (T_at_wall + 50):
    #    print("sosososososos!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

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
                                                 factor_increase=0.8, delta_T_max=1.0, loglevel=1, refine_grid=refine_grid, auto=False)
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
        if names:
            f.set_initial_guess(file_path, group=names[-1])
            print(f"Restor solution: {names[-1]} ")
        #refine_grid = False
        temperature_increment = 0.5 * temperature_increment
        error_count += 1
        print(f"Solver did not converge on iteration {i}. Trying again with "
                       f"dT = {temperature_increment:.2f}")

    if solved:
        solved = False
        #refine_grid = True
        spacing = 0.95 
        name = f"iteration{i}"
        names.append(name)
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
    df.to_csv(csv_path)


    if error_count >= max_error_count:
        print(f'FAILURE! Stopping after {error_count} successive solver '
                       'errors.')
        break

print(f'Stopped after {i} iterations')




df = pd.DataFrame.from_records(data)
df.to_csv(csv_path)

plt.figure()
plt.semilogx(df.strain_rate, df.T_max)
plt.xlabel('Maximum Axial Velocity Gradient [1/s]')
plt.ylabel('Maximum Temperature [K]')
plt.savefig(fig_path)




#Plot                                                                  
fig, ax = plt.subplots(2, 1)                                            
fig.suptitle("H2/O2")                                                   

fig2, ax2 = plt.subplots(3, 1)                                            
fig2.suptitle("H2/O2")                                                   


idx_H2 = f.gas.species_index("H2")                                              
idx_O2 = f.gas.species_index("O2")
idx_OH = f.gas.species_index("OH")

for name in names:                                                      
    if not name == "initial":
        f.restore(file_path, name=name)                                     
        label = name
        # Subplot 1 Temperature
        ax[0].plot(f.mixture_fraction("H"), f.T, label=label)
        # Subplot 2  enthalpy                                                   
        ax[1].plot(f.mixture_fraction("H"), f.h, label=label)
        # Sublots species
        ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label=label)
        ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], label=label)
        ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_OH], label=label)

for a in ax:
    a.grid()                                                            
    a.legend()                                                          
    a.set_xlabel("<- ox z fuel ->")                                     

for a in ax2:
    a.grid()                                                            
    a.legend()                                                          
    a.set_xlabel("<- ox z fuel ->")                                     
                                                                            
ax[0].set_ylabel("T")                                                   
ax[1].set_ylabel("h in  J/kg")                                          

ax2[0].set_ylabel("H2")                                                   
ax2[1].set_ylabel("O2")                                                   
ax2[2].set_ylabel("OH")                                                   

plt.tight_layout()
plt.show()                                 
