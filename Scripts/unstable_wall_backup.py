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
grid = np.linspace(0, width, 250)                                           
f = ct.CounterflowDiffusionFlame(gas, grid=grid)                            
f.P = 1.e5                                                                  
f.fuel_inlet.X = "H2:1"                                                     
f.oxidizer_inlet.X = "O2:1"                                                 
f.fuel_inlet.T = 300                                                        
f.oxidizer_inlet.T = 300                                                    
f.transport_model = "unity-Lewis-number"                            
f.set_refine_criteria(ratio=3.0, slope=0.5, curve=0.5, prune=0.03)          

# Save data in:
file_path = f"Scripts/Data/unstable_utilNo2.h5"                     
csv_path = f"Scripts/Data/unstable_utilNo2.csv"
fig_path = f"Scripts/Data/unstable_utilNo2.png"


# Start from last stable burning flamelet
load_from_file = "Scripts/DataPlot/z06.h5"
h5_file = h5py.File(load_from_file, "r")                                              
last_run = [str(name) for name in h5_file.keys()][-1]
print(last_run)
z_wall = h5_file[last_run]["flame"]["z"].attrs["z_wall"]                     
f.set_initial_guess(data=load_from_file, group=last_run)           

# Wall                                                                      
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
# Names of runs list                                                                   
names = []                                                             
# Initial Solution
#f.solve(loglevel=0, refine_grid=True)                               
#f.save(file_path, name=name, overwrite=True)
#solve_with_wall(f, wall_params, delta_T_max=1.0, loglevel=0)
#ut.save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)


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
unstable_spacing =  0.95

# Amount to adjust temperature at the control point each step [K]
temperature_increment = 20.0
max_increment = 100 

# Try to keep T_max from changing more than this much each step [K]
target_delta_T_max = 20 

# Stop after this many successive errors
max_error_count = 5
error_count = 0

# Stop after any failure if the strain rate has dropped to this fraction of the maximum
strain_rate_tol = 0.10

f.two_point_control_enabled = True

# Prevent two point control from finding solutions with negative inlet velocities
f.flame.set_bounds(spread_rate=(-1e-5, 1e20))
f.max_time_step_count = 1000 
T_max = max(f.T)
a_max = strain_rate = f.strain_rate('max')
data = []  # integral output quantities for each step

for i in range(n_max):
    print(f"i = {i}")
    if strain_rate > 0.98 * a_max:
        spacing = initial_spacing
    else:
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
    print(control_temperature, f.left_control_point_temperature, f.right_control_point_temperature)
    # This decrement is what drives the two-point control. If failure
    # occurs, try decreasing the decrement.
    delta_left =  f.left_control_point_temperature - control_temperature 
    delta_right = f.right_control_point_temperature - control_temperature 
    if delta_left > 0:
        f.left_control_point_temperature -= (temperature_increment*0.5 + delta_left)
    else: 
        f.left_control_point_temperature -= temperature_increment

    if delta_right > 0:
        f.right_control_point_temperature -= (temperature_increment*0.5 + delta_right)
    else: 
        f.right_control_point_temperature -= temperature_increment
    
    print(control_temperature, f.left_control_point_temperature, f.right_control_point_temperature)
   
    f.clear_stats()

    if (f.left_control_point_temperature < f.fuel_inlet.T + 100
        or f.right_control_point_temperature < f.oxidizer_inlet.T + 100
    ):
        print("SUCCESS! Stopping because control point temperature is "
                    "sufficiently close to inlet temperature.")
        break

    try:
        #try:
        factor = 1000
        solved = ut.solve_with_wall(f, wall_params, factor_last_working=factor, delta_T_max=1.0, loglevel=0)
        #except BaseException as err:
        #    print("Try wihtout last_working_factor.=============================================")
        #    factor_last_working = solve_with_wall(f, wall_params, name_fallback=names[-1],factor_last_working=False, delta_T_max=1.0, loglevel=0)


        #print("After solve with wall")
        if abs(max(f.T) - T_max) < 0.8 * target_delta_T_max:
            # Max temperature is changing slowly. Try a larger increment next step
            temperature_increment = min(temperature_increment + 3, max_increment)
        elif abs(max(f.T) - T_max) > target_delta_T_max:
            # Max temperature is changing quickly. Scale down increment for next step
            temperature_increment *= 0.9 * target_delta_T_max / (abs(max(f.T) - T_max))
        error_count = 0
    except ct.CanteraError as err:
        print(err)
        if strain_rate / a_max < strain_rate_tol:
            print('SUCCESS! Traversed unstable branch down to '
                        f'{100 * strain_rate / a_max:.2f}% of the maximum strain rate.')
            break

        # Restore the previous solution and try a smaller temperature increment for the
        # next iteration
        factor = 100
        f.restore(file_path, name=names[-1])
        print(f"Restor solution: {names[-1]} ")
        #f.from_array(backup_state)
        temperature_increment = 0.5 * temperature_increment
        error_count += 1
        print(f"Solver did not converge on iteration {i}. Trying again with "
                       f"dT = {temperature_increment:.2f}")

    if ct.hdf_support():
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
plt.plot(df.strain_rate, df.T_max)
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
