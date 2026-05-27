import numpy as np                                                              
import cantera as ct                                                            
from scipy import special                                                       
from matplotlib import pyplot as plt                                            
import time                                                                     
import h5py
import logging
import sys
import pandas as pd 

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logging.basicConfig(stream=sys.stdout)

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
    hdf5_file[group_z].attrs["mix_frac"] = wall_params["mix_frac"]
    hdf5_file[group_z].attrs["fuel"] = wall_params["fuel"]
    hdf5_file[group_z].attrs["oxidizer"] = wall_params["oxidizer"]
    hdf5_file[group_z].attrs["basis"] = wall_params["basis"]
    hdf5_file[group_z].attrs["chi_st"] = chi_stoich(f, z_stoich)
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
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)          

# Wall                                                                      
wall_params = {                                                             
'Z_wall': 0.5,                                                                
'T_wall': 300.0,                                                            
'factor': 1,                                                                
'mix_frac': 'Bilger',                                                       
'fuel': 'H2',                                                               
'oxidizer': 'O2',                                                           
'basis': 'mass'                                                             
}                                                                           
 
z_stoich = get_z_stoich(gas, wall_params, reaction_mechanism)
file_path = f"Scripts/Data/z_wall_050.h5"                     
csv_path = f"Scripts/Data/z_wall_050.csv"
fig_path = f"Scripts/Data/z_wall_050.png"
# Names                                                                     
name = "initial"                                                            
names = [name,]                                                             

# Initial Solution
logger.info('Creating the initial solution')
f.solve(loglevel=0, refine_grid=True)                               
f.save(file_path, name=name, overwrite=True)
solve_with_wall(f, wall_params, delta_T_max=1.0, loglevel=0)
save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)

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
temperature_increment = 50.0
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

f.max_time_step_count = 500
T_max = max(f.T)
a_max = strain_rate = f.strain_rate('max')
data = []  # integral output quantities for each step
logger.info('Starting two-point control')
factor_last_working = 0

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
    
    logger.debug(f'Iteration {i}: Control temperature = {control_temperature:.2f} K')
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
        logger.info("SUCCESS! Stopping because control point temperature is "
                    "sufficiently close to inlet temperature.")
        break

    try:
        #try:
        factor_last_working = solve_with_wall(f, wall_params, name_fallback=names[-1],factor_last_working=factor_last_working, delta_T_max=1.0, loglevel=0)
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
        logger.debug(err)
        if strain_rate / a_max < strain_rate_tol:
            logger.info('SUCCESS! Traversed unstable branch down to '
                        f'{100 * strain_rate / a_max:.2f}% of the maximum strain rate.')
            break

        # Restore the previous solution and try a smaller temperature increment for the
        # next iteration
        factor_last_working = 100
        f.restore(file_path, name=names[-1])
        print(f"Restor solution: {names[-1]} ")
        #f.from_array(backup_state)
        temperature_increment = 0.5 * temperature_increment
        error_count += 1
        logger.warning(f"Solver did not converge on iteration {i}. Trying again with "
                       f"dT = {temperature_increment:.2f}")

    if ct.hdf_support():
        name = f"iteration{i}"
        names.append(name)
        save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)

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
        logger.warning(f'FAILURE! Stopping after {error_count} successive solver '
                       'errors.')
        break

logger.info(f'Stopped after {i} iterations')




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
