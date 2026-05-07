import numpy as np                                                              
import cantera as ct                                                            
from scipy import special                                                       
from matplotlib import pyplot as plt                                            
import time                                                                     
import h5py

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

def solve_wiht_wall(f, wall_params, name_fallback="initial", delta_T_max=1., factor_increase=2, loglevel=0):
    error_count = 0
    z_wall = wall_params["Z_wall"]
    delta_T_ok = False
    try:                                                                    
        # Factor 1 means residual_T = (T(x,j) - T_wall)^4
        wall_params["factor"] = 1                                           
        f.set_initial_guess(data=file_path, group=names[0])                 
        f.flame.set_non_adiabatic_wall(wall_params)                         
        print(f"Solving for z_wall: {z_wall} with factor = 1")                        
        f.solve(loglevel=loglevel, fefine_grid=True)                               
        print(f"Success with factor = 1 at z_wall: {z_wall}")    
    except:                                                                 
        print(f"Failed initial solution with factor = 1 at z_wall: {z_wall}")
        print(f"Conitune with increasing factor.")
        wall_params["factor"] = 100                                             
    while not delta_T_ok:                                                             
        try:                                                                
            wall_params["factor"] *= factor_increase                        
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=loglevel, refine_grid=True)                           
            delta_T_wall = get_delta_T(f, wall_params)
            if delta_T_wall < delta_T_max:
                delta_T_ok = True
            print(f"Solution found at delta_T_wall: {delta_T_wall}")
        except BaseException as err:      
            print(err)                                                      
            error_count += 1
            if error_count <= 3:
                # Reset factor and reduce factor_increase
                wall_params["factor"] /= factor_increase                       
                if factor_increase > 1.2:                                       
                    error_recovery_code = 0 
                    factor_increase *= 0.9                                      

            if  error_count == 4:
                # Start from initial solution wit wall
                # reset factor and factor_increase
                print("Try with initial solution as ininital guess.")       
                f.set_initial_guess(data=file_path, group=name_fallback)         
                wall_params["factor"] = 10                                  
                factor_increase = 2                                         

            if error_count > 4:
                # Reset factor and reduce factor_increase
                wall_params["factor"] /= factor_increase                       
                if factor_increase > 1.2:                                       
                    error_recovery_code = 0 
                    factor_increase *= 0.9                                      

            if error_recovery_code > 5:
                print("No solution found. error_count: {error_count}")
                break
                

         

# Flame settings                                                            
reaction_mechanism = "h2o2.yaml"                                            
gas = ct.Solution(reaction_mechanism)                                       
width = 18e-3                                                               
grid = np.linspace(0, width, 150)                                           
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
'Z_wall': 0.8,                                                                
'T_wall': 300.0,                                                            
'factor': 1,                                                                
'mix_frac': 'Bilger',                                                       
'fuel': 'H2',                                                               
'oxidizer': 'O2',                                                           
'basis': 'mass'                                                             
}                                                                           
 
z_stoich = get_z_stoich(gas, wall_params, reaction_mechanism)

# Names                                                                     
name = "z_wall_0.80"                                                            
names = ["initial", name,]                                                             
file_path = f"Scripts/Data/test.h5"                     

f.solve(loglevel=0, refine_grid=True)                               
f.save(file_path, name="initial", overwrite=True)                                  

solve_wiht_wall(f, wall_params, delta_T_max=0.1, loglevel=0)
save_with_attributes(f, file_path, name, wall_params, z_stoich, info=True)



#Plot                                                                  
fig, ax = plt.subplots(2, 1)                                            
fig.suptitle("H2/O2")                                                   
for name in names:                                                      
    f.restore(file_path, name=name)                                     
    if name == "initial":
        z = 1.0
    else:
        z = float(name.strip("z_wall_"))                                    

    # Subplot 1 Temperature
    ax[0].plot(f.mixture_fraction("H"), f.T, label=f"{z:.2f} chi_st: {chi_stoich(f, z_stoich):.2f}")
    # Subplot 2  enthalpy                                                   
    ax[1].plot(f.mixture_fraction("H"), f.h, label=f"{z:.2f} chi_st: {chi_stoich(f, z_stoich):.2f}")
                                                                            
    ax[0].grid()                                                            
    ax[0].legend()                                                          
    ax[0].set_ylabel("T")                                                   
    ax[0].set_xlabel("<- ox z fuel ->")                                     
                                                                            
    ax[1].grid()                                                            
    ax[1].legend()                                                          
    ax[1].set_ylabel("h in  J/kg")                                          
    ax[1].set_xlabel("<- ox z fuel ->")                                     
                                                                        
plt.show()                                 
