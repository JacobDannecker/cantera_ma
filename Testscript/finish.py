import numpy as np
import cantera as ct
from scipy import special
from matplotlib import pyplot as plt
import time 
import h5py

# Flame settings
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3  
grid = np.linspace(0, width, 150)
f = ct.CounterflowDiffusionFlame(gas, width=width)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.5  
f.fuel_inlet.X = "H2:1"
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 3.0 
f.oxidizer_inlet.X = "O2:1"
f.oxidizer_inlet.T = 300 
z_stoich = 0.111
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

base_path = "Testscript/Data/precalc.h5" 
file_path = "Testscript/Data/finish.h5"

# Wall 
wall_params = {
    'Z_wall': 0,
    'T_wall': 300.0,
    'factor': 1e3,
    'mix_frac': 'Bilger',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}

# Loop over z
delta_z = 0.1
last_z_working = 1
z = 1
error_counter = 0
max_errors_reached = False
start_time = time.time()
i = 1 

base_file = h5py.File(base_path, "r")
names = list(base_file.keys())[::-1]
print(f"Names: {names}")
for name in names:
    if name == "no_wall":
        continue
    f.set_initial_guess(data=base_path, group=name) 
    z = float(name.strip("z_wall_"))
    wall_params["Z_wall"] = z 
    wall_params["factor"] *= 1e3

    while True:
        try:                                                                
            print(f"Try {name}: ") 
            wall_params["factor"] *= 5                                      
            f.flame.set_non_adiabatic_wall(wall_params)                     
            f.solve(loglevel=0, refine_grid=True)                           
            idx_wall = np.abs(f.mixture_fraction("H") - wall_params["Z_wall"]).argmin()
            delta_T_wall =  f.T[idx_wall] - wall_params["T_wall"]           
            last_z_working = z                                              
            idx_wall = np.abs(f.mixture_fraction("H") - wall_params["Z_wall"]).argmin()
            delta_T_wall =  f.T[idx_wall] - wall_params["T_wall"]           
            print(f"delta T wall: {delta_T_wall}")                          
        except BaseException as err:                                        
            error_counter += 1                                              
            print(err)                                                      
            if error_counter > 2:                                           
                wall_params["factor"] = 10e2                                
                f.set_initial_guess(data=base_path, group="no_wall")         
                print("Try with initial solution as ininital guess.")       
                error_counter = 0                                           
                if max_errors_reached:                                      
                    max_errors_reached = False                              
                    break                                                   
                max_errors_reached = True                                   
                continue                                                    
        if delta_T_wall < 1:                                                
            break          

    print("\n##############################################################################")
    print(f"Solved at z: {z}")
    print("##############################################################################")
    name = "z_wall_" + str(z)
    f.save(file_path, name=name, overwrite=True)

end_time = time.time()
print(f"Total time loop: {end_time-start_time}")

def chi_stoich(f, z_stoich):
        #a = f.strain_rate("mean")
        a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
        chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
        return chi_stoich


# Plot
fig, ax = plt.subplots(2, 1)
fig.suptitle("H2/O2") 
for name in names:
    f.restore(file_path, name=name)
    # Subplot 1 Temp 
    ax[0].plot(f.mixture_fraction("H"), f.T, label=f"With Wall chi_st: {chi_stoich(f, z_stoich):.2f}")
    # Subplot 2  enthalpy
    ax[1].plot(f.mixture_fraction("H"), f.h, label=f"With Wall chi_st: {chi_stoich(f, z_stoich):.2f}")

ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("T")
ax[0].set_xlabel("<- ox z fuel ->")

ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("h in  J/kg")
ax[1].set_xlabel("<- ox z fuel ->")

plt.show()
