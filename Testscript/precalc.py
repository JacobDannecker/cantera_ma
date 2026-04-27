import numpy as np
import cantera as ct
from scipy import special
from matplotlib import pyplot as plt
import time 

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

# Names
name = "no_wall"
names = [name,]
file_path = "Testscript/Data/precalc.h5" 

# Initial sosution no wall

f.transport_model = "unity-Lewis-number"
f.solve(loglevel=0, refine_grid=True)
f.save(file_path, name=name, overwrite=True)

# Wall 
wall_params = {
    'Z_wall': 0,
    'T_wall': 300.0,
    'factor': 1,
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
while True: 
    f.set_initial_guess(data=file_path, group=names[0]) 
    z -= delta_z
    if z < 0.1:
        break
    wall_params["Z_wall"] = z 
    f.flame.set_non_adiabatic_wall(wall_params)
    f.transport_model = "unity-Lewis-number"
    print(f"delta_z: {delta_z}")
    try:
        f.solve(loglevel=0, refine_grid=True)
        print("\n##############################################################################")
        print(f"Solved at z: {z}")
        print("##############################################################################")
    except BaseException as err:
        print(err)
        print(f"\nError at wallpos. z: {z}")
        print(f"Number off errors: {error_counter}")
        print(f"Factor: {wall_params["factor"]}")
        continue

    name = "z_wall_" + str(z)
    f.save(file_path, name=name, overwrite=True)
    names.append(name)

end_time = time.time()
print(f"Total time loop: {end_time-start_time}")
print(names)
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
