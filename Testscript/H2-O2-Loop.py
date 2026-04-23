import numpy as np
import cantera as ct
from scipy import special
from matplotlib import pyplot as plt

# Flame settings
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3  
f = ct.CounterflowDiffusionFlame(gas, width=width)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.5  
f.fuel_inlet.X = "H2:1"
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 3.0 
f.oxidizer_inlet.X = "O2:1"
f.oxidizer_inlet.T = 300 
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

# Names
name = "no_wall"
names = [name,]
file_path = "Testscript/Data/H2-O2-loop.h5" 

# Initial sosution no wall

f.transport_model = "unity-Lewis-number"
f.solve(loglevel=1)
f.save(file_path, name=name, overwrite=True)

# Wall 
wall_params = {
    'Z_wall': 1,
    'T_wall': 300.0,
    'factor': 1,
    'mix_frac': 'Bilger',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}

# Loop over z
z_vals = np.arange(0, 1, 0.1)[::-1]
for z in z_vals:
    f.set_initial_guess(data=file_path, group=name) 
    wall_params["Z_wall"] = z 
    f.flame.set_non_adiabatic_wall(wall_params)
    f.transport_model = "unity-Lewis-number"
    try:
        f.solve(loglevel=1)
        print(f"Solved at z: {z}")
    except:
        print(f"Error at wallpos. z: {z}")
        continue
    name = "z_wall_" + str(z)
    f.save(file_path, name=name, overwrite=True)
    names.append(name)


def chi_stoich(f):
        #a = f.strain_rate("mean")
        a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
        z = 0.055
        chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z))**2)))
        return chi_stoich


# Plot
fig, ax = plt.subplots(2, 1)
fig.suptitle("H2/O2") 
for name in names:
    f.restore(file_path, name=name)
    # Subplot 1 Temp 
    ax[0].plot(f.mixture_fraction("H"), f.T, label=f"With Wall chi_st: {chi_stoich(f):.2f}")
    # Subplot 2  enthalpy
    ax[1].plot(f.mixture_fraction("H"), f.h, label=f"With Wall chi_st: {chi_stoich(f):.2f}")

ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("T")
ax[1].set_xlabel("<- ox z fuel ->")

ax[2].grid()
ax[2].legend()
ax[2].set_ylabel("h in  J/kg")
ax[2].set_xlabel("<- ox z fuel ->")

plt.show()
