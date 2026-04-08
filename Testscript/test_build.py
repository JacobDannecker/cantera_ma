import numpy as np
import cantera as ct
from matplotlib import pyplot as plt

reaction_mechanism = 'h2o2.yaml'
gas = ct.Solution(reaction_mechanism)
width = 18e-3  

grid = np.linspace(0., width, 250)
f = ct.CounterflowDiffusionFlame(gas, grid=grid)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.5  
f.fuel_inlet.X = 'H2:1'
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 3.0
f.oxidizer_inlet.X = 'O2:1'
f.oxidizer_inlet.T = 300 

#f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

f.solve(loglevel=1)
f.save("Testscript/Data/TestNewBuild.h5", name="NewBuild", overwrite=True)
fig, ax = plt.subplots(2, 1)
strain_new_max = f.strain_rate("max")

ax[0].plot(f.grid, f.T, label=f"NewBuild strain_max: {strain_new_max:.2f}")
ax[0].scatter(f.grid[100], f.T[100], color="r")
ax[1].plot(f.mixture_fraction("H"), f.T, label=f"NewBuild strain_max: {strain_new_max:.2f}")

print(f.fuel_inlet.mdot)
print(f.oxidizer_inlet.mdot)

f.restore("Testscript/Data/TestNewBuild.h5", name="Reference")
print(f.fuel_inlet.mdot)
print(f.oxidizer_inlet.mdot)
strain_ref_max = f.strain_rate("max")

ax[0].plot(f.grid, f.T, label=f"Reference strain_max: {strain_ref_max:.2f}")
ax[1].plot(f.mixture_fraction("H"), f.T, label=f"Reference strain_max: {strain_ref_max:.2f}")

ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("T")
ax[1].set_xlabel("<- ox z fuel ->")

ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("T")
ax[0].set_xlabel("<- fuel x ox->")

plt.show()
