# CH4 Air Cunterflow Diffusion Flame with non adiabatic wall.


import numpy as np
import cantera as ct
from matplotlib import pyplot as plt

reaction_mechanism = 'gri30.yaml'
gas = ct.Solution(reaction_mechanism)
width = 18e-3  

grid = np.linspace(0., width, 250)
f = ct.CounterflowDiffusionFlame(gas, grid=grid)
f2 = ct.CounterflowDiffusionFlame(gas, grid=grid)
f.P = 1.e5  
f.fuel_inlet.mdot = 2.0  
f.fuel_inlet.X = 'CH4:1'
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 2.0
f.oxidizer_inlet.X = 'O2:21;N2:79'
f.oxidizer_inlet.T = 300 

#f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

wall_pos = 125 
factor = 1000000 

f.solve(loglevel=1, wall_pos=wall_pos, factor=factor)

f.save("Testscript/Data/CH4-Air.h5", name="Initial", overwrite=True)
f2.restore("Testscript/Data/CH4-Air.h5", name="Base")


# Info
print(f"T at wall: {f.T[wall_pos]}")
print(f"mdot fuel new build: {f.fuel_inlet.mdot}")
print(f"mdot ox new build: {f.oxidizer_inlet.mdot}")
print(f"mdot fuel reference: {f2.fuel_inlet.mdot}")
print(f"mdot ox reference: {f2.oxidizer_inlet.mdot}")



idx_CH4 = f.gas.species_index("CH4")
idx_O2 = f.gas.species_index("O2")

strain_new_max = f.strain_rate("max")
strain_ref_max = f2.strain_rate("max")

# Fig 1 temp subplot 1
fig, ax = plt.subplots(2, 1)
ax[0].plot(f.grid, f.T, label=f"NewBuild strain_max: {strain_new_max:.2f}")
ax[0].vlines(f.grid[wall_pos], 0, 3000, color="b", linestyle="-.")
ax[0].plot(f2.grid, f2.T, label=f"Reference strain_max: {strain_ref_max:.2f}", linestyle="--")
ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("T")
ax[0].set_xlabel("<- fuel x ox->")
# Fig 1 temp subplot 2
ax[1].plot(f.mixture_fraction("C"), f.T, label=f"NewBuild strain_max: {strain_new_max:.2f}")
ax[1].plot(f2.mixture_fraction("C"), f2.T, label=f"Reference strain_max: {strain_ref_max:.2f}", linestyle="--")
ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("T")
ax[1].set_xlabel("<- ox z fuel ->")




plt.show()

