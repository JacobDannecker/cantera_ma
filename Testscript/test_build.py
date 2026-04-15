import numpy as np
import cantera as ct
from matplotlib import pyplot as plt

reaction_mechanism = 'h2o2.yaml'
gas = ct.Solution(reaction_mechanism)
width = 18e-3  

grid = np.linspace(0., width, 250)
f = ct.CounterflowDiffusionFlame(gas, grid=grid)
f2 = ct.CounterflowDiffusionFlame(gas, grid=grid)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.5  
f.fuel_inlet.X = 'H2:1'
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 3.0
f.oxidizer_inlet.X = 'O2:1'
f.oxidizer_inlet.T = 300 

#f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

f.solve(loglevel=1, wall_pos=100, factor=1000000)

f.save("Testscript/Data/TestNewBuild.h5", name="NewBuild", overwrite=True)
f2.restore("Testscript/Data/TestNewBuild.h5", name="Reference")


# Info
print(f"T at wall: {f.T[80]}")
print(f"mdot fuel new build: {f.fuel_inlet.mdot}")
print(f"mdot ox new build: {f.oxidizer_inlet.mdot}")
print(f"mdot fuel reference: {f2.fuel_inlet.mdot}")
print(f"mdot ox reference: {f2.oxidizer_inlet.mdot}")



idx_H2 = f.gas.species_index("H2")
idx_O2 = f.gas.species_index("O2")

strain_new_max = f.strain_rate("max")
strain_ref_max = f2.strain_rate("max")

# Fig 1 temp subplot 1
fig, ax = plt.subplots(2, 1)
ax[0].plot(f.grid, f.T, label=f"NewBuild strain_max: {strain_new_max:.2f}")
ax[0].vlines(f.grid[80], 0, 3000, color="b", linestyle="-.")
ax[0].plot(f2.grid, f2.T, label=f"Reference strain_max: {strain_ref_max:.2f}", linestyle="--")
ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("T")
ax[0].set_xlabel("<- fuel x ox->")
# Fig 1 temp subplot 2
ax[1].plot(f.mixture_fraction("H"), f.T, label=f"NewBuild strain_max: {strain_new_max:.2f}")
ax[1].plot(f2.mixture_fraction("H"), f2.T, label=f"Reference strain_max: {strain_ref_max:.2f}", linestyle="--")
ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("T")
ax[1].set_xlabel("<- ox z fuel ->")



# Fig 2 species subplot 1
idx_H2 = f.gas.species_index("H2")
idx_O2 = f.gas.species_index("O2")
idx_OH = f.gas.species_index("OH")
fig2, ax2 = plt.subplots(3, 1)

 #ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label="H2 NewBuild")
ax2[0].plot(f.grid, f.Y[idx_H2], label="H2 NewBuild")
 #ax2[0].plot(f2.mixture_fraction("H"), f2.Y[idx_H2], label="H2 Reference")
ax2[0].plot(f2.grid, f2.Y[idx_H2], label="H2 Reference")

ax2[0].grid()
ax2[0].legend()
ax2[0].set_ylabel("H2")
ax2[0].set_xlabel("<- fuel z ox ->")


# Fig 2 species suplot 2
 #ax2[1].plot(f.mixture_fraction("h"), f.Y[idx_O2], label="O2 newbuild")
ax2[1].plot(f.grid, f.Y[idx_O2], label="O2 newbuild")
 #ax2[1].plot(f2.mixture_fraction("H"), f2.Y[idx_O2], label="O2 Reference")
ax2[1].plot(f2.grid, f2.Y[idx_O2], label="O2 Reference")

ax2[1].grid()
ax2[1].legend()
ax2[1].set_ylabel("O2")
ax2[1].set_xlabel("<- fuel x ox ->")

# Fig 2 species suplot 3
 #ax2[2].plot(f.mixture_fraction("h"), f.Y[idx_OH], label="OH newbuild")
ax2[2].plot(f.grid, f.Y[idx_OH], label="OH newbuild")
 #ax2[2].plot(f2.mixture_fraction("H"), f2.Y[idx_OH], label="OH Reference")
ax2[2].plot(f2.grid, f2.Y[idx_OH], label="OH Reference")

ax2[2].grid()
ax2[2].legend()
ax2[2].set_ylabel("OH")
ax2[2].set_xlabel("<- fuel x ox ->")

plt.show()
