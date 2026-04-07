import numpy as np
import cantera as ct
from matplotlib import pyplot as plt

reaction_mechanism = 'h2o2.yaml'
gas = ct.Solution(reaction_mechanism)
width = 18e-3  

grid = np.linspace(0., width, 200)
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
f.save("Testscript/Data/TestNewBuild.h5", name="Reference", overwrite=True)

plt.plot(f.grid, f.T)
plt.show()
