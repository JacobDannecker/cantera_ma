# CH4 Air Cunterflow Diffusion Flame with non adiabatic wall.


import numpy as np
import cantera as ct
from matplotlib import pyplot as plt

reaction_mechanism = "gri30.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3  
grid = np.linspace(0, 18e-3, 200)
f = ct.CounterflowDiffusionFlame(gas, width=width)
f2 = ct.CounterflowDiffusionFlame(gas, width=width)
f.P = 1.e5  
f.fuel_inlet.mdot = 1.5  
f.fuel_inlet.X = "CH4:1"
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 1.5 
f.oxidizer_inlet.X = "O2:0.21, N2:0.79"
f.oxidizer_inlet.T = 300 

f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)


f.set_initial_guess(data="Testscript/Data/CH4-Air.h5", group="no_wall")    
# Set up wall                                                                   
params = {                                                                      
    'Z_wall': 0.9,                                                              
    'T_wall': 300.0,                                                            
    'factor': 1,                                                               
    'mix_frac': 'Bilger',                                                       
    'fuel': 'CH4',                                                               
    'oxidizer': 'O2',                                                           
    'basis': 'mass'                                                             
}                                                                               
f.flame.set_non_adiabatic_wall(params)                                          
f.transport_model = "unity-Lewis-number"                                        
f.solve(loglevel=1)                                                             
f.save("Testscript/Data/CH4-Air.h5", name="wall", overwrite=True)          
f2.restore("Testscript/Data/CH4-Air.h5", name="no_wall")

print(f.T)

# Info
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

