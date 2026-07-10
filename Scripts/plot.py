import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy import special
import utils as ut

files = ["stable_Z0.9000.h5"]
path = "Scripts/Data/"
file_list = [path + file for file in files]

plt.style.use("seaborn-v0_8-bright")

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

fig, ax = plt.subplots(2, 1)                                                    
fig.suptitle(f"H2/O2")                                                           
fig2, ax2 = plt.subplots(5, 1)                                                  
fig2.suptitle(f"H2/O2")                                                          
fig3, ax3 = plt.subplots(1, 1)                                                  
fig3.suptitle(f"H2/O2")                                                          
                                                                                
idx_H2 = f.gas.species_index("H2")                                              
idx_O2 = f.gas.species_index("O2")                                              
idx_OH = f.gas.species_index("H2O")                                              
idx_H = f.gas.species_index("H")                                              
idx_O = f.gas.species_index("O")                                              


for file in file_list:
    h5_file = h5py.File(file, "r")                                          
    f = ct.CounterflowDiffusionFlame(gas, grid=grid)                                
    # Get keys to acces h5 groups depending on how they were saved
    if list(h5_file.keys())[0] == "extinction":                                 
        names = ["extinction/" + str(name) for name in h5_file["extinction"].keys()]
    else:                                                                       
        names = [str(name) for name in h5_file.keys()]        


    species = []
    amax = []
    tmax = []


    for  name in names:                                                              
        f.restore(file, name=name)                                             

        amax.append(f.strain_rate("max"))
        tmax.append(np.max(f.T))
        label = name                                                                

        # Plot T
        ax[0].plot(f.mixture_fraction("H"), f.T, label=label)                       
        # Plot h
        ax[1].plot(f.mixture_fraction("H"), f.h, label=label)                       
        # Plot Y
        ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label=label)              
        ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], label=label)              
        ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_OH], label=label)              
        ax2[3].plot(f.mixture_fraction("H"), f.Y[idx_O], label=label)              
        ax2[4].plot(f.mixture_fraction("H"), f.Y[idx_H], label=label)              

    idx = np.argsort(amax)
    amax = np.array(amax)[idx]
    tmax = np.array(tmax)[idx]

    ax3.scatter(amax, tmax, marker="x", label=str(file))
for a in ax:                                                                    
    a.grid()                                                                    
    a.legend()                                                                  
    a.set_xlabel("<- ox z fuel ->")                                             
                                                                                
for a in ax2:                                                                   
    a.grid()                                                                    
    a.legend()                                                                  
    a.set_xlabel("<- ox z fuel ->")                                             

ax3.grid()                                                                    
ax3.legend()                                                                  
                                                                                
ax[0].set_ylabel("T")                                                           
ax[1].set_ylabel("h in  J/kg")                                                  
                                                                                
ax2[0].set_ylabel("H2")                                                         
ax2[1].set_ylabel("O2")                                                         
ax2[2].set_ylabel("H2O")                                                         
ax2[3].set_ylabel("O")                                                         
ax2[4].set_ylabel("H")                                                         

ax3.set_xlabel("a_max")
ax3.set_ylabel("T_max")

plt.tight_layout()                                                              
plt.show()                                                                      


