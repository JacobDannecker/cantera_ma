import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy import special
import utils as ut

file_path = f"Scripts/Data/abcde.h5"

#file_path = f"Scripts/Data/No32_test_final_Z0.4000.h5"

#file_path = f"Scripts/Data/Run2_final_Z0.7500.h5"
#csv_path = f"Scripts/Data/stable_050.csv"


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
fig.suptitle(f"H2/O2 {file_path}")                                                           
                                                                                
fig2, ax2 = plt.subplots(5, 1)                                                  
fig2.suptitle(f"H2/O2 {file_path}")                                                          
                                                                                
fig3, ax3 = plt.subplots(1, 1)                                                  
fig3.suptitle(f"H2/O2 {file_path}")                                                          
                                                                                
idx_H2 = f.gas.species_index("H2")                                              
idx_O2 = f.gas.species_index("O2")                                              
idx_OH = f.gas.species_index("H2O")                                              
idx_H = f.gas.species_index("H")                                              
idx_O = f.gas.species_index("O")                                              
h5_file = h5py.File(file_path, "r")                                          
#names = ["extinction/" + str(name)  for name in h5_file["extinction"].keys()][::]
names = [str(name)  for name in h5_file.keys()][::]
#names.append("initial")
species_idx = f.gas.species_index("OH") 
species_name = "OH"
species = []
amax = []
tmax = []

for  name in names:                                                              
    f.restore(file_path, name=name)                                             
    idx = np.argmin(np.abs(f.mixture_fraction("H") - 0.50))

    if (f.T[idx] - 300) < 5. and (f.mixture_fraction("H")[idx] - 0.50) < 0.01:
        print(f.mixture_fraction("H")[idx])
        amax.append(f.strain_rate("mean"))
        tmax.append(np.max(f.T))

    label = name                                                                
    idx_T_max = np.abs(f.T - np.max(f.T)).argmin()
    species_T_max = f.Y[species_idx][:][idx_T_max]
    species.append(species_T_max)

    # Subplot 1 Temperature                                                     
    ax[0].plot(f.mixture_fraction("H"), f.T, label=label)                       
    # Subplot 2  enthalpy                                                       
    ax[1].plot(f.mixture_fraction("H"), f.h, label=label)                       
#    ax[2].plot(f.mixture_fraction("H"), f.equivalence_ratio, label=label)                       
    # Sublots species                                                           
    ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label=label)              
    ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], label=label)              
    ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_OH], label=label)              
    ax2[3].plot(f.mixture_fraction("H"), f.Y[idx_O], label=label)              
    ax2[4].plot(f.mixture_fraction("H"), f.Y[idx_H], label=label)              
                                                                                
for a in ax:                                                                    
    a.grid()                                                                    
    #a.legend()                                                                  
    a.set_xlabel("<- ox z fuel ->")                                             
                                                                                
for a in ax2:                                                                   
    a.grid()                                                                    
    #a.legend()                                                                  
    a.set_xlabel("<- ox z fuel ->")                                             
                                                                                
ax[0].set_ylabel("T")                                                           
ax[1].set_ylabel("h in  J/kg")                                                  
                                                                                
ax2[0].set_ylabel("H2")                                                         
ax2[1].set_ylabel("O2")                                                         
ax2[2].set_ylabel("H2O")                                                         
ax2[3].set_ylabel("O")                                                         
ax2[4].set_ylabel("H")                                                         

#df = pd.read_csv(csv_path)
#fig3, ax3 = plt.subplots(1, 2)
#ax3[0].plot(df.strain_rate, df.T_max)
#species = np.array(species)
#ax3[1].plot(np.sort(species[1:]), df.T_max)
#
#ax3[0].set_xlabel('Maximum Axial Velocity Gradient [1/s]')
#ax3[0].set_ylabel('Maximum Temperature [K]')
#
#ax3[1].set_xlabel(species_name + " at Tmax")
#ax3[1].set_ylabel('Maximum Temperature [K]')
#

ax3.plot(amax, tmax, marker="x", color="b")
ax3.set_xlabel("a_max")
ax3.set_ylabel("T_max")
plt.tight_layout()                                                              
plt.show()                                                                      


