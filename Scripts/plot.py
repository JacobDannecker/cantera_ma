import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy import special
import utils as ut
#plt.rcParams.update({
#    "font.family": "serif",  # use serif/main font for text elements
#    "text.usetex": True,     # use inline math for ticks
#    "pgf.rcfonts": False     # don't setup fonts from rc parameters
#    })
files = [                                                                       
           "stable_Z0.1600.h5",
          "stable_Z0.1700.h5",
         "stable_Z0.2000.h5",
          "stable_Z0.3000.h5",
          "stable_Z0.4000.h5",
         "stable_Z0.5000.h5",
          "stable_Z0.6000.h5",
         "stable_Z0.7000.h5",
          "stable_Z0.8000.h5",
          "stable_Z0.9000.h5",
          "stable_Z0.9500.h5",
          "stable_Z1.0000.h5",
          "test_stable_Z1.0000.h5",
          "2unstable_Z0.1600.h5",
          "2unstable_Z0.1700.h5",
          "2unstable_Z0.2000.h5",
          "2unstable_Z0.3000.h5",
          "2unstable_Z0.4000.h5",
         "2unstable_Z0.5000.h5",
          "2unstable_Z0.6000.h5",
          "2unstable_Z0.7000.h5",
          "2unstable_Z0.8000.h5",
          "2unstable_Z0.9000.h5",
          "2unstable_Z0.9500.h5",
          "2unstable_Z1.0000.h5",

          "stable_Z0.8500.h5",                                                   
          "stable_Z0.7500.h5",                                                   
          "stable_Z0.6500.h5",                                                   
          "stable_Z0.5500.h5",                                                   
          "stable_Z0.4500.h5",                                                   
          "stable_Z0.5500.h5",                                                   
          "stable_Z0.3500.h5",                                                   
          "stable_Z0.2500.h5",                                                   
         #    "2unstable_Z0.4500.h5"
          ]                                  

files = [
    "2unstable_Z0.1600.h5", 
    "2unstable_Z0.1700.h5", 
    "2unstable_Z0.2000.h5",
    "2unstable_Z0.3000.h5", 
    "2unstable_Z0.3500.h5", 
    "2unstable_Z0.4000.h5",
    "2unstable_Z0.4500.h5", 
    "2unstable_Z0.5000.h5", 
    "2unstable_Z0.5500.h5",
    "2unstable_Z0.6000.h5", 
    "2unstable_Z0.6500.h5", 
    "2unstable_Z0.7000.h5",
    "2unstable_Z0.7500.h5", 
    "2unstable_Z0.8000.h5", 
    "2unstable_Z0.8500.h5",
    "2unstable_Z0.9000.h5", 
    "2unstable_Z0.9500.h5", 
    "2unstable_Z1.0000.h5",
    "stable_Z0.1600.h5", 
    "stable_Z0.1700.h5", 
    "stable_Z0.2000.h5",
    "stable_Z0.2500.h5", 
    "stable_Z0.3000.h5", 
    "stable_Z0.3500.h5",
    "stable_Z0.4000.h5", 
    "stable_Z0.4500.h5", 
    "stable_Z0.5000.h5",
    "stable_Z0.5500.h5", 
    "stable_Z0.6000.h5", 
    "stable_Z0.6500.h5",
    "stable_Z0.7000.h5", 
    "stable_Z0.7500.h5", 
    "stable_Z0.8000.h5",
    "stable_Z0.8500.h5", 
    "stable_Z0.9000.h5", 
    "stable_Z0.9500.h5",
    "stable_Z1.0000.h5",
    "stable_ad.h5",
    "unstable_ad.h5"
]
#files = ["stable_Z1.0000.h5", "stable_Z0.9000.h5", "stable_Z0.9500.h5", "stable_Z0.8500.h5", "test_stable_Z0.9900.h5", ]
#files = ["stable_Z0.3000.h5"]
path = "./Data/"
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

fig4, ax4 = plt.subplots(1, 1)

idx_H2 = f.gas.species_index("H2")                                              
idx_O2 = f.gas.species_index("O2")                                              
idx_H2O = f.gas.species_index("H2O")                                              
idx_H = f.gas.species_index("H")                                              
idx_O = f.gas.species_index("O")                                              
print(f.gas.species_index("OH"))
print(f.gas.species_index("H"))

z_spec = 0.8

for file in file_list:
    h5_file = h5py.File(file, "r")                                          
    f = ct.CounterflowDiffusionFlame(gas, grid=grid)                                
    # Get keys to acces h5 groups depending on how they were saved
    if list(h5_file.keys())[0] == "extinction":                                 
        names = ["extinction/" + str(name) for name in h5_file["extinction"].keys()][::]
    else:                                                                       
        names = [str(name) for name in h5_file.keys()][::]
    
    species = []
    amax = []
    tmax = []
    c_spec = []
    h_spec = []
    T_spec = []

    for  name in names:                                                              
        f.restore(file, name=name)                                             
        rms_b = ut.rms(f.T)
        try:
            if file == "test_stable_Z0.9900.h5":
                ut.correct_enthalpy_flame(f, "H2O", style="line")
            else:
                ut.correct_enthalpy_flame(f, "H2O")
        except:
            continue

        #ut.correct_enthalpy_flame(f, "H2O", style="vshape")
        rms_a = ut.rms(f.T)
        print(f"Diff rms {np.abs(rms_a-rms_b)}")
        print(f.h[0], f.h[-1])
        
        idx_z_spec = np.argmin(np.abs(f.mixture_fraction("H") - z_spec))

        h_spec.append(f.h[idx_z_spec])
        T_spec.append(f.T[idx_z_spec])
        c_spec.append(f.Y[idx_H2O, idx_z_spec])

        if (np.max(f.T) < 400):
            ut.print_r(np.max(f.T))
        if np.max(f.T) > 300:
            amax.append(f.strain_rate("max"))
            #z_stoich = 0.115
            #amax.append(ut.chi_stoich(f, z_stoich))
            tmax.append(np.max(f.T))

        label = name                                                                

        # Plot T
        ax[0].plot(f.mixture_fraction("H"), f.T, label=label)                       
        # Plot h
        ax[1].plot(f.mixture_fraction("H"), f.h, label=label)                       
        # Plot Y
        ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label=label)              
        ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], label=label)              
        ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_H2O], label=label)              
        ax2[3].plot(f.mixture_fraction("H"), f.Y[idx_O], label=label)              
        ax2[4].plot(f.mixture_fraction("H"), f.Y[idx_H], label=label)              
        print(np.max(f.h))
    idx = np.argsort(amax)
    amax = np.array(amax)[idx]
    tmax = np.array(tmax)[idx]

    ax3.scatter(amax, tmax, marker="x", label=str(file))
    ax4.scatter(h_spec, c_spec)
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
plt.savefig("figure.pgf", format="pgf")
plt.show()                                                                      


