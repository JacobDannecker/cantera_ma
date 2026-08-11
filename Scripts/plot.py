import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy import special
import utils as ut

plt.rcParams.update({
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": True     # don't setup fonts from rc parameters
    })

files = [
    "unstable_Z0.1600.h5", 
#    "unstable_Z0.1700.h5", 
#    "unstable_Z0.2000.h5",
#    "unstable_Z0.3000.h5", 
#    "unstable_Z0.3500.h5", 
#    "unstable_Z0.4000.h5",
#    "unstable_Z0.4500.h5", 
#    "unstable_Z0.5000.h5", 
    "unstable_Z0.5500.h5",
#    "unstable_Z0.6000.h5", 
#    "unstable_Z0.6500.h5", 
#    "unstable_Z0.7000.h5",
#    "unstable_Z0.7500.h5", 
    "unstable_Z0.8000.h5", 
#    "unstable_Z0.8500.h5",
#    "unstable_Z0.9000.h5", 
    "unstable_Z0.9500.h5", 
    "unstable_Z1.0000.h5",
    "stable_Z0.1600.h5", 
#    "stable_Z0.1700.h5", 
#    "stable_Z0.2000.h5",
#    "stable_Z0.2500.h5", 
#    "stable_Z0.3000.h5", 
#    "stable_Z0.3500.h5",
#    "stable_Z0.4000.h5", 
#    "stable_Z0.4500.h5", 
    #    "stable_Z0.5000.h5",
    "stable_Z0.5500.h5", 
#    "stable_Z0.6000.h5", 
#    "stable_Z0.6500.h5",
#    "stable_Z0.7000.h5", 
#    "stable_Z0.7500.h5", 
    "stable_Z0.8000.h5",
#    "stable_Z0.8500.h5", 
#    "stable_Z0.9000.h5", 
    "stable_Z0.9500.h5",
    "stable_Z1.0000.h5",
#    "stable_ad.h5",
#    "unstable_ad.h5"
]

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

fig1, ax1 = plt.subplots(1, 2)                                                    
fig2, ax2 = plt.subplots(5, 1)                                                  
fig3, ax3 = plt.subplots(1, 1)                                                  


idx_H2 = f.gas.species_index("H2")                                              
idx_O2 = f.gas.species_index("O2")                                              
idx_H2O = f.gas.species_index("H2O")                                              
idx_H = f.gas.species_index("H")                                              
idx_O = f.gas.species_index("O")                                              
print("O", f.gas.species_index("O"))

markers = ["D", "o", "+", "s", "x", "o", "3", "4"]
styles = ["dotted", "dashed", "dashdot", (0, (3, 5, 1, 5, 1, 5)), (0, (3, 1, 1, 1, 1, 1)), "dotted"]
colors = ["r", "b", "g", "k", "c", "m"]

for i, file in enumerate(file_list):
    h5_file = h5py.File(file, "r")                                          

    # Get keys to acces h5 groups depending on how they were saved              
    if list(h5_file.keys())[0] == "extinction":                                 
        # Stable branch                                                         
        names = ["extinction/" + str(name) for name in h5_file["extinction"].keys()]
        inlet_T = 300.0                                                         
        extinct_thr = inlet_T + 100.0                                           
        extinct = [n for n in names                                             
                   if h5_file[n]["flame"]["T"][:].max() <= extinct_thr]         
        if len(extinct) > 1:                                                    
            # Only use the extict flame that is closest to the last burning flame
            strains = {}                                                        
            for n in extinct:                                                   
                _, _, _, a_max, _, _ = ut.load_data(file, n, "H2O")             
                strains[n] = a_max                                              
            keep = max(extinct, key=lambda n: strains[n])                       
            names = [n for n in names if n not in extinct or n == keep]         
    else:                                                                       
        # Unstable branch                                                       
        names = [str(name) for name in h5_file.keys()] 

    if not "ad" in file:
        z_wall = h5_file[names[0]]["flame/z"].attrs["z_wall"]               
    species = []

    amax = []
    tmax = []

    for  name in names:                                                              
        f.restore(file, name=name)                                             
        
        amax.append(f.strain_rate("max"))
        tmax.append(np.max(f.T))
        
         #   # Plot T
     #   ax1[0].plot(f.mixture_fraction("H"), f.T, linestyle=styles[i], color=colors[i])                       
     #   # Plot h
     #   ax1[1].plot(f.mixture_fraction("H"), f.h, linestyle=styles[i], color=colors[i])                       
     #   # Plot Y
     #   ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], linestyle=styles[i], marker=markers[i], color=colors[i])              
     #   ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], linestyle=styles[i], marker=markers[i], color=colors[i])              
     #   ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_H2O], linestyle=styles[i], marker=markers[i], color=colors[i])              
     #   ax2[3].plot(f.mixture_fraction("H"), f.Y[idx_O], linestyle=styles[i], marker=markers[i], color=colors[i])              
     #   ax2[4].plot(f.mixture_fraction("H"), f.Y[idx_H], linestyle=styles[i], marker=markers[i], color=colors[i])              
        

         #   color = 'r'
     #   style = ':' 
     #   marker = 'x'
     #   labelT = r"Korrigiert bei $z_{wall}$ = " + str(z_wall) + ", $\Delta$ rms(T) = " + f"{d_T:.2f} / K"
     #   labelh = r"Korrigiert bei $z_{wall}$ = " + str(z_wall) + ", $\Delta$ rms(h) = " + f"{d_T:.2f} / J/kg"


    idx = np.argsort(amax)
    amax = np.array(amax)[idx]
    tmax = np.array(tmax)[idx]
    # Plot Tmax over amax
    if "Z0.1600" in file:
        marker = markers[0]
        color = colors[0]
        style = styles[0]
    elif "Z0.5500" in file:
        marker = markers[1]
        color = colors[1]
        style = styles[1]
    elif "Z0.8000" in file:
        marker = markers[2]
        color = colors[2]
        style = styles[2]
    elif "Z0.9500" in file:
        marker = markers[3]
        color = colors[3]
        style = styles[3]
    elif "Z1.0000" in file:
        marker = markers[4]
        color = colors[4]
        style = styles[4]

 

    ax3.scatter(amax, tmax, color=color, marker=marker)


#fig1.suptitle(r"H2/O2 Temperatur und Enthalpie bei Verschiedende $z_{wall}$", size="large")
#fig2.suptitle(r"H2/O2 Spezies bei verschiedenden $z_{wall}$", size="large")
fig3.suptitle(r"H2/O2 $T_{max}$ über $a_{max}$ bei verschiedenen $z_{wall}$", size="large")

for a in ax1:                                                                    
    a.plot([],[],color='k', linestyle='-', label='Jewiliger Originalverlauf')
    a.legend()
    a.grid()                                                                    
    a.legend()                                                                  
    a.set_xlabel("z")                                             
                                                                                
for a in ax2:                                                                   
    a.plot([],[],color='k', linestyle='-', label='Original')
    a.plot([],[],color='r', linestyle=':', label='Korrektur')
    a.legend()
    a.grid()                                                                    
    a.set_xlabel("z")                                             

ax3.scatter([],[],color=colors[0], marker=markers[0], label='$z_{wall} = $ 0.16 ')
ax3.scatter([],[],color=colors[1], marker=markers[1], label='$z_{wall} = $ 0.55 ')
ax3.scatter([],[],color=colors[2], marker=markers[2], label='$z_{wall} = $ 0.80')
ax3.scatter([],[],color=colors[3], marker=markers[3], label='$z_{wall} = $ 0.95')
ax3.scatter([],[],color=colors[4], marker=markers[4], label='$z_{wall} = $ 1.0')

ax3.legend()
ax3.grid()                                                                    
ax3.set_xlabel("z")                                             


#ax3.legend()
#ax3.grid()                                                                    

                                                                                
ax1[0].set_ylabel("$T$ in $K$")                                                           
ax1[1].set_ylabel("$h$ in  $J/kg$")                                                  
                                                                                
ax2[0].set_ylabel("H2")                                                         
ax2[1].set_ylabel("O2")                                                         
ax2[2].set_ylabel("H2O")                                                         
ax2[3].set_ylabel("O")                                                         
ax2[4].set_ylabel("H")                                                         
#
ax3.set_xlabel(r"$a_{max}$ in $1/s$")
ax3.set_ylabel(r"$T_{max}$ in $K$")

plt.tight_layout()                                                              
plt.show()                                                                      


