import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy import special, interpolate
from scipy.ndimage import map_coordinates
import utils as ut

plt.style.use("seaborn-v0_8-bright")
#plt.rcParams.update({                                                           
#     "font.family": "serif",  # use serif/main font for text elements            
#     "text.usetex": True,     # use inline math for ticks                        
#     "pgf.rcfonts": False     # don't setup fonts from rc parameters             
#     })

path = "Final_Data/"
path = "Data/"
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
]   
file_list = [path + file for file in files]

fig, ax = plt.subplots(1, 1)                                                    
fig.suptitle(f"H2/O2")                                                           

z_spec = 0.2

not_added = []
added = 0
# Collect all points at z_spec, correct enthalpy and filter
T_all = []
h_all = []
C_all = []
for file in file_list:
    h5_file = h5py.File(file, "r")                                          
    #f = ct.CounterflowDiffusionFlame(gas, grid=grid)                                
    # Get keys to acces h5 groups depending on how they were saved
    if list(h5_file.keys())[0] == "extinction":                                 
        names = ["extinction/" + str(name) for name in h5_file["extinction"].keys()][::10]
    else:                                                                       
        names = [str(name) for name in h5_file.keys()][::10] 

    z_wall = h5_file[names[0]]["flame/z"].attrs["z_wall"]
    T_file = []
    h_file = []
    C_file = []

    h_prev = 500000
    c_prev = 0.0
    for name in names:
        print(file, name) 
        try:
            # Correct enthalpy
            val_tuple = ut.correct_enthalpy(file, name, "H2O")
            if val_tuple == False:
                # Skip, if enthalpy correcion fails
                continue
        except BaseException as e:
            print(e)
            continue
        T, Y, h, z, a_m, idx_C, max_dh, max_dT = val_tuple
        idx_z_spec = np.argmin(np.abs(z - z_spec))
        idx_wall = np.argmin(np.abs(z - z_wall))
        h_now = h[idx_wall]
        c_now = Y[idx_C, idx_z_spec] 
        if np.isclose(z_wall, z[idx_wall], atol=0.01)  and c_now < c_prev :
            # Check if actual z_wall of the flamelet is close to z_wall of file
            # and if enthalpy is decreasing
            T_file.append(T[idx_z_spec]) 
            h_file.append(h[idx_z_spec])
            C_file.append(Y[idx_C, idx_z_spec])
            added  += 1
        else:
            print(f"This is not added {z_wall} {z[idx_wall]} {T[idx_wall]}, {h_prev}, {h_now}")
            not_added.append((file, name))
        h_prev = h_now
        c_prev = c_now
    T_all.append(T_file)
    h_all.append(h_file)
    C_all.append(C_file)
plt.scatter(np.concat(C_all).flatten(), np.concat(h_all).flatten())
plt.show()
print(f"{len(not_added)} flames were discarded. Using total of {added} Flames")
T_all = np.concat(T_all).flatten()
h_all = np.concat(h_all).flatten()
C_all = np.concat(C_all).flatten()


# Interpolate scattered data to rectlinear grid
x_coord = np.linspace(C_all.min(), C_all.max() , 150)
y_coord = np.linspace(h_all.min(), h_all.max(), 150)
x_grid, y_grid = np.meshgrid(x_coord, y_coord)
T_grid = interpolate.griddata((C_all, h_all), T_all ,(x_grid, y_grid), method="linear")
T_arrays = []                                                                   
C_arrays = []                                                                   
h_array = []                                                                    
                                                                                
for i, row in enumerate(T_grid):                                                
    not_nan = np.where(~np.isnan(row))[0]                                       
    if np.count_nonzero(~np.isnan(row)) >= 2:                                   
        T_arrays.append(row[not_nan])                                           
        C_arrays.append(x_grid[i][not_nan])                                     
        h_array.append(y_grid[i,0])                                             

x_coord_norm = np.linspace(0, 1, 300)                                            
y_coord_norm = np.linspace(0, 1, 300)                                            
x_grid_norm, y_grid_norm = np.meshgrid(x_coord_norm, y_coord_norm)              
                                                                                
T_c = np.zeros((len(T_arrays), x_grid_norm.shape[0]))                           
C = np.zeros((y_grid_norm.shape[0], x_grid_norm.shape[0]))                      
h = np.array(h_array)                                                           
                                                                                
for i in range(len(T_arrays)):                                                  
    # Scale C by row and iterpolate T by row                                    
    C_min = C_arrays[i].min()                                                   
    C_max = C_arrays[i].max()                                                   
    C_norm = (C_arrays[i] - C_min)/(C_max - C_min)                              
    C[i] = x_grid_norm[i]                                                       
    T_c[i] = np.interp(C[i], C_norm, T_arrays[i])                               
                                                                                
sort_h_idc = np.argsort(h)                                                
h = h[sort_h_idc]                                                   
T_c = T_c[sort_h_idc]                                                           
C = C[sort_h_idc]                                                               
h_min = h.min()                                                                 
h_max = h.max()                                                                 
h_norm = (h - h_min)/(h_max - h_min)                                      
T_comp = np.zeros((y_coord_norm.shape[0], x_coord_norm.shape[0]))               
for i in range(y_coord_norm.shape[0]):                                               
    T_comp[:,i] = np.interp(y_grid_norm[:,i], h_norm, T_c[:,i])                      
                                                                                

plt.pcolormesh(x_grid_norm, y_grid_norm, T_comp, cmap=mpl.colormaps["afmhot"])                                
plt.suptitle(f"z: {z_spec}")
plt.xlabel("C norm")
plt.ylabel("h")
plt.colorbar()

plt.show()      



# Plot result
##plt.pcolormesh(x_grid, y_grid, z_gridded, cmap=mpl.colormaps["afmhot"])
plt.savefig("figure.svg", format="svg")
#plt.show()



