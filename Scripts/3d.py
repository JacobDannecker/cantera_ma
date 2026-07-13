import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy import special, interpolate
import utils as ut
plt.style.use("seaborn-v0_8-bright")

path = "Final_Data/"
path = "Data/"
#files = ["stable_Z0.5000.h5", "stable_3_Z0.3000.h5", "stable_Z1.0000.h5", "unstable_Z0.5000.h5", "stable_2_Z0.6000.h5", "stable_unstable_Z0.2000.h5", "stable_Z0.2000.h5"]
#files = ["stable_Z0.5000.h5", "unstable_Z0.5000.h5", "stable_Z1.0000.h5", "stable_2_Z0.6000.h5"]#, "stable_Z1.0000.h5", "stable_2_Z0.6000.h5"] #, "stable_3_Z0.3000.h5", "stable_Z1.0000.h5", "unstable_Z0.5000.h5", "stable_2_Z0.6000.h5", "stable_unstable_Z0.2000.h5", "stable_Z0.2000.h5"]
files = ["stable_Z1.0000.h5"]
file_list = [path + file for file in files]

fig, ax = plt.subplots(1, 1)                                                    
fig.suptitle(f"H2/O2")                                                           

z_spec = 0.8

    
T_all = []
h_all = []
C_all = []
a_max = []
T_max = []

for file in file_list:
    h5_file = h5py.File(file, "r")                                          
    #f = ct.CounterflowDiffusionFlame(gas, grid=grid)                                
    # Get keys to acces h5 groups depending on how they were saved
    if list(h5_file.keys())[0] == "extinction":                                 
        names = ["extinction/" + str(name) for name in h5_file["extinction"].keys()]
    else:                                                                       
        names = [str(name) for name in h5_file.keys()]        
    z_wall = h5_file[names[0]]["flame/z"].attrs["z_wall"]
    names = names[:-4]    
    for name in names:
        
        T, Y, h, z, a_m, idx_C = ut.correct_enthalpy(file, name, "H2O")
        idx_z_spec = np.argmin(np.abs(z - z_spec))
        idx_wall = np.argmin(np.abs(z - z_wall))
        a_max.append(a_m)
        T_max.append(np.max(T))
        if np.isclose(z_wall, z[idx_wall], atol=0.001): #np.isclose(T[idx_wall], 300., atol=1.0) and 
            T_all.append(T[idx_z_spec]) 
            h_all.append(h[idx_z_spec])
            C_all.append(Y[idx_C, idx_z_spec])
        else:
            print(f"This is not added {z_wall} {z[idx_wall]} {T[idx_wall]}")
    
#    plt.scatter(a_max, T_max)
    a_max = []
    T_max = []

T_all = np.array(T_all).flatten()
h_all = np.array(h_all).flatten()
C_all = np.array(C_all).flatten()
C_min = np.min(C_all)
C_max = np.max(C_all)
h_min = np.min(h_all)
h_max = np.max(h_all)

#plt.show()

print(np.max(C_all))
print(np.max(h_all))

print(C_all.shape)
x_coord = np.linspace(C_min, C_max , 800)
y_coord = np.linspace(h_min, h_max, 800)
x_grid, y_grid = np.meshgrid(x_coord, y_coord)

z_gridded = interpolate.griddata((C_all, h_all), T_all ,(x_grid, y_grid), method="linear", fill_value=300.)

# Scale 
x_grid = (x_grid - C_min)/(C_max - C_min)
y_grid = (y_grid - h_min)/(h_max - h_min)

plt.pcolormesh(x_grid, y_grid, z_gridded, cmap=mpl.colormaps["afmhot"])
plt.suptitle(f"z: {z_spec}")
plt.xlabel("C norm")
plt.ylabel("h")
plt.colorbar()
plt.show()



