import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy import special, interpolate
from scipy.ndimage import map_coordinates
import utils as ut
# Collect all points at z_spec, correct enthalpy and filter
path = "Data/"
files = [                                                                       

   "stable_Z1.0000.h5",                                                        
   "unstable_Z1.0000.h5",                                                     
   "unstable_Z0.1600.h5",                                                     
   "unstable_Z0.1700.h5",                                                     
   "unstable_Z0.2000.h5",                                                     
   "unstable_Z0.3000.h5",                                                     
   "unstable_Z0.3500.h5",                                                     
   "unstable_Z0.4000.h5",                                                     
   "unstable_Z0.4500.h5",                                                     
   "unstable_Z0.5000.h5",                                                     
   "unstable_Z0.5500.h5",                                                     
   "unstable_Z0.6000.h5",                                                     
   "unstable_Z0.6500.h5",                                                     
   "unstable_Z0.7000.h5",                                                     
   "unstable_Z0.7500.h5",                                                     
   "unstable_Z0.8000.h5",                                                     
   "unstable_Z0.8500.h5",                                                     
   "unstable_Z0.9000.h5",                                                     
   "unstable_Z0.9500.h5",                                                     
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
   "stable_Z0.8500.h5",                                                        
   "stable_Z0.9000.h5",                                                        
   "stable_Z0.9500.h5",                                                        
   "stable_ad.h5",
   "unstable_ad.h5",
]   
file_list = [path + file for file in files]


z_spec = 0.5
quantity = "T" # "H20", "H", ...
cmax_thr = 0.562

gas = ct.Solution("h2o2.yaml")
quantity_all = []
h_all = []
C_all = []

bound_h_max_h = []
bound_h_max_quantity = []
bound_h_max_C = []

bound_h_min_h = []
bound_h_min_quantity = []
bound_h_min_C = []

bound_C_max_h = []
bound_C_max_quantity = []
bound_C_max_C = []


for file in file_list:
    h5_file = h5py.File(file, "r")                                          
    # Get keys to acces h5 groups depending on how they were saved
    if list(h5_file.keys())[0] == "extinction":                                 
        # Stable branch
        names = ["extinction/" + str(name) for name in h5_file["extinction"].keys()][::]
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
        names = [str(name) for name in h5_file.keys()][::] 


    quantity_file = []
    h_file = []
    C_file = []


    cmax_file_C = []
    cmax_file_h = []
    cmax_file_quantity = []

    for name in names:
        print(file, name) 

        if file == "Data/stable_ad.h5" or file == "Data/unstable_ad.h5":
            val_tuple = ut.correct_enthalpy(file, name, "H2O", style="line", active=True)
        else:
            val_tuple = ut.correct_enthalpy(file, name, "H2O", style="vshape", active=True)
            z_wall = h5_file[names[0]]["flame/z"].attrs["z_wall"]

        T, Y, h, z, a_m, idx_C, max_dh, max_dT = val_tuple

        idx_z_spec = np.argmin(np.abs(z - z_spec))

        if not (file == "Data/stable_ad.h5" or  file == "Data/unstable_ad.h5"): 
            idx_wall = np.argmin(np.abs(z - z_wall))

        if quantity == "T":
            quantity_array = T
            print(f"Collecting T")
        else:
            idx_species = ut.get_species_index(quantity)
            print(f"Collecting {quantity} at index: {idx_species}")
            quantity_array = Y[idx_species]
        
        quantity_file.append(quantity_array[idx_z_spec]) 
        print(quantity_array[idx_z_spec])
        h_file.append(h[idx_z_spec])
        C_file.append(Y[idx_C, idx_z_spec])


        if file == "Data/stable_Z0.1600.h5":
            bound_h_min_h.append(h[idx_z_spec])
            bound_h_min_quantity.append(quantity_array[idx_z_spec])
            bound_h_min_C.append(Y[idx_C, idx_z_spec])


        if "ad" in file and h[idx_z_spec] > 0 and (not "unstable" in file):
            #Adiabatic
            bound_h_max_h.append(h[idx_z_spec])
            bound_h_max_quantity.append(quantity_array[idx_z_spec])
            bound_h_max_C.append(Y[idx_C, idx_z_spec])

        if not "unstable" in file:
            # Cmax 
            cmax_file_C.append(Y[idx_C, idx_z_spec])
            cmax_file_h.append(h[idx_z_spec])
            cmax_file_quantity.append(quantity_array[idx_z_spec])
                
    if not (not cmax_file_C):
        idxcmax = np.argmin(cmax_file_h)
        if cmax_file_C[idxcmax] > cmax_thr:
            bound_C_max_h.append(cmax_file_h[idxcmax])
            bound_C_max_quantity.append(cmax_file_quantity[idxcmax])
            bound_C_max_C.append(cmax_file_C[idxcmax])


    quantity_all.append(quantity_file)
    print(quantity_file)
    h_all.append(h_file)
    C_all.append(C_file)

header = quantity + ",h,C"
np.savetxt(f"{quantity}_Z{z_spec}_Ideal.csv", np.column_stack((np.concat(quantity_all), np.concat(h_all), np.concat(C_all))), delimiter=",", header=header)
np.savetxt(f"{quantity}_Z{z_spec}_hmax.csv", np.column_stack((bound_h_max_quantity, bound_h_max_h, bound_h_max_C)), delimiter=",", header=header)
np.savetxt(f"{quantity}_Z{z_spec}_hmin.csv", np.column_stack((bound_h_min_quantity, bound_h_min_h, bound_h_min_C)), delimiter=",", header=header)
np.savetxt(f"{quantity}_Z{z_spec}_Cmax.csv", np.column_stack((bound_C_max_quantity, bound_C_max_h, bound_C_max_C)), delimiter=",", header=header)



plt.scatter(np.concat(C_all), np.concat(h_all), marker='x')
plt.scatter(bound_h_max_C, bound_h_max_h, marker='o', color='b', label="hmax")
plt.scatter(bound_h_min_C, bound_h_min_h, marker='o', color='r', label="hmin")
plt.plot(bound_C_max_C, bound_C_max_h, marker='o', color='g', label="cmax")
plt.show()
