import numpy as np                                                              
import matplotlib.pyplot as plt                                                 
from scipy import interpolate                                                   
from scipy.interpolate import RegularGridInterpolator                           
import h5py 



def readRawData(file, quantity="T", idx_C=3, z_spec=0.9):
    h5_file = h5py.File(file, "r")                                          
    runs = [str(name) for name in h5_file.keys()]                           
    list_z = []
    list_C = []
    list_quantity = []
    for run in runs:                                                        
        if run == "initial":
            continue
        quantity_all = np.array(h5_file[run]["flame"][quantity])        
        z_all = np.array(h5_file[run]["flame/z"])                       
        C_all = np.array(h5_file[run]["flame/Y"]).T[idx_C]
        idx_z = np.abs(z_all - z_spec).argmin()
        list_z.append(z_all[idx_z])
        list_C.append(C_all[idx_z])
        list_quantity.append(quantity_all[idx_z])
       

     #   if z == z_wall_spec:                                                
     #       print(f"run: {run}")                                            
     #       raw_grid = np.array(h5_file[run]["flame/grid"])                 
     #       raw_z = np.array(h5_file[run]["flame/z"])                       
     #       raw_quantity = np.array(h5_file[run]["flame"][quantity])        
     #                                                                       
     #       # Sort quantity according to z                                  
     #       idc = np.argsort(raw_z)                                         
     #       raw_z = raw_z[idc]                                              
     #       raw_quantity = raw_quantity[idc]                                
     #       interpol_grid_x = np.linspace(0, 1, n_x)                        
     #       interpol_func = interpolate.interp1d(raw_z, raw_quantity)       
     #       quantity_interpolated = interpol_func(interpol_grid_x)          
     #                                                                       
     #       quantity_all_runs.append(quantity_interpolated)                 
     #       chi_st_all_runs.append(float(h5_file[run]["flame/z"].attrs["chi_st"]))

       
    all_C = np.array(list_C)
    all_quantity = np.array(list_quantity)
    z_wall = "to be done"
    h5_file.close()  
    return (z_wall, all_C, all_quantity)                                                      

def sortDataAccordingToC(all_C, all_quantity):
    idx_sort = np.argsort(all_C)
    quantity_all_sorted = all_quantity[idx_sort]
    C_all_sorted = all_C[idx_sort]
    
    return (C_all_sorted, quantity_all_sorted)


def interpolate1dData(C, T, n=400):
    xnew = np.linspace(0, 1, num=n)
    ynew = np.interp(xnew, C, T)
    return (xnew, ynew)

file_list = ["Scripts/Data/z090.h5", "Scripts/Data/z060.h5", "Scripts/Data/z070.h5"]
idx_C = 5
z_spec = 0.5

fig, ax = plt.subplots()
max_Cs = []

for file in file_list:
    z_wall, x, y, = readRawData(file, idx_C=idx_C, z_spec=z_spec)
    max_Cs.append(np.max(x))
max_C = np.max(max_Cs)

for file in file_list:
    z_wall, x, y = readRawData(file, idx_C=idx_C, z_spec=z_spec)
    x, y = sortDataAccordingToC(x, y)
    x, y = interpolate1dData(x, y)
    x = x/max_C
    print(z_wall)
    ax.plot(x, y, label=z_wall)

ax.legend()
plt.show()

