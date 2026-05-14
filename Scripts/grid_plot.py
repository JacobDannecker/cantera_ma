import numpy as np                                                              
import matplotlib.pyplot as plt                                                 
from scipy import interpolate                                                   
from scipy.interpolate import RegularGridInterpolator                           
import h5py 



def readRawData(file, quantity_ident="T", idx_progress_var=3, z_spec=0.9, n=10000):
    h5_file = h5py.File(file, "r")                                          
    runs = [str(name) for name in h5_file.keys()]                           
    list_z = []
    list_progress_var = []
    list_quantity = []
    for run in runs:                                                        
        if run == "initial":
            continue

        if type(quantity_ident) == str:
            quantity = np.array(h5_file[run]["flame"][quantity_ident])        
        elif type(quantity_ident) == int:
            quantity = np.array(h5_file[run]["flame/Y"]).T[quantity_ident]

        grid = np.array(h5_file[run]["flame/grid"])
        z = np.array(h5_file[run]["flame/z"])                       
        progress_var = np.array(h5_file[run]["flame/Y"]).T[idx_progress_var]

        idx_z = np.abs(z - z_spec).argmin()
        print(z[idx_z])
        print(quantity[idx_z])
        xnew = np.linspace(0, np.max(grid), num=n)
        z = np.interp(xnew, grid, z)
        idx_z = np.abs(z - z_spec).argmin()

        quantity = np.interp(xnew, grid, quantity)
        print(z[idx_z])
        print(quantity[idx_z])
        progress_var = np.interp(xnew, grid, progress_var)

        idx_z = np.abs(z - z_spec).argmin()
        
        list_z.append(z[idx_z])
        list_progress_var.append(progress_var[idx_z])
        list_quantity.append(quantity[idx_z])
       
    all_progress_var = np.array(list_progress_var)
    all_quantity = np.array(list_quantity)
    z_wall = "to be done"
    h5_file.close()  
    return (z_wall, all_progress_var, all_quantity)                                                      

def sortDataAccordingToprogress_var(all_progress_var, all_quantity):
    idx_sort = np.argsort(all_progress_var)
    all_quantity_sorted = all_quantity[idx_sort]
    all_progress_var_sorted = all_progress_var[idx_sort]
    return (all_progress_var_sorted, all_quantity_sorted)


def interpolate1dData(progress_var, quantity, norm, n=100):
    progress_var = progress_var/norm
    xnew = np.linspace(0, 1, num=n)
    ynew = np.interp(xnew, progress_var, quantity, left=0.0, right=0.0)
    return (xnew, ynew)

file_list = ["Scripts/Data/z090.h5", "Scripts/Data/z070.h5"]
idx_progress_var = 5
z_spec = 0.4 

fig, ax = plt.subplots()

progress_var_max_list = []
for file in file_list:
    z_wall, x, y = readRawData(file, idx_progress_var=idx_progress_var, z_spec=z_spec)
    progress_var_max_list.append(np.max(x))

norm = np.max(progress_var_max_list)

for file in file_list:
    z_wall, x, y = readRawData(file, idx_progress_var=idx_progress_var, z_spec=z_spec)

    x, y = sortDataAccordingToprogress_var(x, y)
    print(np.min(y), np.max(y))
    #x, y = interpolate1dData(x, y, norm)

    #Overwrite z_wall, since it is not in .h5 yet only in name of file
    z_wall = float(''.join([c for c in file.strip("h5") if c in '1234567890.']))/100

    ax.plot(x, y, label="z_wall = " + str(z_wall))
    ax.scatter(x, y, marker="x")
fig.suptitle("z = " + str(z_spec))
ax.set_xlabel("progress_var OH")
ax.set_ylabel("T")
ax.legend()
plt.show()

