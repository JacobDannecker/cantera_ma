import numpy as np
import cantera as ct
from matplotlib import pyplot as plt
import h5py


file_list = [ "Scripts/Data/extinction_Z0.8000.h5", "Scripts/Data/extinction_Z0.7500.h5", "Scripts/Data/extinction_Z0.7000.h5", "Scripts/Data/extinction_Z0.6500.h5",
             "Scripts/Data/extinction_Z0.6000.h5","Scripts/Data/extinction_Z0.5500.h5", "Scripts/Data/extinction_Z0.5000.h5", "Scripts/Data/extinction_Z0.4500.h5",
"Scripts/Data/extinction_Z0.4000.h5", "Scripts/Data/extinction_Z0.4500.h5", "Scripts/Data/extinction_Z0.3000.h5", "Scripts/Data/extinction_Z0.2500.h5", 
"Scripts/Data/extinction_Z0.2000.h5", "Scripts/Data/extinction_Z0.1500.h5", "Scripts/Data/extinction_Z0.1000.h5", "Scripts/Data/No5extinction_Z0.9000.h5", "Scripts/Data/unstable_No17.h5"]


for file in file_list:
    h5_file = h5py.File(file, "r")
    if file == "Scripts/Data/unstable_No17.h5":
        runs = [str(name)  for name in h5_file.keys()]
    else:
        runs = ["extinction/" + str(name)  for name in h5_file["extinction"].keys()]
    chi_list = []
    Tmax_list = []
    for run in runs:
        T = h5_file[run]["flame"]["T"]
        Tmax = np.max(T) 
        chi_st = h5_file[run]["flame"]["z"].attrs["chi_st"]
        chi_list.append(chi_st)
        Tmax_list.append(Tmax)
        
    plt.semilogx(chi_list, Tmax_list, label=file)
plt.xlabel("chi_st")
plt.ylabel("Tmax")
plt.legend()
plt.show()

