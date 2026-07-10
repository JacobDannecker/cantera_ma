import numpy as np
import cantera as ct
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy import special, interpolate
import utils as ut
plt.style.use("seaborn-v0_8-bright")

def load_data(file_path, name, C):
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
    f.restore(file_path, name)
    return f.grid, f.T, f.Y, f.strain_rate("max"), f.gas.species_index(C), reaction_mechanism


def compute_enthalpy_and_Z(gas, T, Y, fuel_idx=0, oxidizer_idx=3):
    n = len(T)
    h = np.empty(n)
    Z = np.empty(n)
    fuel = np.zeros(gas.n_species)
    oxidizer = np.zeros(gas.n_species)
    fuel[fuel_idx] = 1.0
    oxidizer[oxidizer_idx] = 1.0
    for j in range(n):
        gas.TPY = T[j], gas.P, Y.T[j]
        h[j] = gas.enthalpy_mass
        Z[j] = gas.mixture_fraction(fuel, oxidizer, basis="mass")
    return h, Z


def perfect_v_shape(Z, h, zero_ends=False):
    i_tip = np.argmin(h)
    if zero_ends:
        left_z = np.array([Z[0], Z[i_tip]])
        left_h = np.array([0.0, h[i_tip]])
        right_z = np.array([Z[i_tip], Z[-1]])
        right_h = np.array([h[i_tip], 0.0])
    else:
        left_z = Z[: i_tip + 1]
        left_h = h[: i_tip + 1]
        right_z = Z[i_tip:]
        right_h = h[i_tip:]
    left = np.polyfit(left_z, left_h, 1)
    right = np.polyfit(right_z, right_h, 1)
    h_v = np.where(np.arange(len(Z)) <= i_tip, np.polyval(left, Z),
                   np.polyval(right, Z))
    return h_v, i_tip


def temperature_from_HPY(gas, h, Y, P=None):
   if P is None:
       P = gas.P
   n = len(h)
   T = np.empty(n)
   for j in range(n):
       gas.HPY = h[j], P, Y.T[j]
       T[j] = gas.T
   return T

def compute(file_path, name, C):                                         
   grid, T_orig, Y, a_max, idx_C, mech = load_data(file_path, name, C)                          
   gas_i = ct.Solution(mech)                                                   
   P = gas_i.P                                                                 
   print(grid.shape, T_orig.shape, Y.shape)
   h_orig, Z = compute_enthalpy_and_Z(gas_i, T_orig, Y)                        
   h_v, i_tip = perfect_v_shape(Z, h_orig, zero_ends=False)                     
   T_v = temperature_from_HPY(gas_i, h_v, Y)                                   
   max_dh = np.max(np.abs(h_orig - h_v))                                       
   max_dT = np.max(np.abs(T_orig - T_v))                                       
   print(f"Max dh: {max_dh},Max dt: {max_dT}")
   return T_v, Y, h_v, Z, a_max, idx_C 


path = "Final_Data/"
#files = ["stable_Z0.5000.h5", "stable_3_Z0.3000.h5", "stable_Z1.0000.h5", "unstable_Z0.5000.h5", "stable_2_Z0.6000.h5", "stable_unstable_Z0.2000.h5", "stable_Z0.2000.h5"]
files = ["stable_Z0.5000.h5", "unstable_Z0.5000.h5"]#, "stable_Z1.0000.h5", "stable_2_Z0.6000.h5"] #, "stable_3_Z0.3000.h5", "stable_Z1.0000.h5", "unstable_Z0.5000.h5", "stable_2_Z0.6000.h5", "stable_unstable_Z0.2000.h5", "stable_Z0.2000.h5"]
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
    
    for name in names:
        
        T, Y, h, z, a_m, idx_C = compute(file, name, "H2O")
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
#plt.show()

print(np.max(C_all))
print(np.max(h_all))

print(C_all.shape)
x_coord = np.linspace(0, 1, 800)
y_coord = np.linspace(np.min(h_all), np.max(h_all), 800)
x_grid, y_grid = np.meshgrid(x_coord, y_coord)

z_gridded = interpolate.griddata((C_all/np.max(C_all), h_all), T_all ,(x_grid, y_grid), method="linear")
print(z_gridded.shape)

plt.pcolormesh(x_grid, y_grid, z_gridded)
plt.suptitle(f"z: {z_spec}")
plt.xlabel("C norm")
plt.ylabel("h")
plt.colorbar()
plt.show()



