import numpy as np
import cantera as ct
from scipy import special
from matplotlib import pyplot as plt
import time

def get_delta_T(f, wall_params):
    # Returns True if delta_T smaller than delta_T_max else reurns False
    idx_wall = np.abs(f.mixture_fraction(wall_params["mix_frac"]) - wall_params["Z_wall"]).argmin()
    delta_T_wall =  f.T[idx_wall] - wall_params["T_wall"] 
    return delta_T_wall


def clustered_grid(L, n, cluster_frac=0.32, beta_left=1.5, beta_right=2.6):     
   s = np.linspace(0, 1, n)                                                    
   x = np.empty(n)                                                             
   # Left side: arctanh clusters points near s=cluster_frac                    
   mask = s <= cluster_frac                                                    
   z = s[mask] / cluster_frac                                                  
   y = 1 - np.arctanh(np.tanh(beta_left) * (1 - z)) / beta_left                
   x[mask] = cluster_frac * L * y                                              
   # Right side: arctanh clusters points near s=cluster_frac                   
   z = (s[~mask] - cluster_frac) / (1 - cluster_frac)                          
   y = np.arctanh(np.tanh(beta_right) * z) / beta_right                        
   x[~mask] = cluster_frac * L + (1 - cluster_frac) * L * y                    
   return x      

reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3  
grid = np.linspace(0, width, 200)
f = ct.CounterflowDiffusionFlame(gas, width=width)
gas = ct.Solution("h2o2.yaml")
f2 = ct.CounterflowDiffusionFlame(gas, width=width)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.1  
f.fuel_inlet.X = "H2:1"
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 0.5 
f.oxidizer_inlet.X = "O2:1"
f.oxidizer_inlet.T = 300 
z_stoich = 0.111

#f.set_refine_criteria(ratio=2.0, slope=0.5, curve=0.5, prune=0.03)
file_name = "Scripts/Data/initialTempCurve07.h5"
#f.set_refine_criteria(ratio=2.0, slope=0.5, curve=0.5, prune=0)
# Set up wall 
wall_params = {
    'Z_wall': 0.7,
    'T_wall': 300,
    'factor': 1000,
    'mix_frac': 'Bilger',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
    }
f.transport_model = "unity-Lewis-number"
#tol_ss= [1.0e-5, 1.0e-9]# [rtol atol] for steady-state problem
#tol_ts= [1.0e-5, 1.0e-9]# [rtol atol] for time stepping
#f.flame.set_steady_tolerances(default=tol_ss)
#f.flame.set_transient_tolerances(default=tol_ts)

#f.flame.set_steady_tolerances(default=(5e-3, 5e-3),
#                            T=(3e-6, 0.),
#                            Y=(7e-8, 0.))
#

f.set_refine_criteria(ratio=3, slope=0.5, curve=0.05, prune=0.04,
                      enthalpy=False, enthalpy_curve=0.05)

#f.set_initial_guess(data=file_name, group="") 

f.max_time_step_count = 1000 
delta_T_max = 1.0 
delta_T_ok = False
while not delta_T_ok:
    try:
        f.flame.set_non_adiabatic_wall(wall_params)                     
        f.solve(loglevel=1, refine_grid=True, auto=True)                           
        delta_T_wall = get_delta_T(f, wall_params)
        if delta_T_wall > 10:
            wall_params["factor"] *= 2
        else:
            wall_params["factor"] *= 1.9 
        print(f"DELTA T WALL: {delta_T_wall}")
        if delta_T_wall < delta_T_max:
            print(f"mdot f, o : {f.fuel_inlet.mdot}, {f.oxidizer_inlet.mdot}")
            print(f"Strain max: {f.strain_rate("max")}")
            delta_T_ok = True
            print(f"Solution found at delta_T_wall: {delta_T_wall}")

    except BaseException as err:
        print(err)

f.save(file_name, name="initial", overwrite=True)

#f.save("Scripts/Data/stable_090.h5", name="initial", overwrite=True)
f2.restore(file_name, name="initial")


def chi_stoich(f, z_stoich):
        #a = f.strain_rate("mean")
        a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
        chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
        return chi_stoich

chi_st_new = chi_stoich(f, z_stoich) 
chi_st_ref = chi_stoich(f2, z_stoich) 

# Info
print(f"mdot fuel new build: {f.fuel_inlet.mdot}")
print(f"mdot ox new build: {f.oxidizer_inlet.mdot}")
print("-----------------------------------------------")
print(f"chi_st new build: {chi_st_new}")
print(f"mdot fuel reference: {f2.fuel_inlet.mdot}")
print(f"mdot ox reference: {f2.oxidizer_inlet.mdot}")
print(f"chi_st reference: {chi_st_ref}")


idx_H2 = f.gas.species_index("H2")
idx_O2 = f.gas.species_index("O2")

# Fig 1 subplot 1
fig, ax = plt.subplots(4, 1)
fig.suptitle(" H2/O2") 

#ax[0].plot(f.mixture_fraction("H"), f.density, label=f"new")
#ax[0].plot(f.mixture_fraction("H"), np.gradient(f.mixture_fracion("H"), f.mixture_fraction("H")), label=f"old")
#ax[0].scatter(f.mixture_fraction("H"), f.grid, label=f"new")
#ax[0].scatter(f2.mixture_fraction("H"), f2.grid, label=f"old")
ax[0].scatter(f.mixture_fraction("H"), np.zeros(f.grid.shape[0]), label=f"new")
ax[0].scatter(f2.mixture_fraction("H"), np.ones(f2.grid.shape[0]), label=f"old")


ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("grid")
ax[0].set_xlabel("<- fuel x ox->")
# Fig 1  subplot 2
ax[1].plot(f.mixture_fraction("H"), f.T, label=f"new")
ax[1].plot(f2.mixture_fraction("H"), f2.T, label=f"old", linestyle="--")

ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("T")
ax[1].set_xlabel("<- ox z fuel ->")
# Fig1  subplot 3
ax[2].plot(f.mixture_fraction("H"), f.enthalpy_mass, label=f"new")
ax[2].plot(f2.mixture_fraction("H"), f2.h, label=f"old", linestyle="--")

ax[2].grid()
ax[2].legend()
ax[2].set_ylabel("h in  J/kg")
ax[2].set_xlabel("<- ox z fuel ->")

# Fig1  subplot 3
ax[3].plot(f.grid, f.enthalpy_mass, label=f"new", marker="x")
ax[3].plot(f2.grid, f2.h, label=f"old", linestyle="--", marker="x")
ax[3].grid()
ax[3].legend()
ax[3].set_ylabel("h in  J/kg")
ax[3].set_xlabel("grid")



idx_H2 = f.gas.species_index("H2")
idx_O2 = f.gas.species_index("O2")
idx_OH = f.gas.species_index("OH")
fig2, ax2 = plt.subplots(3, 1)

# Fig 2 species subplot 1
ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label="new")
ax2[0].plot(f2.mixture_fraction("H"), f2.Y[idx_H2], label="old")

ax2[0].grid()
ax2[0].legend()
ax2[0].set_ylabel("H2")
ax2[0].set_xlabel("<- fuel z ox ->")

# Fig 2 species suplot 2
ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], label="new")
ax2[1].plot(f2.mixture_fraction("H"), f2.Y[idx_O2], label="old")

ax2[1].grid()
ax2[1].legend()
ax2[1].set_ylabel("O2")
ax2[1].set_xlabel("<- fuel x ox ->")

# Fig 2 species suplot 3
ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_OH], label="new")
ax2[2].plot(f2.mixture_fraction("H"), f2.Y[idx_OH], label="old")

ax2[2].grid()
ax2[2].legend()
ax2[2].set_ylabel("OH")
ax2[2].set_xlabel("<- fuel x ox ->")


plt.show()
