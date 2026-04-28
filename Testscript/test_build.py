import numpy as np
import cantera as ct
from scipy import special
from matplotlib import pyplot as plt
import time

reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3  

f = ct.CounterflowDiffusionFlame(gas, width=width)
f2 = ct.CounterflowDiffusionFlame(gas, width=width)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.5  
f.fuel_inlet.X = "H2:1"
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 3.0 
f.oxidizer_inlet.X = "O2:1"
f.oxidizer_inlet.T = 300 
z_stoich = 0.111

f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
file_name = "Testscript/Data/test.h5"
f.set_initial_guess(data=file_name, group="0400") 
# Set up wall 
params = {
    'Z_wall': 0.4,
    'T_wall': 300,
    'factor': 1550000,
    'mix_frac': 'Bilger',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}
f.flame.set_non_adiabatic_wall(params)
f.transport_model = "unity-Lewis-number"


start = time.time()
f.solve(loglevel=1, refine_grid=True)
end = time.time()
print(end - start)

print(f.transport_model)
f.save(file_name, name="0400", overwrite=True)
f2.restore(file_name, name="no_wall")


def chi_stoich(f, z_stoich):
        #a = f.strain_rate("mean")
        a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
        chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
        return chi_stoich

chi_st_new = chi_stoich(f, z_stoich) 
chi_st_ref = chi_stoich(f2, z_stoich) 

# Info
idx_wall = np.abs(f.mixture_fraction("H") - params["Z_wall"]).argmin()
print(f"delta T wall: {f.T[idx_wall] - params["T_wall"]}")
print(f"mdot fuel new build: {f.fuel_inlet.mdot}")
print(f"mdot ox new build: {f.oxidizer_inlet.mdot}")
print("-----------------------------------------------")
print(f"chi_st new build: {chi_st_new}")
print(f"mdot fuel reference: {f2.fuel_inlet.mdot}")
print(f"mdot ox reference: {f2.oxidizer_inlet.mdot}")
print(f"chi_st reference: {chi_st_ref}")


idx_H2 = f.gas.species_index("H2")
idx_O2 = f.gas.species_index("O2")

# Fig 1 temp subplot 1
fig, ax = plt.subplots(3, 1)
fig.suptitle(" H2/O2") 

ax[0].plot(f.grid, f.T, label=f"With Wall chi_st: {chi_st_new:.2f}")
#ax[0].vlines(f.grid[wall_pos], 0, 3000, color="b", linestyle="-.")
ax[0].plot(f2.grid, f2.T, label=f"No Wall chi_st: {chi_st_ref:.2f}", linestyle="--")
ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("T")
ax[0].set_xlabel("<- fuel x ox->")
# Fig 1 temp subplot 2
ax[1].plot(f.mixture_fraction("H"), f.T, label=f"With Wall chi_st: {chi_st_new:.2f}")
ax[1].plot(f2.mixture_fraction("H"), f2.T, label=f"No Wall chi_st: {chi_st_ref:.2f}", linestyle="--")
ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("T")
ax[1].set_xlabel("<- ox z fuel ->")
# Fig1 enthalpy subplot 3
ax[2].plot(f.mixture_fraction("H"), f.h, label=f"With Wall chi_st: {chi_st_new:.2f}")
ax[2].plot(f2.mixture_fraction("H"), f2.h, label=f"No Wall chi_st: {chi_st_ref:.2f}", linestyle="--")
#ax[2].scatter(f.mixture_fraction("H"), f.h)
#f2.mixture_fraction("H"),
ax[2].grid()
ax[2].legend()
ax[2].set_ylabel("h in  J/kg")
ax[2].set_xlabel("<- ox z fuel ->")




# Fig 2 species subplot 1
idx_H2 = f.gas.species_index("H2")
idx_O2 = f.gas.species_index("O2")
idx_OH = f.gas.species_index("OH")
fig2, ax2 = plt.subplots(3, 1)

 #ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label="H2 With Wall")
ax2[0].plot(f.mixture_fraction("H"), f.Y[idx_H2], label="H2 With Wall")
 #ax2[0].plot(f2.mixture_fraction("H"), f2.Y[idx_H2], label="H2 No Wall")
ax2[0].plot(f2.mixture_fraction("H"), f2.Y[idx_H2], label="H2 No Wall")

ax2[0].grid()
ax2[0].legend()
ax2[0].set_ylabel("H2")
ax2[0].set_xlabel("<- fuel z ox ->")


# Fig 2 species suplot 2
 #ax2[1].plot(f.mixture_fraction("h"), f.Y[idx_O2], label="O2 newbuild")
ax2[1].plot(f.mixture_fraction("H"), f.Y[idx_O2], label="O2 newbuild")
 #ax2[1].plot(f2.mixture_fraction("H"), f2.Y[idx_O2], label="O2 No Wall")
ax2[1].plot(f2.mixture_fraction("H"), f2.Y[idx_O2], label="O2 No Wall")

ax2[1].grid()
ax2[1].legend()
ax2[1].set_ylabel("O2")
ax2[1].set_xlabel("<- fuel x ox ->")

# Fig 2 species suplot 3
 #ax2[2].plot(f.mixture_fraction("h"), f.Y[idx_OH], label="OH newbuild")
ax2[2].plot(f.mixture_fraction("H"), f.Y[idx_OH], label="OH newbuild")
 #ax2[2].plot(f2.mixture_fraction("H"), f2.Y[idx_OH], label="OH No Wall")
ax2[2].plot(f2.mixture_fraction("H"), f2.Y[idx_OH], label="OH No Wall")

ax2[2].grid()
ax2[2].legend()
ax2[2].set_ylabel("OH")
ax2[2].set_xlabel("<- fuel x ox ->")

plt.show()
