import numpy as np
import cantera as ct
from scipy import special
from matplotlib import pyplot as plt
import time 
import h5py

def chi_stoich(f, z_stoich):
        #a = f.strain_rate("mean")
        a = np.mean(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))
        chi_stoich = a*np.pi*(np.exp(-2*((special.erfinv(1-2*z_stoich))**2)))
        return chi_stoich


def add_attributes(f, file_path, wall_params, z_stoich):
    z_array = f.mixture_fraction(wall_params["mix_frac"])
    h_mass_array = f.enthalpy_mass
    h_mole_array = f.enthalpy_mole
    hdf5_file = h5py.File(file_path, "a")
    group_z = name + "/flame/z"
    group_h_mass = name + "/flame/h_mass"
    group_h_mole = name + "/flame/h_mole"
    # Data
    hdf5_file.create_dataset(name=group_z, data=z_array)
    hdf5_file.create_dataset(name=group_h_mass, data=h_mass_array)
    hdf5_file.create_dataset(name=group_h_mole, data=h_mole_array)
    # Attributes
    hdf5_file[group_z].add_attr = "mix_frac"
    hdf5_file[group_z].add_attr = "fuel"
    hdf5_file[group_z].add_attr = "oxidizer"
    hdf5_file[group_z].add_attr = "basis"
    hdf5_file[group_z].add_attr = "chi_st"
    hdf5_file[group_z].attrs["mix_frac"] = wall_params["mix_frac"]
    hdf5_file[group_z].attrs["fuel"] = wall_params["fuel"]
    hdf5_file[group_z].attrs["oxidizer"] = wall_params["oxidizer"]
    hdf5_file[group_z].attrs["basis"] = wall_params["basis"]
    hdf5_file[group_z].attrs["chi_st"] = chi_stoich(f, z_stoich)
    hdf5_file.close()
 
# Flame settings
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3  
grid = np.linspace(0, width, 150)
f = ct.CounterflowDiffusionFlame(gas, grid=grid)
f.P = 1.e5  
f.fuel_inlet.mdot = 0.5 
f.fuel_inlet.X = "H2:1"
f.fuel_inlet.T = 300 
f.oxidizer_inlet.mdot = 3.0 
f.oxidizer_inlet.X = "O2:1"
f.oxidizer_inlet.T = 300 
z_stoich = 0.111
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

# Wall 
wall_params = {
    'Z_wall': 1,
    'T_wall': 300.0,
    'factor': 1,
    'mix_frac': 'Bilger',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}

# Names
name = "no_wall"
names = [name,]
file_path = "Testscript/Data/RunEnth1.h5" 

# Initial solution no wall
f.transport_model = "unity-Lewis-number"
f.solve(loglevel=0, refine_grid=True)
f.save(file_path, name=name, overwrite=True)
add_attributes(f, file_path, wall_params, z_stoich)

# Loop settings 
z = 1
delta_z = 0.1
last_z_working = 1
error_counter = 0
max_errors = 3
failed_once = False
max_errors_reached = False
flame_is_extinct = False
factor_increase = 2 
failed_z = []

start_time = time.time()
while True: 
    z -= delta_z
    if z < 0.1: # Add check for extinct flame
        break
    if flame_is_extinct:
        break
    # Try calculating initial guess with factor = 1
    try:
        wall_params["Z_wall"] = z 
        wall_params["factor"] = 1 
        f.set_initial_guess(data=file_path, group=names[0]) 
        f.flame.set_non_adiabatic_wall(wall_params)
        f.transport_model = "unity-Lewis-number"
        print(f"Solving for z: {z} with factor = 1") 
        f.solve(loglevel=0, fefine_grid=True)
        print(f"Success for initial solution with factor = 1 at z: {z}")
    except:
        print(f"Failed initial solution with factor = 1 at z: {z}\n Conitune with increasing factor")
        pass

    
    wall_params["factor"] = 100
    while True:
        try:
            wall_params["factor"] *= factor_increase
            f.flame.set_non_adiabatic_wall(wall_params)
            f.solve(loglevel=0, refine_grid=True)
            idx_wall = np.abs(f.mixture_fraction("H") - wall_params["Z_wall"]).argmin()
            delta_T_wall =  f.T[idx_wall] - wall_params["T_wall"]
            print(f"Succes in factor increase loop. Delta T wall: {delta_T_wall}")
            if (np.max(f.T) < 310):
                print("Flame is extinct!")
                flame_is_extinct=True
        except BaseException as err:
            error_counter += 1
            print(f"Errors ins z = {z}: {error_counter}")
            print(err)
            #wall_params["factor"] /= factor_increase
            if factor_increase > 1.2:
                factor_increase *= 0.9
            if error_counter > max_errors:
                print(f"max_errors_reached true, error_counter: {error_counter}")
              #  # Reset factor 
              #  wall_params["factor"] /= factor_increase
              #  # Decrease factor_increase
              #  if factor_increase > 3
              #      factor_increase -= 1 
                # Try starting from no_wall initial solution
                f.set_initial_guess(data=file_path, group=names[0]) 
                wall_params["factor"] = 10
                factor_increase = 2 
                error_counter = 0
                if failed_once:
                    # Abort when starting form no_wall solution also fails
                    print(f"No solution found at z{z}")
                    failed_z.append(z)
                    failed_once = False
                    break
                print("Try with initial solution as ininital guess.")
                failed_once  = True
                continue
        if delta_T_wall < 1:
            # Solution is close enough. Try next z.
            error_counter = 0
            print("\n##############################################################################")
            print(f"Solved at z: {z}")
            print("##############################################################################")
            name = "z_wall_" + str(z)
            f.save(file_path, name=name, overwrite=True)
            add_attributes(f, file_path, wall_params, z_stoich)
            names.append(name)
            break

end_time = time.time()
print(f"Total time loop: {end_time-start_time}")
print(f"Failed at z: {failed_z}")




# Plot
fig, ax = plt.subplots(2, 1)
fig.suptitle("H2/O2") 
for name in names:
    f.restore(file_path, name=name)
    # Subplot 1 Temp 
    if name == "no_wall":
        z = 1
    else:
        z = float(name.strip("z_wall_"))
    ax[0].plot(f.mixture_fraction("H"), f.T, label=f"{z:.2f} chi_st: {chi_stoich(f, z_stoich):.2f}")
    # Subplot 2  enthalpy
    ax[1].plot(f.mixture_fraction("H"), f.h, label=f"{z:.2f} chi_st: {chi_stoich(f, z_stoich):.2f}")

ax[0].grid()
ax[0].legend()
ax[0].set_ylabel("T")
ax[0].set_xlabel("<- ox z fuel ->")

ax[1].grid()
ax[1].legend()
ax[1].set_ylabel("h in  J/kg")
ax[1].set_xlabel("<- ox z fuel ->")

plt.show()
