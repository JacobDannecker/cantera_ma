import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import h5py


file_list = ["Testscript/Data/Run1.h5", "Testscript/Data/Run2.h5", "Testscript/Data/Run3.h5"]
quantity = ["T"]
z_wall_spec = 0.6
z_stoich = 0.111

n_x = 200 
n_y = 20 



def interpolateQuantity(file_list, quantity="T", z_wall_spec="no_wall", z_stoich="0.111", n_x=200, n_y=50):
    n_files = len(file_list)
    z_all_runs = np.zeros((n_files, n_x))
    quantity_all_runs = np.zeros((n_files, n_x))
    chi_st_all_runs = np.zeros(n_files)
    
    for i, file in enumerate(file_list):
        h5_file = h5py.File(file, "r") 
        runs = [str(name) for name in h5_file.keys()]
        for run in runs:
            if "no_wall" in run:
                z = "no_wall"
            else:
                z = round(float(run.strip("z_wall_")), 2)
            if z == z_wall_spec:
                print(f"run: {run}")
                raw_grid = h5_file[run]["flame/grid"]
                raw_z = h5_file[run]["flame/z"]
                raw_data = h5_file[run]["flame"][quantity]

                interpol_grid_x = np.linspace(0, raw_grid[-1], n_x)
                
                interpol_func_z = interpolate.interp1d(raw_grid, raw_z)
                z_interpolated = interpol_func_z(interpol_grid_x)

                interpol_func_quantity = interpolate.interp1d(raw_grid, raw_data)
                quantity_interpolated = interpol_func_quantity(interpol_grid_x)
                
                z_all_runs[i] = z_interpolated
                quantity_all_runs[i] = quantity_interpolated 
                chi_st_all_runs[i] = h5_file[run]["flame/z"].attrs["chi_st"]
        h5_file.close()
    
    # Sort arrays according to asccending phi_stoich
    idc_sorted = np.argsort(chi_st_all_runs)

    chi_st_sorted = chi_st_all_runs[idc_sorted]
    quantity_sorted = quantity_all_runs[idc_sorted, :]
    z_sorted = z_all_runs[idc_sorted, :]

    print(quantity_sorted.shape)
    print(z_sorted.shape)
    # Inerpolate on grid
    x_raw = z_interpolated[0] 
    y_raw = chi_st_sorted
    # Flipp ggf. enfernen.
    interpol_func_grid = RegularGridInterpolator((interpol_grid_x, y_raw), np.flip(quantity_sorted.T), bounds_error=False, fill_value=None)
    interpol_grid_y = np.linspace(np.min(chi_st_sorted), np.max(chi_st_sorted), n_y)
    X, Y = np.meshgrid(np.linspace(0, 1, 100), interpol_grid_y)
    quantity_interpol_grid = interpol_func_grid((X, Y))

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel("z")
    ax.set_ylabel(r"chi_st / $s^{-1}$")
    mesh = ax.pcolormesh(X, Y, quantity_interpol_grid)
    fig.colorbar(mesh, ax=ax, label=f"{str(quantity)}")
    ax.set_title("Title")
    for i, chi_st in enumerate(chi_st_sorted):
        print("Scatttter")
        ax.scatter(z_sorted[i, idc_sorted[i]], chi_st , marker="x", color="r", label=f"chi_stoich{chi_st_sorted[i]}")
        ax.legend()
    plt.show()
    # Save as csv

interpolateQuantity(file_list, z_wall_spec=0.5, n_x=300, n_y=20)
