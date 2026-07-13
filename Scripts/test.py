# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
Diffusion flame extinction strain rate
======================================

This example computes the extinction point of a counterflow diffusion flame.
A hydrogen-oxygen diffusion flame at 1 bar is studied.

The tutorial makes use of the scaling rules derived by Fiala and Sattelmayer
(doi:10.1155/2014/484372). Please refer to this publication for a detailed
explanation. Also, please don't forget to cite it if you make use of it.

Requires: cantera >= 3.2, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, diffusion flame, strained flame, extinction,
          saving output, plotting
"""

from pathlib import Path
import numpy as np
import h5py
from scipy import special
import cantera as ct
import utils as ut
from matplotlib import pyplot as plt


# Flame settings (shared across all Z_wall runs)
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3
output_path = Path("./Data")
output_path.mkdir(parents=True, exist_ok=True)

delta_T_min = 1.


# Create fresh flame object
grid = np.linspace(0, width, 250)
f = ct.CounterflowDiffusionFlame(gas, grid=grid)
f.P = 1.e5
f.fuel_inlet.mdot = 0.1
f.oxidizer_inlet.mdot = 0.5
f.fuel_inlet.X = "H2:1"
f.oxidizer_inlet.X = "O2:1"
f.fuel_inlet.T = 300
f.oxidizer_inlet.T = 300
f.transport_model = "unity-Lewis-number"
z_wall = 0.98
file_tag = z_wall
file_path = str(output_path / f"test_{file_tag}.h5")

# Wall parameters for this run
wall_params = {
    'Z_wall': z_wall,
    'T_wall': 300.0,
    'factor': 1000,
    'mix_frac': 'H',
    'fuel': 'H2',
    'oxidizer': 'O2',
    'basis': 'mass'
}

z_stoich_val = ut.get_z_stoich(gas, wall_params, reaction_mechanism)
# Output files with Z_wall in the name
initial_file = "Data/stable_Z1.0000.h5"

h5_file = h5py.File(initial_file, "r")
runs = list(h5_file["extinction"].keys())
runs = ["0039"]

for run in runs:
    print(initial_file, run)
    f.restore(initial_file, "extinction/" + run)
    plt.plot(f.mixture_fraction("H"), f.enthalpy_mass)
    T, Y, h, z, a_m, idx_C = ut.correct_enthalpy(initial_file, "extinction/" + run, "H2O", style="line")         
    f.flame.set_values("T", T)
    plt.plot(z,h)
    plt.show()
    factor_last_working = h5_file["extinction"][run]["flame"]["z"].attrs["factor"]

try:

    f.set_refine_criteria(ratio=3, slope=0.05, curve=0.005, prune=0.003,
                          enthalpy=True, enthalpy_curve=0.005)
    runtime, _ =ut.solve_with_wall(f, wall_params,
                       factor_last_working=factor_last_working,
                       delta_T_max=1., loglevel=1, factor_increase=1.5, refine_grid=True, auto=True)
    ut.save_with_attributes(f, file_path, run, wall_params, z_stoich_val, info=True, runtime=runtime)
except (ut.SolveFailure, ct.CanteraError) as e:
    ut.print_r("Failed!")
    print(e)


h5_file.close()
