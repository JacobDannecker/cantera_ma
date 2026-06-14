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



# Flame settings (shared across all Z_wall runs)
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3
output_path = Path("Scripts/Data")
output_path.mkdir(parents=True, exist_ok=True)

# Exponents for scaling rules (Fiala and Sattelmayer 2014)
exp_d_a = -1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

delta_alpha_factor = 50.
delta_alpha_min = 0.001

# Z_wall values to scan: from 1.0 down to 0.1115 in steps of 0.5
Z_wall_values = np.arange(0.90, 0.e2, -0.02)

for Z_wall_val in Z_wall_values:
    print(f"\n{'='*60}")
    print(f"Computing for Z_wall = {Z_wall_val}")
    print('='*60)

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
    f.set_refine_criteria(ratio=3, slope=0.5, curve=0.05, prune=0.04,
                          enthalpy=False, enthalpy_curve=0.05)

    # Wall parameters for this run
    wall_params = {
        'Z_wall': Z_wall_val,
        'T_wall': 300.0,
        'factor': 1000,
        'mix_frac': 'H',
        'fuel': 'H2',
        'oxidizer': 'O2',
        'basis': 'mass'
    }

    z_stoich_val = ut.get_z_stoich(gas, wall_params, reaction_mechanism)

    # Output files with Z_wall in the name
    file_tag = f"Z{Z_wall_val:.4f}"
    file_path = str(output_path / f"Xextinction_{file_tag}.h5")

    name = "initial"
    names = [name]

    temperature_limit_extinction = max(f.oxidizer_inlet.T, f.fuel_inlet.T)

    f.set_initial_guess(data="Scripts/Data/enthalpy.h5", group="initial")
    ut.save_with_attributes(f, file_path, name, wall_params, z_stoich_val, info=True)

    # --- Extinction loop ---
    alpha = [1.]
    delta_alpha = 5
    n = 0
    n_last_burning = 0
    T_max = [np.max(f.T)]
    a_max = [np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))]

    while True:
        n += 1
        alpha.append(alpha[n_last_burning] + delta_alpha)
        strain_factor = alpha[-1] / alpha[n_last_burning]

        f.flame.grid *= strain_factor ** exp_d_a
        f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
        f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
        f.flame.set_values("velocity", f.flame.velocity * strain_factor ** exp_u_a)
        f.flame.set_values("spreadRate", f.flame.spread_rate * strain_factor ** exp_V_a)
        f.flame.set_values(
            "Lambda", f.flame.radial_pressure_gradient * strain_factor ** exp_lam_a)

        solved = False
        failure_type = None
        try:
            factor_last_working = 1000
            ut.solve_with_wall(f, wall_params, name_fallback=names[-1],
                               factor_last_working=factor_last_working,
                               delta_T_max=1.0, loglevel=0)
            solved = True
        except (ut.SolveFailure, ct.CanteraError) as e:
            failure_type = getattr(e, 'failure_type', ut.classify_failure(str(e)))
            print(f"Error: Did not converge at n = {n}, type={failure_type}")

        if solved:
            t_max_val = float(np.max(f.T))
            a_max_val = float(np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid))))
            T_max.append(t_max_val)
            a_max.append(a_max_val)
            print(f"MAX Temp {t_max_val}")

            if not np.isclose(t_max_val, temperature_limit_extinction):
                n_last_burning = n
                name = f"extinction/{n:04d}"
                names.append(name)
                ut.save_with_attributes(f, file_path, name, wall_params, z_stoich_val, info=True)
                print(f'Flame burning at alpha = {alpha[n_last_burning]:8.4f}. Proceeding, '
                      f'delta_alpha = {delta_alpha}')
                continue
            else:
                print(f"Flame extinguished (solved, T~ambient) at alpha = {alpha[-1]:8.4f}")
        else:
            print(f"Flame extinguished (solver failed: {failure_type}) at alpha = {alpha[-1]:8.4f}")

        if delta_alpha < delta_alpha_min:
            name = f"extinction/{n:04d}"
            names.append(name)
            if solved:
                ut.save_with_attributes(f, file_path, name, wall_params, z_stoich_val, info=True)
            else:
                print("  (not saving — solver failed, no valid solution)")
            print(f'Flame extinguished at alpha = {alpha[-1]:8.4f}. '
                  'Abortion criterion satisfied.')
            break
        else:
            delta_alpha = delta_alpha / delta_alpha_factor
            alpha.pop()
            n = n - 1
            print(f'Flame extinguished at alpha = {alpha[-1]:8.4f} (discarded). '
                  f'Restoring alpha = {alpha[n_last_burning]:8.4f} and '
                  f'trying delta_alpha = {delta_alpha}')
            print(f"Setting initial guess {names[-1]}")
            f.set_initial_guess(data=file_path, group=names[-1])

    # --- Results for this Z_wall ---
    name = f"extinction/{n_last_burning:04d}"
    f.restore(file_path, name)
    print(f"\n--- Results for Z_wall = {Z_wall_val} ---")
    print('----------------------------------------------------------------------')
    print('Parameters at the extinction point:')
    print('Pressure p={0} bar'.format(f.P / 1e5))
    print('Peak temperature T={0:4.0f} K'.format(np.max(f.T)))
    print('Mean axial strain rate a_mean={0:.2e} 1/s'.format(f.strain_rate('mean')))
    print('Maximum axial strain rate a_max={0:.2e} 1/s'.format(f.strain_rate('max')))
    print('Fuel inlet potential flow axial strain rate a_fuel={0:.2e} 1/s'.format(
          f.strain_rate('potential_flow_fuel')))
    print('Oxidizer inlet potential flow axial strain rate a_ox={0:.2e} 1/s'.format(
          f.strain_rate('potential_flow_oxidizer')))
    print('Axial strain rate at stoichiometric surface a_stoich={0:.2e} 1/s'.format(
          f.strain_rate('stoichiometric', fuel='H2')))

print("\nAll Z_wall runs completed.")

