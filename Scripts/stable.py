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
import cantera as ct
import utils as ut


# Flame settings (shared across all Z_wall runs)
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3
output_path = Path("./Data")
output_path.mkdir(parents=True, exist_ok=True)

# Exponents for scaling rules (Fiala and Sattelmayer 2014)
exp_d_a = -1. / 2.
exp_u_a = 1. / 2.
exp_V_a = 1.
exp_lam_a = 2.
exp_mdot_a = 1. / 2.

delta_alpha_min = 0.001
delta_T_min = 1.
Z_wall_values = [0.90]

print(Z_wall_values)

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
    enthalpy_refinement = True
    enthalpy_curve = 0.005
    f.set_refine_criteria(ratio=3, slope=0.5, curve=0.005, prune=0.003,
                          enthalpy=enthalpy_refinement, enthalpy_curve=enthalpy_curve)

    # Wall parameters for this run. Z_wall_width is the half-width of the
    # compact mixture-fraction window over which the wall sink is smoothly
    # blended in (see Flow1D::evalEnergy); factor is managed by
    # utils.solve_with_wall below.
    wall_params = {
        'Z_wall': Z_wall_val,
        'Z_wall_width': 0.001,
        'T_wall': 300.0,
        'mix_frac': 'H',
        'fuel': 'H2',
        'oxidizer': 'O2',
        'basis': 'mass'
    }

    z_stoich_val = ut.get_z_stoich(gas, wall_params, reaction_mechanism)

    # Output files with Z_wall in the name
    file_tag = f"Z{Z_wall_val:.4f}"
    file_path = str(output_path / f"test_stable_{file_tag}.h5")

    # Tracks the group name (within file_path) of the last-known-good state
    # for each step index n, so a failed step can restore the right one.
    name = "initial"
    names = [name]

    temperature_limit_extinction = max(f.oxidizer_inlet.T, f.fuel_inlet.T)

    f.set_initial_guess(data="Scripts/initial.h5", group="initial")

    # --- Establish the wall at the base (unstrained) state via factor search. ---
    factor_last_working = False
    try:
        runtime, factor_last_working = ut.solve_with_wall(
            f, wall_params, factor_last_working=factor_last_working,
            delta_T_max=0.1, loglevel=0, factor_increase=2., refine_grid=True, auto=True)
    except (ut.SolveFailure, ct.CanteraError) as e:
        failure_type = getattr(e, 'failure_type', ut.classify_failure(str(e)))
        ut.print_r(f"Could not establish the wall at Z_wall={Z_wall_val} "
                   f"(type={failure_type}): {e}")
        continue
    ut.print_g(f"Wall established at Z_wall={Z_wall_val}, factor={factor_last_working:.3g}. "
               f"T_max={np.max(f.T):.1f} K, n_grid={len(f.grid)}")
    ut.save_with_attributes(f, file_path, "initial", wall_params, z_stoich_val, enthalpy_refinement=enthalpy_refinement, enthalpy_curve=enthalpy_curve, runtime=runtime, info=True)

    # --- Extinction loop ---
    # last_good_alpha is the alpha value of the last converged state (the
    # baseline the next strain step is taken from). alpha[] below is a log
    # of every *attempted* value (successful or not) and must not be
    # indexed by n_last_burning: it grows by one entry every iteration
    # regardless of outcome, while n_last_burning only advances on success,
    # so alpha[n_last_burning] silently points at a stale entry as soon as
    # any failure has occurred.
    last_good_alpha = 1.
    alpha = [last_good_alpha]
    # delta_alpha only shrinks on failure below; grown back gently (rather
    # than reset outright) after each success so an early rough patch
    # doesn't force tiny steps for the rest of the run.
    delta_alpha_start = 20
    delta_alpha = delta_alpha_start
    growth_factor = 1.3
    n = 0
    n_last_burning = 0
    T_max = [np.max(f.T)]
    a_max = [np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid)))]
    refine_grid = True
    auto = False
    while True:
        n += 1
        alpha.append(last_good_alpha + delta_alpha)
        print(f'Proceeding, delta_alpha = {delta_alpha}')
        strain_factor = alpha[-1] / last_good_alpha
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
            runtime, factor_last_working = ut.solve_with_wall(f, wall_params,
                               factor_last_working=factor_last_working,
                               delta_T_max=0.1, loglevel=0, factor_increase=2., refine_grid=refine_grid, auto=auto)
            solved = True

        except (ut.SolveFailure, ct.CanteraError) as e:
            failure_type = getattr(e, 'failure_type', ut.classify_failure(str(e)))
            ut.print_r(f"Error: Did not converge at n = {n}, type={failure_type}")

        if solved:
            refine_grid = True
            t_max_val = float(np.max(f.T))
            a_max_val = float(np.max(np.abs(np.gradient(f.velocity) / np.gradient(f.grid))))
            T_max.append(t_max_val)
            a_max.append(a_max_val)
            if not np.isclose(t_max_val, temperature_limit_extinction):
                n_last_burning = n
                last_good_alpha = alpha[-1]
                name = f"extinction/{n:04d}"
                names.append(name)
                ut.save_with_attributes(f, file_path, name, wall_params, z_stoich_val, enthalpy_refinement=enthalpy_refinement, enthalpy_curve=enthalpy_curve, runtime=runtime, info=True)
                ut.print_y(f"Flame burning at alpha = {last_good_alpha:8.4f}.")
                delta_alpha = min(delta_alpha * growth_factor, delta_alpha_start)
                continue
            else:
                print(f"Flame extinguished (solved, T~ambient) at alpha = {alpha[-1]:8.4f}")
        else:
            ut.print_r(f"Flame extinguished (solver failed: {failure_type}) at alpha = {alpha[-1]:8.4f}")

        if (delta_alpha < delta_alpha_min):
            name = f"extinction/{n:04d}"
            names.append(name)
            if solved and (T_max[-2] - T_max[-1] < delta_T_min):
                ut.save_with_attributes(f, file_path, name, wall_params, z_stoich_val, enthalpy_refinement=enthalpy_refinement, enthalpy_curve=enthalpy_curve, runtime=runtime, info=True)
            else:
                print("  (not saving — solver failed, no valid solution)")
            print(f'Flame extinguished at alpha = {alpha[-1]:8.4f}. '
                  'Abortion criterion satisfied.')
            break
        else:
            auto = False
            delta_alpha *= 0.5
            refine_grid = False
            n = n - 1
            print(f'Flame extinguished at alpha = {alpha[-1]:8.4f} (discarded). '
                  f'Restoring alpha = {last_good_alpha:8.4f} and '
                  f'trying delta_alpha = {delta_alpha}')
            print(f"Setting initial guess {names[-1]}")
            # names[n_last_burning] is the actual saved group for that step
            # (names[0] == "initial", not "extinction/0000").
            name = names[n_last_burning]
            print(f"Before f:{f.fuel_inlet.mdot} o:{f.oxidizer_inlet.mdot}")
            f.restore(file_path, name)
            print(f"After f:{f.fuel_inlet.mdot} o:{f.oxidizer_inlet.mdot}")

    # --- Results for this Z_wall ---
    name = names[n_last_burning]
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
