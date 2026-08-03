# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
Unstable branch continuation
=============================

Traces the unstable branch of the S-curve via two-point control, starting
from stable.py's own extinction-point output so it continues from exactly
where the strain-rate sweep left off.

Requires: cantera >= 3.2, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, diffusion flame, strained flame,
          unstable branch, two-point control
"""

import numpy as np
import cantera as ct
import h5py
import pandas as pd
import utils as ut
from matplotlib import pyplot as plt

# Flame settings
reaction_mechanism = "h2o2.yaml"
gas = ct.Solution(reaction_mechanism)
width = 18e-3
grid = np.linspace(0, width, 300)
f = ct.CounterflowDiffusionFlame(gas, grid=grid)
f.P = 1.e5
f.fuel_inlet.X = "H2:1"
f.oxidizer_inlet.X = "O2:1"
f.fuel_inlet.T = 300
f.oxidizer_inlet.T = 300

 
# Initialize from stable.py's own extinction-point output, so the unstable
# branch continues from exactly where the strain-rate sweep left off.

start_path = "Data/"
save_path = "Data/"

files = [
       #  "stable_Z0.4500.h5", 
        "stable_Z0.5000.h5",
         #"stable_Z0.2500.h5",  
         ]
start_files = [start_path + file for file in files]    
for start_file in start_files:
    h5_file = h5py.File(start_file, "r")
    keys = ["extinction/" + key for key in sorted(h5_file["extinction"].keys())]
    
    Z_wall_val = h5_file[keys[0]]["flame/z"].attrs["z_wall"]
    file_tag = f"Z{Z_wall_val:.4f}"                                             
    save_file = save_path + f"3unstable_{file_tag}.h5"


    # stable.py's np.isclose extinction check has a tight tolerance, so the last
    # few saved keys can already be extinguished (T_max at ambient) before the
    # abort criterion actually triggers -- using keys[-1] directly would start
    # the two-point control from a cold, non-burning state. Search backward for
    # the last key that's still meaningfully burning instead.
    inlet_T = max(f.fuel_inlet.T, f.oxidizer_inlet.T)
    last_run = keys[-1]
    for key in reversed(keys):
        if h5_file[key]["flame"]["T"][:].max() > inlet_T + 100:
            last_run = key
            break
    else:
        raise RuntimeError(f"No burning state found in {start_file}; all saved "
                            "states are near ambient temperature.")
    first_run = keys[0]

    # Restore first solution
    f.restore(start_file, first_run)

    a_max = strain_rate = f.strain_rate("max")

    wall_params = {
        'Z_wall': h5_file[last_run]["flame"]["z"].attrs["z_wall"],
        'Z_wall_width': 0.001,
        'T_wall': 300.0,
        'factor': h5_file[last_run]["flame"]["z"].attrs["factor"],
        'mix_frac': 'H',
        'fuel': 'H2',
        'oxidizer': 'O2',
        'basis': 'mass'
    }

    z_stoich = ut.get_z_stoich(gas, wall_params, reaction_mechanism)

    # Restore last (most-strained, closest to extinction) solution
    f.restore(start_file, last_run)

    f.transport_model = "unity-Lewis-number"


    ##############################################################################
    # Flame Continuation
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    # Maximum number of steps to take
    n_max = 500

    # Relative temperature defining control point locations, with 1 being the peak
    # temperature and 0 being the inlet temperature. Lower values tend to avoid solver
    # failures early on, while using higher values on the unstable branch tend to help
    # with finding solutions where the peak temperature is very low.
    initial_spacing = 0.6
    unstable_spacing = 0.95

    # Amount to adjust temperature at the control point each step [K]
    temperature_increment = 20
    max_increment = 100

    # Try to keep T_max from changing more than this much each step [K]
    target_delta_T_max = 10.

    # Stop after this many successive errors
    max_error_count = 5
    error_count = 0

    # Stop after any failure if the strain rate has dropped to this fraction of the maximum
    strain_rate_tol = 0.1

    f.two_point_control_enabled = True
    # Prevent two point control from finding solutions with negative inlet velocities
    f.flame.set_bounds(spread_rate=(-1e-5, 1e20))
    f.max_time_step_count = 1000
    T_max = max(f.T)
    print("tmax: ", T_max)
    data = []  # integral output quantities for each step
    solved = False

    factor_last_working = wall_params["factor"]
    refine_grid = True
    auto = False
    enthalpy_refinement = True
    enthalpy_curve = 0.05
    f.set_refine_criteria(ratio=3.0, slope=0.5, curve=0.05, prune=0.01, enthalpy=enthalpy_refinement, enthalpy_curve=enthalpy_curve)
    for i in range(n_max):
        print(f"i = {i}")
        if strain_rate > 0.98 * a_max:
            spacing = initial_spacing
        else:
            spacing = unstable_spacing

        control_temperature = np.min(f.T) + spacing * (np.max(f.T) - np.min(f.T))

        # Store the flame state in case the iteration fails and we need to roll back
        backup_state = f.to_array()

        print(f"T_max: {np.max(f.T)}")
        print(f'Iteration {i}: Control temperature = {control_temperature:.2f} K')

        f.set_left_control_point(control_temperature)
        f.set_right_control_point(control_temperature)

        print(f"Temperature increment = {temperature_increment}")

        ## This decrement is what drives the two-point control. If failure
        ## occurs, try decreasing the decrement.
        if i == 1:
            f.right_control_point_temperature -= 1.
            f.left_control_point_temperature -= 1.
        else:
            f.right_control_point_temperature -= temperature_increment
            f.left_control_point_temperature -= temperature_increment
        f.clear_stats()

        if (f.left_control_point_temperature < f.fuel_inlet.T + 100
            or f.right_control_point_temperature < f.oxidizer_inlet.T + 100
        ):
            print("SUCCESS! Stopping because control point temperature is "
                        "sufficiently close to inlet temperature.")
            break

        try:
            runtime, factor_last_working = ut.solve_with_wall(f, wall_params, factor_last_working=factor_last_working,
                                                     factor_increase=1.1, factor_decrease=0.9, delta_T_max=1.0, loglevel=1, refine_grid=refine_grid, auto=auto)

            print(f"Factor last working: {factor_last_working}")
            if abs(max(f.T) - T_max) < 0.8 * target_delta_T_max:
                # Max temperature is changing slowly. Try a larger increment next step
                temperature_increment = min(temperature_increment + 3, max_increment)
            elif abs(max(f.T) - T_max) > target_delta_T_max:
                # Max temperature is changing quickly. Scale down increment for next step
                temperature_increment *= 0.9 * target_delta_T_max / (abs(max(f.T) - T_max))
            solved = True
            error_count = 0
        except (ut.SolveFailure, ct.CanteraError) as err:
            solved = False
            print(err)
            if strain_rate / a_max < strain_rate_tol:
                print('SUCCESS! Traversed unstable branch down to '
                            f'{100 * strain_rate / a_max:.2f}% of the maximum strain rate.')
                break
            # Restore the previous solution and try a smaller temperature increment for the
            # next iteration
            f.from_array(backup_state)
            temperature_increment = 0.7 * temperature_increment
            error_count += 1
            print(f"Did not converge on iteration {i}. New dT = {temperature_increment:.2f}")

        if solved:
            print(f"Saving on iteration{i}")
            name = f"iteration{i}"
            ut.save_with_attributes(f, save_file, name, wall_params, z_stoich, enthalpy_refinement=enthalpy_refinement, enthalpy_curve=enthalpy_curve, info=True)

            # Collect output stats
            T_max = max(f.T)
            T_mid = 0.5 * (min(f.T) + max(f.T))
            s = np.where(f.T > T_mid)[0][[0, -1]]
            width = f.grid[s[1]] - f.grid[s[0]]
            strain_rate = f.strain_rate('max')
            a_max = max(strain_rate, a_max)

            data.append({
                'T_max': max(f.T),
                'strain_rate': strain_rate,
                'heat_release_rate': trapezoid(f.heat_release_rate, f.grid),
                'n_points': len(f.grid),
                'flame_width': width,
                'Tc_increment': temperature_increment,
                'time_steps': sum(f.time_step_stats),
                'eval_count': sum(f.eval_count_stats),
                'cpu_time': sum(f.jacobian_time_stats + f.eval_time_stats),
                'errors': error_count
            })

            df = pd.DataFrame.from_records(data)

        if error_count >= max_error_count:
            print(f'FAILURE! Stopping after {error_count} successive solver '
                           'errors.')
            break

    print(f'Stopped after {i} iterations')

if data:
    plt.figure()
    plt.plot(df.strain_rate, df.T_max)
    plt.xlabel('a_max')
    plt.ylabel('T_max')
    plt.show()
else:
    print("No successful iterations to plot.")
