#!/usr/bin/env python3
"""
settings.py
Print ALL settings used in a Cantera CounterflowDiffusionFlame simulation
for a H2/O2 flame with unity Lewis number.
"""

import cantera as ct
import numpy as np

# ===========================================================================
# 1.  GAS PHASE / MECHANISM
# ===========================================================================
gas = ct.Solution('h2o2.yaml')

gas.TP = 300, ct.one_atm

gas.transport_model = 'unity-Lewis-number'

print('=' * 72)
print('Cantera CounterflowDiffusionFlame  —  Settings')
print('=' * 72)
print()
print('--- Gas Phase / Mechanism ------------------------------------------')
print(f'  Mechanism file:        {gas.source}')
print(f'  Cantera version:       {ct.__version__}')
print(f'  Equation of state:     {gas.phase_of_matter}  (ideal gas)')
print(f'  Elements:              {", ".join(gas.element_names)}')
print(f'  Species ({gas.n_species}):        {", ".join(gas.species_names)}')
print(f'  Reactions ({gas.n_reactions}):')
for i in range(gas.n_reactions):
    r = gas.reaction(i)
    print(f'    {i:2d}:  {r.equation:40s}  type={r.reaction_type}')
print(f'  Initial T:             {gas.T} K')
print(f'  Initial P:             {gas.P / 1e5:.5f} bar  ({gas.P:.0f} Pa)')
print(f'  Thermodynamics:        NASA-7 coefficient polynomials (2 ranges)')
print(f'    Species  |  T_min [K]  |  T_max [K]  |  T_mid [K]')
for sp in gas.species():
    thermo = sp.thermo
    print(f'    {sp.name:8s}  |  {thermo.min_temp:9.1f}  |  {thermo.max_temp:9.1f}  |  {thermo.coeffs[0]:9.1f}')

print()
print('--- Transport ------------------------------------------------------')
print(f'  Model:                 {gas.transport_model}')
print(f'  Lewis numbers:         Le_k = 1  for all species')
print(f'  Schmidt numbers:       Sc_k = 1  for all species')
print(f'  Diffusivity:           D_k = λ / (ρ·c_p)  (thermal diffusivity)')

# ===========================================================================
# 2.  DOMAIN / GRID
# ===========================================================================
width = 0.02          # [m]
n_init = 12
grid = np.linspace(0, width, n_init)

print()
print('--- Domain (Grid) --------------------------------------------------')
print(f'  Domain width:          {width*1e3:.2f} mm')
print(f'  Initial grid points:   {n_init}')
print(f'  Initial spacing:       {np.diff(grid).mean()*1e6:.1f} µm  (uniform)')
print(f'  Initial grid:')
for i, x in enumerate(grid):
    print(f'    {i:2d}:  x = {x*1e3:8.4f} mm')

# ===========================================================================
# 3.  BOUNDARY CONDITIONS
# ===========================================================================
fuel_T = 300.0            # [K]
oxidizer_T = 300.0         # [K]
fuel_X = 'H2:1'
oxidizer_X = 'O2:1'
fuel_mdot = 0.5            # [kg/m²·s]
oxidizer_mdot = 0.5        # [kg/m²·s]

fuel_gas = ct.Solution('h2o2.yaml')
fuel_gas.TPX = fuel_T, ct.one_atm, fuel_X
oxidizer_gas = ct.Solution('h2o2.yaml')
oxidizer_gas.TPX = oxidizer_T, ct.one_atm, oxidizer_X

print()
print('--- Boundary Conditions --------------------------------------------')
print(f'  Fuel inlet  (x = 0 mm):')
print(f'    Temperature:         {fuel_T} K')
print(f'    Composition:         {fuel_X}')
print(f'    Mass fractions:')
for i, sp in enumerate(gas.species_names):
    if fuel_gas.Y[i] > 0:
        print(f'      {sp:6s}  =  {fuel_gas.Y[i]:.4f}')
print(f'    Density:             {fuel_gas.density:.4f} kg/m³')
print(f'    Mass flux:           {fuel_mdot} kg/m²·s')
print()
print(f'  Oxidizer inlet  (x = {width*1e3:.2f} mm):')
print(f'    Temperature:         {oxidizer_T} K')
print(f'    Composition:         {oxidizer_X}')
print(f'    Mass fractions:')
for i, sp in enumerate(gas.species_names):
    if oxidizer_gas.Y[i] > 0:
        print(f'      {sp:6s}  =  {oxidizer_gas.Y[i]:.4f}')
print(f'    Density:             {oxidizer_gas.density:.4f} kg/m³')
print(f'    Mass flux:           {oxidizer_mdot} kg/m²·s')

print()
print(f'  Pressure:              {ct.one_atm / 1e5:.5f} bar  ({ct.one_atm:.0f} Pa)')

# ===========================================================================
# 4.  FLAME OBJECT & SOLVER SETTINGS
# ===========================================================================
f = ct.CounterflowDiffusionFlame(gas, grid=grid)

f.fuel_inlet.X = fuel_X
f.fuel_inlet.T = fuel_T
f.oxidizer_inlet.X = oxidizer_X
f.oxidizer_inlet.T = oxidizer_T

f.fuel_inlet.mdot = fuel_mdot
f.oxidizer_inlet.mdot = oxidizer_mdot

f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1, prune=0.05)

f.flame.set_steady_tolerances(default=(1e-5, 1e-12))
f.flame.set_transient_tolerances(default=(1e-5, 1e-12))

f.fuel_inlet.set_steady_tolerances(default=(1e-5, 1e-8))
f.fuel_inlet.set_transient_tolerances(default=(1e-5, 1e-8))
f.oxidizer_inlet.set_steady_tolerances(default=(1e-5, 1e-8))
f.oxidizer_inlet.set_transient_tolerances(default=(1e-5, 1e-8))

f.energy_enabled = True
f.soret_enabled = False
f.radiation_enabled = False
f.boundary_emissivities = (0.0, 0.0)

refine = f.get_refine_criteria()
flame_settings = f.flame.get_settings3()

print()
print('--- Solver Settings -------------------------------------------------')
print(f'  Energy equation:       {f.energy_enabled}')
print(f'  Soret effect:          {f.soret_enabled}')
print(f'  Radiation:             {f.radiation_enabled}')
print(f'  Boundary emissivity:   fuel={f.boundary_emissivities[0]},'
      f' oxidizer={f.boundary_emissivities[1]}')
print(f'  Max grid points:       {flame_settings["refine-criteria"]["max-points"]}')
print(f'  Max time steps:        {f.max_time_step_count}')
print()
print(f'  Refinement criteria:')
print(f'    ratio:               {refine["ratio"]}')
print(f'    slope:               {refine["slope"]}')
print(f'    curve:               {refine["curve"]}')
print(f'    prune:               {refine["prune"]}')
print(f'    grid-min:            {flame_settings["refine-criteria"]["grid-min"]:.1e} m')
print()
print(f'  Tolerances  (rel, abs):')
tols = f.flame.get_settings3()["tolerances"]
print(f'    Flame  steady:       rel={tols["steady-reltol"]:.0e},'
      f' abs={tols["steady-abstol"]:.0e}')
print(f'    Flame  transient:    rel={tols["transient-reltol"]:.0e},'
      f' abs={tols["transient-abstol"]:.0e}')
print(f'    Inlets steady:       rel=1e-5, abs=1e-8')
print(f'    Inlets transient:    rel=1e-5, abs=1e-8')
print()
print(f'  Flux-gradient basis:   {flame_settings["flux-gradient-basis"]}')
print(f'  Solution components    ({f.flame.n_components}):')
for c in f.flame.component_names:
    print(f'    {c}')

# ===========================================================================
# 5.  SOLVE
# ===========================================================================
print()
print('--- Solving --------------------------------------------------------')
print('  Calling f.solve(loglevel=0, auto=True) ...')
f.solve(loglevel=0, auto=True)
print('  Done.')

# ===========================================================================
# 6.  RESULTS
# ===========================================================================
print()
print('--- Results --------------------------------------------------------')
print(f'  Peak temperature:            {max(f.T):.1f} K')
print(f'  Heat release rate (max):    {max(f.heat_release_rate) / 1e9:.3f} GW/m³')
print(f'  Final grid points:           {f.flame.n_points}')
print(f'  Final grid spacing (min):   {np.min(np.diff(f.grid))*1e6:.2f} µm')
print(f'  Final grid spacing (max):   {np.max(np.diff(f.grid))*1e6:.2f} µm')

# Strain rate
print(f'  Strain rate at extinction:  (run continuation to find)')
print(f'  Strain rate (max):           {f.strain_rate("max"):.1f} 1/s')
print(f'  Strain rate (mean):          {f.strain_rate("mean"):.1f} 1/s')

# Thermal flame thickness
dx = f.grid
T = f.T
grad = np.gradient(T, dx)
delta_T = max(T) - min(T)
dTdx_max = max(abs(grad))
thickness = delta_T / dTdx_max if dTdx_max > 0 else 0
print(f'  Flame thickness (thermal):   {thickness*1e3:.3f} mm')

# Mixture fraction
print(f'  Mixture fraction at flame:   {f.mixture_fraction("Bilger")[np.argmax(f.T)]:.4f}')

# Species peak values
print()
print(f'  Species peak mass fractions:')
for i, name in enumerate(gas.species_names):
    peak = max(f.Y[i, :])
    pos = f.grid[np.argmax(f.Y[i, :])]
    print(f'    {name:6s}:  {peak:.4e}  at x={pos*1e3:.2f} mm')

# Key consumption/production info
print()
print(f'  Integrated production rates:')
for i, name in enumerate(gas.species_names):
    net = np.trapezoid(f.net_production_rates[i, :] * gas.molecular_weights[i], f.grid)
    if abs(net) > 1e-20:
        print(f'    {name:6s}:  {net:.4e}  kg/m²·s')

print()
print('=' * 72)
print('End of settings.')
print('=' * 72)
