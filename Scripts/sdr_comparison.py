import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv

p = ct.one_atm
tin_f = 300.0
tin_o = 300.0
comp_f = 'H2:1'
comp_o = 'O2:0.21, N2:0.78, AR:0.01'
width = 0.03

gas = ct.Solution('gri30.yaml')

# --- Sweep over strain rates ---
mdot_f_values = [0.2, 0.3, 0.4, 0.5, 0.6]
mdot_ratio = 2.5  # mdot_o / mdot_f

results = []
last_f = None

for mdot_f in mdot_f_values:
    mdot_o = mdot_f * mdot_ratio
    f = ct.CounterflowDiffusionFlame(gas, width=width)
    f.fuel_inlet.mdot = mdot_f
    f.fuel_inlet.X = comp_f
    f.fuel_inlet.T = tin_f
    f.oxidizer_inlet.mdot = mdot_o
    f.oxidizer_inlet.X = comp_o
    f.oxidizer_inlet.T = tin_o
    f.boundary_emissivities = 0.0, 0.0
    f.radiation_enabled = False
    f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)
    f.solve(loglevel=1, auto=True)
    last_f = f

    x = f.grid
    T = f.T
    rho = f.density
    cp = f.cp_mass
    lam = f.thermal_conductivity

    # Mixture fraction from H element (Z = Y_H since Y_H_fuel=1, Y_H_ox=0)
    Y_H = np.zeros_like(x)
    M_H = 1.00794
    for i, name in enumerate(gas.species_names):
        nH = gas.n_atoms(name, 'H')
        if nH > 0:
            Y_H += f.Y[i] * (nH * M_H) / gas.molecular_weights[i]
    Z = Y_H

    # Stoichiometric point from peak temperature
    i_st = np.argmax(T)
    Z_st = Z[i_st]

    # Direct chi = 2 * alpha * |dZ/dx|^2
    alpha = lam / (rho * cp)
    dZ_dx = np.gradient(Z, x)
    chi_direct = 2.0 * alpha * dZ_dx**2
    chi_st_direct = chi_direct[i_st]

    # Strain rates from Cantera
    a_ox = f.strain_rate('potential_flow_oxidizer')
    a_fuel = f.strain_rate('potential_flow_fuel')
    a_max = f.strain_rate('max')
    a_stoich = f.strain_rate('stoichiometric', fuel='H2', oxidizer='O2', stoich=0.5)

    # erf model: chi(Z) = (a/pi) * exp(-2 * erfinv(2Z-1)^2)
    # at Z_st: chi_st = (a/pi) * exp(-2 * erfinv(2*Z_st-1)^2)
    exp_factor = np.exp(-2.0 * erfinv(2.0 * Z_st - 1.0)**2)

    chi_st_erf_ox = (a_ox / np.pi) * exp_factor
    chi_st_erf_fuel = (a_fuel / np.pi) * exp_factor

    results.append({
        'mdot_f': mdot_f,
        'chi_st_direct': chi_st_direct,
        'chi_st_erf_ox': chi_st_erf_ox,
        'chi_st_erf_fuel': chi_st_erf_fuel,
        'a_ox': a_ox,
        'a_fuel': a_fuel,
        'a_max': a_max,
        'a_stoich': a_stoich,
        'Z_st': Z_st,
        'T_peak': T[i_st],
    })

# --- Print comparison table ---
print(f"{'mdot_f':>6s} | {'Z_st':>6s} | {'T_peak':>7s} | {'a_ox':>8s} | {'a_fuel':>8s} | "
      f"{'chi_st_dir':>10s} | {'chi_st_erf_ox':>12s} | {'chi_st_erf_fuel':>13s} | "
      f"{'ratio_ox':>8s} | {'ratio_fuel':>9s}")
print("-" * 110)
for r in results:
    ratio_ox = r['chi_st_direct'] / r['chi_st_erf_ox']
    ratio_fuel = r['chi_st_direct'] / r['chi_st_erf_fuel']
    print(f"{r['mdot_f']:6.1f} | {r['Z_st']:6.4f} | {r['T_peak']:7.0f} | "
          f"{r['a_ox']:8.1f} | {r['a_fuel']:8.1f} | "
          f"{r['chi_st_direct']:10.2f} | {r['chi_st_erf_ox']:12.2f} | "
          f"{r['chi_st_erf_fuel']:13.2f} | {ratio_ox:8.2f} | {ratio_fuel:9.2f}")

# --- Detailed profile for last flame ---
f = last_f
x = f.grid
T = f.T
rho = f.density
cp = f.cp_mass
lam = f.thermal_conductivity

Y_H = np.zeros_like(x)
M_H = 1.00794
for i, name in enumerate(gas.species_names):
    nH = gas.n_atoms(name, 'H')
    if nH > 0:
        Y_H += f.Y[i] * (nH * M_H) / gas.molecular_weights[i]
Z = Y_H

i_st = np.argmax(T)
Z_st = Z[i_st]
alpha = lam / (rho * cp)
dZ_dx = np.gradient(Z, x)
chi_direct = 2.0 * alpha * dZ_dx**2
chi_st_direct = chi_direct[i_st]

a_ox = f.strain_rate('potential_flow_oxidizer')
a_fuel = f.strain_rate('potential_flow_fuel')
a_max = f.strain_rate('max')

# Full erf profiles
mask = (Z > 1e-6) & (Z < 1.0 - 1e-6)
chi_erf_ox = np.full_like(Z, np.nan)
chi_erf_fuel = np.full_like(Z, np.nan)
chi_erf_ox[mask] = (a_ox / np.pi) * np.exp(-2.0 * erfinv(2.0 * Z[mask] - 1.0)**2)
chi_erf_fuel[mask] = (a_fuel / np.pi) * np.exp(-2.0 * erfinv(2.0 * Z[mask] - 1.0)**2)

# --- Plots ---
fig, axes = plt.subplots(2, 2, figsize=(14, 9))

# 1. Temperature and mixture fraction
ax = axes[0, 0]
ax_twin = ax.twinx()
ax.plot(x * 1000, T, 'k-', label='T')
ax_twin.plot(x * 1000, Z, 'b--', label='Z')
ax.axvline(x[i_st] * 1000, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('x [mm]')
ax.set_ylabel('T [K]', color='k')
ax_twin.set_ylabel('Z [-]', color='b')
ax.set_title(f'Temperature and Z (mdot_f = {f.fuel_inlet.mdot}, a_ox = {a_ox:.0f} 1/s)')

# 2. chi vs x
ax = axes[0, 1]
ax.semilogy(x * 1000, chi_direct, 'r-', label=r'Direct: $2\alpha|\nabla Z|^2$')
ax.semilogy(x * 1000, chi_erf_ox, 'b--', label=rf'erf (a_ox = {a_ox:.0f})')
ax.semilogy(x * 1000, chi_erf_fuel, 'g--', label=rf'erf (a_fuel = {a_fuel:.0f})')
ax.axvline(x[i_st] * 1000, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('x [mm]')
ax.set_ylabel(r'$\chi$ [1/s]')
ax.set_title(r'$\chi$ in physical space')
ax.legend()

# 3. chi vs Z
ax = axes[1, 0]
ax.semilogy(Z, chi_direct, 'r-', label='Direct', linewidth=2)
ax.semilogy(Z, chi_erf_ox, 'b--', label=rf'erf (a_ox = {a_ox:.0f})')
ax.semilogy(Z, chi_erf_fuel, 'g--', label=rf'erf (a_fuel = {a_fuel:.0f})')
ax.axvline(Z_st, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel('Z [-]')
ax.set_ylabel(r'$\chi$ [1/s]')
ax.set_title(r'$\chi$ vs Z')
ax.legend()
ax.set_xlim(0, 0.15)

# 4. Extinction curve: chi_st vs strain rate
ax = axes[1, 1]
chi_dirs = [r['chi_st_direct'] for r in results]
a_oxs = [r['a_ox'] for r in results]
a_fuels = [r['a_fuel'] for r in results]
chi_erf_oxs = [r['chi_st_erf_ox'] for r in results]
chi_erf_fuels = [r['chi_st_erf_fuel'] for r in results]

ax.plot(a_oxs, chi_dirs, 'ro-', label=r'$\chi_{st}$ direct')
ax.plot(a_oxs, chi_erf_oxs, 'bs--', label=r'$\chi_{st}$ erf (a_ox)')
ax.plot(a_fuels, chi_erf_fuels, 'g^--', label=r'$\chi_{st}$ erf (a_fuel)')
ax.set_xlabel('Strain rate a [1/s]')
ax.set_ylabel(r'$\chi_{st}$ [1/s]')
ax.set_title(r'$\chi_{st}$ vs strain rate')
ax.legend()

plt.tight_layout()
plt.savefig('sdr_comparison.png', dpi=150)
print(f"\nDetailed plot saved to sdr_comparison.png")
