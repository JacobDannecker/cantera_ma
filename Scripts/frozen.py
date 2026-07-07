#!/usr/bin/env python3
"""
frozen.py

Simulate a counterflow diffusion flame (H2/O2) under three chemistry regimes:
  1. Normal           — finite-rate chemistry (full mechanism)
  2. Frozen           — no reactions (pure mixing)
  3. Infinitely fast  — equilibrium (Burke-Schumann limit)

Boundary conditions:
  Fuel:       H2   at T=300 K, mdot=0.5  kg/m²/s
  Oxidizer:   O2   at T=300 K, mdot=3.0  kg/m²/s
  Pressure:   1 bar
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
MECH = 'h2o2.yaml'
WIDTH = 0.02          # 20 mm
N_GRID = 600
P = 1e5               # 1 bar
FUEL_MDOT = 0.5       # kg/m²/s
OX_MDOT = 3.0         # kg/m²/s
T_IN = 300.0          # K


def build_flame(gas, grid):
    f = ct.CounterflowDiffusionFlame(gas, grid=grid)
    f.P = P
    f.fuel_inlet.mdot = FUEL_MDOT
    f.fuel_inlet.X = 'H2:1'
    f.fuel_inlet.T = T_IN
    f.oxidizer_inlet.mdot = OX_MDOT
    f.oxidizer_inlet.X = 'O2:1'
    f.oxidizer_inlet.T = T_IN
    f.transport_model = 'unity-Lewis-number'
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1, prune=0.05)
    return f


def solve_flame(label, gas, grid):
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    f = build_flame(gas, grid)
    f.solve(loglevel=0, auto=True)
    print(f"  Peak temperature:  {max(f.T):.1f} K")
    print(f"  Grid points:       {f.flame.n_points}")
    return f


def equilibrium_limit(n_pts, gas, h_fuel, h_ox, p, z_min=1e-8, z_max=1-1e-8):
    """Burke-Schumann equilibrium limit on a refined mixture fraction grid."""
    iH2 = gas.species_index('H2')
    iO2 = gas.species_index('O2')
    z_eq = np.linspace(z_min, z_max, n_pts)
    T_eq = np.empty(n_pts)
    h_mix_arr = np.empty(n_pts)
    Y_eq = np.empty((gas.n_species, n_pts))
    Y = np.zeros(gas.n_species)
    for i in range(n_pts):
        Y[:] = 0.0
        Y[iH2] = z_eq[i]
        Y[iO2] = 1.0 - z_eq[i]
        h_mix = z_eq[i] * h_fuel + (1.0 - z_eq[i]) * h_ox
        h_mix_arr[i] = h_mix
        gas.HPY = h_mix, p, Y
        gas.equilibrate('HP')
        T_eq[i] = gas.T
        Y_eq[:, i] = gas.Y
    return z_eq, T_eq, Y_eq, h_mix_arr


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    grid = np.linspace(0, WIDTH, N_GRID)

    # ---- 1. Normal (finite-rate chemistry) ----
    gas = ct.Solution(MECH)
    f_norm = solve_flame("1. Normal (finite-rate chemistry)", gas, grid)

    # ---- 2. Frozen chemistry (no reactions) ----
    gas_frz = ct.Solution(MECH)
    for i in range(gas_frz.n_reactions):
        gas_frz.set_multiplier(0.0, i)
    f_frz = solve_flame("2. Frozen chemistry (no reactions)", gas_frz, grid)

    # ---- 3. Infinitely fast chemistry (equilibrium limit) ----
    print(f"\n{'='*60}")
    print("  3. Infinitely fast chemistry (equilibrium limit)")
    print(f"{'='*60}")
    gas_eq = ct.Solution(MECH)
    gas_eq.TP = T_IN, P
    gas_eq.X = 'H2:1'
    h_fuel = gas_eq.enthalpy_mass
    gas_eq.X = 'O2:1'
    h_ox = gas_eq.enthalpy_mass
    N_EQ = 2000
    z_eq, T_eq, Y_eq, h_eq = equilibrium_limit(N_EQ, gas_eq, h_fuel, h_ox, P)
    print(f"  Peak temperature:  {max(T_eq):.1f} K")
    print(f"  Mixture fraction grid: {N_EQ} points from {z_eq[0]:.1e} to {z_eq[-1]:.2f}")
    gas_st = ct.Solution(MECH)
    gas_st.set_equivalence_ratio(1.0, 'H2', 'O2')
    z_stoich = gas_st.mixture_fraction('H2', 'O2', basis='mass', element='H')
    print(f"  Stoichiometric mixture fraction: {z_stoich:.4f}")
    idx_st = np.argmin(np.abs(z_eq - z_stoich))
    print(f"  Equilibrium T at stoichiometric: {T_eq[idx_st]:.1f} K")
    z_frz = f_frz.mixture_fraction('H')

    # ---- Summary ----
    print(f"\n{'='*60}")
    print("  SUMMARY")
    print(f"{'='*60}")
    print(f"  1. Normal:              T_max = {max(f_norm.T):8.1f} K")
    print(f"  2. Frozen:              T_max = {max(f_frz.T):8.1f} K")
    print(f"  3. Equilibrium:         T_max = {max(T_eq):8.1f} K")
    print(f"{'='*60}")

    # ---- Save data ----
    iH2 = gas.species_index('H2')
    iO2 = gas.species_index('O2')
    iH2O = gas.species_index('H2O')

    np.savetxt('Scripts/Data/flame_normal.csv',
               np.column_stack([f_norm.grid, f_norm.T,
                                f_norm.mixture_fraction('H')]),
               header='x [m], T [K], z_H', delimiter=',')
    np.savetxt('Scripts/Data/flame_frozen.csv',
               np.column_stack([f_frz.grid, f_frz.T,
                                f_frz.mixture_fraction('H')]),
               header='x [m], T [K], z_H', delimiter=',')
    np.savetxt('Scripts/Data/flame_equilibrium.csv',
               np.column_stack([z_eq, T_eq, h_eq, Y_eq[iH2, :], Y_eq[iO2, :], Y_eq[iH2O, :]]),
               header='z_H [-], T_eq [K], h_eq [J/kg], Y_H2 [-], Y_O2 [-], Y_H2O [-]', delimiter=',')

    print("\n  Data saved to Scripts/Data/")

    # ---- Plot ----
    z_norm = f_norm.mixture_fraction('H')
    z_frz = f_frz.mixture_fraction('H')

    fig, axs = plt.subplots(3, 2, figsize=(10, 11))

    # Temperature
    axs[0, 0].plot(z_norm, f_norm.T, label='Normal')
    axs[0, 0].plot(z_frz, f_frz.T, label='Frozen')
    axs[0, 0].plot(z_eq, T_eq, '--', label='Equilibrium')
    axs[0, 0].set_xlabel('Mixture fraction $z_H$')
    axs[0, 0].set_ylabel('$T$ [K]')
    axs[0, 0].legend()
    axs[0, 0].grid(True)

    # H2
    axs[0, 1].plot(z_norm, f_norm.Y[iH2, :], label='Normal')
    axs[0, 1].plot(z_frz, f_frz.Y[iH2, :], label='Frozen')
    axs[0, 1].plot(z_eq, Y_eq[iH2, :], '--', label='Equilibrium')
    axs[0, 1].set_xlabel('Mixture fraction $z_H$')
    axs[0, 1].set_ylabel('$Y_{H_2}$')
    axs[0, 1].legend()
    axs[0, 1].grid(True)

    # O2
    axs[1, 0].plot(z_norm, f_norm.Y[iO2, :], label='Normal')
    axs[1, 0].plot(z_frz, f_frz.Y[iO2, :], label='Frozen')
    axs[1, 0].plot(z_eq, Y_eq[iO2, :], '--', label='Equilibrium')
    axs[1, 0].set_xlabel('Mixture fraction $z_H$')
    axs[1, 0].set_ylabel('$Y_{O_2}$')
    axs[1, 0].legend()
    axs[1, 0].grid(True)

    # H2O
    axs[1, 1].plot(z_norm, f_norm.Y[iH2O, :], label='Normal')
    axs[1, 1].plot(z_frz, f_frz.Y[iH2O, :], label='Frozen')
    axs[1, 1].plot(z_eq, Y_eq[iH2O, :], '--', label='Equilibrium')
    axs[1, 1].set_xlabel('Mixture fraction $z_H$')
    axs[1, 1].set_ylabel('$Y_{H_2O}$')
    axs[1, 1].legend()
    axs[1, 1].grid(True)

    # Enthalpy
    axs[2, 0].plot(z_norm, f_norm.enthalpy_mass, label='Normal')
    axs[2, 0].plot(z_frz, f_frz.enthalpy_mass, label='Frozen')
    axs[2, 0].plot(z_eq, h_eq, '--', label='Equilibrium')
    axs[2, 0].set_xlabel('Mixture fraction $z_H$')
    axs[2, 0].set_ylabel('$h$ [J/kg]')
    axs[2, 0].legend()
    axs[2, 0].grid(True)

    axs[2, 1].axis('off')

    plt.tight_layout()
    plt.savefig('Scripts/Data/flame_comparison.png', dpi=150)
    plt.close()
    print("  Plot saved to Scripts/Data/flame_comparison.png")


if __name__ == '__main__':
    main()
