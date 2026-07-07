#!/usr/bin/env python3
"""
Print the boundary conditions and residual equations for a
CounterflowDiffusionFlame, referencing the C++ source code.

Domain layout:

  [Inlet1D (fuel)]  <->  [AxisymmetricFlow]  <->  [Inlet1D (oxidizer)]

Residual evaluation is two-pass:
  1. Flow domain (AxisymmetricFlow) sets DEFAULT residuals at j=0 and j=N-1
  2. Boundary domains (Inlet1D) MODIFY them by subtracting the desired value
"""

import cantera as ct
import numpy as np

MECH = 'h2o2.yaml'
FUEL_MDOT = 0.5
OX_MDOT = 3.0
FUEL_T = 300.0
OX_T = 300.0
P = 1e5

gas = ct.Solution(MECH)
f = ct.CounterflowDiffusionFlame(gas, grid=np.linspace(0, 0.02, 600))
f.P = P
f.fuel_inlet.mdot = FUEL_MDOT
f.fuel_inlet.X = 'H2:1'
f.fuel_inlet.T = FUEL_T
f.oxidizer_inlet.mdot = OX_MDOT
f.oxidizer_inlet.X = 'O2:1'
f.oxidizer_inlet.T = OX_T
f.transport_model = 'unity-Lewis-number'
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1, prune=0.05)

# Source code references
SRC = {
    'Inlet1D::eval': 'src/oneD/Boundary1D.cpp:200-286',
    'Flow1D::evalContinuity': 'src/oneD/Flow1D.cpp:512-567',
    'Flow1D::evalMomentum': 'src/oneD/Flow1D.cpp:569-603',
    'Flow1D::evalLambda': 'src/oneD/Flow1D.cpp:605-641',
    'Flow1D::evalEnergy': 'src/oneD/Flow1D.cpp:643-702',
    'Flow1D::evalSpecies': 'src/oneD/Flow1D.cpp:704-780',
}

# Component offsets (from include/cantera/oneD/Flow1D.h)
# c_offset_U = 0  (axial velocity / mass flow)
# c_offset_V = 1  (radial spread rate)
# c_offset_T = 2  (temperature)
# c_offset_L = 3  (radial pressure gradient, Lambda)
# c_offset_Y = 4  (first species)

c_U, c_V, c_T, c_L = 0, 1, 2, 3
c_Y = 4

print("=" * 72)
print("  Boundary Conditions for CounterflowDiffusionFlame")
print("=" * 72)

print(f"\n{'─'*72}")
print("  DOMAIN STACK")
print(f"{'─'*72}")
print(f"  Index  Type         Name")
print(f"  {0:5d}  {f.fuel_inlet.domain_type:<12s} {f.fuel_inlet.name}")
print(f"  {1:5d}  {f.flame.domain_type:<12s} {f.flame.name}")
print(f"  {2:5d}  {f.oxidizer_inlet.domain_type:<12s} {f.oxidizer_inlet.name}")

print(f"\n{'─'*72}")
print("  FUEL INLET (left)  —  Inlet1D::eval()  [src/oneD/Boundary1D.cpp:200-286]")
print(f"{'─'*72}")
print(f"  Set temperature:     {f.fuel_inlet.T:.1f} K")
print(f"  Set composition:     {f.fuel_inlet.X}")
print(f"  Set mass flux:       {f.fuel_inlet.mdot:.4f} kg/m²/s")
print(f"  Set spread rate V₀:  {f.fuel_inlet.spread_rate:.2f} 1/s")
print()

print("  Variable   Flow1D default residual        Inlet modification            Combined residual (r=0)")
print("  ─────────  ──────────────────────────────  ───────────────────────────  ─────────────────────────────")
for i, label in enumerate(['V', 'T']):
    if label == 'V':
        default = f"rsd[c_V] = V(0)"
        mod = f"rsd[c_V] -= V₀  = 0"
        combined = f"V(0) = V₀  = 0"
    else:
        default = f"rsd[c_T] = T(0)"
        mod = f"rsd[c_T] -= T_in  = {f.fuel_inlet.T:.1f}"
        combined = f"T(0) = T_in  = {f.fuel_inlet.T:.1f} K"
    print(f"  {label:<8s}  {default:<40s}  {mod:<35s}  {combined}")

print(f"  Yₖ        rsd[c_Y+k] = -Jₖ(0) − ρu·Yₖ(0)   rsd[c_Y+k] += ṁ·Yₖ,in         convective flux: ṁ·Yₖ = ṁ·Yₖ,in")
print(f"           (diff flux + conv flux left)                                     (excess species → sum Yₖ = 1)")
print(f"  Λ         rsd[c_L] = −ρu(0)              rsd[c_L] += ṁ                  ρu = ṁ  = {f.fuel_inlet.mdot:.1f} kg/m²/s")
print()
print(f"  → Dirichlet: V=0, T={f.fuel_inlet.T:.1f}K")
print(f"  → Convective flux: ṁ·Yₖ,in for each species (excess species = sum-to-unity)")
print(f"  → Mass flow: ρu = ṁ = {f.fuel_inlet.mdot:.1f} kg/m²/s  (through Λ equation)")

print(f"\n{'─'*72}")
print("  OXIDIZER INLET (right)  —  Inlet1D::eval()  [src/oneD/Boundary1D.cpp:257-285]")
print(f"{'─'*72}")
print(f"  Set temperature:     {f.oxidizer_inlet.T:.1f} K")
print(f"  Set composition:     {f.oxidizer_inlet.X}")
print(f"  Set mass flux:       {f.oxidizer_inlet.mdot:.4f} kg/m²/s")
print(f"  Set spread rate V₀:  {f.oxidizer_inlet.spread_rate:.2f} 1/s")
print()

print("  Variable   Flow1D default residual        Inlet modification            Combined residual (r=0)")
print("  ─────────  ──────────────────────────────  ───────────────────────────  ─────────────────────────────")
print(f"  V          rsd[c_V] = V(N-1)              rsd[c_V] -= V₀  = 0           V(N-1) = V₀  = 0")
print(f"  T          rsd[c_T] = T(N-1)              rsd[c_T] -= T_in = {f.oxidizer_inlet.T:.1f}   T(N-1) = T_in = {f.oxidizer_inlet.T:.1f} K")
print(f"  Yₖ         rsd = Jₖ(N-2) + ρu·Yₖ(N-1)     rsd[c_Y+k] += ṁ·Yₖ,in         convective flux: ṁ·Yₖ = ṁ·Yₖ,in")
print(f"  u/cont.    rsd[c_U] = ρu(N-1) (zero flux) rsd[c_U] += ṁ                  ρu = ṁ  = {f.oxidizer_inlet.mdot:.1f} kg/m²/s")
print()
print(f"  → All residuals are **subtracted** from flow domain defaults")
print(f"  → Excess species is the one with largest mass fraction at each boundary")

print(f"\n{'─'*72}")
print("  FLOW DOMAIN DEFAULT BOUNDARY RESIDUALS  (before Inlet1D modification)")
print("  These are set by AxisymmetricFlow::eval()  [src/oneD/Flow1D.cpp]")
print(f"{'─'*72}")
print()
print("  Continuity equation  —  Flow1D::evalContinuity()  [src/oneD/Flow1D.cpp:512-567]")
print("  ───────────────────────────────────────────────────────────────────────────────")
print(f"   Left  (j=0):")
print(f"     rsd[c_U] = -(ρu₁ − ρu₀)/Δz − (ρ₁V₁ + ρ₀V₀)")
print(f"     diag = 0   (algebraic, no time derivative)")
print(f"   Right (j=N-1):")
print(f"     rsd[c_U] = ρu(N-1)   (zero mass flux through right boundary)")
print(f"     diag = 0")
print(f"   Interior:")
print("     rsd[c_U] = -(ρu_{j+1} - ρu_j)/dz - (ρ_{j+1}V_{j+1} + ρ_jV_j)")
print()
print("  Momentum (V / spread rate)  —  Flow1D::evalMomentum()  [src/oneD/Flow1D.cpp:569-603]")
print("  ───────────────────────────────────────────────────────────────────────────────")
print(f"   Both boundaries:")
print(f"     rsd[c_V] = V(0)        left")
print(f"     rsd[c_V] = V(N-1)      right")
print(f"     → Inlet1D modifies: rsd[c_V] -= V₀")
print(f"   Interior:")
print(f"     rsd[c_V] = (shear − Λ − ρu·dV/dz − ρV²) / ρ")
print()
print("  Lambda (radial pressure gradient)  —  Flow1D::evalLambda()  [src/oneD/Flow1D.cpp:605-641]")
print("  ───────────────────────────────────────────────────────────────────────────────")
print(f"   Left  (j=0):")
print(f"     rsd[c_L] = −ρu(0)   → Inlet1D: rsd[c_L] += ṁ  ⇒  ρu = ṁ")
print(f"   Right (j=N-1):")
print(f"     rsd[c_L] = Λ(N-1) − Λ(N-2)   (zero gradient)")
print(f"   Interior:")
print("     rsd[c_L] = Λ_j - Λ_{j-1}   (Λ is spatially constant)")
print()
print("  Energy (T)  —  Flow1D::evalEnergy()  [src/oneD/Flow1D.cpp:643-702]")
print("  ───────────────────────────────────────────────────────────────────────────────")
print(f"   Both boundaries:")
print(f"     rsd[c_T] = T(0)       left")
print(f"     rsd[c_T] = T(N-1)     right")
print(f"     → Inlet1D modifies: rsd[c_T] -= T_inlet")
print(f"   Interior:")
print(f"     rsd[c_T] = ρu·dT/dz − (λ/cp)·d²T/dz² − Σ(hₖ·ω̇ₖ·Wₖ)/cp − radiation")
print()
print("  Species (Yₖ)  —  Flow1D::evalSpecies()  [src/oneD/Flow1D.cpp:704-780]")
print("  ───────────────────────────────────────────────────────────────────────────────")
print(f"   Left  (j=0):")
print(f"     rsd[c_Y+k] = −(fluxₖ₀)₋ − ρu·Yₖ   (diffusive + convective flux)")
print(f"     → Inlet1D: rsd[c_Y+k] += ṁ·Yₖ,in")
print(f"     Excess species: rsd = 1 − ΣYₖ   (sum-to-unity)")
print(f"   Right (j=N-1):")
print(f"     rsd[c_Y+k] = fluxₖ(N-2) + ρu·Yₖ")
print(f"     → Inlet1D: rsd[c_Y+k] += ṁ·Yₖ,in")
print(f"     Excess species skipped")
print(f"   Interior:")
print(f"     rsd[c_Y+k] = ρu·dYₖ/dz − dJₖ/dz − ω̇ₖ·Wₖ")
print()

print(f"{'─'*72}")
print("  DEFINITION OF RESIDUAL ZERO")
print(f"{'─'*72}")
print("  The solver converges when ‖r‖ < tolerance for all grid points.")
print("  At the boundaries, this means the boundary conditions are enforced")
print("  exactly (to within the solver tolerance).")
print()
print("  Key design: The flow domain sets default residuals that would give")
print("  'neutral' BCs (V=0, T=0, etc.). The Inlet1D then modifies these")
print("  by subtracting the desired value:")
print()
print("    r[T] = T(0)          (flow domain default, r=0 means T(0)=0)")
print("    r[T] -= T_inlet      (inlet modifies: r=0 means T(0)=T_inlet)")
print()
print("  When Newton converges to r=0, the BC is satisfied.")
