import cantera as ct
gas = ct.Solution("h2o2.yaml")
fuel = "H2:1"
oxidizer = "O2:1"
gas.set_equivalence_ratio(1.0, fuel, oxidizer)
Z_st = gas.mixture_fraction(fuel, oxidizer, element="Bilger")
print(f"mix_fraction stoich: {Z_st}")
