import h5py
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


def load_data(file_path, name):
    with h5py.File(file_path, "r") as f:
        grp = f[f"{name}/flame"]
        T = grp["T"][()]
        Y = grp["Y"][()]
        grid = grp["grid"][()]
        source = grp["phase"].attrs["source"]
        mech = source.split("/")[-1]
    return grid, T, Y, mech


def compute_enthalpy_and_Z(gas, T, Y, fuel_idx=0, oxidizer_idx=3):
    n = len(T)
    h = np.empty(n)
    Z = np.empty(n)
    fuel = np.zeros(gas.n_species)
    oxidizer = np.zeros(gas.n_species)
    fuel[fuel_idx] = 1.0
    oxidizer[oxidizer_idx] = 1.0
    for j in range(n):
        gas.TPY = T[j], gas.P, Y[j]
        h[j] = gas.enthalpy_mass
        Z[j] = gas.mixture_fraction(fuel, oxidizer, basis="mass")
    return h, Z


def perfect_v_shape(Z, h, zero_ends=False):
    i_tip = np.argmin(h)
    if zero_ends:
        left_z = np.array([Z[0], Z[i_tip]])
        left_h = np.array([0.0, h[i_tip]])
        right_z = np.array([Z[i_tip], Z[-1]])
        right_h = np.array([h[i_tip], 0.0])
    else:
        left_z = Z[: i_tip + 1]
        left_h = h[: i_tip + 1]
        right_z = Z[i_tip:]
        right_h = h[i_tip:]
    left = np.polyfit(left_z, left_h, 1)
    right = np.polyfit(right_z, right_h, 1)
    h_v = np.where(np.arange(len(Z)) <= i_tip, np.polyval(left, Z),
                   np.polyval(right, Z))
    return h_v, i_tip


def temperature_from_HPY(gas, h, Y, P=None):
    if P is None:
        P = gas.P
    n = len(h)
    T = np.empty(n)
    for j in range(n):
        gas.HPY = h[j], P, Y[j]
        T[j] = gas.T
    return T


def process_name(file_path, name, gas):
    grid, T_orig, Y, mech = load_data(file_path, name)
    gas_i = ct.Solution(mech)
    P = gas_i.P
    h_orig, Z = compute_enthalpy_and_Z(gas_i, T_orig, Y)
    h_v, i_tip = perfect_v_shape(Z, h_orig, zero_ends=True)
    T_v = temperature_from_HPY(gas_i, h_v, Y)
    max_dh = np.max(np.abs(h_orig - h_v))
    max_dT = np.max(np.abs(T_orig - T_v))
    return {
        "name": name,
        "grid": grid,
        "Z": Z,
        "n": len(grid),
        "h_orig": h_orig,
        "h_v": h_v,
        "T_orig": T_orig,
        "T_v": T_v,
        "i_tip": i_tip,
        "max_dh": max_dh,
        "max_dT": max_dT,
        "h_min": h_orig.min(),
        "h_max": h_orig.max(),
        "Z_tip": Z[i_tip],
        "h_tip": h_orig[i_tip],
    }


if __name__ == "__main__":
    file_path = "Scripts/Data/enthalpy.h5"
    names = ["old", "new"]
    save_path = "Scripts/Data/complete.png"

    gas_ref = ct.Solution("h2o2.yaml")
    results = [process_name(file_path, n, gas_ref) for n in names]

    print(f"{'':>30s}  {'old':>12s}  {'new':>12s}")
    print("-" * 58)
    print(f"{'grid points':>30s}  {results[0]['n']:>12d}  {results[1]['n']:>12d}")
    print(f"{'h range [MJ/kg]':>30s}  {results[0]['h_min']/1e6:>8.3f}..{results[0]['h_max']/1e6:>5.2f}"
          f"  {results[1]['h_min']/1e6:>8.3f}..{results[1]['h_max']/1e6:>5.2f}")
    print(f"{'Z at tip':>30s}  {results[0]['Z_tip']:>12.4f}  {results[1]['Z_tip']:>12.4f}")
    print(f"{'h at tip [MJ/kg]':>30s}  {results[0]['h_tip']/1e6:>12.3f}  {results[1]['h_tip']/1e6:>12.3f}")
    print(f"{'tip index':>30s}  {results[0]['i_tip']:>12d}  {results[1]['i_tip']:>12d}")
    print(f"{'max |h - h_V| [MJ/kg]':>30s}  {results[0]['max_dh']/1e6:>12.3e}  {results[1]['max_dh']/1e6:>12.3e}")
    print(f"{'max |T - T_V| [K]':>30s}  {results[0]['max_dT']:>12.3f}  {results[1]['max_dT']:>12.3f}")

    print()
    for key, unit, div in [("h", "MJ/kg", 1e6), ("T", "K", 1)]:
        for res, lbl in zip(results, ["old", "new"]):
            rms_orig = np.sqrt(np.mean(res[f"{key}_orig"]**2))
            rms_v = np.sqrt(np.mean(res[f"{key}_v"]**2))
            red = (1 - rms_v / rms_orig) * 100
            print(f"RMS({lbl})  {key}: {rms_orig/div:.4e} {unit},"
                  f"  {key}_V: {rms_v/div:.4e} {unit},"
                  f"  reduction: {red:.1f}%")

    old, new = results
    T_v_old_on_new = np.interp(new["Z"][::-1], old["Z"][::-1],
                                old["T_v"][::-1])[::-1]

    styles = [dict(color="C0", label="no ref"), dict(color="C1", label="ref")]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for res, sty in zip(results, styles):
        c = sty["color"]
        lbl = sty["label"]

        ax = axes[0, 0]
        ax.plot(res["Z"], res["h_orig"] / 1e6, color=c, linestyle="-",
                label=f"{lbl} orig", linewidth=1)
        ax.plot(res["Z"], res["h_v"] / 1e6, color=c, linestyle="--",
                label=f"{lbl} calc", linewidth=2)
        ax.plot(res["Z_tip"], res["h_tip"] / 1e6, "o", color=c, markersize=6)

        ax = axes[1, 0]
        ax.plot(res["Z"], res["T_orig"], color=c, linestyle="-",
                label=f"{lbl} orig", linewidth=1)
        ax.plot(res["Z"], res["T_v"], color=c, linestyle="--",
                label=f"{lbl} calc", linewidth=2)

        ax = axes[1, 1]
        ax.plot(res["Z"], res["T_v"] - res["T_orig"], color=c,
                linestyle="-", label=lbl)
        ax.axhline(0, color="gray", linewidth=0.5)

    ax = axes[0, 1]
    ax.plot(new["Z"], new["T_v"] - T_v_old_on_new, "k-",
            label="ref calc - no ref calc", linewidth=1)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_xlabel("Z")
    ax.set_ylabel("ΔT in K")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_title("ΔT calculated ref - no ref")

    axes[0, 0].set_xlabel("Z")
    axes[0, 0].set_ylabel("h in MJ/kg")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True)
    axes[0, 0].set_title("Enthalpy")

    axes[1, 0].set_xlabel("Z")
    axes[1, 0].set_ylabel("T in K")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True)
    axes[1, 0].set_title("Temperature")

    axes[1, 1].set_xlabel("Z")
    axes[1, 1].set_ylabel("ΔT T in K")
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(True)
    axes[1, 1].set_title("ΔT  original vs. calculated")

    fig.suptitle(f"no ref ({results[0]['n']} pts) and ref ({results[1]['n']} pts)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"\nSaved {save_path}")
