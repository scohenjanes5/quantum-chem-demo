import pyscf
import numpy as np
import os
from glob import glob
from pyscf.tools import cubegen
from pyscf import gto, scf
from rich.progress import track


def calculate_mo_energies(num_nh3):
    """Calculates MO energies for Co(NH3)_n^(3+) complexes."""

    mols = glob(f"CoNH3{num_nh3}*opt.xyz")
    mo_energies_series = []
    names = []
    for geometry in mols:
        mol = gto.M(
            atom=geometry,
            basis="6-31g",
            charge=3,
        )
        mol.build()
        mf = scf.RHF(mol).x2c()
        if "tet" in geometry:
            mf.chkfile = f"CoNH3{num_nh3}-tet.chk"
            suffix = "tet"
        else:
            mf.chkfile = f"CoNH3{num_nh3}-square.chk"
            suffix = "square"

        names.append(suffix)
        mf.init_guess = "chkfile"
        print(f"Beginning calculation for {geometry}...")
        mf.kernel()
        pyscf.mp.MP2(mf).kernel()

        mo_energies_series.append(mf.mo_energy)
        homo_index = np.where(mf.mo_occ > 0)[-1][-1]
        print(homo_index)

        # Generate cube files for visualization
        for orb_idx in track(range(20, 30), description="Generating .cube files"):
            name = f"CoNH3{num_nh3}_mo_{orb_idx}_{suffix}.cube"
            if not os.path.exists(name):
                cubegen.orbital(
                    mol,
                    f"CoNH3{num_nh3}_mo_{orb_idx}_{suffix}.cube",
                    mf.mo_coeff[:, orb_idx - 1],
                )

    return mo_energies_series, names


if __name__ == "__main__":
    # Calculate MO energies for each complex
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.title(r"MO Energies for Co(NH$_3$)$_4$$^{3+}$, Tetrahedral vs. Square Planar")
    plt.xlabel("Molecular Orbital Index (0-indexed)")
    plt.ylabel("Energy (Hartree)")
    plt.grid(True)

    num_nh3 = 4
    mo_energies_series, names = calculate_mo_energies(num_nh3)
    shift = 20

    real_names = {"tet": "tetrahedral", "square": "square planar"}
    for mo_energies, name in zip(mo_energies_series, names):
        # Plotting the MO energies
        plt.plot(
            mo_energies[shift:], marker="o", linestyle="--", label=f"{real_names[name]}"
        )
        print(f"Energy of MO 30 for {real_names[name]}: {mo_energies[30]:.6f}")

    xticks = np.arange(len(mo_energies[shift:]))
    plt.xticks(xticks, [f"{i}" for i in xticks + shift])

    # # Add a vertical line at the HOMO
    # homo_index = 31
    # plt.axvline(
    #     x=homo_index - shift - 1,
    #     color="gray",
    #     linestyle="--",
    #     label="HOMO",
    # )  # Shift back to plot index

    plt.legend()
    plt.show()
