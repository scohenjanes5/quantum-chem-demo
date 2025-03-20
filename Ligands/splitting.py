import pyscf
import numpy as np
import os
from pyscf.tools import cubegen
from pyscf import gto, scf, grad, mcscf


def calculate_mo_energies(num_nh3):
    """Calculates MO energies for Co(NH3)_n^(3+) complexes with geometry optimization."""

    # Define the complex geometry (simplified, adjust as needed)
    geometry = f"CoNH3{num_nh3}-opt.xyz"
    if not os.path.exists(geometry):
        geometry = f"CoNH3{num_nh3}.xyz"
        if not os.path.exists(geometry):
            raise FileNotFoundError(f"File {geometry} not found.")

    # Create the PySCF molecule object
    try:
        mol = gto.M(
            atom=geometry,
            basis="sto-3g",  # Minimal basis for speed
            charge=3,  # +3 charge for the complex
        )
        mol.build()
    except:
        mol = gto.M(
            atom=geometry,
            basis="sto-3g",
            charge=3,
            spin=1,  # Open-shell system
        )
        mol.build()

    # Geometry optimization
    mf = scf.RHF(mol)
    mf.chkfile = f"CoNH3{num_nh3}.chk"
    mf.init_guess = "chkfile"
    mf.kernel()
    pyscf.mp.MP2(mf).kernel()

    if "opt" not in geometry:
        # Perform geometry optimization
        from pyscf.geomopt import berny_solver

        mol_eq = berny_solver.optimize(mf, maxcycles=50)  # Increase maxcycles if needed
        mol_eq.tofile(f"CoNH3{num_nh3}-opt.xyz")
        mf = scf.RHF(mol_eq)
        mf.kernel()
        pyscf.mp.MP2(mf).kernel()
    else:
        mol_eq = mol
    # Get MO energies from the optimized geometry
    mo_energies = mf.mo_energy

    # Generate cube files for orbitals 8-10
    for orb_idx in range(30, 41):
        name = f"CoNH3{num_nh3}_mo_{orb_idx}_optimized.cube"
        if not os.path.exists(name):
            cubegen.orbital(
                mol_eq,
                f"CoNH3{num_nh3}_mo_{orb_idx}_optimized.cube",
                mf.mo_coeff[:, orb_idx - 1],
            )

    return mo_energies


if __name__ == "__main__":
    # Calculate MO energies for each complex
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.title(f"MO Energies for Co(NH3)_i^(3+), i = 4, 5, 6 (Geometry Optimized)")
    plt.xlabel("(Shifted) Molecular Orbital Index")
    plt.ylabel("Energy (Hartree)")
    plt.grid(True)

    for i, num_nh3 in enumerate([4, 5, 6]):
        mo_energies = calculate_mo_energies(num_nh3)

        # Plotting the MO energies
        plt.plot(
            mo_energies[20:], marker="o", linestyle="--", label=f"Co(NH3)_{num_nh3}^3+"
        )

    plt.legend()
    plt.show()
