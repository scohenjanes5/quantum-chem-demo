import pyscf
import numpy as np
import os
from glob import glob
from pyscf.tools import cubegen
from pyscf import gto, scf, grad, mcscf
from rich.progress import track


def calculate_mo_energies():
    """Calculates MO energies for Co3+ ion."""

    # Define the Co3+ ion
    mol = gto.M(
        atom="Co 0 0 0",  # Cobalt at the origin
        basis="6-31g",
        charge=3,  # +3 charge for the ion
    )
    mol.build()

    mf = scf.RHF(mol).x2c()
    mf.kernel()
    pyscf.mp.MP2(mf).kernel()

    # Get MO energies
    mo_energies = mf.mo_energy

    return mo_energies


if __name__ == "__main__":
    # Calculate MO energies for the Co3+ ion
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.title(r"MO Energies for Co$^{3+}$ Ion")
    plt.xlabel("Molecular Orbital Index (0-indexed)")
    plt.ylabel("Energy (Hartree)")
    plt.grid(True)

    mo_energies = calculate_mo_energies()
    shift = 0  # No shift needed for a single atom

    # Plotting the MO energies
    plt.plot(mo_energies[shift:], marker="o", linestyle="--", label="Co$^{3+}$")

    xticks = np.arange(len(mo_energies[shift:]))
    plt.xticks(xticks, [f"{i}" for i in xticks + shift])

    plt.legend()
    plt.show()
