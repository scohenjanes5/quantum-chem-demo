import pyscf
import numpy as np
import matplotlib.pyplot as plt


def calculate_mo_energies():
    """Calculates and plots MO energies for argon."""

    # Define the argon atom
    argon_mol = pyscf.gto.M(
        atom="Ar 0 0 0",  # Argon at the origin
        basis="6-31G",
        spin=0,  # Singlet state (closed-shell)
    )
    argon_mol.build()

    # Run the Restricted Hartree-Fock calculation
    mf = argon_mol.RHF().run()
    mf.MP2().run()

    # Get MO energies
    mo_energies = mf.mo_energy

    # Get occupation numbers
    mo_occ = mf.mo_occ

    # Determine occupied orbital indices
    occupied_indices = np.where(mo_occ > 0)[0]

    return mo_energies[1:], occupied_indices[1:]


def plot_mo_energies(mo_energies, occupied_indices):
    """Plots MO energies for argon."""

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot MO energies
    plt.plot(mo_energies, marker="o", linestyle="--", label="Argon", color="blue")
    plt.plot(
        occupied_indices - 1,
        mo_energies[occupied_indices - 1],
        marker="*",
        color="blue",
        linestyle="None",
        label="Occupied",
        markersize=12,
    )

    # Add labels and title
    plt.xlabel("# Molecular Orbitals past 1s")
    plt.ylabel("Energy (Hartree)")
    plt.title("Molecular Orbital Energies of Argon")
    plt.legend()

    # Customize the plot
    plt.grid(True)
    plt.xticks(np.arange(0, len(mo_energies), 1))  # Show all MO indices

    # Show the plot
    plt.savefig("argon_mo_energies.png")
    plt.show()


if __name__ == "__main__":
    # Calculate MO energies for argon
    mo_energies, occupied_indices = calculate_mo_energies()

    # Plot the MO energies
    plot_mo_energies(mo_energies, occupied_indices)
