import pyscf
import numpy as np
import matplotlib.pyplot as plt


def calculate_mo_energies():
    """Calculates and plots MO energies for nitrogen in singlet form."""

    # Define the nitrogen molecule
    nitrogen_mol = pyscf.gto.M(
        atom="N 0 0 0; N 0 0 1.1",  # Example geometry
        basis="6-31G",
        spin=0,  # Singlet state
    )
    nitrogen_mol.build()

    # Run the Restricted Hartree-Fock calculation
    mf = nitrogen_mol.RHF().run()
    mf.MP2().run()

    # Get MO energies
    mo_energies = mf.mo_energy

    # Get occupation numbers
    mo_occ = mf.mo_occ

    # Determine occupied orbital indices
    occupied_indices = np.where(mo_occ > 0)[0]

    return mo_energies[2:], occupied_indices[2:]


def plot_mo_energies(mo_energies, occupied_indices):
    """Plots MO energies for nitrogen."""

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot MO energies
    plt.plot(mo_energies, marker="o", linestyle="--", label="Nitrogen", color="blue")
    plt.plot(
        occupied_indices - 2,
        mo_energies[occupied_indices - 2],
        marker="*",
        color="blue",
        linestyle="None",
        label="Occupied",
        markersize=12,
    )

    # Add labels and title
    plt.xlabel("Molecular Orbital Index")
    plt.ylabel("Energy (Hartree)")
    plt.title("Molecular Orbital Energies of Nitrogen (Singlet)")
    plt.legend()

    # Customize the plot
    plt.grid(True)
    plt.xticks(np.arange(0, len(mo_energies), 1))  # Show all MO indices

    # Show the plot
    plt.savefig("nitrogen_mo_energies.png")
    plt.show()


if __name__ == "__main__":
    # Calculate MO energies for nitrogen
    mo_energies, occupied_indices = calculate_mo_energies()

    # Plot the MO energies
    plot_mo_energies(mo_energies, occupied_indices)
