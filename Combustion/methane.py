import pyscf
import numpy as np
import matplotlib.pyplot as plt


def calculate_mo_energies():
    """Calculates and plots MO energies for methane in singlet form."""

    # Define the methane molecule
    methane_mol = pyscf.gto.M(
        atom="C 0 0 0; H 0.634 0.634 0.634; H -0.634 -0.634 0.634; H 0.634 -0.634 -0.634; H -0.634 0.634 -0.634",  # Example geometry
        basis="6-31G",
        spin=0,  # Singlet state
    )
    methane_mol.build()

    # Run the Restricted Hartree-Fock calculation
    mf = methane_mol.RHF().run()
    mf.MP2().run()

    # Get MO energies
    mo_energies = mf.mo_energy

    # Get occupation numbers
    mo_occ = mf.mo_occ

    # Determine occupied orbital indices
    occupied_indices = np.where(mo_occ > 0)[0]

    # hcore
    # hcore = mf.get_hcore()
    # print(f"{hcore = }")

    return mo_energies[1:], occupied_indices[1:]


def plot_mo_energies(mo_energies, occupied_indices):
    """Plots MO energies for methane."""

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot MO energies
    plt.plot(mo_energies, marker="o", linestyle="--", label="Methane", color="red")
    plt.plot(
        occupied_indices - 1,
        mo_energies[occupied_indices - 1],
        marker="*",
        color="red",
        linestyle="None",
        label="Occupied",
        markersize=12,
    )

    # Add labels and title
    plt.xlabel("Molecular Orbital Index - 1")
    plt.ylabel("Energy (Hartree)")
    plt.title("Molecular Orbital Energies of Methane")
    plt.legend()

    # Customize the plot
    plt.grid(True)
    plt.xticks(np.arange(0, len(mo_energies), 1))  # Show all MO indices

    # Show the plot
    plt.savefig("methane_mo_energies.png")
    plt.show()


if __name__ == "__main__":
    # Calculate MO energies for methane
    mo_energies, occupied_indices = calculate_mo_energies()

    # Plot the MO energies
    plot_mo_energies(mo_energies, occupied_indices)
