import pyscf
import numpy as np
import matplotlib.pyplot as plt

def calculate_mo_energies():
    """Calculates and plots MO energies for CO2 in singlet form."""

    # Define the CO2 molecule
    co2_mol = pyscf.gto.M(
        atom="C 0 0 0; O 0 0 1.16; O 0 0 -1.16",  # Example geometry
        basis="6-31G",
        spin=0,  # Singlet state
    )
    co2_mol.build()

    # Run the Restricted Hartree-Fock calculation
    mf = co2_mol.RHF().run()
    mf.MP2().run()

    # Get MO energies
    mo_energies = mf.mo_energy

    # Get occupation numbers
    mo_occ = mf.mo_occ

    # Determine occupied orbital indices
    occupied_indices = np.where(mo_occ > 0)[0]

    return mo_energies[3:], occupied_indices[3:]


def plot_mo_energies(mo_energies, occupied_indices):
    """Plots MO energies for CO2."""

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot MO energies
    plt.plot(mo_energies, marker='o', linestyle='--', label='CO2', color='blue')
    plt.plot(occupied_indices-3, mo_energies[occupied_indices-3], marker='*', color='blue', linestyle='None', label='Occupied', markersize=12)

    # Add labels and title
    plt.xlabel('Molecular Orbital Index')
    plt.ylabel('Energy (Hartree)')
    plt.title('Molecular Orbital Energies of CO2 (Singlet)')
    plt.legend()

    # Customize the plot
    plt.grid(True)
    plt.xticks(np.arange(0, len(mo_energies), 1))  # Show all MO indices

    # Show the plot
    plt.savefig("co2_mo_energies.png")
    plt.show()


if __name__ == "__main__":
    # Calculate MO energies for CO2
    mo_energies, occupied_indices = calculate_mo_energies()

    # Plot the MO energies
    plot_mo_energies(mo_energies, occupied_indices)