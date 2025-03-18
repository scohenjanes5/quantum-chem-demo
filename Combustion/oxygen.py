import pyscf
import numpy as np
import matplotlib.pyplot as plt


def calculate_mo_energies(spin):
    """Calculates and plots MO energies for O2 in a given spin state."""

    # Define the O2 molecule
    o2_mol = pyscf.gto.M(
        atom="O 0 0 0; O 0 0 1.2",  # Example geometry
        basis="6-31G",
        spin=spin,  # Set the spin state (0 for singlet, 2 for triplet)
    )
    o2_mol.build()

    # Run the Unrestricted Hartree-Fock calculation
    mf = o2_mol.UHF().run()
    # mf = o2_mol.UKS().set(xc="b3lyp").run()
    mf.MP2().run()

    print(f"Total energy: {mf.e_tot}")

    # Get MO energies (alpha and beta)
    mo_energies_alpha = mf.mo_energy[0]
    mo_energies_beta = mf.mo_energy[1]

    # Get occupation numbers (alpha and beta)
    mo_occ_alpha = mf.mo_occ[0]
    mo_occ_beta = mf.mo_occ[1]

    # print(f"{mf.mo_occ = }")

    # Determine occupied orbital indices (alpha and beta)
    occupied_indices_alpha = np.where(mo_occ_alpha > 0)[0]
    occupied_indices_beta = np.where(mo_occ_beta > 0)[0]

    occupied_indices_alpha = occupied_indices_alpha[2:]
    occupied_indices_beta = occupied_indices_beta[2:]
    mo_energies_alpha = mo_energies_alpha[2:]
    mo_energies_beta = mo_energies_beta[2:]

    print(f"Occupied indices (alpha): {occupied_indices_alpha}")
    print(f"Occupied indices (beta): {occupied_indices_beta}")

    return (
        mo_energies_alpha,
        occupied_indices_alpha,
        mo_energies_beta,
        occupied_indices_beta,
    )


def plot_mo_energies(
    mo_energies_singlet_alpha,
    occupied_indices_singlet_alpha,
    mo_energies_singlet_beta,
    occupied_indices_singlet_beta,
    mo_energies_triplet_alpha,
    occupied_indices_triplet_alpha,
    mo_energies_triplet_beta,
    occupied_indices_triplet_beta,
):
    """Plots MO energies for singlet and triplet O2."""

    # Plot singlet MO energies (alpha)
    plt.plot(
        mo_energies_singlet_alpha,
        marker="o",
        linestyle="--",
        label="Singlet O2",
        color="red",
    )
    plt.plot(
        occupied_indices_singlet_alpha - 2,
        mo_energies_singlet_alpha[: len(occupied_indices_singlet_alpha)],
        marker="*",
        color="red",
        markersize=12,
        linestyle="None",
        label="Occupied (Singlet alpha)",
    )

    # Plot triplet MO energies (alpha)
    plt.plot(
        mo_energies_triplet_alpha,
        marker="x",
        linestyle="--",
        label="Triplet O2 (alpha)",
        color="purple",
    )
    plt.plot(
        occupied_indices_triplet_alpha - 2,
        mo_energies_triplet_alpha[: len(occupied_indices_triplet_alpha)],
        marker="*",
        color="purple",
        markersize=12,
        linestyle="None",
        label="Occupied (Triplet alpha)",
    )

    # Plot triplet MO energies (beta)
    plt.plot(
        mo_energies_triplet_beta,
        marker="x",
        linestyle="--",
        label="Triplet O2 (beta)",
        color="blue",
    )
    plt.plot(
        occupied_indices_triplet_beta - 2,
        mo_energies_triplet_beta[: len(occupied_indices_triplet_beta)],
        marker="*",
        color="blue",
        markersize=12,
        linestyle="None",
        label="Occupied (Triplet beta)",
    )

    plt.title("O2 Electronic Structure")
    plt.xlabel("Molecular Orbital Index - 2")
    plt.ylabel("Energy (Hartree)")
    plt.legend()
    plt.grid(True)

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig("O2_MO_energies_separate.png")
    plt.show()


if __name__ == "__main__":
    # Calculate MO energies for singlet O2 (spin=0)
    (
        mo_energies_singlet_alpha,
        occupied_indices_singlet_alpha,
        mo_energies_singlet_beta,
        occupied_indices_singlet_beta,
    ) = calculate_mo_energies(spin=0)

    # Calculate MO energies for triplet O2 (spin=2)
    (
        mo_energies_triplet_alpha,
        occupied_indices_triplet_alpha,
        mo_energies_triplet_beta,
        occupied_indices_triplet_beta,
    ) = calculate_mo_energies(spin=2)

    # Plot the MO energies
    plot_mo_energies(
        mo_energies_singlet_alpha,
        occupied_indices_singlet_alpha,
        mo_energies_singlet_beta,
        occupied_indices_singlet_beta,
        mo_energies_triplet_alpha,
        occupied_indices_triplet_alpha,
        mo_energies_triplet_beta,
        occupied_indices_triplet_beta,
    )
