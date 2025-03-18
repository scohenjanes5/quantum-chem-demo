import pyscf

BASIS = "6-31G"
METHOD = "HF"  # Options: "HF", "B3LYP", "CCSD"


def get_cation_charge(metal_sym):
    # Determine charge of metal
    if metal_sym in ["Li", "Na", "K", "Rb", "Cs"]:
        metal_charge = 1
    elif metal_sym in ["Be", "Mg", "Ca", "Sr", "Ba"]:
        metal_charge = 2
    else:
        raise ValueError(f"Metal {metal_sym} not supported.")
    return metal_charge


def OH_energy():
    # Calculate energy of hydroxide ion
    oh_mol = pyscf.gto.M(atom="O 0 0 0; H 0 0 0.9", basis=BASIS, charge=-1)
    oh_mol.build()
    if METHOD == "B3LYP":
        oh_mf = oh_mol.KS().set(xc="b3lyp").run()
    elif METHOD == "CCSD":
        oh_mf = oh_mol.CCSD().run()
    else:
        oh_mf = oh_mol.HF().run()
    e_hydroxide_ion = oh_mf.e_tot
    return e_hydroxide_ion


def metal_energy(metal_sym, charge=0):
    """Calculates the energy of a metal atom or cation with a given charge."""
    try:
        metal_atom = pyscf.gto.M(atom=metal_sym + " 0 0 0", basis=BASIS, charge=charge)
    except:
        metal_atom = pyscf.gto.M(
            atom=metal_sym + " 0 0 0", basis=BASIS, spin=1, charge=charge
        )
    metal_atom.build()
    if METHOD == "B3LYP":
        metal_mf = metal_atom.KS().set(xc="b3lyp").run()
    elif METHOD == "CCSD":
        metal_mf = metal_atom.CCSD().run()
    else:
        metal_mf = metal_atom.HF().run()
    e_metal = metal_mf.e_tot
    return e_metal


def dissociation_energy(mol_string, metal_cation_energies, full=True):
    """
    Calculates the dissociation energy of a metal hydroxide molecule.
    If full=True, it calculates the energy of dissociation into metal ion and all hydroxide ions.
    If full=False and the metal is divalent, it calculates the energy of dissociation into metal+OH and remaining hydroxide ion.
    """
    try:
        mol = pyscf.gto.M(atom=mol_string, basis=BASIS)
    except:
        raise ValueError(f"Basis set not found for molecule {mol_string}.")
    mol.build()
    if METHOD == "B3LYP":
        mf = mol.KS().set(xc="b3lyp").run()
    elif METHOD == "CCSD":
        mf = mol.CCSD().run()
    else:
        mf = mol.HF().run()
    e_mol = mf.e_tot

    # Determine metal and number of hydroxide ions
    atoms = mol_string.split(";")
    metal_sym = atoms[0].split()[0]  # e.g., Li, Na, Mg, Ba
    num_hydroxides = len(
        [a for a in atoms if "O" in a]
    )  # Count 'O' atoms to determine number of OH groups

    metal_charge = get_cation_charge(metal_sym)

    if full or metal_charge != 2:
        # Calculate dissociation energy
        e_products = (
            metal_cation_energies[(metal_sym, metal_charge)]
            + num_hydroxides * HYDROXIDE_ENERGY
        )
    else:
        # Calculate energy of metal+OH ion
        metal_oh_string = metal_sym + " 0 0 0; O 0 0 1.6; H 0 0 2.5"
        try:
            metal_oh = pyscf.gto.M(atom=metal_oh_string, basis=BASIS, charge=1)
        except:
            metal_oh = pyscf.gto.M(atom=metal_oh_string, basis=BASIS, charge=1, spin=1)
        metal_oh.build()
        if METHOD == "B3LYP":
            metal_oh_mf = metal_oh.KS().set(xc="b3lyp").run()
        elif METHOD == "CCSD":
            metal_oh_mf = metal_oh.CCSD().run()
        else:
            metal_oh_mf = metal_oh.HF().run()
        e_metal_oh = metal_oh_mf.e_tot

        # Calculate dissociation energy
        e_products = e_metal_oh + (num_hydroxides - 1) * HYDROXIDE_ENERGY
    diss_energy = e_products - e_mol
    return diss_energy


def ionization_energy(
    metal_sym, metal_atom_energies, metal_cation_energies, first=False
):
    """
    Calculates the ionization energy of an atom.
    """
    e_atom = metal_atom_energies[metal_sym]

    if not first:
        metal_charge = get_cation_charge(metal_sym)
    else:
        metal_charge = 1

    e_atom_ion = metal_cation_energies[(metal_sym, metal_charge)]

    ionization_energy = e_atom_ion - e_atom
    return ionization_energy


# Define molecules
LIOH = "Li 0 0 0; O 0 0 1.6; H 0 0 2.5"
NAOH = "Na 0 0 0; O 0 0 1.6; H 0 0 2.5"
KOH = "K 0 0 0; O 0 0 1.6; H 0 0 2.5"
BEOH2 = "Be 0 0 0; O 0 0 1.6; H 0 0 2.5; O 0 0 -1.6; H 0 0 -2.5"
MGOH2 = "Mg 0 0 0; O 0 0 1.6; H 0 0 2.5; O 0 0 -1.6; H 0 0 -2.5"
CAOH2 = "Ca 0 0 0; O 0 0 1.6; H 0 0 2.5; O 0 0 -1.6; H 0 0 -2.5"

# Calculate hydroxide energy
HYDROXIDE_ENERGY = OH_energy()

# Calculate metal cation energies
metal_symbols = ["Li", "Na", "K", "Be", "Mg", "Ca"]
metal_cation_energies = {}
for metal_sym in metal_symbols:
    metal_charge = get_cation_charge(metal_sym)
    if metal_charge == 1:
        metal_cation_energies[(metal_sym, 1)] = metal_energy(metal_sym, 1)
    elif metal_charge == 2:
        metal_cation_energies[(metal_sym, 1)] = metal_energy(metal_sym, 1)
        metal_cation_energies[(metal_sym, 2)] = metal_energy(metal_sym, 2)

# Calculate metal atom energies
metal_atom_energies = {}
for metal_sym in metal_symbols:
    metal_atom_energies[metal_sym] = metal_energy(metal_sym)

# Calculate dissociation energies
diss_energy_LiOH = dissociation_energy(LIOH, metal_cation_energies)
diss_energy_NaOH = dissociation_energy(NAOH, metal_cation_energies)
diss_energy_KOH = dissociation_energy(KOH, metal_cation_energies)
# diss_energy_BeOH2 = dissociation_energy(
#     BEOH2, metal_cation_energies, metal_atom_energies
# )
# diss_energy_MgOH2 = dissociation_energy(
#     MGOH2, metal_cation_energies, metal_atom_energies
# )
# diss_energy_CaOH2 = dissociation_energy(
#     CAOH2, metal_cation_energies, metal_atom_energies
# )

print("Dissociation Energies (Hartree):")
print(f"  Li+, OH-: {diss_energy_LiOH:.6f}")
print(f"  Na+, OH-: {diss_energy_NaOH:.6f}")
print(f"  K+, OH-: {diss_energy_KOH:.6f}")
# print(f"  BeOH2: {diss_energy_BeOH2:.6f}, {diss_energy_BeOH2/2:.6f} per OH group")
# print(f"  MgOH2: {diss_energy_MgOH2:.6f}, {diss_energy_MgOH2/2:.6f} per OH group")
# print(f"  CaOH2: {diss_energy_CaOH2:.6f}, {diss_energy_CaOH2/2:.6f} per OH group")

diss_energy_BeOH = dissociation_energy(BEOH2, metal_cation_energies, full=False)
diss_energy_MgOH = dissociation_energy(MGOH2, metal_cation_energies, full=False)
diss_energy_CaOH = dissociation_energy(CAOH2, metal_cation_energies, full=False)

print(f"  BeOH+, OH-: {diss_energy_BeOH:.6f}")
print(f"  MgOH+, OH-: {diss_energy_MgOH:.6f}")
print(f"  CaOH+, OH-: {diss_energy_CaOH:.6f}")

# Calculate ionization energies
ion_energy_Li = ionization_energy("Li", metal_atom_energies, metal_cation_energies)
ion_energy_Na = ionization_energy("Na", metal_atom_energies, metal_cation_energies)
ion_energy_K = ionization_energy("K", metal_atom_energies, metal_cation_energies)
# ion_energy_Be = ionization_energy("Be", metal_atom_energies, metal_cation_energies)
# ion_energy_Mg = ionization_energy("Mg", metal_atom_energies, metal_cation_energies)
# ion_energy_Ca = ionization_energy("Ca", metal_atom_energies, metal_cation_energies)

print("\nIonization Energies (Hartree):")
print(f"  Li: {ion_energy_Li:.6f}")
print(f"  Na: {ion_energy_Na:.6f}")
print(f"  K: {ion_energy_K:.6f}")
# print(f"  Be: {ion_energy_Be:.6f}")
# print(f"  Mg: {ion_energy_Mg:.6f}")
# print(f"  Ca: {ion_energy_Ca:.6f}")

ion_energy_Be_2 = ionization_energy(
    "Be", metal_atom_energies, metal_cation_energies, first=True
)
ion_energy_Mg_2 = ionization_energy(
    "Mg", metal_atom_energies, metal_cation_energies, first=True
)
ion_energy_Ca_2 = ionization_energy(
    "Ca", metal_atom_energies, metal_cation_energies, first=True
)

print(f"  Be: {ion_energy_Be_2:.6f}")
print(f"  Mg: {ion_energy_Mg_2:.6f}")
print(f"  Ca: {ion_energy_Ca_2:.6f}")

import matplotlib.pyplot as plt
import numpy as np

# Data for plotting
group1_metals = ["Li", "Na", "K"]
group2_metals = ["Be", "Mg", "Ca"]

# Period of the element
group1_periods = [2, 3, 4]
group2_periods = [2, 3, 4]

# Conversion factor from Hartree to kJ/mol
HARTREE_TO_KJ_PER_MOL = 2625.5

group1_diss_energies = (
    np.array([diss_energy_LiOH, diss_energy_NaOH, diss_energy_KOH])
    * HARTREE_TO_KJ_PER_MOL
)
group2_diss_energies = (
    np.array([diss_energy_BeOH, diss_energy_MgOH, diss_energy_CaOH])
    * HARTREE_TO_KJ_PER_MOL
)

group1_ion_energies = (
    np.array([ion_energy_Li, ion_energy_Na, ion_energy_K]) * HARTREE_TO_KJ_PER_MOL
)
group2_ion_energies = (
    np.array([ion_energy_Be_2, ion_energy_Mg_2, ion_energy_Ca_2])
    * HARTREE_TO_KJ_PER_MOL
)


# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 10))

# Plot dissociation energies
axs[0].plot(group1_periods, group1_diss_energies, marker="o", label="Group 1")
axs[0].plot(group2_periods, group2_diss_energies, marker="o", label="Group 2")
axs[0].set_xlabel("Period")
axs[0].set_ylabel("Dissociation Energy (kJ/mol)")
axs[0].set_title("Dissociation Energy vs Period")
axs[0].set_xticks(group1_periods)
axs[0].legend()

# Plot ionization energies
axs[1].plot(group1_periods, group1_ion_energies, marker="o", label="Group 1")
axs[1].plot(group2_periods, group2_ion_energies, marker="o", label="Group 2")
axs[1].set_xlabel("Period")
axs[1].set_ylabel("Ionization Energy (kJ/mol)")
axs[1].set_title("Ionization Energy vs Period")
axs[1].set_xticks(group1_periods)
axs[1].legend()

# Adjust layout and display the plot
plt.tight_layout()
plt.show()
