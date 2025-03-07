import pyscf
import numpy as np

BASIS = '6-31G'
def get_cation_charge(metal_sym):
    # Determine charge of metal
    if metal_sym in ['Li', 'Na', 'K']:
        metal_charge = 1
    elif metal_sym in ['Be', 'Mg', 'Ca']:
        metal_charge = 2
    else:
        raise ValueError(f"Metal {metal_sym} not supported.")
    return metal_charge

def OH_energy():
    # Calculate energy of hydroxide ion
    oh_mol = pyscf.gto.M(atom='O 0 0 0; H 0 0 0.9', basis=BASIS, charge=-1)
    oh_mol.build()
    oh_mf = oh_mol.HF().run()
    e_hydroxide_ion = oh_mf.e_tot
    return e_hydroxide_ion

def dissociation_energy(mol_string):
    """
    Calculates the dissociation energy of a metal hydroxide molecule into metal ion and hydroxide ion(s).
    """
    try:
        mol = pyscf.gto.M(atom=mol_string, basis=BASIS)
    except:
        raise ValueError(f"Basis set not found for molecule {mol_string}.")
    mol.build()
    mf = mol.HF().run()
    e_mol = mf.e_tot

    # Determine metal and number of hydroxide ions
    atoms = mol_string.split(';')
    metal_sym = atoms[0].split()[0]  # e.g., Li, Na, Mg, Ba
    num_hydroxides = len([a for a in atoms if 'O' in a])  # Count 'O' atoms to determine number of OH groups

    metal_charge = get_cation_charge(metal_sym)
    
    # Calculate energy of metal ion
    metal_atom = pyscf.gto.M(atom=metal_sym, basis=BASIS, charge=metal_charge)
    metal_atom.build()
    metal_mf = metal_atom.HF().run()
    e_metal_ion = metal_mf.e_tot

    # Calculate dissociation energy
    e_products = e_metal_ion + num_hydroxides * HYDROXIDE_ENERGY
    diss_energy = e_products - e_mol
    return diss_energy

def ionization_energy(atom_string):
    """
    Calculates the ionization energy of an atom.
    """
    try:
        atom = pyscf.gto.M(atom=atom_string, basis=BASIS)
    except:
        atom = pyscf.gto.M(atom=atom_string, basis=BASIS, spin=1)
    atom.build()
    mf = atom.HF().run()
    e_atom = mf.e_tot

    metal_sym = atom_string.split()[0]
    metal_charge = get_cation_charge(metal_sym)
    # Ionized atom
    try:
        atom_ion = pyscf.gto.M(atom=atom_string, basis=BASIS, charge=metal_charge)
    except:
        atom_ion = pyscf.gto.M(atom=atom_string, basis=BASIS, charge=metal_charge, spin=1)
    atom_ion.build()
    mf_ion = atom_ion.HF().run()
    e_atom_ion = mf_ion.e_tot

    ionization_energy = e_atom_ion - e_atom
    return ionization_energy

#ionization energies as orbital energies

# Dissociation Energies
HYDROXIDE_ENERGY = OH_energy()
diss_energy_LiOH = dissociation_energy('Li 0 0 0; O 0 0 1.6; H 0 0 2.5')
diss_energy_NaOH = dissociation_energy('Na 0 0 0; O 0 0 1.6; H 0 0 2.5')
diss_energy_KOH = dissociation_energy('K 0 0 0; O 0 0 1.6; H 0 0 2.5')
diss_energy_BeOH2 = dissociation_energy('Be 0 0 0; O 0 0 1.6; H 0 0 2.5; O 0 0 -1.6; H 0 0 -2.5')
diss_energy_MgOH2 = dissociation_energy('Mg 0 0 0; O 0 0 1.6; H 0 0 2.5; O 0 0 -1.6; H 0 0 -2.5')
diss_energy_CaOH2 = dissociation_energy('Ca 0 0 0; O 0 0 1.6; H 0 0 2.5; O 0 0 -1.6; H 0 0 -2.5')


print("Dissociation Energies (Hartree):")
print(f"  LiOH: {diss_energy_LiOH:.6f}")
print(f"  NaOH: {diss_energy_NaOH:.6f}")
print(f"  KOH: {diss_energy_KOH:.6f}")
print(f"  BeOH2: {diss_energy_BeOH2:.6f}, {diss_energy_BeOH2/2:.6f} per OH group")
print(f"  MgOH2: {diss_energy_MgOH2:.6f}, {diss_energy_MgOH2/2:.6f} per OH group")
print(f"  CaOH2: {diss_energy_CaOH2:.6f}, {diss_energy_CaOH2/2:.6f} per OH group")

# Ionization Energies
ion_energy_Li = ionization_energy('Li 0 0 0')
ion_energy_Na = ionization_energy('Na 0 0 0')
ion_energy_K = ionization_energy('K 0 0 0')
ion_energy_Be = ionization_energy('Be 0 0 0')
ion_energy_Mg = ionization_energy('Mg 0 0 0')
ion_energy_Ca = ionization_energy('Ca 0 0 0')

print("\nIonization Energies (Hartree):")
print(f"  Li: {ion_energy_Li:.6f}")
print(f"  Na: {ion_energy_Na:.6f}")
print(f"  K: {ion_energy_K:.6f}")
print(f"  Be: {ion_energy_Be:.6f}")
print(f"  Mg: {ion_energy_Mg:.6f}")
print(f"  Ca: {ion_energy_Ca:.6f}")

"""
converged SCF energy = -82.8927462231428
converged SCF energy = -7.23548002444817
converged SCF energy = -75.3055242495255
converged SCF energy = -237.214111279513
converged SCF energy = -161.659276566549
converged SCF energy = -75.3055242495255
converged SCF energy = -350.492291881676
converged SCF energy = -198.811709459527
converged SCF energy = -75.3055242495254
converged SCF energy = -827.243545697224
converged SCF energy = -676.10389345016
converged SCF energy = -75.3055242495255
Dissociation Energies (Hartree):
  LiOH: 0.351742
  NaOH: 0.249310
  MgOH2: 1.069534
  CaOH2: 0.528604
converged SCF energy = -7.43123581108367  <S^2> = 0.75000068  2S+1 = 2.0000007
converged SCF energy = -7.23548002444817
converged SCF energy = -161.841425092229  <S^2> = 0.75004965  2S+1 = 2.0000497
converged SCF energy = -161.659276566549
converged SCF energy = -199.595219247262
converged SCF energy = -198.811709459527
converged SCF energy = -676.707922922738
converged SCF energy = -676.10389345016

Ionization Energies (Hartree):
  Li: 0.195756
  Na: 0.182149
  Mg: 0.783510
  Ca: 0.604029

Results vary with basis set and method used, but here we see that the dissociation energies are mich higher 
for the divalent cations than for the monovalent cations. The ionization energies are also higher for the 
divalent cations than for the monovalent cations. This contributes to the strength of the mono-hydroxide bases
relative to the di-hydroxide bases. We also see that the larger divalent cation has a lower dissociation energy
than the smaller divalent cation, which is consistent with the increased sheilding of the nucleus by the larger additonal
electron shell.
"""