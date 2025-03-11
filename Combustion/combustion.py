import pyscf

BASIS = "6-31G"  # Options: "sto-3g", "6-31G", "cc-pVDZ"
METHOD = "B3LYP"  # Options: "HF", "B3LYP", "CCSD"


def atom_energy(atom_string):
    """Calculates the energy of an atom."""
    try:
        atom = pyscf.gto.M(atom=atom_string, basis=BASIS)
    except:
        print("Using spin=1 for open shell atoms: ", atom_string.split()[0])
        atom = pyscf.gto.M(atom=atom_string, basis=BASIS, spin=1)
    atom.build()
    if METHOD == "B3LYP":
        mf = atom.KS().set(xc="b3lyp")
    elif METHOD == "CCSD":
        mf = atom.CCSD()
    else:
        mf = atom.HF()
    mf.chkfile = f"{atom_string.split()[0]}.chk"
    mf.init_guess = 'chkfile'
    mf.kernel()
    e_atom = mf.e_tot
    #convert to kJ/mol
    e_atom = e_atom * 2625.5
    return e_atom


def molecule_energy(mol_string):
    """Calculates the energy of a molecule."""
    mol = pyscf.gto.M(atom=mol_string, basis=BASIS)
    mol.build()
    if METHOD == "B3LYP":
        mf = mol.KS().set(xc="b3lyp")
    elif METHOD == "CCSD":
        mf = mol.CCSD()
    else:
        mf = mol.HF()

    mf.chkfile = f"{mol_string.split(".")[0]}.chk"
    mf.init_guess = 'chkfile'
    mf.kernel()
    e_mol = mf.e_tot
    #convert to kJ/mol
    e_mol = e_mol * 2625.5
    return e_mol


# Define molecules
CH4 = "methane.xyz"
O2 = "oxygen.xyz" #"O 0 0 0; O 0 0 1.21"
CO2 = "CO2.xyz"
H2O = "water.xyz"

# Calculate atom energies
e_C = atom_energy("C 0 0 0")
e_H = atom_energy("H 0 0 0")
e_O = atom_energy("O 0 0 0")

# Calculate molecule energies
e_CH4 = molecule_energy(CH4)
e_O2 = molecule_energy(O2)
e_CO2 = molecule_energy(CO2)
e_H2O = molecule_energy(H2O)

# Calculate heats of formation
hf_CH4 = e_CH4 - (e_C + 4 * e_H)
hf_O2 = e_O2 - (2 * e_O)
hf_CO2 = e_CO2 - (e_C + 2 * e_O)
hf_H2O = e_H2O - (e_O + 2 * e_H)

# Combustion reaction: CH4 + 2O2 -> CO2 + 2H2O
e_reactants = hf_CH4 + 2 * hf_O2 + (e_C + 4 * e_H) + 2 * (2 * e_O)
e_products = hf_CO2 + 2 * hf_H2O + (e_C + 2 * e_O) + 2 * (e_O + 2 * e_H)

# Calculate reaction energy
reaction_energy = e_products - e_reactants

print("Heats of Formation from atoms (kJ/mol):")
print(f"  CH4: {hf_CH4:.0f}")
print(f"  O2:  {hf_O2:.0f}")
print(f"  CO2: {hf_CO2:.0f}")
print(f"  H2O: {hf_H2O:.0f}")

print("Energy per bond/bond order (kJ/mol):")
print(f"  C-H: {hf_CH4/4:.0f}")
print(f"  O=O: {hf_O2/2:.0f}")
print(f"  C=O: {hf_CO2/4:.0f}")
print(f"  O-H: {hf_H2O/2:.0f}")

print("\nReaction Energy (kJ/mol):")
print(f"  CH4 + 2O2 -> CO2 + 2H2O: {reaction_energy:.0f}")

lit_hf_CH4 = -74.873
lit_hf_O2 = 0
lit_hf_CO2 = -393.522
lit_hf_H2O =  	-285.830
lit_reaction_energy = 2 * lit_hf_H2O + lit_hf_CO2 - lit_hf_CH4 - 2 * lit_hf_O2
print("\nHeats of formation from std states (Lit: https://janaf.nist.gov/) (kJ/mol):")
print(f"  CH4: {lit_hf_CH4:.0f}")
print(f"  O2:  {lit_hf_O2:.0f}")
print(f"  CO2: {lit_hf_CO2:.0f}")
print(f"  H2O: {lit_hf_H2O:.0f}")
print(f"  Reaction: {lit_reaction_energy:.0f}")