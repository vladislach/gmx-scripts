import re


def extract_sections(itp_path):
    """Extracts relevant sections (atomtypes, atoms, bonds, etc.) from a given .itp file."""
    sections = {
        "atomtypes": [],
        "atoms": [],
        "bonds": [],
        "pairs": [],
        "angles": [],
        "dihedrals": []
    }

    section_pattern = re.compile(r'\[\s*(\w+)\s*\]')  # Pattern to match section headers (e.g., [ atoms ])
    current_section = None  # Keeps track of the current section being parsed

    with open(itp_path, 'r') as file:
        for line in file:
            # Check if the line matches a section header
            match = section_pattern.match(line)
            if match:
                section_name = match.group(1).lower()
                # Set the current section if it matches one of the sections of interest
                current_section = section_name if section_name in sections else None
            elif current_section and line.strip() and not line.startswith(";"):
                # Add non-comment, non-empty lines to the current section
                sections[current_section].append(line.split(';')[0].strip())

    return sections


def process_atom(nr, typeA, resnr, res, atom, cgnr, chargeA, massA, typeB='', chargeB='', massB=''):
    """Helper function to process atoms and add them to the merged list."""
    resnr = '1'
    res = 'LIG'  # Use 'LIG' as residue name for merged ligand
    return [nr, typeA, resnr, res, atom, cgnr, chargeA, massA, typeB, chargeB, massB]


def merge_atomtypes(atomtypes_A, atomtypes_B):
    """
    Merges and deduplicates atomtypes from ligand A and ligand B (lists). Returns a dictionary
    containing merged atomtypes.
    """
    atomtypes_merged = {
        name: (sigma, epsilon)
        for atomtype in set(atomtypes_A + atomtypes_B)
        for name, _, _, _, _, sigma, epsilon in [atomtype.split()]
    }
    return atomtypes_merged


def merge_atoms(sections_A, sections_B, unique_atoms_A, unique_atoms_B, mcs_atoms_A, A_to_B):
    """
    Merges atoms from ligand A and ligand B, identifying common and unique atoms.
    
    Parameters:
    - sections_A: Dictionary containing atom data from ligand A.
    - sections_B: Dictionary containing atom data from ligand B.
    - unique_atoms_A: List of atom indices unique to ligand A.
    - unique_atoms_B: List of atom indices unique to ligand B.
    - mcs_atoms_A: List of atom indices for the Maximum Common Substructure (MCS) in ligand A.
    - A_to_B: Dictionary mapping MCS atom indices from ligand A to ligand B.
    
    Returns:
    - A sorted list of merged atoms.
    - A merged dictionary of atom types.
    """
    atomtypes_A = sections_A['atomtypes']
    atomtypes_B = sections_B['atomtypes']

    # Merge atomtypes from both ligands
    atomtypes_merged = merge_atomtypes(atomtypes_A, atomtypes_B)

    atoms_merged = []

    # Process atoms unique to ligand A
    for i in unique_atoms_A:
        atomA = sections_A['atoms'][i]
        nr, typeA, resnr, res, atom, cgnr, chargeA, massA = atomA.split()
        typeB = f'dum_{typeA}'  # Use dummy atom type for ligand B
        if typeB not in atomtypes_merged:
            atomtypes_merged[typeB] = ('0.000000000', '0.000000000')  # Add dummy atomtype if not present
        atoms_merged.append(process_atom(nr, typeA, resnr, res, atom, cgnr, chargeA, massA, typeB, '0.000000', massA))

    # Process atoms common to both ligands (common substructure)
    for i in mcs_atoms_A:
        atomA = sections_A['atoms'][i]
        atomB = sections_B['atoms'][A_to_B[i]]
        nr, typeA, resnr, res, atom, cgnr, chargeA, massA = atomA.split()
        _, typeB, _, _, _, _, chargeB, massB = atomB.split()
        
        # Handle atoms with different charges between ligands A and B
        if chargeA != chargeB:
            atoms_merged.append(process_atom(nr, typeA, resnr, res, atom, cgnr, chargeA, massA, typeB, chargeB, massB))
        else:
            atoms_merged.append(process_atom(nr, typeA, resnr, res, atom, cgnr, chargeA, massA))

    # Process atoms unique to ligand B
    for i in unique_atoms_B:
        atomB = sections_B['atoms'][i]
        nr, typeB, resnr, res, atom, cgnr, chargeB, massB = atomB.split()
        nr = str(len(atoms_merged) + 1)  # Reindex the atom for the merged topology
        cgnr = str(len(atoms_merged) + 1)  # Use the new index for the charge group
        atom = 'D' + atom  # Prefix the atom name with 'D' to indicate dummy status
        typeA = f'dum_{typeB}'  # Dummy atom type for ligand A
        if typeA not in atomtypes_merged:
            atomtypes_merged[typeA] = ('0.000000000', '0.000000000')  # Add dummy atomtype if not present
        atoms_merged.append(process_atom(nr, typeA, resnr, res, atom, cgnr, '0.000000', massB, typeB, chargeB, massB))

    # Sort the merged atoms by their index
    atoms_merged = sorted(atoms_merged, key=lambda x: int(x[0]))

    return atoms_merged, atomtypes_merged


def process_bonds(sections_A, sections_B, B_to_merged_reindexed):
    """Processes and merges bonds from ligand A and ligand B."""
    # Process bonds of ligand A
    bonds_A = {
        (int(i), int(j)): (b0, kb)  # (i, j) are atom indices, b0 is bond length, kb is force constant
        for line in sections_A['bonds']
        for i, j, _, b0, kb in [line.split()]
    }

    # Process bonds of ligand B with reindexed atom indices for the merged molecule
    bonds_B = {
        (B_to_merged_reindexed[int(i)], B_to_merged_reindexed[int(j)]): (b0, kb)
        for line in sections_B['bonds']
        for i, j, _, b0, kb in [line.split()]
    }

    # Combine and deduplicate bonds
    bonds_merged = bonds_A | {bond: bonds_B[bond] for bond in bonds_B if bond not in bonds_A}

    return bonds_merged


def process_pairs(sections_A, sections_B, B_to_merged_reindexed):
    """Processes and merges pairs from ligand A and ligand B."""
    # Process pairs of ligand A
    pairs_A = [
        (int(i), int(j))  # (i, j) are atom indices for the interaction pairs
        for line in sections_A['pairs']
        for i, j, _ in [line.split()]
    ]

    # Process pairs of ligand B with reindexed atom indices for the merged molecule
    pairs_B = [
        (B_to_merged_reindexed[int(i)], B_to_merged_reindexed[int(j)])
        for line in sections_B['pairs']
        for i, j, _ in [line.split()]
    ]

    # Combine and deduplicate pairs from both ligands
    pairs_merged = pairs_A + sorted(list(set(pairs_B) - set(pairs_A)))

    return pairs_merged


def process_angles(sections_A, sections_B, B_to_merged_reindexed):
    """Processes and merges angles from ligand A and ligand B."""
    # Process angles of ligand A
    angles_A = {
        (int(i), int(j), int(ak)): (th0, cth)  # (i, j, ak) are atom indices, th0 is the angle, cth is the force constant
        for line in sections_A['angles']
        for i, j, ak, _, th0, cth in [line.split()]
    }

    # Process angles of ligand B with reindexed atom indices for the merged molecule
    angles_B = {
        (B_to_merged_reindexed[int(i)], B_to_merged_reindexed[int(j)], B_to_merged_reindexed[int(ak)]): (th0, cth)
        for line in sections_B['angles']
        for i, j, ak, _, th0, cth in [line.split()]
    }

    # Combine and deduplicate angles from both ligands
    angles_merged = angles_A | {angle: angles_B[angle] for angle in angles_B if angle not in angles_A}

    return angles_merged


def process_dihedrals(sections_A, sections_B, B_to_merged_reindexed):
    """Processes and merges dihedrals from ligand A and ligand B."""
    # Process dihedrals of ligand A
    dihedrals_A = {
        (int(i), int(j), int(k), int(l)): (funct, phi0, kphi, n)  # (i, j, k, l) are atom indices
        for line in sections_A['dihedrals']                      # phi0 is the phase angle, kphi is the force constant, n is multiplicity
        for i, j, k, l, funct, phi0, kphi, n in [line.split()]
    }

    # Process dihedrals of ligand B with reindexed atom indices for the merged molecule
    dihedrals_B = {
        (B_to_merged_reindexed[int(i)], B_to_merged_reindexed[int(j)], B_to_merged_reindexed[int(k)],
         B_to_merged_reindexed[int(l)]): (funct, phi0, kphi, n)
        for line in sections_B['dihedrals']
        for i, j, k, l, funct, phi0, kphi, n in [line.split()]
    }

    # Combine and deduplicate dihedrals from both ligands
    dihedrals_merged = dihedrals_A | {dihedral: dihedrals_B[dihedral] for dihedral in dihedrals_B if dihedral not in dihedrals_A}

    return dihedrals_merged


def format_atomtype_entry(name, sigma, epsilon):
    """
    Formats a single atomtype entry for the force field file."""
    return f"  {name:<7}  0.00000  0.00000    A    {sigma:>12}  {epsilon:>12}\n"


def format_atom_entry(nr, type, resnr, residue, atom, cgnr, charge, mass, typeB='', chargeB='', massB=''):
    """Formats a single atom entry for the merged .itp file."""
    return (
        f"  {nr:<4}  "
        f"{type:<7} "
        f"{resnr:>3}      "
        f"{residue:>3}     "
        f"{atom:<5}  "
        f"{cgnr:<3}   "
        f"{charge:>9}   "
        f"{mass:>9}    "
        f"{typeB:<7} "
        f"{chargeB:>9}   "
        f"{massB:>9}"
    )


def format_bond_entry(i, j, b0, kb):
    """Formats a single bond entry for the merged .itp file."""
    return f"  {i:<5}  {j:<5}   1   {b0:>12}   {kb:>12}"


def format_pair_entry(i, j):
    """Formats a single pair entry for the merged .itp file."""
    return f"  {i:<5}  {j:<5}   1"


def format_angle_entry(ai, aj, ak, th0, cth):
    """Formats a single angle entry for the merged .itp file."""
    return f"  {ai:<5}  {aj:<5}  {ak:<5}   1   {th0:>12}   {cth:>12}"


def format_dihedral_entry(ai, aj, ak, al, funct, phi0, kphi, n):
    """Formats a single dihedral entry for the merged .itp file."""
    return f"  {ai:<5}  {aj:<5}  {ak:<5}  {al:<5}   {funct}  {phi0:>10}  {kphi:>10}    {n:>2}"


def write_merged_itp(sections_merged, output_path):
    """Writes the merged ligand topology to a .itp file."""
    with open(output_path, 'w') as file:
        # Loop through each section and write the data to the file
        for section_name in ['moleculetype', 'atoms', 'bonds', 'pairs', 'angles', 'dihedrals']:
            file.write(f"[ {section_name} ]\n")
            for line in sections_merged[section_name]:
                file.write(f"{line}\n")
            file.write("\n")


def write_ffmerged_itp(atomtypes_merged, output_path):
    """Writes the ligand force field parameters to a .itp file."""
    # Format the ligand force field parameters (atomtypes)
    atomtypes_formatted = [
        format_atomtype_entry(name, sigma, epsilon)
        for name, (sigma, epsilon) in atomtypes_merged.items()
    ]

    # Write the formatted atomtypes to the file
    with open(output_path, "w") as file:
        file.write('[ atomtypes ]\n')
        file.write('; name       mass   charge   ptype      sigma        epsilon\n')
        file.writelines(atomtypes_formatted)


def write_merged_top(ffmerged_itp, merged_itp, output_top):
    """
    Generates a topology file for the merged ligand by combining the force field parameters 
    and ligand topology with water topology and system details.
    """
    with open(ffmerged_itp, 'r') as f:
        ffmerged_content = f.read()
    with open(merged_itp, 'r') as f:
        merged_content = f.read()

    with open(output_top, 'w') as file:
        # Include the force field parameters
        file.write(
            '''; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"\n
'''
        )
        
        # Add ligand force field parameters and ligand topology
        file.write(ffmerged_content)
        file.write("\n")
        file.write(merged_content)

        # Include the water topology
        file.write(
            '''; Include water topology
#include "amber99sb-ildn.ff/tip3p.itp"\n
'''
        )

        # Define the system and molecules
        file.write(
            '''[ system ]
Merged Ligand in Water

[ molecules ]
; Compound      #mols
Ligand            1
'''
        )


def merge_topol_files(
    protein_top_file='output/data/protein.top',
    ffmerged_itp_file='output/data/ffmerged.itp',
    merged_itp_file='output/data/merged.itp',
    output_file="output/data/complex.top"
):
    """Merges the protein, ligand force field, and ligand topology files into a single complex topology file."""

    # Read the protein topology file
    with open(protein_top_file, "r") as file:
        lines = file.readlines()

    # Find and keep only the lines after the force field inclusion comment
    start_idx = next(
        i for i, line in enumerate(lines) if line.strip().startswith('; Include forcefield parameters')
    )
    clean_lines = lines[start_idx:]

    # Read the ligand force field and topology files
    with open(ffmerged_itp_file, "r") as f:
        ffmerged_content = f.read()
    with open(merged_itp_file, "r") as f:
        lig_content = f.read()

    # Insert ligand parameters before the [ moleculetype ] section
    for i, line in enumerate(clean_lines):
        if line.strip().startswith('[ moleculetype ]'):
            clean_lines.insert(i, f'; Include ligand parameters\n{ffmerged_content}\n')
            break

    # Insert ligand topology before the water topology
    for i, line in enumerate(clean_lines):
        if line.strip().startswith('; Include water topology'):
            clean_lines.insert(i, f'; Include ligand topology\n{lig_content}\n')
            break

    # Rename the system in the [ system ] section
    for i, line in enumerate(clean_lines):
        if line.strip().startswith('[ system ]'):
            clean_lines[i + 2] = "Protein-Ligand Complex\n"

    # Add Ligand to the [ molecules ] section
    for i, line in enumerate(clean_lines):
        if line.strip().startswith('[ molecules ]'):
            clean_lines.insert(i + 3, "Ligand              1\n")
            break

    # Write the merged content to the output complex topology file
    with open(output_file, "w") as file:
        file.writelines(clean_lines)

