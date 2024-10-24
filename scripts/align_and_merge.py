from rdkit import Chem
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem.rdFMCS import FindMCS


def align_ligands(ligA_file, ligB_file, output_file, B_to_A, unique_atoms_B):
    """
    Aligns ligand B to ligand A, applying greater weight to atoms near the attachment points, 
    and saves the aligned ligand B as a PDB file.

    Parameters:
    - ligA_file: Path to the MOL2 file for ligand A.
    - ligB_file: Path to the MOL2 file for ligand B.
    - output_file: Path to save the aligned ligand B as a PDB file.
    - B_to_A: Dictionary mapping atoms in ligand B to ligand A.
    - unique_atoms_B: Set of unique atoms in ligand B.
    """
    
    # Load ligA and ligB
    ligA = Chem.MolFromMol2File(ligA_file, removeHs=False)
    ligB = Chem.MolFromMol2File(ligB_file, removeHs=False)

    # Identify attachment points: atoms in ligB directly connected to unique atoms
    attachment_atoms_B = {
        neighbor.GetIdx()
        for idx in unique_atoms_B
        for neighbor in ligB.GetAtomWithIdx(idx).GetNeighbors()
        if neighbor.GetIdx() not in unique_atoms_B
    }

    # Extend important atoms by including neighbors of the attachment points
    important_atoms_B = attachment_atoms_B | {
        neighbor.GetIdx()
        for idx in attachment_atoms_B
        for neighbor in ligB.GetAtomWithIdx(idx).GetNeighbors()
        if neighbor.GetIdx() not in unique_atoms_B
    }

    # Extend important atoms by including neighbors to a depth of 2
    important_atoms_B = important_atoms_B | {
        neighbor.GetIdx()
        for idx in important_atoms_B
        for neighbor in ligB.GetAtomWithIdx(idx).GetNeighbors()
        if neighbor.GetIdx() not in unique_atoms_B
    }

    # Map B_to_A with adjusted weights for important atoms
    B_to_A_map = list(B_to_A.items())
    weights = [
        10000.0 if pair[0] in important_atoms_B or pair[1] in important_atoms_B else 1.0
        for pair in B_to_A.items()
    ]

    # Align ligB to ligA using atom map and weights
    AlignMol(ligB, ligA, atomMap=B_to_A_map, maxIters=100, weights=weights)
    
    # Save the aligned ligand B as a PDB file
    Chem.MolToPDBFile(ligB, output_file)


def merge_ligands(ligA, ligB, atoms_merged, unique_atoms_B, output_file="output/data/merged.pdb"):
    """
    Saves the merged ligand (ligA + unique atoms from ligB) coordinates to a PDB file.

    Parameters:
    - ligA: RDKit molecule object for ligand A.
    - ligB: RDKit molecule object for ligand B.
    - atoms_merged: List of merged atom data.
    - unique_atoms_B: Set of unique atoms from ligand B to be included in the merged structure.
    - output_file: Path to save the merged PDB file (default: "output/data/merged.pdb").
    """
    
    # Merge coordinates: ligand A + unique atoms from ligand B
    coords_merged = (
        ligA.GetConformer().GetPositions().tolist() +
        ligB.GetConformer().GetPositions()[unique_atoms_B].tolist()
    )

    # Open the output PDB file and write the merged coordinates
    with open(output_file, "w") as pdb_file:
        pdb_file.write("TITLE    Merged Ligand\n")
        pdb_file.write("MODEL    1\n")

        # Write each atom and its coordinates to the PDB file
        for i, (atom_row, coords) in enumerate(zip(atoms_merged, coords_merged), start=1):
            atom_name = atom_row[4]
            x, y, z = coords
            pdb_file.write(
                f"ATOM  {i:5d} {atom_name:<4} LIG     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
            )
        
        pdb_file.write("ENDMDL\n")


def read_gro_file(filename):
    """Reads a GROMACS .gro file and extracts the title, atom count, atom information, and box vector."""
    with open(filename, "r") as file:
        lines = file.readlines()
        atom_count = int(lines[1].strip())
        atoms = lines[2:-1]
        box_vector = lines[-1].strip()
    return atom_count, atoms, box_vector


def merge_gro_files(protein_gro, ligand_gro, output_gro):
    """Merges the protein and ligand .gro files into a single complex and saves it as a new .gro file."""

    # Read the protein and ligand .gro files
    protein_count, protein_atoms, protein_box = read_gro_file(protein_gro)
    ligand_count, ligand_atoms, _ = read_gro_file(ligand_gro)

    # Calculate the total number of atoms
    total_atoms = protein_count + ligand_count

    # Write the merged .gro file
    with open(output_gro, "w") as output_file:
        output_file.write("Protein-Ligand Complex\n")
        output_file.write(f"{total_atoms:5d}\n")
        output_file.writelines(protein_atoms)
        output_file.writelines(ligand_atoms)
        output_file.write(protein_box + "\n")

