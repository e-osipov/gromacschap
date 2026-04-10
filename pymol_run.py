"""
Run PyMOL in console mode: load test.pdb, process with RDKit, save session.
"""

import sys
import os

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# PyMOL headless/console mode
import pymol
from pymol import cmd


def load_and_inspect_with_rdkit(pdb_path: str):
    """Use RDKit to inspect the PDB before loading into PyMOL."""
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
    if mol is None:
        print(f"[RDKit] Warning: could not parse {pdb_path}")
        return

    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    mw = Descriptors.MolWt(mol)
    print(f"[RDKit] Atoms: {num_atoms}, Bonds: {num_bonds}, MW: {mw:.2f} Da")
    return mol


def run_pymol_session(pdb_path: str, session_out: str = "output.pse"):
    """Launch PyMOL in console mode, run commands, and save session."""

    # Initialize PyMOL without GUI
    pymol.finish_launching(["pymol", "-cq"])  # -c: no GUI, -q: quiet

    # Load the structure
    mol_name = os.path.splitext(os.path.basename(pdb_path))[0]
    cmd.load(pdb_path, mol_name)
    print(f"[PyMOL] Loaded '{pdb_path}' as object '{mol_name}'")

    # --- Basic visualization commands ---

    # Show cartoon for protein backbone
    cmd.show("cartoon", mol_name)

    # Show sticks for ligands / HETATM residues
    cmd.show("sticks", f"{mol_name} and hetatm")

    # Color by secondary structure
    cmd.color("cyan", f"{mol_name} and ss h")   # helices
    cmd.color("yellow", f"{mol_name} and ss s")  # sheets
    cmd.color("white", f"{mol_name} and ss l+")  # loops

    # Color heteroatoms by element
    cmd.util.cbag(f"{mol_name} and hetatm")

    # Remove solvent for a cleaner view
    cmd.remove("solvent")

    # Zoom to fit the whole structure
    cmd.zoom(mol_name)

    # Set a nicer background
    cmd.bg_color("black")

    # Orient along principal axes
    cmd.orient(mol_name)

    # Print some info
    n_atoms = cmd.count_atoms(mol_name)
    print(f"[PyMOL] Atoms in session (after solvent removal): {n_atoms}")

    # Save the session
    cmd.save(session_out)
    print(f"[PyMOL] Session saved to '{session_out}'")

    # Also export a PNG for quick inspection (optional)
    png_out = session_out.replace(".pse", ".png")
    cmd.png(png_out, width=1920, height=1080, dpi=150, ray=1)
    print(f"[PyMOL] Ray-traced image saved to '{png_out}'")

    cmd.quit()


def main():
    pdb_path = sys.argv[1] if len(sys.argv) > 1 else "test.pdb"
    session_out = sys.argv[2] if len(sys.argv) > 2 else "output.pse"

    if not os.path.isfile(pdb_path):
        print(f"Error: '{pdb_path}' not found.")
        sys.exit(1)

    print(f"=== RDKit inspection ===")
    load_and_inspect_with_rdkit(pdb_path)

    print(f"\n=== PyMOL session ===")
    run_pymol_session(pdb_path, session_out)


if __name__ == "__main__":
    main()
