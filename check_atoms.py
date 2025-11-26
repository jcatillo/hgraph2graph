import sys
from rdkit import Chem

def get_atom_types(filename):
    atom_types = set()
    with open(filename, 'r') as f:
        for line in f:
            smiles = line.strip().split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                for atom in mol.GetAtoms():
                    atom_types.add((atom.GetSymbol(), atom.GetFormalCharge()))
    return atom_types

if __name__ == "__main__":
    atoms = get_atom_types('data/antifungal/all_single.txt')
    print("Found atoms:")
    for symbol, charge in sorted(atoms):
        print(f"('{symbol}', {charge})")
