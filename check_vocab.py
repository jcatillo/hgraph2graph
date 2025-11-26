import sys
from rdkit import Chem

def check_vocab(filename):
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) < 2:
                print(f"Line {i+1}: Invalid format (less than 2 parts): {line.strip()}")
                continue
            
            smiles1, smiles2 = parts[0], parts[1]
            
            mol1 = Chem.MolFromSmiles(smiles1)
            if mol1 is None:
                print(f"Line {i+1}: Invalid SMILES 1: {smiles1}")
            
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol2 is None:
                print(f"Line {i+1}: Invalid SMILES 2: {smiles2}")

if __name__ == "__main__":
    check_vocab('data/antifungal/vocab.txt')
