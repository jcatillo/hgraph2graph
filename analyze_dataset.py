import sys
from rdkit import Chem

# Read all molecules
with open('data/antifungal/all.txt', 'r') as f:
    molecules = [line.strip() for line in f if line.strip()]

print(f"Total molecules: {len(molecules)}")

# Analyze the dataset
multi_component = []
single_component = []
invalid_smiles = []

for smiles in molecules:
    if '.' in smiles:
        multi_component.append(smiles)
    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_smiles.append(smiles)
        else:
            single_component.append(smiles)

print(f"\nBreakdown:")
print(f"  Multi-component (salts/mixtures): {len(multi_component)}")
print(f"  Single-component valid: {len(single_component)}")
print(f"  Invalid SMILES: {len(invalid_smiles)}")

print(f"\nPercentage of single-component: {len(single_component)/len(molecules)*100:.1f}%")

# Save single-component molecules
with open('data/antifungal/all_single.txt', 'w') as f:
    for smiles in single_component:
        f.write(smiles + '\n')

print(f"\nSaved {len(single_component)} single-component molecules to data/antifungal/all_single.txt")
