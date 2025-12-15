import time
import pandas as pd
from chembl_webresource_client.new_client import new_client

molecule = new_client.molecule

print("Fetching data from ChEMBL...")
res_large = molecule.filter(
    molecule_properties__full_mwt__range=[451, 600],
    molecule_properties__qed_weighted__gt=0.7
).only(['molecule_chembl_id', 'molecule_structures'])[:2000]

for attempt in range(5):
    try:
        all_results = list(res_large)  
        time.sleep(0.2)              
        break
    except Exception as e:
        if attempt == 4:
            raise
        time.sleep(0.2 * (attempt + 1)) 

print(f"Fetched {len(all_results)} large molecules.")

data = []
for res in all_results:
    structures = res.get('molecule_structures')
    if structures:
        smiles = structures.get('canonical_smiles')
        chembl_id = res.get('molecule_chembl_id')

        if smiles and chembl_id:
            data.append({'ChEMBL_ID': chembl_id, 'SMILES': smiles})

df = pd.DataFrame(data)

print(f"Collected {len(df)} molecules.")
print(df.head())

df.to_csv('../data/large_molecules.csv', index=False)
