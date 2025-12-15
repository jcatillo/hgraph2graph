import pandas as pd
from chembl_webresource_client.new_client import new_client

molecule = new_client.molecule


# 1. Fetch data (Your provided filters)
print("Fetching data from ChEMBL...")
res_small = molecule.filter(molecule_properties__full_mwt__range=[250, 350], molecule_properties__qed_weighted__gt=0.7).only(['molecule_chembl_id', 'molecule_structures'])[:700]

# 2. Combine results into a single list
# We convert the QuerySets to lists to force execution and merge them
all_results = list(res_small)

print(f"Fetched {len(all_results)} small molecules.")


# 3. Extract IDs and SMILES
data = []
for res in all_results:
    # safely get structure dict and smiles
    structures = res.get('molecule_structures')
    if structures:
        smiles = structures.get('canonical_smiles')
        chembl_id = res.get('molecule_chembl_id')
        if smiles and chembl_id:
            data.append({'ChEMBL_ID': chembl_id, 'SMILES': smiles})

# 4. Create DataFrame
df = pd.DataFrame(data)

print(f"Collected {len(df)} molecules.")
print(df.head())

df.to_csv('../data/small_molecule.csv', index=False)