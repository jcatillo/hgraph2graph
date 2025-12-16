import torch
import argparse
import random
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm

from hgraph import *
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import MolToFile

lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)

parser = argparse.ArgumentParser(description="Evaluate reconstruction performance")
parser.add_argument("--vocab", required=True)
parser.add_argument("--atom_vocab", default=common_atom_vocab)
parser.add_argument("--model", required=True)
parser.add_argument("--test", required=True, help="Path to test.txt")
parser.add_argument("--output_dir", type=str, default="./evaluation_results")
parser.add_argument("--seed", type=int, default=7)

parser.add_argument("--rnn_type", type=str, default="LSTM")
parser.add_argument("--hidden_size", type=int, default=250)
parser.add_argument("--embed_size", type=int, default=250)
parser.add_argument("--batch_size", type=int, default=50)
parser.add_argument("--latent_size", type=int, default=32)
parser.add_argument("--depthT", type=int, default=15)
parser.add_argument("--depthG", type=int, default=15)
parser.add_argument("--diterT", type=int, default=1)
parser.add_argument("--diterG", type=int, default=3)
parser.add_argument("--dropout", type=float, default=0.0)

args = parser.parse_args()

with open(args.vocab, 'r', encoding='utf-8-sig') as vf:
    vocab = [x.strip("\r\n ").split() for x in vf if x.strip()]
args.vocab = PairVocab(vocab, cuda=torch.cuda.is_available())

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = HierVAE(args).to(device)
ckpt = torch.load(args.model, map_location=device)
model.load_state_dict(ckpt[0])
model.eval()

torch.manual_seed(args.seed)
random.seed(args.seed)

out_dir = Path(args.output_dir)
out_dir.mkdir(parents=True, exist_ok=True)

# Helpers

# canonicalize SMILES
def canonicalize(sm):
    try:
        m = Chem.MolFromSmiles(sm)
        return Chem.MolToSmiles(m) if m else None
    except:
        return None

# compute fingerprint
def fp(sm):
    m = Chem.MolFromSmiles(sm)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None

# compute tanimoto similarity
def tanimoto(a,b):
    if a is None or b is None:
        return np.nan
    return Chem.DataStructs.TanimotoSimilarity(a,b)

# count atoms and bonds
def count_atoms(sm):
    m = Chem.MolFromSmiles(sm)
    return m.GetNumAtoms() if m else None

# count bonds
def count_bonds(sm):
    m = Chem.MolFromSmiles(sm)
    return m.GetNumBonds() if m else None

def count_rings(sm):
    """Count rings in molecule (error diagnostic)"""
    m = Chem.MolFromSmiles(sm)
    if m is None:
        return None
    return len(Chem.GetSSSR(m))

def count_aromatic_atoms(sm):
    """Count aromatic atoms (error diagnostic)"""
    m = Chem.MolFromSmiles(sm)
    if m is None:
        return None
    return sum(1 for atom in m.GetAtoms() if atom.GetIsAromatic())

def count_heteroatoms(sm):
    """Count non-carbon atoms (error diagnostic)"""
    m = Chem.MolFromSmiles(sm)
    if m is None:
        return None
    return sum(1 for atom in m.GetAtoms() if atom.GetAtomicNum() != 6)


# main function to reconstruct one molecule (returns "INVALID" if fails) [smiles = list of one]
def reconstruct_one(smiles):
    try:
        graphs, tensors, orders = MolGraph.tensorize([smiles], args.vocab, args.atom_vocab)
        batch = (graphs, tensors, orders)
        with torch.no_grad():
            out_list = model.reconstruct(batch)
        return out_list[0] if out_list else "INVALID"
    except Exception:
        return "INVALID"

# load test set
with open(args.test, 'r') as tf:
    test_smiles = [l.strip() for l in tf if l.strip()]

# evaluate
records = []
for i, sm in tqdm(enumerate(test_smiles), total=len(test_smiles)):
    sm_in = sm
    sm_out = reconstruct_one(sm_in)
    valid_out = 1 if canonicalize(sm_out) is not None else 0
    cin = canonicalize(sm_in)
    cout = canonicalize(sm_out) if valid_out else None
    exact = 1 if (cin and cout and cin == cout) else 0
    tin = fp(cin) if cin else None
    tout = fp(cout) if cout else None
    tani = tanimoto(tin, tout)
    ai = count_atoms(cin) if cin else None
    ao = count_atoms(cout) if cout else None
    bi = count_bonds(cin) if cin else None
    bo = count_bonds(cout) if cout else None
    d_atoms = (ao - ai) if (ai is not None and ao is not None) else np.nan
    d_bonds = (bo - bi) if (bi is not None and bo is not None) else np.nan

    # Additional error diagnostics
    ri = count_rings(cin) if cin else None
    ro = count_rings(cout) if cout else None
    d_rings = (ro - ri) if (ri is not None and ro is not None) else np.nan

    arom_i = count_aromatic_atoms(cin) if cin else None
    arom_o = count_aromatic_atoms(cout) if cout else None
    d_aromatic = (arom_o - arom_i) if (arom_i is not None and arom_o is not None) else np.nan

    hetero_i = count_heteroatoms(cin) if cin else None
    hetero_o = count_heteroatoms(cout) if cout else None
    d_heteroatoms = (hetero_o - hetero_i) if (hetero_i is not None and hetero_o is not None) else np.nan

    records.append({
        'id': i,
        'smiles_in': sm_in,
        'smiles_out': cout if valid_out else '',
        'valid_out': valid_out,
        'exact_match': exact,
        'tanimoto': tani,
        'delta_atoms': d_atoms,
        'delta_bonds': d_bonds,
        'delta_rings': d_rings,
        'delta_aromatic': d_aromatic,
        'delta_heteroatoms': d_heteroatoms,
        'atoms_in': ai,
        'atoms_out': ao,
    })

# Save results
df = pd.DataFrame(records)

# Replace NaN with "NA" for invalid cases
df.loc[df['valid_out'] == 0, 'tanimoto'] = 'NA'
df.loc[df['valid_out'] == 0, 'delta_atoms'] = 'NA'
df.loc[df['valid_out'] == 0, 'delta_bonds'] = 'NA'
df.loc[df['valid_out'] == 0, 'delta_rings'] = 'NA'
df.loc[df['valid_out'] == 0, 'delta_aromatic'] = 'NA'
df.loc[df['valid_out'] == 0, 'delta_heteroatoms'] = 'NA'
df.loc[df['valid_out'] == 0, 'atoms_out'] = 'NA'

df.to_csv(out_dir / 'results.csv', index=False)

# Summary stats
vm = df['valid_out'] == 1
summary_data = {
    'Metric': [
        'Total Molecules',
        'Valid Reconstructions',
        'Validity Rate (%)',
        'Exact Matches (valid)',
        'Exact Match Rate (%)',
        'Mean Tanimoto (valid)',
        'Std Tanimoto (valid)',
        'Min Tanimoto (valid)',
        'Max Tanimoto (valid)',
        'Mean Δ Atoms (valid)',
        'Mean Δ Bonds (valid)',
        'Mean Δ Rings (valid)',
        'Mean Δ Aromatic (valid)',
        'Mean Δ Heteroatoms (valid)',
    ],
    'Value': [
        len(df),
        vm.sum(),
        f"{vm.mean() * 100:.2f}",
        (df.loc[vm, 'exact_match'] == 1).sum() if vm.sum() > 0 else 0,
        f"{(df.loc[vm, 'exact_match'] == 1).mean() * 100:.2f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'tanimoto'], errors='coerce').mean():.4f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'tanimoto'], errors='coerce').std():.4f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'tanimoto'], errors='coerce').min():.4f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'tanimoto'], errors='coerce').max():.4f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'delta_atoms'], errors='coerce').mean():.2f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'delta_bonds'], errors='coerce').mean():.2f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'delta_rings'], errors='coerce').mean():.2f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'delta_aromatic'], errors='coerce').mean():.2f}" if vm.sum() > 0 else "N/A",
        f"{pd.to_numeric(df.loc[vm, 'delta_heteroatoms'], errors='coerce').mean():.2f}" if vm.sum() > 0 else "N/A",
    ]
}
pd.DataFrame(summary_data).to_csv(out_dir / 'summary_stats.csv', index=False)
print("\n=== Summary Statistics ===")
print(pd.DataFrame(summary_data).to_string(index=False))

print(f"\n✓ Evaluation complete! Results saved to {out_dir}")
print(f"  - results.csv: Full metric table with {len(df)} molecules")
print(f"  - summary_stats.csv: Aggregate statistics")
print(f"\nTo generate visualizations, run:")
print(f"  python visualize_results.py --input {out_dir / 'results.csv'}")
