# Part_B/evaluate_part_b.py - CHECKPOINT DYNAMICS EVALUATION
import torch
import argparse
import random
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
import glob

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from hgraph import *
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)

parser = argparse.ArgumentParser(description="Evaluate reconstruction across training checkpoints (Part B)")
parser.add_argument("--vocab", required=True)
parser.add_argument("--atom_vocab", default=common_atom_vocab)
parser.add_argument("--checkpoint_dir", required=True, help="Directory containing checkpoints")
parser.add_argument("--test", required=True, help="Path to test.txt")
parser.add_argument("--output_dir", type=str, default="./Part_B")
parser.add_argument("--num_checkpoints", type=int, default=6, help="Number of evenly-spaced checkpoints (minimum 5 required)")
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

torch.manual_seed(args.seed)
random.seed(args.seed)

out_dir = Path(args.output_dir)
out_dir.mkdir(parents=True, exist_ok=True)

# Load vocabulary
with open(args.vocab, 'r', encoding='utf-8-sig') as vf:
    vocab = [x.strip("\r\n ").split() for x in vf if x.strip()]
args.vocab = PairVocab(vocab, cuda=torch.cuda.is_available())

# Load test set
with open(args.test, 'r') as tf:
    test_smiles = [l.strip() for l in tf if l.strip()]

print(f"Loaded {len(test_smiles)} test molecules")

# ============ HELPER FUNCTIONS ============

def canonicalize(sm):
    try:
        m = Chem.MolFromSmiles(sm)
        return Chem.MolToSmiles(m) if m else None
    except:
        return None

def fp(sm):
    m = Chem.MolFromSmiles(sm)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=4096) if m else None

def tanimoto(a, b):
    if a is None or b is None:
        return np.nan
    return Chem.DataStructs.TanimotoSimilarity(a, b)

def reconstruct_one(smiles, model):
    try:
        graphs, tensors, orders = MolGraph.tensorize([smiles], args.vocab, args.atom_vocab)
        batch = (graphs, tensors, orders)
        with torch.no_grad():
            out_list = model.reconstruct(batch)
        return out_list[0] if out_list else "INVALID"
    except Exception:
        return "INVALID"

# ============ SELECT CHECKPOINTS ============

checkpoint_dir = Path(args.checkpoint_dir)
all_checkpoints = sorted(glob.glob(str(checkpoint_dir / "model.ckpt.*")))

checkpoint_numbers = []
for ckpt_path in all_checkpoints:
    try:
        num = int(Path(ckpt_path).name.split('.')[-1])
        checkpoint_numbers.append(num)
    except:
        pass

checkpoint_numbers = sorted(set(checkpoint_numbers))
print(f"Found {len(checkpoint_numbers)} checkpoints total")

# Ensure at least 5 checkpoints
if len(checkpoint_numbers) < 5:
    print(f"Warning: Only {len(checkpoint_numbers)} checkpoints found, but at least 5 required for Part B")
    selected_ckpts = checkpoint_numbers
else:
    # Auto-select evenly spaced checkpoints (minimum 5)
    if len(checkpoint_numbers) >= args.num_checkpoints:
        indices = np.linspace(0, len(checkpoint_numbers) - 1, args.num_checkpoints, dtype=int)
        selected_ckpts = [checkpoint_numbers[i] for i in indices]
    else:
        selected_ckpts = checkpoint_numbers

selected_ckpts = sorted(selected_ckpts)
print(f"\n{'='*70}")
print(f"Evaluating {len(selected_ckpts)} TRAINING CHECKPOINTS")
print(f"Checkpoints: {selected_ckpts}")
print(f"Test set size: {len(test_smiles)} molecules")
print(f"{'='*70}\n")

# ============ EVALUATE ACROSS CHECKPOINTS ============

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
checkpoint_results = []

for idx, ckpt_num in enumerate(selected_ckpts, 1):
    print(f"\n{'='*70}")
    print(f"[{idx}/{len(selected_ckpts)}] CHECKPOINT {ckpt_num} EVALUATION")
    print(f"{'='*70}")
    
    ckpt_path = checkpoint_dir / f"model.ckpt.{ckpt_num}"
    
    if not ckpt_path.exists():
        print(f"❌ Checkpoint {ckpt_path} not found, skipping...")
        continue
    
    # Load model
    try:
        model = HierVAE(args).to(device)
        ckpt = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ckpt[0])
        model.eval()
        print(f"✓ Loaded model from {ckpt_path}")
    except Exception as e:
        print(f"❌ Error loading checkpoint {ckpt_num}: {e}")
        continue
    
    # Evaluate on test set
    exact_matches = []
    tanimoto_scores = []
    valid_count = 0
    invalid_smiles = []
    
    # Track results by stage
    early_exact = []
    early_tanimoto = []
    mid_exact = []
    mid_tanimoto = []
    late_exact = []
    late_tanimoto = []
    
    # Define stage boundaries
    third = len(test_smiles) // 3
    early_end = third
    mid_end = 2 * third
    
    print(f"Processing {len(test_smiles)} test molecules in 3 stages...")
    print(f"  Early: molecules 0-{early_end}")
    print(f"  Mid:   molecules {early_end}-{mid_end}")
    print(f"  Late:  molecules {mid_end}-{len(test_smiles)}\n")
    
    for i, sm_in in tqdm(enumerate(test_smiles), total=len(test_smiles), 
                         desc=f"  Evaluating", leave=True):
        sm_out = reconstruct_one(sm_in, model)
        
        valid_out = 1 if canonicalize(sm_out) is not None else 0
        cin = canonicalize(sm_in)
        cout = canonicalize(sm_out) if valid_out else None
        exact = 1 if (cin and cout and cin == cout) else 0

        tin = fp(cin) if cin else None
        tout = fp(cout) if cout else None
        tani = tanimoto(tin, tout)
        
        if valid_out:
            valid_count += 1
        else:
            invalid_smiles.append((sm_in, sm_out))
            
        exact_matches.append(exact)
        if not np.isnan(tani):
            tanimoto_scores.append(tani)
        
        # Categorize by stage
        if i < early_end:
            early_exact.append(exact)
            if not np.isnan(tani):
                early_tanimoto.append(tani)
        elif i < mid_end:
            mid_exact.append(exact)
            if not np.isnan(tani):
                mid_tanimoto.append(tani)
        else:
            late_exact.append(exact)
            if not np.isnan(tani):
                late_tanimoto.append(tani)
    
    # Calculate statistics by stage
    early_exact_rate = np.mean(early_exact) if early_exact else 0.0
    early_tanimoto_mean = np.mean(early_tanimoto) if early_tanimoto else 0.0
    
    mid_exact_rate = np.mean(mid_exact) if mid_exact else 0.0
    mid_tanimoto_mean = np.mean(mid_tanimoto) if mid_tanimoto else 0.0
    
    late_exact_rate = np.mean(late_exact) if late_exact else 0.0
    late_tanimoto_mean = np.mean(late_tanimoto) if late_tanimoto else 0.0
    
    # Calculate overall statistics
    exact_match_rate = np.mean(exact_matches) if exact_matches else 0.0
    mean_tanimoto = np.mean(tanimoto_scores) if tanimoto_scores else 0.0
    median_tanimoto = np.median(tanimoto_scores) if tanimoto_scores else 0.0
    validity_rate = valid_count / len(test_smiles)
    
    # Print results in structured format
    print(f"\n{'─'*70}")
    print(f"CHECKPOINT {ckpt_num} - RECONSTRUCTION METRICS BY STAGE")
    print(f"{'─'*70}")
    
    print(f"\n--- EARLY STAGE (Molecules 0-{early_end}) ---")
    print(f"Samples: {len(early_exact)}")
    print(f"Exact Match Rate:        {early_exact_rate*100:6.2f}%")
    print(f"Mean Tanimoto:           {early_tanimoto_mean:6.4f}")
    
    print(f"\n--- MID STAGE (Molecules {early_end}-{mid_end}) ---")
    print(f"Samples: {len(mid_exact)}")
    print(f"Exact Match Rate:        {mid_exact_rate*100:6.2f}%")
    print(f"Mean Tanimoto:           {mid_tanimoto_mean:6.4f}")
    
    print(f"\n--- LATE STAGE (Molecules {mid_end}-{len(test_smiles)}) ---")
    print(f"Samples: {len(late_exact)}")
    print(f"Exact Match Rate:        {late_exact_rate*100:6.2f}%")
    print(f"Mean Tanimoto:           {late_tanimoto_mean:6.4f}")
    
    print(f"\n{'─'*70}")
    print(f"CHECKPOINT {ckpt_num} - OVERALL METRICS")
    print(f"{'─'*70}")
    print(f"Total Molecules Evaluated: {len(test_smiles)}")
    print(f"Valid Reconstructions:     {valid_count}/{len(test_smiles)}")
    print(f"Validity Rate:             {validity_rate*100:6.2f}% (Chemical validity)")
    print(f"")
    print(f"Exact Match Rate:          {exact_match_rate*100:6.2f}% (Perfect reconstruction)")
    print(f"Mean Tanimoto Similarity:  {mean_tanimoto:6.4f} (Avg fingerprint similarity)")
    print(f"Median Tanimoto:           {median_tanimoto:6.4f} (Median fingerprint similarity)")
    print(f"")
    print(f"Invalid Outputs:           {len(invalid_smiles)}/{len(test_smiles)} molecules")
    print(f"{'─'*70}\n")
    
    checkpoint_results.append({
        'checkpoint': ckpt_num,
        'exact_match_rate': exact_match_rate,
        'mean_tanimoto': mean_tanimoto,
        'median_tanimoto': median_tanimoto,
        'validity_rate': validity_rate,
        'valid_molecules': valid_count,
        'total_molecules': len(test_smiles),
        'early_exact_rate': early_exact_rate,
        'early_tanimoto': early_tanimoto_mean,
        'mid_exact_rate': mid_exact_rate,
        'mid_tanimoto': mid_tanimoto_mean,
        'late_exact_rate': late_exact_rate,
        'late_tanimoto': late_tanimoto_mean,
    })

# ============ SAVE RESULTS ============

df_dynamics = pd.DataFrame(checkpoint_results)
df_dynamics.to_csv(out_dir / 'checkpoint_dynamics.csv', index=False)

print("\n" + "="*70)
print("CHECKPOINT DYNAMICS - FINAL SUMMARY")
print("="*70)
print(df_dynamics.to_string(index=False))
print("="*70)

# Print training progression analysis
print("\n" + "="*70)
print("TRAINING PROGRESSION ANALYSIS")
print("="*70)
print(f"\nExact Match Accuracy Evolution:")
print(f"  Early (Ckpt {selected_ckpts[0]}):  {df_dynamics.iloc[0]['exact_match_rate']*100:6.2f}%")
print(f"  Mid   (Ckpt {selected_ckpts[len(selected_ckpts)//2]}):  {df_dynamics.iloc[len(df_dynamics)//2]['exact_match_rate']*100:6.2f}%")
print(f"  Late  (Ckpt {selected_ckpts[-1]}): {df_dynamics.iloc[-1]['exact_match_rate']*100:6.2f}%")

print(f"\nTanimoto Similarity Evolution:")
print(f"  Early (Ckpt {selected_ckpts[0]}):  {df_dynamics.iloc[0]['mean_tanimoto']:6.4f}")
print(f"  Mid   (Ckpt {selected_ckpts[len(selected_ckpts)//2]}):  {df_dynamics.iloc[len(df_dynamics)//2]['mean_tanimoto']:6.4f}")
print(f"  Late  (Ckpt {selected_ckpts[-1]}): {df_dynamics.iloc[-1]['mean_tanimoto']:6.4f}")

print(f"\nValidity Rate Evolution:")
print(f"  Early (Ckpt {selected_ckpts[0]}):  {df_dynamics.iloc[0]['validity_rate']*100:6.2f}%")
print(f"  Mid   (Ckpt {selected_ckpts[len(selected_ckpts)//2]}):  {df_dynamics.iloc[len(df_dynamics)//2]['validity_rate']*100:6.2f}%")
print(f"  Late  (Ckpt {selected_ckpts[-1]}): {df_dynamics.iloc[-1]['validity_rate']*100:6.2f}%")

print("\n" + "="*70)
print(f"✓ Checkpoint dynamics saved to: {out_dir / 'checkpoint_dynamics.csv'}")
print(f"\nNext step - Generate visualizations:")
print(f"  python Part_B/plot_part_b.py --input {out_dir / 'checkpoint_dynamics.csv'}")
print("="*70)