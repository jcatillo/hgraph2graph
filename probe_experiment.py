import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import rdkit.Chem as Chem
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.rdMolDescriptors as MolDescriptors
import os, sys, argparse
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from hgraph import * 
from tqdm import tqdm
import csv

def get_mol_props(smiles):
    """Compute basic molecular properties from SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'mw': None, 'rings': None, 'atoms': None, 'logp': None}
    return {
        'mw': Descriptors.MolWt(mol),
        'rings': MolDescriptors.CalcNumRings(mol),
        'atoms': mol.GetNumAtoms(),
        'logp': Descriptors.MolLogP(mol),
    }

def make_cuda(tensors):
    tree_tensors, graph_tensors = tensors
    make_tensor = lambda x: x if type(x) is torch.Tensor else torch.tensor(x)
    tree_tensors = [make_tensor(x).cuda().long() for x in tree_tensors[:-1]] + [tree_tensors[-1]]
    graph_tensors = [make_tensor(x).cuda().long() for x in graph_tensors[:-1]] + [graph_tensors[-1]]
    return tree_tensors, graph_tensors

def run_probing(args):
    # 1. Load Vocab
    vocab = [x.strip("\r\n ").split() for x in open(args.vocab)]
    args.vocab = PairVocab(vocab)
    
    # 2. Load Model
    model = HierVAE(args).cuda()
    model.load_state_dict(torch.load(args.model)[0])
    model.eval()
    print(f"Loaded model from {args.model}")

    # 3. Load Data
    with open(args.test_data) as f:
        data = [line.strip() for line in f]
    
    dataset = MoleculeDataset(data, args.vocab, args.atom_vocab, args.batch_size)
    print(f"Loaded {len(data)} test molecules. Processing in {len(dataset)} batches.")

    # Collect hidden states with their decode steps
    all_data = []  # list of (step, vec_np, props)
    
    count_skipped = 0

    print("Collecting decoder hidden states...")
    with torch.no_grad():
        for i in tqdm(range(len(dataset))):
            batch_smiles = dataset.batches[i]
            batch_tensors = dataset[i]
            
            # Ground truth properties per SMILES
            batch_props = [get_mol_props(s) for s in batch_smiles]

            # 2. Run the Model (The Subject)
            graphs, tensors, orders = batch_tensors
            tree_tensors, graph_tensors = tensors = make_cuda(tensors)
            
            try:
                root_vecs, tree_vecs, _, graph_vecs = model.encoder(tree_tensors, graph_tensors)
                root_vecs = model.R_mean(root_vecs) 
                
                # Prepare inputs for the decoder
                decoder_inputs = (root_vecs, root_vecs, root_vecs)
                
                # Capture internal decoder states
                hidden_states_log = model.decoder(decoder_inputs, graphs, tensors, orders, probe_mode=True)
                
                # hidden_states_log is a list of tuples: (step_number, batch_index, hidden_vector)
                
                for step, b_idx, vec in hidden_states_log:
                    b_idx = int(b_idx)
                    props = batch_props[b_idx]
                    if props['mw'] is None: continue
                    
                    vec_np = vec.numpy()
                    
                    all_data.append((step, vec_np, props))
                        
            except Exception as e:
                # Some batches might fail due to graph issues or incompatible sizes in corner cases
                # print(f"Skipping batch due to error: {e}")
                count_skipped += 1
                continue

    print(f"Collected. Failed batches: {count_skipped}")

    if not all_data:
        print("No hidden states collected; aborting.")
        return

    # Dynamic stage thresholds by step (quantiles: 33%, 66%)
    steps = np.array([s for s, _, _ in all_data])
    smin, smax = steps.min(), steps.max()
    t1, t2 = np.quantile(steps, [1/3, 2/3])
    print(f"Step range: {smin}..{smax} | Quantiles -> Early<= {t1:.2f}, Mid<= {t2:.2f}, Late> {t2:.2f}")

    # Bucket by dynamic stages
    X_early, y_mw_early, y_ring_early, y_atoms_early, y_logp_early = [], [], [], [], []
    X_mid, y_mw_mid, y_ring_mid, y_atoms_mid, y_logp_mid = [], [], [], [], []
    X_late, y_mw_late, y_ring_late, y_atoms_late, y_logp_late = [], [], [], [], []

    for step, vec_np, props in all_data:
        if props['mw'] is None:
            continue
        if step <= t1:
            X_early.append(vec_np); y_mw_early.append(props['mw']); y_ring_early.append(props['rings']); y_atoms_early.append(props['atoms']); y_logp_early.append(props['logp'])
        elif step <= t2:
            X_mid.append(vec_np); y_mw_mid.append(props['mw']); y_ring_mid.append(props['rings']); y_atoms_mid.append(props['atoms']); y_logp_mid.append(props['logp'])
        else:
            X_late.append(vec_np); y_mw_late.append(props['mw']); y_ring_late.append(props['rings']); y_atoms_late.append(props['atoms']); y_logp_late.append(props['logp'])

    results = []

    def train_and_eval(name, X, y_mw, y_ring, y_atoms, y_logp):
        if len(X) < 100:
            print(f"[{name}] Not enough data points ({len(X)}). Skipping.")
            return None
            
        X = np.array(X)
        y_mw = np.array(y_mw)
        y_ring = np.array(y_ring)
        y_atoms = np.array(y_atoms)
        y_logp = np.array(y_logp)
        
        # Molecular weight (regression)
        reg = LinearRegression()
        reg.fit(X, y_mw)
        pred_mw = reg.predict(X)
        score_mw = r2_score(y_mw, pred_mw)
        
        # Ring count (regression)
        reg_ring = LinearRegression()
        reg_ring.fit(X, y_ring)
        pred_ring = reg_ring.predict(X)
        score_ring = r2_score(y_ring, pred_ring)
        
        # Atom count (regression)
        reg_atoms = LinearRegression()
        reg_atoms.fit(X, y_atoms)
        pred_atoms = reg_atoms.predict(X)
        score_atoms = r2_score(y_atoms, pred_atoms)

        # logP (regression)
        reg_logp = LinearRegression()
        reg_logp.fit(X, y_logp)
        pred_logp = reg_logp.predict(X)
        score_logp = r2_score(y_logp, pred_logp)

        print(f"--- {name} Results ---")
        print(f"Samples: {len(X)}")
        print(f"MW R2:   {score_mw:.4f}")
        print(f"Rings R2:{score_ring:.4f}")
        print(f"Atoms R2:{score_atoms:.4f}")
        print(f"logP R2: {score_logp:.4f}")
        print("-----------------------")
        
        res = {
            'Stage': name,
            'Samples': len(X),
            'MW_R2': score_mw,
            'Ring_R2': score_ring,
            'Atoms_R2': score_atoms,
            'LogP_R2': score_logp
        }
        return res

    result_early = train_and_eval("Early", X_early, y_mw_early, y_ring_early, y_atoms_early, y_logp_early)
    if result_early:
        results.append(result_early)
    
    result_mid = train_and_eval("Mid", X_mid, y_mw_mid, y_ring_mid, y_atoms_mid, y_logp_mid)
    if result_mid:
        results.append(result_mid)
    
    result_late = train_and_eval("Late", X_late, y_mw_late, y_ring_late, y_atoms_late, y_logp_late)
    if result_late:
        results.append(result_late)
    
    # Write results to CSV if output file is specified
    if args.output_csv and results:
        # Build dynamic header
        all_keys = set()
        for r in results:
            all_keys.update(r.keys())
        preferred = ['Stage', 'Samples', 'MW_R2', 'Ring_R2', 'Atoms_R2', 'LogP_R2']
        fieldnames = [k for k in preferred if k in all_keys] + [k for k in sorted(all_keys) if k not in preferred]
        with open(args.output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"Results saved to {args.output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vocab", required=True)
    parser.add_argument("--atom_vocab", default=common_atom_vocab)
    parser.add_argument("--model", required=True)
    parser.add_argument("--test_data", required=True, help="Path to test.txt")
    parser.add_argument("--output_csv", default=None, help="Path to output CSV file for results")
    
    # Model Args (Must match generate.py defaults)
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
    
    run_probing(args)
