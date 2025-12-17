import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import rdkit.Chem as Chem
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.rdMolDescriptors as MolDescriptors
import os, sys, argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from tqdm import tqdm

# Add parent directory to path to import hgraph
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from hgraph import *

def get_mol_props(smiles):
    """Extract molecular properties from SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None, None, None, None
    mw = Descriptors.MolWt(mol)
    rings = MolDescriptors.CalcNumRings(mol)
    atoms = mol.GetNumAtoms()
    logp = Descriptors.MolLogP(mol)
    return mw, rings, atoms, logp

def make_cuda(tensors):
    tree_tensors, graph_tensors = tensors
    make_tensor = lambda x: x if type(x) is torch.Tensor else torch.tensor(x)
    tree_tensors = [make_tensor(x).cuda().long() for x in tree_tensors[:-1]] + [tree_tensors[-1]]
    graph_tensors = [make_tensor(x).cuda().long() for x in graph_tensors[:-1]] + [graph_tensors[-1]]
    return tree_tensors, graph_tensors

def train_and_eval(name, X, y_mw, y_ring, y_atoms, y_logp):
    """Train linear probes and evaluate R2 scores"""
    if len(X) < 100:
        print(f"[{name}] Not enough data points ({len(X)}). Skipping.")
        return None
        
    X = np.array(X)
    y_mw = np.array(y_mw)
    y_ring = np.array(y_ring)
    y_atoms = np.array(y_atoms)
    y_logp = np.array(y_logp)
    
    results = {'Stage': name, 'Samples': len(X)}
    
    # PROBE 1: Molecular Weight
    reg = LinearRegression()
    reg.fit(X, y_mw)
    pred_mw = reg.predict(X)
    score_mw = r2_score(y_mw, pred_mw)
    results['MW_R2'] = score_mw
    
    # PROBE 2: Ring Count
    reg_ring = LinearRegression()
    reg_ring.fit(X, y_ring)
    pred_ring = reg_ring.predict(X)
    score_ring = r2_score(y_ring, pred_ring)
    results['Ring_R2'] = score_ring
    
    # PROBE 3: Atom Count
    reg_atoms = LinearRegression()
    reg_atoms.fit(X, y_atoms)
    pred_atoms = reg_atoms.predict(X)
    score_atoms = r2_score(y_atoms, pred_atoms)
    results['Atoms_R2'] = score_atoms
    
    # PROBE 4: LogP
    reg_logp = LinearRegression()
    reg_logp.fit(X, y_logp)
    pred_logp = reg_logp.predict(X)
    score_logp = r2_score(y_logp, pred_logp)
    results['LogP_R2'] = score_logp
    
    print(f"\n--- {name} Results ---")
    print(f"Samples: {len(X)}")
    print(f"Molecular Weight R2: {score_mw:.4f}")
    print(f"Ring Count R2:       {score_ring:.4f}")
    print(f"Atom Count R2:       {score_atoms:.4f}")
    print(f"LogP R2:             {score_logp:.4f}")
    print("-" * 50)
    
    return results

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
    
    # Data collectors for all properties
    X_early, y_mw_early, y_ring_early, y_atoms_early, y_logp_early = [], [], [], [], []
    X_mid, y_mw_mid, y_ring_mid, y_atoms_mid, y_logp_mid = [], [], [], [], []
    X_late, y_mw_late, y_ring_late, y_atoms_late, y_logp_late = [], [], [], [], []
    
    count_skipped = 0

    print("Starting Brain Scan (Data Collection)...")
    with torch.no_grad():
        for i in tqdm(range(len(dataset))):
            batch_smiles = dataset.batches[i]
            batch_tensors = dataset[i]
            
            # 1. Establish Reality (Ground Truth)
            batch_props = []
            for s in batch_smiles:
                mw, rings, atoms, logp = get_mol_props(s)
                batch_props.append({'mw': mw, 'rings': rings, 'atoms': atoms, 'logp': logp})

            # 2. Run the Model
            graphs, tensors, orders = batch_tensors
            tree_tensors, graph_tensors = tensors = make_cuda(tensors)
            
            try:
                root_vecs, tree_vecs, _, graph_vecs = model.encoder(tree_tensors, graph_tensors)
                root_vecs = model.R_mean(root_vecs) 
                
                # Prepare inputs for the decoder
                decoder_inputs = (root_vecs, root_vecs, root_vecs)
                
                # 3. The Scan - capture internal states
                hidden_states_log = model.decoder(decoder_inputs, graphs, tensors, orders, probe_mode=True)
                
                for step, b_idx, vec in hidden_states_log:
                    b_idx = int(b_idx)
                    props = batch_props[b_idx]
                    if props['mw'] is None: continue
                    
                    vec_np = vec.cpu().numpy()
                    
                    # 4. Bucketing by stage
                    if step <= 2:
                        X_early.append(vec_np)
                        y_mw_early.append(props['mw'])
                        y_ring_early.append(props['rings'])
                        y_atoms_early.append(props['atoms'])
                        y_logp_early.append(props['logp'])
                    elif 5 <= step <= 7:
                        X_mid.append(vec_np)
                        y_mw_mid.append(props['mw'])
                        y_ring_mid.append(props['rings'])
                        y_atoms_mid.append(props['atoms'])
                        y_logp_mid.append(props['logp'])
                    elif 10 <= step <= 15:
                        X_late.append(vec_np)
                        y_mw_late.append(props['mw'])
                        y_ring_late.append(props['rings'])
                        y_atoms_late.append(props['atoms'])
                        y_logp_late.append(props['logp'])
                        
            except Exception as e:
                count_skipped += 1
    
    print(f"Skipped {count_skipped} batches due to errors\n")
    
    # Train and evaluate probes
    results_list = []
    results_list.append(train_and_eval("Early", X_early, y_mw_early, y_ring_early, y_atoms_early, y_logp_early))
    results_list.append(train_and_eval("Mid", X_mid, y_mw_mid, y_ring_mid, y_atoms_mid, y_logp_mid))
    results_list.append(train_and_eval("Late", X_late, y_mw_late, y_ring_late, y_atoms_late, y_logp_late))
    
    # Save results to CSV
    results_df = pd.DataFrame([r for r in results_list if r is not None])
    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(args.output_dir, 'probe_results.csv')
    results_df.to_csv(output_file, index=False)
    print(f"\n✓ Results saved to {output_file}")
    print(results_df.to_string(index=False))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vocab", required=True)
    parser.add_argument("--atom_vocab", default=common_atom_vocab)
    parser.add_argument("--model", required=True)
    parser.add_argument("--test_data", required=True, help="Path to test.txt")
    parser.add_argument("--output_dir", type=str, default="./Part_C", help="Directory to save results")
    
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
