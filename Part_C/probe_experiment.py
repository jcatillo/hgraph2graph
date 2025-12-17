import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import rdkit.Chem as Chem
import rdkit.Chem.Descriptors as Descriptors
import rdkit.Chem.rdMolDescriptors as MolDescriptors
import os, sys, argparse
import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import r2_score, accuracy_score
from hgraph import * 
from tqdm import tqdm

def get_mol_props(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None, None
    mw = Descriptors.MolWt(mol)
    rings = MolDescriptors.CalcNumRings(mol)
    return mw, rings

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
    
    X_early, y_mw_early, y_ring_early = [], [], []
    X_mid, y_mw_mid, y_ring_mid = [], [], []
    X_late, y_mw_late, y_ring_late = [], [], []
    
    count_skipped = 0

    print("Starting Brain Scan (Data Collection)...")
    with torch.no_grad():
        for i in tqdm(range(len(dataset))):
            batch_smiles = dataset.batches[i]
            batch_tensors = dataset[i]
            
            # 1. Establish Reality (Ground Truth)
            # We calculate the ACTUAL properties of the molecule we are about to build.
            batch_props = []
            for s in batch_smiles:
                mw, rings = get_mol_props(s)
                batch_props.append({'mw': mw, 'rings': rings})

            # 2. Run the Model (The Subject)
            graphs, tensors, orders = batch_tensors
            tree_tensors, graph_tensors = tensors = make_cuda(tensors)
            
            try:
                root_vecs, tree_vecs, _, graph_vecs = model.encoder(tree_tensors, graph_tensors)
                root_vecs = model.R_mean(root_vecs) 
                
                # Prepare inputs for the decoder
                decoder_inputs = (root_vecs, root_vecs, root_vecs)
                
                # 3. The Scan 
                # Run the decoder with probe_mode=True to capture internal states
                hidden_states_log = model.decoder(decoder_inputs, graphs, tensors, orders, probe_mode=True)
                
                # hidden_states_log is a list of tuples: (step_number, batch_index, hidden_vector)
                
                for step, b_idx, vec in hidden_states_log:
                    b_idx = int(b_idx)
                    props = batch_props[b_idx]
                    if props['mw'] is None: continue
                    
                    vec_np = vec.numpy()
                    
                    # 4. Bucketing the "Thoughts"
                    # We sort the captured vectors into buckets based on when they occurred.
                    if step <= 2:
                        X_early.append(vec_np)
                        y_mw_early.append(props['mw'])
                        y_ring_early.append(props['rings'])
                    elif 5 <= step <= 7:
                        X_mid.append(vec_np)
                        y_mw_mid.append(props['mw'])
                        y_ring_mid.append(props['rings'])
                    elif 10 <= step <= 15:
                        X_late.append(vec_np)
                        y_mw_late.append(props['mw'])
                        y_ring_late.append(props['rings'])
                        
            except Exception as e:
                # Some batches might fail due to graph issues or incompatible sizes in corner cases
                # print(f"Skipping batch due to error: {e}")
                count_skipped += 1
                continue

    print(f"Probing Complete. Training Probes... (Skipped {count_skipped} batches)")

    def train_and_eval(name, X, y_mw, y_ring):
        if len(X) < 100:
            print(f"[{name}] Not enough data points ({len(X)}). Skipping.")
            return
            
        X = np.array(X)
        y_mw = np.array(y_mw)
        y_ring = np.array(y_ring)
        
        # PROBE 1: Molecular Weight (Regression)
        # Question: Can a simple line predict the final weight from this thought vector?
        reg = LinearRegression()
        reg.fit(X, y_mw)
        pred_mw = reg.predict(X)
        score_mw = r2_score(y_mw, pred_mw)
        
        # PROBE 2: Ring Count (Regression)
        # Question: Can we predict the final number of rings?
        reg_ring = LinearRegression()
        reg_ring.fit(X, y_ring)
        pred_ring = reg_ring.predict(X)
        score_ring = r2_score(y_ring, pred_ring)
        
        print(f"--- {name} Results ---")
        print(f"Samples: {len(X)}")
        print(f"Molecular Weight R2: {score_mw:.4f} (Ability to predict weight)")
        print(f"Ring Count R2:       {score_ring:.4f} (Ability to predict structure)")
        print("-----------------------")

    train_and_eval("Early (Step 0-2)", X_early, y_mw_early, y_ring_early)
    train_and_eval("Mid   (Step 5-7)", X_mid, y_mw_mid, y_ring_mid)
    train_and_eval("Late  (Step 10-15)", X_late, y_mw_late, y_ring_late)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vocab", required=True)
    parser.add_argument("--atom_vocab", default=common_atom_vocab)
    parser.add_argument("--model", required=True)
    parser.add_argument("--test_data", required=True, help="Path to test.txt")
    
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
