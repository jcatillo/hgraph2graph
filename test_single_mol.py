from hgraph import MolGraph

# Test with the first molecule from all.txt
test_smiles = "Clc1ccccc1C(c1ccccc1)(c1ccccc1)n1ccnc1"

try:
    print(f"Testing molecule: {test_smiles}")
    hmol = MolGraph(test_smiles)
    print(f"Success! Clusters: {hmol.clusters}")
    print(f"Tree nodes: {len(hmol.mol_tree.nodes)}")
    
    vocab = set()
    for node, attr in hmol.mol_tree.nodes(data=True):
        smiles = attr["smiles"]
        vocab.add(attr["label"])
        for i, s in attr["inter_label"]:
            vocab.add((smiles, s))
    
    print(f"Vocab size: {len(vocab)}")
    print(f"Sample vocab items: {list(vocab)[:5]}")
    
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
