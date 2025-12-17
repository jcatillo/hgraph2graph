# Experiment Instructions: Parts A, B, C

This document provides step-by-step instructions for running the three core experiments analyzing the **Hierarchical Generation of Molecular Graphs using Structural Motifs (HGraph)** model on the **Finals dataset**.

---

## Table of Contents
1. [Environment Setup](#environment-setup)
2. [Dependencies](#dependencies)
3. [Part A: Autoregressive Decoding (Inference Gap)](#part-a-autoregressive-decoding-inference-gap)
4. [Part B: Checkpoint Dynamics](#part-b-checkpoint-dynamics)
5. [Part C: Decoder Localization (Linear Probing)](#part-c-decoder-localization-linear-probing)

---

## Environment Setup

### Prerequisites
- Python >= 3.6
- CUDA-capable GPU (recommended)
- Conda, pip, or uv package manager

### Create Virtual Environment

**Option 1: Using Conda (recommended)**
```bash
conda create -n hgraph_env python=3.8
conda activate hgraph_env
```

**Option 2: Using venv**
```bash
python -m venv hgraph_env
source hgraph_env/bin/activate  # On Windows: hgraph_env\Scripts\activate
```

### Install Dependencies

Install all required packages:
```bash
pip install -r requirements.txt
```

Or manually install core dependencies with CUDA support:

**For CUDA 11.8:**
```bash
pip install torch>=1.0.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install networkx rdkit numpy scikit-learn pandas tqdm
```

**For CUDA 12.1:**
```bash
pip install torch>=1.0.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
pip install networkx rdkit numpy scikit-learn pandas tqdm
```

**For CPU only (not recommended):**
```bash
pip install torch>=1.0.0 networkx rdkit numpy scikit-learn pandas tqdm
```

### Install Package
```bash
pip install -e .
```

---

## Dependencies

| Package | Version | Purpose | CUDA |
|---------|---------|---------|------|
| PyTorch | >= 1.0.0 | Deep learning framework | **Required** (11.8 or 12.1) |
| RDKit | >= 2019.03 | Molecular property calculations | — |
| NetworkX | Latest | Graph operations | — |
| NumPy | Latest | Numerical computing | — |
| Pandas | Latest | Data manipulation | — |
| scikit-learn | Latest | Machine learning (linear regression for probing) | — |
| tqdm | Latest | Progress bars | — |

**Optional:**
- Matplotlib / Seaborn: For visualization

---

## Part A: Autoregressive Decoding (Inference Gap)

**Objective:** Quantify the performance gap between **teacher-forced** decoding (training) and **autoregressive** decoding (inference).

### Quick Start

```bash
python Part_A/evaluate.py \
    --vocab data/Finals/vocab.txt \
    --model ckpt/Finals-ckpt/model.ckpt.2000 \
    --test data/Finals/test.txt \
    --output_dir ./Part_A
```

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--vocab` | *required* | Path to vocabulary file (Finals/vocab.txt) |
| `--model` | *required* | Path to trained model checkpoint |
| `--test` | *required* | Path to test SMILES file |
| `--output_dir` | ./Part_A | Directory to save results |
| `--batch_size` | 50 | Batch size for inference |
| `--hidden_size` | 250 | Hidden layer dimension |
| `--latent_size` | 32 | Latent space dimension |
| `--seed` | 7 | Random seed for reproducibility |

### Outputs

The script generates:
- `summary_stats.csv`: Reconstruction metrics (validity, Tanimoto similarity, structural deltas)
- `results.csv`: Per-molecule predictions with SMILES, validity, and similarity scores
- `comparisons/`: Qualitative examples of original vs. reconstructed molecules

### Key Metrics

- **Validity Rate (%)**: Fraction of chemically valid reconstructions
- **Exact Match Rate (%)**: Exact SMILES reconstruction (typically 0% due to inference)
- **Mean Tanimoto**: Average structural similarity (0–1, 1=identical)
- **Mean Δ Atoms**: Average difference in atom count between original and reconstructed
- **Mean Δ Bonds**: Average difference in bond count

### Example Output
```
Validity Rate: 94.61%
Mean Tanimoto: 0.1716
Mean Δ Atoms: 1.23
Mean Δ Bonds: 1.45
```

---

## Part B: Checkpoint Dynamics

**Objective:** Track how reconstruction performance evolves across training checkpoints.

### Quick Start

```bash
python Part_B/evaluate_part_b.py \
    --vocab data/Finals/vocab.txt \
    --checkpoint_dir ckpt/Finals-ckpt \
    --test data/Finals/test.txt \
    --output_dir ./Part_B
```

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--vocab` | *required* | Path to vocabulary file |
| `--checkpoint_dir` | *required* | Directory containing model.ckpt.* files |
| `--test` | *required* | Path to test SMILES file |
| `--output_dir` | ./Part_B | Directory to save results |
| `--batch_size` | 50 | Batch size for inference |
| `--skip_checkpoints` | 100 | Interval for checkpoint selection (e.g., every 100th) |

### Outputs

- `checkpoint_dynamics.csv`: Metrics vs. checkpoint iteration
  - Columns: `Checkpoint`, `Validity_Rate`, `Mean_Tanimoto`, `Mean_Atoms_Delta`, etc.
- `plot_checkpoint_dynamics.png`: Visualization of metric trends

### Analysis Focus

- **Training convergence**: Validity rate and similarity metrics vs. checkpoint
- **Peak performance**: Identify optimal checkpoint for downstream tasks
- **Overfitting detection**: Monitor test-set divergence

### Example Visualization Commands

```bash
python Part_B/plot_checkpoint_dynamics.py \
    --csv Part_B/checkpoint_dynamics.csv \
    --outdir ./Part_B
```

---

## Part C: Decoder Localization (Linear Probing)

**Objective:** Determine **where** molecular properties are encoded in the decoder's hidden states across training stages (early, mid, late).

### Setup

Ensure the decoder has `probe_mode=True` support in `hgraph/decoder.py`. This flag logs hidden states at specified steps during decoding.

### Quick Start

```bash
python Part_C/probe_experiment.py \
    --vocab data/Finals/vocab.txt \
    --model ckpt/Finals-ckpt/model.ckpt.2000 \
    --test_data data/Finals/test.txt \
    --output_dir ./Part_C
```

### Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--vocab` | *required* | Path to vocabulary file |
| `--atom_vocab` | common_atom_vocab | Atom vocabulary list |
| `--model` | *required* | Path to trained model checkpoint |
| `--test_data` | *required* | Path to test SMILES file |
| `--output_dir` | ./Part_C | Directory to save probe_results.csv |
| `--batch_size` | 50 | Batch size for data collection |
| `--hidden_size` | 250 | Hidden layer dimension |
| `--latent_size` | 32 | Latent space dimension |

### Execution Flow

1. **Initialize model** with loaded checkpoint
2. **Load test molecules** and extract ground-truth properties:
   - Molecular Weight (MW)
   - Ring Count
   - Atom Count
   - LogP (lipophilicity)
3. **Capture decoder states** at three stages:
   - **Early** (steps 0–2): Initial motif selection
   - **Mid** (steps 5–7): Mid-generation structure
   - **Late** (steps 10–15): Final refinement
4. **Linear probing**: Train separate regression probes for each property
   - Model: `sklearn.linear_model.LinearRegression`
   - Metric: R² score (0–1, 1=perfect prediction)
5. **Save results** to `Part_C/probe_results.csv`

### Outputs

**probe_results.csv:**
```
Stage,Samples,MW_R2,Ring_R2,Atoms_R2,LogP_R2
Early,5000,0.6234,0.5812,0.6021,0.5432
Mid,5000,0.8534,0.8102,0.8456,0.7934
Late,5000,0.9904,0.9969,0.9903,0.9712
```

### Visualization

Generate probe results plots:
```bash
python Part_C/visualize_probe.py \
    --csv Part_C/probe_results.csv \
    --outdir ./Part_C
```

This creates:
- `probe_results_line.png`: R² scores vs. stage
- `probe_results_grouped.png`: Grouped bar chart (property × stage)
- `probe_results_heatmap.png`: Heatmap of R² scores

### Interpretation

- **High R² (late)**: Property information becomes *more* localized in decoder states as generation progresses
- **Low R² (early)**: Initial states contain less structured property information
- **Gradual increase**: Evidence of hierarchical information encoding
- **Property ranking**: LogP typically harder to predict than MW (higher complexity)

---

## Running All Experiments

### Sequential Execution

```bash
# Part A: Inference evaluation
python Part_A/evaluate.py \
    --vocab data/Finals/vocab.txt \
    --model ckpt/Finals-ckpt/model.ckpt.2000 \
    --test data/Finals/test.txt

# Part B: Checkpoint dynamics
python Part_B/evaluate_part_b.py \
    --vocab data/Finals/vocab.txt \
    --checkpoint_dir ckpt/Finals-ckpt \
    --test data/Finals/test.txt

# Part C: Decoder localization
python Part_C/probe_experiment.py \
    --vocab data/Finals/vocab.txt \
    --model ckpt/Finals-ckpt/model.ckpt.2000 \
    --test_data data/Finals/test.txt

# Visualize Part C results
python Part_C/visualize_probe.py \
    --csv Part_C/probe_results.csv \
    --outdir ./Part_C
```

### Batch Script (Linux/Mac)

```bash
#!/bin/bash
VOCAB="data/Finals/vocab.txt"
MODEL="ckpt/Finals-ckpt/model.ckpt.2000"
TEST="data/Finals/test.txt"
CKPT_DIR="ckpt/Finals-ckpt"

echo "=== Part A ==="
python Part_A/evaluate.py --vocab $VOCAB --model $MODEL --test $TEST

echo "=== Part B ==="
python Part_B/evaluate_part_b.py --vocab $VOCAB --checkpoint_dir $CKPT_DIR --test $TEST

echo "=== Part C ==="
python Part_C/probe_experiment.py --vocab $VOCAB --model $MODEL --test_data $TEST
python Part_C/visualize_probe.py --csv Part_C/probe_results.csv --outdir ./Part_C

echo "=== All experiments complete ==="
```

---

## Troubleshooting

### Import Errors
```
ModuleNotFoundError: No module named 'hgraph'
```
**Solution:** Install the package with `pip install -e .` from the root directory.

### CUDA Out of Memory
Reduce `--batch_size`:
```bash
python Part_A/evaluate.py ... --batch_size 32
```

### Missing Checkpoints
Verify checkpoint files exist:
```bash
ls -la ckpt/Finals-ckpt/model.ckpt.*
```

### Vocabulary Mismatch
Ensure `--vocab` matches the model's training vocabulary.

---

## References

- **Original Paper:** Hierarchical Generation of Molecular Graphs using Structural Motifs  
  https://arxiv.org/pdf/2002.03230.pdf
- **Repository:** https://github.com/wengong-jin/hgraph2graph
- **Dataset:** Finals (subset of ChEMBL for evaluation)

---

## Contact & Citation

For detailed methodology, see `Part_D/training_explanation.md` and the corresponding flowcharts.
