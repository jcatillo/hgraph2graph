import argparse
import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem

# Try to import MolToFile, fallback to alternative if not available
try:
    from rdkit.Chem.Draw import MolToFile
    MOLTOFILE_AVAILABLE = True
except ImportError:
    MOLTOFILE_AVAILABLE = False
    print("Warning: MolToFile not available, using alternative rendering method")
    from rdkit.Chem import Draw

from PIL import Image

# Suppress RDKit warnings
import rdkit
lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)

# Set environment variable for OpenMP
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

parser = argparse.ArgumentParser(description="Generate visualizations from evaluation results")
parser.add_argument("--input", "-i", required=True, help="Path to results.csv file")
parser.add_argument("--output_dir", "-o", default=None, help="Output directory for visualizations (default: same as input file directory)")
parser.add_argument("--num_examples", type=int, default=10, help="Number of top examples to visualize (default: 10)")
args = parser.parse_args()

results_file = Path(args.input)

if not results_file.exists():
    print(f"Error: {results_file} not found!")
    exit(1)

# Determine output directory
if args.output_dir:
    results_dir = Path(args.output_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
else:
    results_dir = results_file.parent

# Load results
print(f"Loading results from {results_file}...")
df = pd.read_csv(results_file)

# Filter valid molecules
vm = df['valid_out'] == 1
tan_valid = pd.to_numeric(df.loc[vm, 'tanimoto'], errors='coerce').dropna()

print("\n=== Generating Distribution Plots ===")

# 1. Tanimoto histogram
print("  Creating Tanimoto distribution...")
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(tan_valid, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
ax.axvline(tan_valid.mean(), color='r', linestyle='--', linewidth=2, label=f'Mean: {tan_valid.mean():.3f}')
ax.set_xlabel('Tanimoto Similarity', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('Distribution of Tanimoto Similarity (Valid Reconstructions)', fontsize=13)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(results_dir / 'tanimoto_distribution.png', dpi=150, bbox_inches='tight')
print(f"  ✓ Saved: tanimoto_distribution.png")
plt.close()

# 2. Δ Atoms distribution
print("  Creating Δ Atoms distribution...")
d_atoms_valid = pd.to_numeric(df.loc[vm, 'delta_atoms'], errors='coerce').dropna()
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(d_atoms_valid, bins=30, edgecolor='black', alpha=0.7, color='coral')
ax.set_xlabel('Δ Atoms', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('Distribution of Atom Count Differences', fontsize=13)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(results_dir / 'atom_diff_distribution.png', dpi=150, bbox_inches='tight')
print(f"  ✓ Saved: atom_diff_distribution.png")
plt.close()

# 3. Validity Rate pie chart
print("  Creating validity rate chart...")
fig, ax = plt.subplots(figsize=(8, 6))
labels_pie = ['Invalid', 'Valid']
ax.pie([1 - vm.mean(), vm.mean()], labels=labels_pie, autopct='%1.1f%%', 
       colors=['#ff9999', '#66b3ff'], textprops={'fontsize': 12})
ax.set_title(f'Validity Rate: {vm.mean()*100:.1f}%', fontsize=13)
plt.tight_layout()
plt.savefig(results_dir / 'validity_rate.png', dpi=150, bbox_inches='tight')
print(f"  ✓ Saved: validity_rate.png")
plt.close()

# 4. Exact Match rate
print("  Creating exact match rate chart...")
if vm.sum() > 0:
    exact_valid = (df.loc[vm, 'exact_match'] == 1).sum()
    not_exact = vm.sum() - exact_valid
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(['Exact Match', 'Different'], [exact_valid, not_exact], 
           color=['#90EE90', '#FFB6C6'], width=0.6)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'Exact Match Rate (Valid): {exact_valid/vm.sum()*100:.1f}%', fontsize=13)
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(results_dir / 'exact_match_rate.png', dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: exact_match_rate.png")
    plt.close()

# 5. Generate diverse molecular structure comparisons
print(f"\n=== Generating 10 Diverse Molecular Structure Comparisons ===")
print("  Distribution: 2 Perfect Match, 2 High, 2 Medium, 4 Low similarity")

# Select diverse examples by Tanimoto score ranges
valid_df = df[df['valid_out'] == 1].copy()
valid_df['tanimoto_numeric'] = pd.to_numeric(valid_df['tanimoto'], errors='coerce')

selected_examples = []

# 1. Perfect matches (Tanimoto = 1.0 or exact_match = 1)
perfect_mask = (valid_df['tanimoto_numeric'] == 1.0) | (valid_df['exact_match'] == 1)
perfect_examples = valid_df[perfect_mask].nlargest(2, 'tanimoto_numeric')
selected_examples.append(perfect_examples)
print(f"  - Selected {len(perfect_examples)} perfect match examples")

# 2. High similarity (0.7 < Tanimoto < 1.0)
high_mask = (valid_df['tanimoto_numeric'] > 0.7) & (valid_df['tanimoto_numeric'] < 1.0)
high_examples = valid_df[high_mask].nlargest(2, 'tanimoto_numeric')
selected_examples.append(high_examples)
print(f"  - Selected {len(high_examples)} high similarity examples")

# 3. Medium similarity (0.4 <= Tanimoto <= 0.7)
medium_mask = (valid_df['tanimoto_numeric'] >= 0.4) & (valid_df['tanimoto_numeric'] <= 0.7)
medium_examples = valid_df[medium_mask].sample(n=min(2, medium_mask.sum()), random_state=42)
selected_examples.append(medium_examples)
print(f"  - Selected {len(medium_examples)} medium similarity examples")

# 4. Low similarity (Tanimoto < 0.4)
low_mask = valid_df['tanimoto_numeric'] < 0.4
low_examples = valid_df[low_mask].sample(n=min(4, low_mask.sum()), random_state=42)
selected_examples.append(low_examples)
print(f"  - Selected {len(low_examples)} low similarity examples")

# Combine all examples
top_examples = pd.concat(selected_examples, ignore_index=False)

examples_list = []
for idx, row in top_examples.iterrows():
    smiles_in = row['smiles_in']
    smiles_out = row['smiles_out']
    tani = row['tanimoto']
    exact = row['exact_match']
    
    print(f"  Processing example {idx} (Tanimoto: {tani:.4f})...")
    
    mol_in = Chem.MolFromSmiles(smiles_in)
    mol_out = Chem.MolFromSmiles(smiles_out)
    
    # Generate interpretation
    if tani == 1.0:
        interp = f"Perfect fingerprint match ({tani:.3f}): Structurally identical or only stereochemistry differs."
    elif exact == 1:
        interp = "Perfect: Exact reconstruction (canonical SMILES match)."
    elif tani > 0.7:
        interp = f"High similarity ({tani:.3f}): Minor structural rearrangements."
    elif tani > 0.5:
        interp = f"Moderate similarity ({tani:.3f}): Significant rearrangements present."
    else:
        interp = f"Low similarity ({tani:.3f}): Major reconstruction errors."
    
    examples_list.append({
        'ID': idx,
        'Input_SMILES': smiles_in,
        'Output_SMILES': smiles_out,
        'Tanimoto': f"{tani:.4f}",
        'Interpretation': interp
    })
    
    if mol_in and mol_out:
        try:
            # Compute 2D coordinates
            AllChem.Compute2DCoords(mol_in)
            AllChem.Compute2DCoords(mol_out)
            
            # Calculate metrics
            tani_val = float(tani)
            atoms_in = mol_in.GetNumAtoms()
            atoms_out = mol_out.GetNumAtoms()
            bonds_in = mol_in.GetNumBonds()
            bonds_out = mol_out.GetNumBonds()
            delta_atoms = atoms_out - atoms_in
            delta_bonds = bonds_out - bonds_in
            
            # Find maximum common substructure to highlight differences
            from rdkit.Chem import rdFMCS
            
            # Find MCS (Maximum Common Substructure)
            mcs_result = rdFMCS.FindMCS([mol_in, mol_out], 
                                       timeout=2, 
                                       completeRingsOnly=True,
                                       matchValences=True,
                                       bondCompare=rdFMCS.BondCompare.CompareAny,
                                       atomCompare=rdFMCS.AtomCompare.CompareAny)
            
            # Get MCS as a molecule
            if mcs_result.numAtoms > 0:
                mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
                
                # Find matching atoms/bonds in input molecule
                match_in = mol_in.GetSubstructMatch(mcs_mol)
                match_out = mol_out.GetSubstructMatch(mcs_mol)
                
                # Highlight atoms NOT in the MCS (differences)
                highlight_atoms_in = [i for i in range(atoms_in) if i not in match_in]
                highlight_atoms_out = [i for i in range(atoms_out) if i not in match_out]
                
                # Highlight bonds NOT in the MCS
                highlight_bonds_in = []
                for bond in mol_in.GetBonds():
                    if bond.GetBeginAtomIdx() not in match_in or bond.GetEndAtomIdx() not in match_in:
                        highlight_bonds_in.append(bond.GetIdx())
                
                highlight_bonds_out = []
                for bond in mol_out.GetBonds():
                    if bond.GetBeginAtomIdx() not in match_out or bond.GetEndAtomIdx() not in match_out:
                        highlight_bonds_out.append(bond.GetIdx())
            else:
                # No common substructure, highlight everything
                highlight_atoms_in = list(range(atoms_in))
                highlight_atoms_out = list(range(atoms_out))
                highlight_bonds_in = [bond.GetIdx() for bond in mol_in.GetBonds()]
                highlight_bonds_out = [bond.GetIdx() for bond in mol_out.GetBonds()]
            
            # Generate temporary images with highlighted differences
            temp_in = results_dir / f'temp_in_{idx}.png'
            temp_out = results_dir / f'temp_out_{idx}.png'
            
            # Red color for differences
            highlight_color = (1, 0, 0)  # RGB red
            
            if MOLTOFILE_AVAILABLE:
                MolToFile(mol_in, str(temp_in), size=(600, 600), kekulize=True, wedgeBonds=True,
                         highlightAtoms=highlight_atoms_in, highlightBonds=highlight_bonds_in,
                         highlightColor=highlight_color)
                MolToFile(mol_out, str(temp_out), size=(600, 600), kekulize=True, wedgeBonds=True,
                         highlightAtoms=highlight_atoms_out, highlightBonds=highlight_bonds_out,
                         highlightColor=highlight_color)
            else:
                # Fallback: use Draw.MolToImage
                img_in_pil = Draw.MolToImage(mol_in, size=(600, 600), kekulize=True, wedgeBonds=True,
                                            highlightAtoms=highlight_atoms_in, 
                                            highlightBonds=highlight_bonds_in)
                img_in_pil.save(temp_in)
                img_out_pil = Draw.MolToImage(mol_out, size=(600, 600), kekulize=True, wedgeBonds=True,
                                             highlightAtoms=highlight_atoms_out,
                                             highlightBonds=highlight_bonds_out)
                img_out_pil.save(temp_out)
            
            # Load images
            img_in = Image.open(temp_in)
            img_out = Image.open(temp_out)
            
            # Create figure with structure and legend
            fig = plt.figure(figsize=(18, 10))
            
            # Create grid: 2 rows, 3 columns (structures span 2 cols each, legend in middle)
            gs = fig.add_gridspec(2, 3, height_ratios=[3, 1], width_ratios=[1, 0.3, 1], hspace=0.3, wspace=0.1)
            
            ax_in = fig.add_subplot(gs[0, 0])
            ax_legend = fig.add_subplot(gs[0, 1])
            ax_out = fig.add_subplot(gs[0, 2])
            ax_smiles_in = fig.add_subplot(gs[1, 0])
            ax_smiles_out = fig.add_subplot(gs[1, 2])
            
            # Display molecular structures
            ax_in.imshow(img_in)
            ax_in.set_title('Input Molecule', fontsize=14, fontweight='bold', pad=10)
            ax_in.axis('off')
            
            ax_out.imshow(img_out)
            ax_out.set_title('Output Molecule', fontsize=14, fontweight='bold', pad=10)
            ax_out.axis('off')
            
            # Display SMILES strings
            ax_smiles_in.axis('off')
            smiles_in_wrapped = '\n'.join([smiles_in[i:i+45] for i in range(0, len(smiles_in), 45)])
            ax_smiles_in.text(0.5, 0.5, f'SMILES:\n{smiles_in_wrapped}', 
                             ha='center', va='center', fontsize=8, family='monospace',
                             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7, pad=0.5))
            
            ax_smiles_out.axis('off')
            smiles_out_wrapped = '\n'.join([smiles_out[i:i+45] for i in range(0, len(smiles_out), 45)])
            ax_smiles_out.text(0.5, 0.5, f'SMILES:\n{smiles_out_wrapped}', 
                              ha='center', va='center', fontsize=8, family='monospace',
                              bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7, pad=0.5))
            
            # Create legend with metrics
            ax_legend.axis('off')
            legend_text = f"""COMPARISON METRICS

Tanimoto Similarity:
{tani_val:.4f}

Input Molecule:
  Atoms: {atoms_in}
  Bonds: {bonds_in}

Output Molecule:
  Atoms: {atoms_out}
  Bonds: {bonds_out}

Delta:
  Δ Atoms: {delta_atoms:+d}
  Δ Bonds: {delta_bonds:+d}

Red = Differences
(atoms/bonds not in
common substructure)
"""
            
            ax_legend.text(0.5, 0.5, legend_text, ha='center', va='center', 
                          fontsize=10, family='monospace', fontweight='bold',
                          bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, pad=1))
            
            # Add main title
            fig.suptitle(f'Example {idx}: Reconstruction Comparison', 
                        fontsize=16, fontweight='bold', y=0.98)
            
            save_path = results_dir / f'example_{idx:03d}_comparison.png'
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  ✓ Saved: example_{idx:03d}_comparison.png")
            plt.close()
            
            # Clean up temp files
            temp_in.unlink()
            temp_out.unlink()
            
        except Exception as e:
            print(f"  ✗ Warning: Could not generate image for example {idx}: {e}")

# Save qualitative examples table
examples_df = pd.DataFrame(examples_list)
examples_df.to_csv(results_dir / 'qualitative_examples.csv', index=False)
print(f"\n✓ Saved qualitative examples table: qualitative_examples.csv")

print(f"\n✓ Visualization complete! All outputs saved to {results_dir}")
print(f"  Generated files:")
print(f"  - tanimoto_distribution.png")
print(f"  - atom_diff_distribution.png")
print(f"  - validity_rate.png")
print(f"  - exact_match_rate.png")
print(f"  - qualitative_examples.csv")
print(f"  - example_XXX_comparison.png ({len(examples_list)} files)")
