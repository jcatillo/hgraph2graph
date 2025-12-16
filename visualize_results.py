import argparse
import os
from pathlib import Path
import pandas as pd
import numpy as np # Import numpy for category comparisons

# --- FIX START ---
# Delete the conflicting 'MPLBACKEND' environment variable set by interactive environments (like Colab),
# which causes a ValueError when trying to import matplotlib before the script can set 'Agg'.
if 'MPLBACKEND' in os.environ:
    del os.environ['MPLBACKEND']
    
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# --- FIX END ---

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
# Removed --num_examples argument as it's no longer used for sampling
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

# --- MODIFIED SECTION: Generate ALL molecular structure comparisons into folders ---
print(f"\n=== Generating ALL Molecular Structure Comparisons by Category ===")

valid_df = df[df['valid_out'] == 1].copy()
valid_df['tanimoto_numeric'] = pd.to_numeric(valid_df['tanimoto'], errors='coerce')

# Define categories and their criteria
CATEGORIES = {
    'Perfect_Match_1_0': (lambda x: np.isclose(x, 1.0)), # Check for Tanimoto exactly 1.0
    'High_Similarity_0_7_to_1_0': (lambda x: (x > 0.7) & (x < 1.0)),
    'Medium_Similarity_0_4_to_0_7': (lambda x: (x >= 0.4) & (x <= 0.7)),
    'Low_Similarity_Below_0_4': (lambda x: (x < 0.4))
}

# Add a category column to the DataFrame
def assign_category(tani):
    if np.isclose(tani, 1.0):
        return 'Perfect_Match_1_0'
    elif tani > 0.7:
        return 'High_Similarity_0_7_to_1_0'
    elif tani >= 0.4:
        return 'Medium_Similarity_0_4_to_0_7'
    else:
        return 'Low_Similarity_Below_0_4'

valid_df['category'] = valid_df['tanimoto_numeric'].apply(assign_category)

examples_list = []
total_generated = 0

# Create category subdirectories and iterate through all valid molecules
for category_name, category_df in valid_df.groupby('category'):
    category_dir = results_dir / 'comparisons' / category_name
    category_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n  Processing category: {category_name} ({len(category_df)} examples)")
    
    for idx, row in category_df.iterrows():
        smiles_in = row['smiles_in']
        smiles_out = row['smiles_out']
        tani = row['tanimoto_numeric']
        exact = row['exact_match']
        
        mol_in = Chem.MolFromSmiles(smiles_in)
        mol_out = Chem.MolFromSmiles(smiles_out)
        
        # Generate interpretation
        if np.isclose(tani, 1.0):
            interp = f"Perfect fingerprint match ({tani:.3f}): Structurally identical or only stereochemistry/hashing differs."
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
            'Category': category_name,
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
                atoms_in = mol_in.GetNumAtoms()
                atoms_out = mol_out.GetNumAtoms()
                bonds_in = mol_in.GetNumBonds()
                bonds_out = mol_out.GetNumBonds()
                delta_atoms = atoms_out - atoms_in
                delta_bonds = bonds_out - bonds_in
                
                # Find maximum common substructure to highlight differences
                from rdkit.Chem import rdFMCS
                
                # Find MCS (Maximum Common Substructure)
                # Setting timeout to 1s to prevent stalls on very large/complex structures
                mcs_result = rdFMCS.FindMCS([mol_in, mol_out], 
                                            timeout=1, 
                                            completeRingsOnly=True,
                                            matchValences=True,
                                            bondCompare=rdFMCS.BondCompare.CompareAny,
                                            atomCompare=rdFMCS.AtomCompare.CompareAny)
                
                highlight_atoms_in, highlight_atoms_out = [], []
                highlight_bonds_in, highlight_bonds_out = [], []
                
                if mcs_result.numAtoms > 0 and mcs_result.numAtoms < mol_in.GetNumAtoms():
                    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
                    
                    # Find matching atoms/bonds
                    match_in = mol_in.GetSubstructMatch(mcs_mol)
                    match_out = mol_out.GetSubstructMatch(mcs_mol)
                    
                    # Highlight atoms NOT in the MCS (differences)
                    highlight_atoms_in = [i for i in range(atoms_in) if i not in match_in]
                    highlight_atoms_out = [i for i in range(atoms_out) if i not in match_out]
                    
                    # Highlight bonds NOT in the MCS
                    highlight_bonds_in = [bond.GetIdx() for bond in mol_in.GetBonds() 
                                          if bond.GetBeginAtomIdx() not in match_in or bond.GetEndAtomIdx() not in match_in]
                    highlight_bonds_out = [bond.GetIdx() for bond in mol_out.GetBonds()
                                           if bond.GetBeginAtomIdx() not in match_out or bond.GetEndAtomIdx() not in match_out]
                else:
                    # If MCS is empty or covers the whole molecule (exact match), highlight nothing.
                    pass 

                # Generate temporary images with highlighted differences
                temp_in = category_dir / f'temp_in_{idx}.png'
                temp_out = category_dir / f'temp_out_{idx}.png'
                
                # Red color for differences
                highlight_color = (1, 0, 0) # RGB red
                
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
                                                 highlightBonds=highlight_bonds_in, highlightColor=highlight_color)
                    img_in_pil.save(temp_in)
                    img_out_pil = Draw.MolToImage(mol_out, size=(600, 600), kekulize=True, wedgeBonds=True,
                                                  highlightAtoms=highlight_atoms_out,
                                                  highlightBonds=highlight_bonds_out, highlightColor=highlight_color)
                    img_out_pil.save(temp_out)
                
                # Load images
                img_in = Image.open(temp_in)
                img_out = Image.open(temp_out)
                
                # Create figure with structure and legend
                fig = plt.figure(figsize=(18, 10))
                
                # Create grid
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
{tani:.4f}

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
                fig.suptitle(f'Example {idx}: Reconstruction Comparison (Category: {category_name.replace("_", " ")})', 
                             fontsize=16, fontweight='bold', y=0.98)
                
                save_path = category_dir / f'example_{idx:03d}_comparison.png'
                plt.savefig(save_path, dpi=150, bbox_inches='tight')
                plt.close(fig) # Close specific figure
                
                # Clean up temp files
                temp_in.unlink()
                temp_out.unlink()
                total_generated += 1
                
            except Exception as e:
                print(f"  ✗ Warning: Could not generate image for example {idx} in {category_name}: {e}")

print(f"\n  ✓ Successfully generated and categorized {total_generated} comparison images.")

# Save qualitative examples table
examples_df = pd.DataFrame(examples_list)
examples_df.to_csv(results_dir / 'qualitative_examples_all.csv', index=False)
print(f"\n✓ Saved qualitative examples table: qualitative_examples_all.csv")

print(f"\n✓ Visualization complete! All outputs saved to {results_dir}")
print(f"  - Comparison images are organized in subfolders under: {results_dir / 'comparisons'}")