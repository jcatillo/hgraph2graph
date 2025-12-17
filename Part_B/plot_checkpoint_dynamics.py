# Part_B/plot_part_b.py - CHECKPOINT DYNAMICS VISUALIZATION
import argparse
import os
from pathlib import Path
import pandas as pd
import numpy as np

if 'MPLBACKEND' in os.environ:
    del os.environ['MPLBACKEND']
    
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

parser = argparse.ArgumentParser(description="Generate visualizations for Part B checkpoint dynamics")
parser.add_argument("--input", "-i", required=True, help="Path to checkpoint_dynamics.csv")
parser.add_argument("--output_dir", "-o", default=None, help="Output directory for visualizations")
args = parser.parse_args()

results_file = Path(args.input)

if not results_file.exists():
    print(f"Error: {results_file} not found!")
    exit(1)

if args.output_dir:
    results_dir = Path(args.output_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
else:
    results_dir = results_file.parent

print(f"Loading checkpoint dynamics from {results_file}...")
df = pd.read_csv(results_file)

print("\n=== Generating Per-Checkpoint Plots ===\n")

# Create plots for each checkpoint showing early/mid/late stages
for idx, row in df.iterrows():
    ckpt_num = int(row['checkpoint'])
    print(f"  Creating plot for Checkpoint {ckpt_num}...")
    
    # Data for this checkpoint
    stages = ['Early', 'Mid', 'Late']
    exact_rates = [
        row['early_exact_rate'] * 100,
        row['mid_exact_rate'] * 100,
        row['late_exact_rate'] * 100
    ]
    tanimoto_vals = [
        row['early_tanimoto'],
        row['mid_tanimoto'],
        row['late_tanimoto']
    ]
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Exact Match Rate by Stage
    colors_exact = ['#06A77D', '#FFB703', '#E76F51']
    bars1 = ax1.bar(stages, exact_rates, color=colors_exact, edgecolor='black', linewidth=2, alpha=0.8)
    ax1.set_ylabel('Exact Match Rate (%)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Reconstruction Stage', fontsize=12, fontweight='bold')
    ax1.set_title(f'Checkpoint {ckpt_num}: Exact Match by Stage', fontsize=13, fontweight='bold')
    ax1.set_ylim([0, 105])
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels on bars
    for bar, val in zip(bars1, exact_rates):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Plot 2: Tanimoto Similarity by Stage
    colors_tani = ['#06A77D', '#FFB703', '#E76F51']
    bars2 = ax2.bar(stages, tanimoto_vals, color=colors_tani, edgecolor='black', linewidth=2, alpha=0.8)
    ax2.set_ylabel('Mean Tanimoto Similarity', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Reconstruction Stage', fontsize=12, fontweight='bold')
    ax2.set_title(f'Checkpoint {ckpt_num}: Tanimoto by Stage', fontsize=13, fontweight='bold')
    ax2.set_ylim([0, 1.05])
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels on bars
    for bar, val in zip(bars2, tanimoto_vals):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.4f}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.suptitle(f'Checkpoint {ckpt_num} - Stage Breakdown', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(results_dir / f'checkpoint_{ckpt_num}_stages.png', dpi=150, bbox_inches='tight')
    print(f"     ✓ Saved: checkpoint_{ckpt_num}_stages.png")
    plt.close()

print(f"\n{'='*70}")
fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(df['checkpoint'], df['exact_match_rate'] * 100, 'o-', linewidth=3, 
        markersize=12, color='#2E86AB', label='Exact Match Rate', markerfacecolor='white', markeredgewidth=2)
ax.fill_between(df['checkpoint'], df['exact_match_rate'] * 100, alpha=0.2, color='#2E86AB')
ax.set_xlabel('Training Checkpoint (Step)', fontsize=13, fontweight='bold')
ax.set_ylabel('Exact Match Accuracy (%)', fontsize=13, fontweight='bold')
ax.set_title('Reconstruction Accuracy During Training', fontsize=14, fontweight='bold')
ax.grid(alpha=0.3, linestyle='--')
ax.legend(fontsize=12, loc='lower right')
ax.set_ylim([0, 105])
plt.tight_layout()
plt.savefig(results_dir / 'exact_match_vs_checkpoint.png', dpi=150, bbox_inches='tight')
print(f"     ✓ Saved: exact_match_vs_checkpoint.png")
plt.close()

# 2. Mean and Median Tanimoto vs Training Step
print("  2. Creating Tanimoto similarity plot...")
fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(df['checkpoint'], df['mean_tanimoto'], 'o-', linewidth=3, markersize=11, 
        color='#A23B72', label='Mean Tanimoto', markerfacecolor='white', markeredgewidth=2)
ax.plot(df['checkpoint'], df['median_tanimoto'], 's--', linewidth=2.5, markersize=9, 
        color='#F18F01', label='Median Tanimoto', markerfacecolor='white', markeredgewidth=2)
ax.fill_between(df['checkpoint'], df['mean_tanimoto'], alpha=0.15, color='#A23B72')
ax.set_xlabel('Training Checkpoint (Step)', fontsize=13, fontweight='bold')
ax.set_ylabel('Tanimoto Similarity', fontsize=13, fontweight='bold')
ax.set_title('Molecular Fingerprint Similarity During Training', fontsize=14, fontweight='bold')
ax.grid(alpha=0.3, linestyle='--')
ax.legend(fontsize=12, loc='lower right')
ax.set_ylim([0, 1.05])
plt.tight_layout()
plt.savefig(results_dir / 'tanimoto_vs_checkpoint.png', dpi=150, bbox_inches='tight')
print(f"     ✓ Saved: tanimoto_vs_checkpoint.png")
plt.close()

# 3. Stage-based progression (Early/Mid/Late)
print("  3. Creating stage-based progression plot...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Exact match by stage
ax1.plot(df['checkpoint'], df['early_exact_rate'] * 100, 'o-', linewidth=2.5, markersize=10,
         color='#06A77D', label='Early Stage', markerfacecolor='white', markeredgewidth=2)
ax1.plot(df['checkpoint'], df['mid_exact_rate'] * 100, 's-', linewidth=2.5, markersize=9,
         color='#FFB703', label='Mid Stage', markerfacecolor='white', markeredgewidth=2)
ax1.plot(df['checkpoint'], df['late_exact_rate'] * 100, '^-', linewidth=2.5, markersize=10,
         color='#E76F51', label='Late Stage', markerfacecolor='white', markeredgewidth=2)
ax1.set_ylabel('Exact Match Rate (%)', fontsize=12, fontweight='bold')
ax1.set_title('Exact Match Accuracy by Molecule Position (Early/Mid/Late)', fontsize=13, fontweight='bold')
ax1.grid(alpha=0.3, linestyle='--')
ax1.legend(fontsize=11, loc='lower right')
ax1.set_ylim([0, 105])

# Tanimoto by stage
ax2.plot(df['checkpoint'], df['early_tanimoto'], 'o-', linewidth=2.5, markersize=10,
         color='#06A77D', label='Early Stage', markerfacecolor='white', markeredgewidth=2)
ax2.plot(df['checkpoint'], df['mid_tanimoto'], 's-', linewidth=2.5, markersize=9,
         color='#FFB703', label='Mid Stage', markerfacecolor='white', markeredgewidth=2)
ax2.plot(df['checkpoint'], df['late_tanimoto'], '^-', linewidth=2.5, markersize=10,
         color='#E76F51', label='Late Stage', markerfacecolor='white', markeredgewidth=2)
ax2.set_xlabel('Training Checkpoint (Step)', fontsize=12, fontweight='bold')
ax2.set_ylabel('Mean Tanimoto Similarity', fontsize=12, fontweight='bold')
ax2.set_title('Fingerprint Similarity by Molecule Position', fontsize=13, fontweight='bold')
ax2.grid(alpha=0.3, linestyle='--')
ax2.legend(fontsize=11, loc='lower right')
ax2.set_ylim([0, 1.05])

plt.tight_layout()
plt.savefig(results_dir / 'stage_progression.png', dpi=150, bbox_inches='tight')
print(f"     ✓ Saved: stage_progression.png")
plt.close()

# 4. Combined metrics overview
print("  4. Creating combined metrics overview...")
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

# Exact match
ax1.plot(df['checkpoint'], df['exact_match_rate'] * 100, 'o-', linewidth=2.5, markersize=10, color='#2E86AB')
ax1.fill_between(df['checkpoint'], df['exact_match_rate'] * 100, alpha=0.2, color='#2E86AB')
ax1.set_ylabel('Rate (%)', fontsize=11, fontweight='bold')
ax1.set_title('Exact Match Accuracy', fontsize=12, fontweight='bold')
ax1.grid(alpha=0.3)
ax1.set_ylim([0, 105])

# Validity
ax2.plot(df['checkpoint'], df['validity_rate'] * 100, 's-', linewidth=2.5, markersize=10, color='#06A77D')
ax2.fill_between(df['checkpoint'], df['validity_rate'] * 100, alpha=0.2, color='#06A77D')
ax2.set_ylabel('Rate (%)', fontsize=11, fontweight='bold')
ax2.set_title('Validity Rate', fontsize=12, fontweight='bold')
ax2.grid(alpha=0.3)
ax2.set_ylim([0, 105])

# Mean Tanimoto
ax3.plot(df['checkpoint'], df['mean_tanimoto'], '^-', linewidth=2.5, markersize=10, color='#A23B72')
ax3.fill_between(df['checkpoint'], df['mean_tanimoto'], alpha=0.2, color='#A23B72')
ax3.set_xlabel('Training Checkpoint', fontsize=11, fontweight='bold')
ax3.set_ylabel('Similarity', fontsize=11, fontweight='bold')
ax3.set_title('Mean Tanimoto Similarity', fontsize=12, fontweight='bold')
ax3.grid(alpha=0.3)
ax3.set_ylim([0, 1.05])

# Median Tanimoto
ax4.plot(df['checkpoint'], df['median_tanimoto'], 'd-', linewidth=2.5, markersize=10, color='#F18F01')
ax4.fill_between(df['checkpoint'], df['median_tanimoto'], alpha=0.2, color='#F18F01')
ax4.set_xlabel('Training Checkpoint', fontsize=11, fontweight='bold')
ax4.set_ylabel('Similarity', fontsize=11, fontweight='bold')
ax4.set_title('Median Tanimoto Similarity', fontsize=12, fontweight='bold')
ax4.grid(alpha=0.3)
ax4.set_ylim([0, 1.05])

plt.suptitle('Part B: Checkpoint Dynamics - All Metrics', fontsize=14, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig(results_dir / 'combined_metrics.png', dpi=150, bbox_inches='tight')
print(f"     ✓ Saved: combined_metrics.png")
plt.close()

# 5. Learning curves comparison (normalized)
print("  5. Creating learning curves comparison...")
fig, ax = plt.subplots(figsize=(12, 7))

norm_exact = df['exact_match_rate'] * 100
norm_tanimoto = df['mean_tanimoto'] * 100
norm_validity = df['validity_rate'] * 100

ax.plot(df['checkpoint'], norm_exact, 'o-', linewidth=3, markersize=11, 
        color='#2E86AB', label='Exact Match Rate (%)', markerfacecolor='white', markeredgewidth=2)
ax.plot(df['checkpoint'], norm_tanimoto, 's-', linewidth=3, markersize=10, 
        color='#A23B72', label='Mean Tanimoto (scaled %)', markerfacecolor='white', markeredgewidth=2)
ax.plot(df['checkpoint'], norm_validity, '^-', linewidth=3, markersize=10, 
        color='#06A77D', label='Validity Rate (%)', markerfacecolor='white', markeredgewidth=2)

ax.set_xlabel('Training Checkpoint (Step)', fontsize=13, fontweight='bold')
ax.set_ylabel('Score (%)', fontsize=13, fontweight='bold')
ax.set_title('Training Progress: All Metrics Comparison (Normalized)', fontsize=14, fontweight='bold')
ax.grid(alpha=0.3, linestyle='--')
ax.legend(fontsize=12, loc='lower right')
ax.set_ylim([0, 105])
plt.tight_layout()
plt.savefig(results_dir / 'learning_curves.png', dpi=150, bbox_inches='tight')
print(f"     ✓ Saved: learning_curves.png")
plt.close()

# 6. Stage comparison heatmap
print("  6. Creating stage comparison heatmap...")
fig, ax = plt.subplots(figsize=(10, 6))

stages_data = pd.DataFrame({
    'Early Exact': df['early_exact_rate'] * 100,
    'Mid Exact': df['mid_exact_rate'] * 100,
    'Late Exact': df['late_exact_rate'] * 100,
})

im = ax.imshow(stages_data.T, cmap='YlGn', aspect='auto')

ax.set_xticks(range(len(df)))
ax.set_xticklabels([f"Ckpt {x}" for x in df['checkpoint']], rotation=45)
ax.set_yticks(range(3))
ax.set_yticklabels(['Early', 'Mid', 'Late'])
ax.set_xlabel('Training Checkpoint', fontsize=12, fontweight='bold')
ax.set_ylabel('Stage', fontsize=12, fontweight='bold')
ax.set_title('Exact Match Rate Heatmap: Progress by Stage', fontsize=13, fontweight='bold')

for i in range(len(stages_data.T)):
    for j in range(len(df)):
        text = ax.text(j, i, f'{stages_data.iloc[j, i]:.1f}%',
                      ha="center", va="center", color="black", fontsize=10, fontweight='bold')

cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Accuracy (%)', fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig(results_dir / 'stage_heatmap.png', dpi=150, bbox_inches='tight')
print(f"     ✓ Saved: stage_heatmap.png")
plt.close()

print(f"\n{'='*70}")
print("CHECKPOINT DYNAMICS VISUALIZATION SUMMARY")
print(f"{'='*70}\n")

print(f"📊 Per-Checkpoint Plots Generated: {len(df)} plots")
for _, row in df.iterrows():
    ckpt_num = int(row['checkpoint'])
    print(f"  ✓ checkpoint_{ckpt_num}_stages.png")

print(f"\n📈 Summary Plots:")
print(f"  1. exact_match_vs_checkpoint.png - Accuracy curve (all checkpoints)")
print(f"  2. tanimoto_vs_checkpoint.png    - Similarity curve (all checkpoints)")
print(f"  3. stage_progression.png         - Early/Mid/Late stage comparison")
print(f"  4. combined_metrics.png          - All 4 metrics in one view")
print(f"  5. learning_curves.png           - Normalized comparison")
print(f"  6. stage_heatmap.png             - Heatmap of stage progression")

print(f"\n📈 Key Metrics:")
print(f"  • Total checkpoints plotted: {len(df)}")
print(f"  • Exact Match improvement: {df['exact_match_rate'].iloc[0]*100:.2f}% → {df['exact_match_rate'].iloc[-1]*100:.2f}%")
print(f"  • Tanimoto improvement: {df['mean_tanimoto'].iloc[0]:.4f} → {df['mean_tanimoto'].iloc[-1]:.4f}")
print(f"  • Validity improvement: {df['validity_rate'].iloc[0]*100:.2f}% → {df['validity_rate'].iloc[-1]*100:.2f}%")

print(f"\n✓ All {6 + len(df)} visualizations saved to: {results_dir}")
print(f"{'='*70}")