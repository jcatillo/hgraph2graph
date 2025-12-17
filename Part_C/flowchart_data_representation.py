"""
Generate a flowchart for the data representation pipeline
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

def create_flowchart():
    fig, ax = plt.subplots(1, 1, figsize=(14, 18))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 24)
    ax.axis('off')
    
    # Define colors
    color_input = '#E8F4F8'
    color_process = '#B3D9E6'
    color_data = '#7FC8DC'
    color_output = '#A8E6CF'
    color_vocab = '#FFD3B6'
    
    def add_box(ax, x, y, width, height, text, color, fontsize=10):
        """Add a rectangular box with text"""
        box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                             boxstyle="round,pad=0.1", 
                             edgecolor='black', facecolor=color, linewidth=2)
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center', fontsize=fontsize, 
                weight='bold', wrap=True)
    
    def add_arrow(ax, x1, y1, x2, y2, label=''):
        """Add arrow between boxes"""
        arrow = FancyArrowPatch((x1, y1), (x2, y2),
                               arrowstyle='->', mutation_scale=20, 
                               linewidth=2, color='black')
        ax.add_patch(arrow)
        if label:
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(mid_x + 0.3, mid_y, label, fontsize=9, 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Title
    ax.text(5, 23.5, 'Data Representation Pipeline', 
            ha='center', va='center', fontsize=16, weight='bold')
    ax.text(5, 23, 'From SMILES to Numerical Tensors', 
            ha='center', va='center', fontsize=12, style='italic', color='gray')
    
    # ===== STAGE 1: VOCABULARY =====
    y_pos = 21.5
    ax.text(5, y_pos + 0.5, 'STAGE 1: VOCABULARY BUILDING', 
            ha='center', fontsize=12, weight='bold', 
            bbox=dict(boxstyle='round', facecolor='#FFE6CC', alpha=0.7))
    
    add_box(ax, 5, y_pos, 3, 0.8, 'train.txt\n(raw SMILES)', color_input)
    add_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
    
    add_box(ax, 5, y_pos - 1.6, 3.5, 0.8, 'get_vocab.py\nparallel processing', color_process)
    add_arrow(ax, 5, y_pos - 2.0, 5, y_pos - 2.7)
    
    add_box(ax, 5, y_pos - 3.2, 3.5, 0.8, 'MolGraph(smiles)\nfor each molecule', color_process)
    add_arrow(ax, 5, y_pos - 3.6, 5, y_pos - 4.3)
    
    add_box(ax, 5, y_pos - 4.8, 3.5, 0.8, 'Extract:\n- Clusters (rings + bonds)\n- Attachments', color_data, fontsize=9)
    add_arrow(ax, 5, y_pos - 5.2, 5, y_pos - 5.9)
    
    add_box(ax, 5, y_pos - 6.4, 3.5, 0.8, 'vocab.txt\n(unique motifs)', color_vocab)
    
    # ===== STAGE 2: VOCAB OBJECT =====
    y_pos = 14.5
    ax.text(5, y_pos + 0.5, 'STAGE 2: VOCABULARY OBJECT', 
            ha='center', fontsize=12, weight='bold',
            bbox=dict(boxstyle='round', facecolor='#FFE6CC', alpha=0.7))
    
    add_box(ax, 5, y_pos, 3.5, 0.8, 'PairVocab(vocab.txt)', color_vocab)
    add_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
    
    # Split into three components
    add_box(ax, 2, y_pos - 1.8, 2.5, 1.2, 'hvocab\n(cluster list)\nC, CC, C1CC...', color_data, fontsize=9)
    add_box(ax, 5, y_pos - 1.8, 2.5, 1.2, 'vocab\n(pairs)\n(full, inter)\n(C, C*), ...', color_data, fontsize=9)
    add_box(ax, 8, y_pos - 1.8, 2.5, 1.2, 'mask\n(valid combos)\n0 or -1000', color_data, fontsize=9)
    
    # ===== STAGE 3: PREPROCESSING =====
    y_pos = 10
    ax.text(5, y_pos + 1.2, 'STAGE 3: PREPROCESSING', 
            ha='center', fontsize=12, weight='bold',
            bbox=dict(boxstyle='round', facecolor='#E6F3FF', alpha=0.7))
    
    add_box(ax, 5, y_pos, 3, 0.8, 'preprocess.py', color_process)
    add_arrow(ax, 5, y_pos - 0.4, 5, y_pos - 1.1)
    
    add_box(ax, 5, y_pos - 1.6, 3.5, 0.8, 'Load SMILES pairs\n(source, target)', color_data, fontsize=9)
    add_arrow(ax, 5, y_pos - 2.0, 5, y_pos - 2.7)
    
    add_box(ax, 5, y_pos - 3.2, 3.5, 0.8, 'Create batches\nparallel processing', color_process, fontsize=9)
    add_arrow(ax, 5, y_pos - 3.6, 5, y_pos - 4.3)
    
    add_box(ax, 5, y_pos - 4.8, 3, 0.8, 'For each batch:\ntensorize_pair()', color_process, fontsize=9)
    
    # ===== STAGE 4: MOLCRAPH DECOMPOSITION =====
    y_pos = 4.5
    ax.text(5, y_pos + 1.2, 'STAGE 4: MOLCRAPH DECOMPOSITION', 
            ha='center', fontsize=12, weight='bold',
            bbox=dict(boxstyle='round', facecolor='#E6FFE6', alpha=0.7))
    
    add_box(ax, 1.5, y_pos, 2.5, 0.8, 'SMILES\n↓\nbuild_mol_graph()', color_data, fontsize=9)
    add_box(ax, 1.5, y_pos - 1.2, 2.5, 0.8, 'Atom Layer\n(all atoms+bonds)', color_output, fontsize=9)
    
    add_box(ax, 4, y_pos, 2.5, 0.8, 'SMILES\n↓\nfind_clusters()', color_data, fontsize=9)
    add_box(ax, 4, y_pos - 1.2, 2.5, 0.8, 'Motif Layer\n(rings+bonds)', color_output, fontsize=9)
    
    add_box(ax, 6.5, y_pos, 2.5, 0.8, 'SMILES\n↓\ntree_decomp()', color_data, fontsize=9)
    add_box(ax, 6.5, y_pos - 1.2, 2.5, 0.8, 'Tree Structure\n(hierarchy)', color_output, fontsize=9)
    
    add_box(ax, 9, y_pos, 1.5, 0.8, 'SMILES\n↓\nlabel_tree()', color_data, fontsize=9)
    add_box(ax, 9, y_pos - 1.2, 1.5, 0.8, 'Orders\n(DFS)', color_output, fontsize=9)
    
    # Arrows down
    for x_pos in [1.5, 4, 6.5, 9]:
        add_arrow(ax, x_pos, y_pos - 0.4, x_pos, y_pos - 0.8)
    
    # ===== STAGE 5: TENSORIZATION =====
    y_pos = 1.5
    ax.text(5, y_pos + 0.8, 'STAGE 5: TENSORIZATION OUTPUT', 
            ha='center', fontsize=12, weight='bold',
            bbox=dict(boxstyle='round', facecolor='#F0E6FF', alpha=0.7))
    
    add_box(ax, 2.5, y_pos - 0.8, 2.5, 1.2, 'Tree Tensors\nfnode, fmess\nagraph, cgraph', color_output, fontsize=9)
    add_box(ax, 5, y_pos - 0.8, 2.5, 1.2, 'Graph Tensors\nfnode, fmess\nagraph, bgraph', color_output, fontsize=9)
    add_box(ax, 7.5, y_pos - 0.8, 2.5, 1.2, 'Orders\ngeneration\nsequences', color_output, fontsize=9)
    
    plt.tight_layout()
    plt.savefig('flowchart_data_representation.png', dpi=300, bbox_inches='tight')
    print("Flowchart saved: flowchart_data_representation.png")
    plt.close()

if __name__ == '__main__':
    create_flowchart()
