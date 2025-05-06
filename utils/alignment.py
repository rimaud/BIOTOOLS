from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import numpy as np

def plot_alignment(seq1, seq2):
    """Generate and display a sequence alignment"""
    alignments = pairwise2.align.globalxx(seq1, seq2)
    
    # Create alignment visualization
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis('off')
    
    # Display top alignment
    alignment_str = format_alignment(*alignments[0])
    ax.text(0.1, 0.5, alignment_str, 
           fontfamily='monospace',
           va='center',
           ha='left')
    
    return fig