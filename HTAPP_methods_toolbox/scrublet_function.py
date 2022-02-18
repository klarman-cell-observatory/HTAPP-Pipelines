"""
Date: 2019/10/04
Author: Orr Ashenberg

This script has a function to run scrublet. It is meant to be called from within R.
"""

import matplotlib
#matplotlib.use('TkAgg')  # backend may need to be changed to this if getting backend errors
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
plt.close("all")
import numpy as np
import os

def run_scrublet(sampleid, sparse_dir, figures_dir):
    """Run doublet detection program scrublet on a sample and its corresponding counts matrix.

    Args:
        *sampleid* (char): Sampleid name.
        *sparse_dir* (char): Directory containing 10x count matrix stored in MM format. The directory corresponds
        to the sample in *sampeid*.
        *figures_dir* (char): Directory to write figures and doublet score text file.

    Returns:
        bool: The return value. True for success, False otherwise.
    """

    print(sampleid)

    # Read the count matrix data and store as cells x genes in CSC format.
    counts_matrix = scipy.io.mmread(sparse_dir + '/matrix.mtx').T.tocsc()
    genes = np.array(scr.load_genes(sparse_dir + '/genes.tsv', delimiter='\t', column=1))
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))

    # Scrublet filtering.
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,  min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)

    # Change threshold if the doublet threshold score was set so high that no cells are called as doublets.
    threshold = max(doublet_scores[~predicted_doublets])  # maximum doublet score for non-doublet, ie threshold used
    if threshold > 0.3:
        predicted_doublets = scrub.call_doublets(threshold=0.3)
        print("changed doublet threshold to 0.3 for %s" % sampleid)

    # Plots of scores histogram and of scores mapped onto UMAP embedding.
#     scrub.plot_histogram()
#     plt.suptitle('%s' % sampleid)
#     plt.savefig('%s/%s_scrublet_histogram.png' % (figures_dir, sampleid))
#     plt.show()
#     scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
#     scrub.plot_embedding('UMAP', order_points=True)
#     plt.suptitle('%s' % sampleid)
#     plt.savefig('%s/%s_scrublet_embedding.png' % (figures_dir, sampleid))
#     plt.show()

    # Write scrublet doublet results to file: predicted doublet state and doublet score.
    scrubletfile = "%s/%s_scrublet_scores.txt" % (figures_dir, sampleid)
    with open(scrubletfile, "w") as f:
        f.write("doublet\tdoubletscore\n")
        z = list(zip(predicted_doublets, doublet_scores))
        for doublet, score in z:
            f.write("%s\t%s\n" % (doublet, score))
