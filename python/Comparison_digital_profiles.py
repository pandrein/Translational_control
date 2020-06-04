
# This code compares the matrices 01 pairwise (row by row) to quantify the similarities between Ribo-seq profiles coming from different datasets and referring to the same ORF.
# We assign a similarity score to each comparison: The score ranging from 0.5 (random matches) to 1 (high degree of similarity).

#Uploads the needed modules --------------------------------

import numpy as np
import pandas as pd
import pyranges as pr
import os
import sys
from tqdm import trange
import time
from joblib import Parallel, delayed
import multiprocessing
import random
import numpy.ma as ma
from natsort import natsorted, ns
from itertools import combinations
from pathlib import Path

# COMPARING DIGITAL PROFILES ##=========================================================================================
# Compares the matrices 01 pairwise, row by row. Assigns a similarity score to each comparison.

def compare_pair(pair, genes_lengths, gene_list):
    pair_names = [pair[0]['file_name'], pair[1]['file_name']]
    pair_matrix = [pair[0]['matrix'], pair[1]['matrix']]
    # uses the GeneID as index
    pair_matrix = [f.set_index("GeneID") for f in pair_matrix]
    # selects for each matrix only the genes in common
    pair_matrix_cleaned = [f.loc[gene_list, :] for f in pair_matrix]
    # binarize (True for each element equal to 1) and convert to numpy
    pair_matrix_binarized = [(f == 1).to_numpy() for f in pair_matrix_cleaned]
    # stacks the pair in a multidimensional array
    pair_matrix_binarized_stack = np.stack(pair_matrix_binarized)
    # generates matrix_comparison and match_scores
    pair_comparison = np.all(pair_matrix_binarized_stack, axis=0)
    match_score = np.sum(pair_comparison, axis=1)
    # selects only the lengths of genes in common
    genes_lengths = genes_lengths.loc[gene_list, :]["GeneLength"]
    # divides for the gene length
    match_score = 2 * match_score / genes_lengths   #calculates the relative number of matches --> the ratio between the number of matches and the length of the ORF
    return match_score, pair_names

#SETUP =================================================================================================================

matrix_data_comparison_dir = "output_matrix/"
output_data_dir = "Output_match_scores_ecoli/"

#DATA UPLOAD ===========================================================================================================

def main():
    # finds the common genes
    working_dir = os.getcwd()
    matrix_comparison_dir = os.path.join(working_dir, matrix_data_comparison_dir)
    output_dir = os.path.join(working_dir, output_data_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(matrix_comparison_dir):
        os.makedirs(matrix_comparison_dir)
    genes_lengths_path = os.path.join(working_dir, "gene_lengths.csv")

    genes_lengths = pd.read_csv(genes_lengths_path)
    genes_lengths = genes_lengths.drop(genes_lengths.columns[0], axis=1)
    genes_lengths = genes_lengths.rename(columns={genes_lengths.columns[0]: "GeneID", genes_lengths.columns[1]: "GeneLength"}).set_index("GeneID")

    # lists the files in the input directory
    input_matrix_type = tuple([".csv"])
    matrix_file_list = [os.path.abspath(os.path.join(matrix_comparison_dir, f)) for f in os.listdir(matrix_comparison_dir) if f.endswith(input_matrix_type) and "01" in os.path.basename(f)]

    # opens all the matrices and rename the first column (each matrix is associated with its name)
    matrices = [{'matrix': pd.read_csv(f, delimiter= ","), 'file_name': f} for f in matrix_file_list]
    matrices = [{'matrix': f['matrix'].rename(columns={f['matrix'].columns[0]: "GeneID"}), 'file_name': os.path.basename(f['file_name'])} for f in matrices]

    # gets the gene list for each matrix
    gene_lists = [pd.DataFrame(f['matrix']["GeneID"]).set_index("GeneID") for f in matrices]
    # gets the genes in common for each matrix
    gene_list = gene_lists[0].join(gene_lists[1:], how="inner")
    gene_list['GeneID'] = gene_list.index
    gene_list = gene_list.to_numpy().squeeze()

    # lists of pair of matrices
    pairs = [list(f) for f in combinations(matrices, 2)]

    # saves match scores

    match_scores = []
    for pair in pairs:
        match_score, pair_names = compare_pair(pair, genes_lengths, gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem + ".csv"
        match_scores.append(match_score)
        match_score.to_csv(os.path.join(output_dir, pair_names), index=True, header=True, decimal='.', sep=',', float_format='%.6f')


if __name__ == '__main__':
    main()
