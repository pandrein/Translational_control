# This code compares the matrices 01 pairwise (row by row) to quantify the similarities between Ribo-seq profiles coming from different datasets and referring to the same ORF.
# We assign a similarity score to each comparison: The score ranging from 0.5 (random matches) to 1 (high degree of similarity).

# Uploads the needed modules --------------------------------

import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


# COMPARING DIGITAL PROFILES ##=========================================================================================
# Compares the matrices 01 pairwise, row by row. Assigns a similarity score to each comparison.

def compare_pair(pair, genes_lengths, gene_list):
    pair_names = [pair[0]['file_name'], pair[1]['file_name']]
    pair_matrix = [pair[0]['matrix'], pair[1]['matrix']]
    # selects for each matrix only the genes in common
    pair_matrix_cleaned = [f.loc[gene_list, :] for f in pair_matrix]
    # binarize (True for each element equal to 1) and convert to numpy
    pair_matrix_with_same_col_dim = [f.to_numpy()[:, ~np.isnan(f).all(0)] for f in pair_matrix_cleaned]
    pair_matrix_binarized = [f == 1 for f in pair_matrix_with_same_col_dim]
    # pair_matrix_binarized = [(f == 1).to_numpy() for f in pair_matrix_cleaned]
    # stacks the pair in a multidimensional array
    pair_matrix_binarized_stack = np.stack(pair_matrix_binarized)
    # generates matrix_comparison and match_scores
    pair_comparison = np.all(pair_matrix_binarized_stack, axis=0)
    match_score = np.sum(pair_comparison, axis=1)
    # selects only the lengths of genes in common
    genes_lengths = genes_lengths.loc[gene_list, :]["GeneLength"]
    # divides for the gene length
    match_score = 2 * match_score / genes_lengths  # calculates the relative number of matches --> the ratio between the number of matches and the length of the ORF
    return match_score, pair_names


# SETUP =================================================================================================================

matrix_data_comparison_dir = "matrix_python/check_reduced/"
output_data_dir = "match_scores_python/"


# DATA UPLOAD ===========================================================================================================

def main():
    # finds the common genes
    working_dir = os.getcwd()
    matrix_comparison_dir = os.path.join(working_dir, matrix_data_comparison_dir)
    output_dir = os.path.join(working_dir, output_data_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(matrix_comparison_dir):
        os.makedirs(matrix_comparison_dir)
    genes_lengths_path = os.path.join(working_dir, "CDShumanGenesLengths.txt")

    genes_lengths = pd.read_csv(genes_lengths_path, sep="\t", header=None)
    genes_lengths = genes_lengths.rename(columns={genes_lengths.columns[0]: "GeneID", genes_lengths.columns[1]: "GeneLength"}).set_index("GeneID")

    # lists the files in the input directory
    # input_matrix_type = tuple([".csv"])
    input_matrix_type = tuple([".parquet"])
    matrix_file_list = [os.path.abspath(os.path.join(matrix_comparison_dir, f)) for f in os.listdir(matrix_comparison_dir) if f.endswith(input_matrix_type) and "01" in os.path.basename(f)]

    matrices = [{'matrix': pd.read_parquet(f), 'file_name': f} for f in matrix_file_list]

    # gets the gene list for each matrix
    gene_lists = [pd.DataFrame(f['matrix'].index, columns={"GeneID"}) for f in matrices]
    # gets the genes in common for each matrix
    gene_list = gene_lists[0]
    if len(gene_lists) > 1:
        for i in range(1, len(gene_lists)):
            l = gene_lists[i]
            gene_list = gene_list.merge(l, on="GeneID")

    # gene_list['GeneID'] = gene_list.index
    gene_list = gene_list.to_numpy().squeeze()
    # gene_list = gene_list[:3]

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
