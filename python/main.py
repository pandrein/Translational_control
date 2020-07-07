from utils import InputFileManager
from utils import create_dir_if_not_exist
from Signal_digitalisation_human import MatricesExtractor
from Comparison_digital_profiles_human import compare_pair
from itertools import combinations
import os
import pandas as pd
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

# I/O directories
input_dir = os.path.join(os.getcwd(), "check_reduced/")  # get the path to the data input directory
output_dir = os.path.join(os.getcwd(), "matrix_python/match_scores/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "CDShumanGenesLengths.txt")  # path to upload the file containing each gene's ID and the correspondent gene length

create_dir_if_not_exist([input_dir, output_dir])


def signal_digitalisation(genes, bed_files_dicts, areReadsRandomized):
    matrix_01_list = []
    for bed_files_dict in bed_files_dicts:
        bed_file = bed_files_dict["bed_file"]
        bed_file_name = bed_files_dict["bed_file_name"]
        me = MatricesExtractor(bed_file, genes)
        # extract the matrices
        pd_matrix_coverage, matrix_01 = me.extract_matrices(areReadsRandomized=areReadsRandomized)
        matrix_01_list.append({'matrix': matrix_01, 'file_name': bed_file_name})
        # coverage_matrix_csv_path = os.path.join(output_dir, bed_file_name + "_matrix_coverage.csv")
        # matrix_01_csv_path = os.path.join(output_dir, bed_file_name + "_matrix_01.csv")

        # # Exports the dataFrames into CSV files
        # matrix_01.to_csv(matrix_01_csv_path, index=True)
        # pd_matrix_coverage.to_csv(coverage_matrix_csv_path, index=True)
    return matrix_01_list


def compute_real_match_scores(genes, bed_files_dicts):
    matrix_01_list = signal_digitalisation(genes, bed_files_dicts, areReadsRandomized=False)

    # gets the gene list for each matrix
    gene_lists = [pd.DataFrame(f['matrix'].index, columns={"GeneID"}) for f in matrix_01_list]
    # gets the genes in common for each matrix
    gene_list = gene_lists[0]
    if len(gene_lists) > 1:
        for i in range(1, len(gene_lists)):
            l = gene_lists[i]
            gene_list = gene_list.merge(l, on="GeneID")

    # gene_list['GeneID'] = gene_list.index
    gene_list = gene_list.to_numpy().squeeze()

    # lists of pair of matrices
    pairs = [list(f) for f in combinations(matrix_01_list, 2)]

    # saves match scores
    match_scores = []
    for pair in pairs:
        match_score, pair_names = compare_pair(pair, genes.set_index('GeneID'), gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem + ".csv"
        match_scores.append(match_score)
        match_score.to_csv(os.path.join(output_dir, pair_names), index=True, header=True, decimal='.', sep=',', float_format='%.6f')

    return gene_list, match_scores


def compute_fake_match_scores(genes, bed_files_dicts, gene_list):
    matrix_01_list = signal_digitalisation(genes, bed_files_dicts, areReadsRandomized=True)

    # lists of pair of matrices
    pairs = [list(f) for f in combinations(matrix_01_list, 2)]

    # saves match scores
    match_scores = []
    for pair in pairs:
        match_score, pair_names = compare_pair(pair, genes.set_index('GeneID'), gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
        match_scores.append({'pair_name': pair_names, 'match_score': match_score})

    return match_scores


def main():
    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes_human()
    bed_files_dicts = ifm.get_bed_files()
    gene_list, match_scores = compute_real_match_scores(genes, bed_files_dicts)
    print("generated real match scores")
    # compute fake match scores
    fake_match_scores = compute_fake_match_scores(genes, bed_files_dicts, gene_list)

    match_scores_hist = {}
    for fake_match_score in fake_match_scores:
        match_scores_hist[fake_match_score['pair_name']] = fake_match_score['match_score'].to_list()

    for i in range(0, 10):
        fake_match_scores = compute_fake_match_scores(genes, bed_files_dicts, gene_list)
        for fake_match_score in fake_match_scores:
            match_scores_hist[fake_match_score['pair_name']] = match_scores_hist[fake_match_score['pair_name']] + fake_match_score['match_score'].to_list()

    for key in match_scores_hist:
        plt.hist(match_scores_hist[key], bins=100)
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        plt.title(key)
        plt.show()


if __name__ == '__main__':
    main()
