import multiprocessing
import os
import time
from itertools import combinations
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from Comparison_digital_profiles_human import compare_pair
from Signal_digitalisation_human import MatricesExtractor
from utils import InputFileManager
from utils import create_dir_if_not_exist
import sys
import numpy as np
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

# from multiprocessing import set_start_method
# set_start_method("spawn")


num_cores = multiprocessing.cpu_count()
num_task = num_cores - 1

# I/O directories
input_dir = os.path.join(os.getcwd(), "check_reduced/")  # get the path to the data input directory
match_scores_output_dir = os.path.join(os.getcwd(), "matrix_python/match_scores/")  # Sets the directory where all the saved outputs will be stored
reproducible_sequence_output_dir = os.path.join(os.getcwd(), "matrix_python/reproducible_sequence/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "gene_lengths.csv")  # path to upload the file containing each gene's ID and the correspondent gene length
histogram_plot_path = os.path.join(os.getcwd(), "genes_histograms/")  # path to upload the file containing each gene's ID and the correspondent gene length
intermediate_results = os.path.join(os.getcwd(), "intermediate_results/")

create_dir_if_not_exist([input_dir, match_scores_output_dir, histogram_plot_path, reproducible_sequence_output_dir,intermediate_results])

num_comparison = 1  # NOTA: numero di confronti random da eseguire per ogni coppia di file bed


def signal_digitalisation(genes, bed_files_dicts, areReadsRandomized):
    matrix_01_list = []
    for bed_files_dict in bed_files_dicts:
        bed_file = bed_files_dict["bed_file"]
        bed_file_name = bed_files_dict["bed_file_name"]
        me = MatricesExtractor(bed_file, genes)
        # extract the matrices
        pd_matrix_coverage, matrix_01 = me.extract_matrices(areReadsRandomized=areReadsRandomized)
        matrix_01_list.append({'matrix': matrix_01, 'file_name': bed_file_name})

        # NOTA: per velocizzare le operazioni per ora ho tolto il salvataggio
        # coverage_matrix_csv_path = os.path.join(match_scores_output_dir, bed_file_name + "_matrix_coverage.csv")
        # matrix_01_csv_path = os.path.join(match_scores_output_dir, bed_file_name + "_matrix_01.csv")

        # # Exports the dataFrames into CSV files
        # matrix_01.to_csv(matrix_01_csv_path, index=True)
        # pd_matrix_coverage.to_csv(coverage_matrix_csv_path, index=True)
    return matrix_01_list


# FIX ME al momento esegue la parte modificata comparison_digital_profiles_human, verificare che sia corretta. Per e.coli provare la comparison_digital_profiles.
def compute_real_match_scores(genes, bed_files_dicts, save_results=True):
    print("start real matrix digitalization...")
    matrix_01_list = signal_digitalisation(genes, bed_files_dicts, areReadsRandomized=False)
    print("digitalization complete...")

    # gets the gene list for each matrix
    gene_lists = [pd.DataFrame(f['matrix'].index, columns={"GeneID"}) for f in matrix_01_list]
    # gets the genes in common for each matrix
    gene_list = gene_lists[0]  # FIX ME sarebbe possibile fin dall'inizio selezionare solo i geni a comune
    if len(gene_lists) > 1:
        for i in range(1, len(gene_lists)):
            l = gene_lists[i]
            gene_list = gene_list.merge(l, on="GeneID")

    # gene_list['GeneID'] = gene_list.index
    gene_list = gene_list.to_numpy().squeeze()

    # lists of pair of matrices
    pairs = [list(f) for f in combinations(matrix_01_list, 2)]

    print("start pair comparison...")
    # saves match scores
    match_scores = []
    pair_names_list = []
    for pair in pairs:
        # print ("compare " + pair[0]['file_name'] + " and " + pair[1]['file_name'])
        match_score, pair_names = compare_pair(pair, genes.set_index('GeneID'), gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
        match_scores.append({"match_score": match_score, "pair_name": pair_names})
        pair_names_list.append(pair_names)
        if save_results:
            pair_names = pair_names + ".csv"
            match_score.to_csv(os.path.join(match_scores_output_dir, pair_names), index=True, header=True, decimal='.', sep=',', float_format='%.6f')
    print("real comparison complete")

    return gene_list, match_scores, pair_names_list, matrix_01_list


def random_comparison(arguments):
    start = time.time()
    matrices_extractors, genes, gene_list = arguments
    matrix_01_pair = []
    for matrices_extractor in matrices_extractors:
        bed_file_name = matrices_extractor["bed_file_name"]
        me = matrices_extractor["me"]
        # extract the matrices
        pd_matrix_coverage, matrix_01 = me.extract_matrices(areReadsRandomized=True)
        matrix_01_pair.append({'matrix': matrix_01, 'file_name': bed_file_name})

    match_score, pair_names = compare_pair(matrix_01_pair, genes.set_index('GeneID'), gene_list)
    pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
    end = time.time()
    print("end_comparison in " + str(end - start) + "seconds")
    return {'pair_name': pair_names, 'match_score': match_score}
    # match_scores.append({'pair_name': pair_names, 'match_score': match_score})


def compare_pair_n_times(bed_files_pair, genes, gene_list, n):
    # extract a pair of bed files
    match_scores = []

    matrices_extractors = []
    for bed_files_dict in bed_files_pair:  # FIX ME creazione delle classi estrattori a monte, verificare correttezza
        bed_file = bed_files_dict["bed_file"]
        bed_file_name = bed_files_dict["bed_file_name"]
        me = MatricesExtractor(bed_file, genes)
        matrices_extractors.append({"me": me, "bed_file_name": bed_file_name})

    arguments = matrices_extractors, genes, gene_list
    pool = multiprocessing.Pool(processes=num_task)
    res = []
    for i in range(n):
        res.append(pool.apply_async(random_comparison, [arguments]))
    pool.close()
    pool.join()
    for i in res:
        match_scores.append(i.get())

    return match_scores


def main(num):
    print("num of core available: " + str(num_cores) + " used: " + str(num_task))

    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes()
    bed_files_dicts = ifm.get_bed_files()
    gene_list, _, _, _ = compute_real_match_scores(genes, bed_files_dicts, save_results=False)  # FIX ME è totalmente inutile fare tutti i calcoli per calcolare il gene_list... andrebbe messa una funzione a parte

    # create pairs of bed files
    bed_files_pairs = [list(f) for f in combinations(bed_files_dicts, 2)]
    print("start fake matrix comparisons...")
    start = time.time()
    match_scores_list = []

    for bed_files_pair in bed_files_pairs:
        match_scores_list.append(compare_pair_n_times(bed_files_pair, genes, gene_list, num_comparison))
    end = time.time()
    print('fake matrix comparisons completed in ' + str(end - start) + " sec(s)")

    with open(intermediate_results+'/match_scores_list_' + str(num) + '_aggregating_' + str(num_comparison) + 'comparisons.npy', 'wb') as f:
        np.save(f, match_scores_list)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='riboseq reproducibility')
    parser.add_argument('--num', metavar='n', required=False, help='experiment identifier')
    args = parser.parse_args()
    main(num=args.num)
