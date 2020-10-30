import multiprocessing
import os
import time
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from Comparison_digital_profiles_human import compare_pair
from Signal_digitalisation_human import MatricesExtractor
from utils import InputFileManager
from utils import create_dir_if_not_exist

np.random.seed(10)

add_small_random_value_to_random_comparison = False

num_cores = multiprocessing.cpu_count()
num_task = num_cores - 1

# I/O directories
input_dir = os.path.join(os.getcwd(), "check_reduced/")  # get the path to the data input directory
match_scores_output_dir = os.path.join(os.getcwd(), "matrix_python/match_scores/")  # Sets the directory where all the saved outputs will be stored
reproducible_sequence_output_dir = os.path.join(os.getcwd(), "matrix_python/reproducible_sequence/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "gene_lengths.csv")  # path to upload the file containing each gene's ID and the correspondent gene length
histogram_plot_path = os.path.join(os.getcwd(), "genes_histograms/")  # path to upload the file containing each gene's ID and the correspondent gene length
intermediate_results = os.path.join(os.getcwd(), "intermediate_results/")
random_comparisons_folder_match_scores = os.path.join(os.getcwd(), "random_comparisons_folder_match_scores/")
random_comparisons_folder_matrix_01 = os.path.join(os.getcwd(), "random_comparisons_folder_matrix_01/")

create_dir_if_not_exist([input_dir, match_scores_output_dir, histogram_plot_path, reproducible_sequence_output_dir, intermediate_results])

num_comparison = 10  # NOTA: numero di confronti random da eseguire per ogni coppia di file bed

save_random_match_scores = True
save_matrix_01 = True


def extract_gene_list(genes, bed_files_dicts):
    gene_lists = []
    for bed_files_dict in bed_files_dicts:
        bed_file = bed_files_dict["bed_file"]
        table_FP = bed_file["Chromosome"].value_counts().sort_index().rename_axis('GeneID').reset_index(name='ReadsCounts')
        table_FP_Geneslengths = pd.merge(table_FP, genes, on="GeneID")

        bed_table_FP_reduced = bed_file.join(table_FP_Geneslengths.set_index('GeneID'), on='Chromosome').set_index("Chromosome")
        gene_name_list = []
        for name_of_the_group, group in bed_table_FP_reduced.groupby("Chromosome"):
            gene_name_list.append(name_of_the_group)
        gene_lists.append(gene_name_list)

    # gets the gene list for each matrix
    gene_lists = [pd.DataFrame(f, columns={"GeneID"}) for f in gene_lists]
    # gets the genes in common for each matrix
    gene_list = gene_lists[0]
    if len(gene_lists) > 1:
        for i in range(1, len(gene_lists)):
            l = gene_lists[i]
            gene_list = gene_list.merge(l, on="GeneID")

    gene_list = gene_list.to_numpy().squeeze()
    return gene_list


def random_comparison(arguments, rep_num=0): # esegue i confronti random
    start = time.time()
    matrices_extractors, genes, gene_list = arguments
    matrix_01_pair = []
    for matrices_extractor in matrices_extractors:
        bed_file_name = matrices_extractor["bed_file_name"]
        me = matrices_extractor["me"]
        # extract the matrices
        pd_matrix_coverage, matrix_01 = me.extract_matrices(areReadsRandomized=True, add_small_random_value=add_small_random_value_to_random_comparison, rep_num=rep_num)  # estrae le matrici 01 e coverage per ogni file bed
        matrix_01_pair.append({'matrix': matrix_01, 'file_name': bed_file_name})

        if save_matrix_01:
            create_dir_if_not_exist([random_comparisons_folder_matrix_01])
            save_dir = os.path.join(random_comparisons_folder_matrix_01, bed_file_name)
            create_dir_if_not_exist([save_dir])
            matrix_01.to_csv(os.path.join(save_dir, str(rep_num) + ".csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')

    match_score, pair_names = compare_pair(matrix_01_pair, genes.set_index('GeneID'), gene_list)  # compara coppie di file bed  (di matrici 0-1)

    pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
    end = time.time()
    print("end_comparison in " + str(end - start) + "seconds")
    return {'pair_name': pair_names, 'match_score': match_score}


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

    if save_random_match_scores:
        create_dir_if_not_exist([random_comparisons_folder_match_scores])
        for i in range(len(match_scores)):
            save_dir = os.path.join(random_comparisons_folder_match_scores, match_scores[i]["pair_name"])
            create_dir_if_not_exist([save_dir])
            match_scores[i]["match_score"].to_csv(os.path.join(save_dir, str(i) + ".csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')
    return match_scores


def compare_pair_n_times_serial(bed_files_pair, genes, gene_list, n): # per ogni coppia di file bed esegue n confronti
    # extract a pair of bed files
    match_scores = []

    matrices_extractors = []
    for bed_files_dict in bed_files_pair:
        bed_file = bed_files_dict["bed_file"]
        bed_file_name = bed_files_dict["bed_file_name"]
        me = MatricesExtractor(bed_file, genes, bed_file_name)
        matrices_extractors.append({"me": me, "bed_file_name": bed_file_name})

    arguments = matrices_extractors, genes, gene_list

    for i in range(n):
        match_scores.append(random_comparison(arguments, i))  # effettua n confronti random

    if save_random_match_scores:
        create_dir_if_not_exist([random_comparisons_folder_match_scores])
        for i in range(len(match_scores)):
            save_dir = os.path.join(random_comparisons_folder_match_scores, match_scores[i]["pair_name"])
            create_dir_if_not_exist([save_dir])
            match_scores[i]["match_score"].to_csv(os.path.join(save_dir, str(i) + ".csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')

    return match_scores


def main(num):
    print("num of core available: " + str(num_cores) + " used: " + str(num_task))

    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes()
    bed_files_dicts = ifm.get_bed_files()
    gene_list = extract_gene_list(genes, bed_files_dicts)

    # create pairs of bed files
    bed_files_pairs = [list(f) for f in combinations(bed_files_dicts, 2)]
    print("start fake matrix comparisons...")

    match_scores_list = []

    # start = time.time()
    # for bed_files_pair in bed_files_pairs:
    #     match_scores_list.append(compare_pair_n_times(bed_files_pair, genes, gene_list, num_comparison))
    # end = time.time()
    # print('fake matrix comparisons completed in ' + str(end - start) + " sec(s)")

    for bed_files_pair in bed_files_pairs:
        match_scores_list.append(compare_pair_n_times_serial(bed_files_pair, genes, gene_list, num_comparison))

    with open(intermediate_results + '/match_scores_list_' + str(num) + '_aggregating_' + str(num_comparison) + 'comparisons.npy', 'wb') as f:
        np.save(f, match_scores_list)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='riboseq reproducibility')
    parser.add_argument('--num', metavar='n', required=False, help='experiment identifier')
    args = parser.parse_args()
    main(num=args.num)
