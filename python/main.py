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

# num_cores = multiprocessing.cpu_count()

# I/O directories
input_dir = os.path.join(os.getcwd(), "check_reduced/")  # get the path to the data input directory
output_dir = os.path.join(os.getcwd(), "matrix_python/match_scores/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "CDShumanGenesLengths.txt")  # path to upload the file containing each gene's ID and the correspondent gene length
histogram_plot_path = os.path.join(os.getcwd(), "genes_histograms/")  # path to upload the file containing each gene's ID and the correspondent gene length

create_dir_if_not_exist([input_dir, output_dir, histogram_plot_path])

num_comparison = 5  # NOTA: numero di confronti random da eseguire per ogni coppia di file bed


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


def compare_pair_n_times(bed_files_pair, genes, gene_list, n):
    # extract a pair of bed files
    match_scores = []
    for i in range(n):
        matrix_01_pair = []
        for bed_files_dict in bed_files_pair:
            bed_file = bed_files_dict["bed_file"]
            bed_file_name = bed_files_dict["bed_file_name"]
            me = MatricesExtractor(bed_file, genes)
            # extract the matrices
            pd_matrix_coverage, matrix_01 = me.extract_matrices(areReadsRandomized=True)
            matrix_01_pair.append({'matrix': matrix_01, 'file_name': bed_file_name})

        match_score, pair_names = compare_pair(matrix_01_pair, genes.set_index('GeneID'), gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
        match_scores.append({'pair_name': pair_names, 'match_score': match_score})
    return match_scores


def main():
    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes_human()
    bed_files_dicts = ifm.get_bed_files()
    gene_list, match_scores = compute_real_match_scores(genes, bed_files_dicts)

    # create pairs of bed files
    bed_files_pairs = [list(f) for f in combinations(bed_files_dicts, 2)]

    match_scores_list = []
    # # randomize and digitalise each pair
    p = multiprocessing.Pool(processes=1)
    # NOTA: processes rappresenta il numero di processi in parallelo che eseguono i calcoli.
    # Idealmente ce ne vorrebbe uno per ogni coppia di file bed (es. con 4 file l'ideale sarebbero 6 processi)
    # La cosa migliore è usare il numero maggiore di processi che possono stare in memoria
    start = time.time()
    res = []
    for bed_files_pair in bed_files_pairs:
        res.append(p.apply_async(compare_pair_n_times, [bed_files_pair, genes, gene_list, num_comparison]))

    p.close()
    p.join()

    for i in res:
        match_scores_list.append(i.get())
    end = time.time()
    print("Complete")
    print('total time (s)= ' + str(end - start))

    print(len(match_scores_list[0]))

    match_scores_hist = {}
    for fake_match_scores in match_scores_list:
        for fake_match_score in fake_match_scores:
            # fake_match_score contains the scores of one pair
            pair_name = fake_match_score['pair_name']
            match_scores = fake_match_score['match_score']
            gene_hist = {}
            for gene, match_score in match_scores.items():
                gene_hist[gene] = [match_score]

            if pair_name in match_scores_hist:
                for gene, match_score in match_scores.items():
                    match_scores_hist[pair_name][gene].append(match_score)  # = [match_score]
            else:
                match_scores_hist[pair_name] = gene_hist
    i = 0
    for pair_name in match_scores_hist:
        for gene in match_scores_hist[pair_name]:
            gene_hist = pd.Series(match_scores_hist[pair_name][gene])
            gene_hist.plot.hist(grid=True, bins=20, rwidth=0.9, color='#607c8e')
            plt.savefig(histogram_plot_path+"dataset_pair:" + pair_name + " gene:" + gene + ".png")
            plt.figure(i) #NOTA: fatto solo per evitare di sovrascrivere i plot, credo ci siano modi migliori...
            i += 1

    # NOTA: ci sarebbe da fare i confronti tra i reali e gli istogrammi creati però non ho ben capito come fare


if __name__ == '__main__':
    main()
