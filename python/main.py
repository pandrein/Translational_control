"""
Lo script inizialmente calcola gli score reali esattamente come prima, poi quando ha finito fa i confronti tra matrici random.
Per le matrici random ho cambiato l’ordine delle operazioni in maniera che possono essere fatti più confronti contemporaneamente per velocizzare il calcolo.
Per prima cosa raggruppo i file bed in coppie, in questo modo tutto il calcolo può essere fatto separatamente per ognuna delle coppie.
Data una coppia di file bed parte un processo separato che applica la funzione compare_pair_n_times.
La funzione calcola le matrici 01 (posizionando le reads random) per ognuno dei due file bed della coppia, poi confronta le due matrici 01 ottenute e trova i match score.
Questa operazione viene ripetuta un numero di volte dato dalla variabile num_comparison.
Quindi per esempio se num_comparison è 10^4, la funzione compare_pair_n_times farà 10^4 confronti tra coppie di matrici 01 con reads random.
Questa operazione viene ripetuta per ogni coppia di file bed.

La cosa interessante è questa. Dato che i calcoli per ogni coppia di file bed sono indipendenti, possono essere mandati in esecuzione su un processore diverso.
Quindi se hai 10 file bed, le possibili coppie sono 45, in teoria ognuno dei 45 processi (ognuno dei processi esegue num_comparison confronti in sequenza) possono essere eseguiti contemporaneamente ma su processori diversi.
In questo modo ci mette poco più che se la coppia fosse solo una. Questo se non finisce la memoria, cosa che sulla maggior parte dei computer è molto probabile, spero non per il computer di Cambridge...
Su un server come quello di Cambridge potresti provare a mettere  un numero di task = al numero di coppie (es 45 task se confronti 10 file) e num_comparison = 10^4 se non va out of memory.
"""

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

num_cores = multiprocessing.cpu_count()
num_task = num_cores

# I/O directories
input_dir = os.path.join(os.getcwd(), "check_reduced/")  # get the path to the data input directory
match_scores_output_dir = os.path.join(os.getcwd(), "matrix_python/match_scores/")  # Sets the directory where all the saved outputs will be stored
reproducible_sequence_output_dir = os.path.join(os.getcwd(), "matrix_python/reproducible_sequence/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "gene_lengths.csv")  # path to upload the file containing each gene's ID and the correspondent gene length
histogram_plot_path = os.path.join(os.getcwd(), "genes_histograms/")  # path to upload the file containing each gene's ID and the correspondent gene length

create_dir_if_not_exist([input_dir, match_scores_output_dir, histogram_plot_path,reproducible_sequence_output_dir])

num_comparison = 5  # NOTA: numero di confronti random da eseguire per ogni coppia di file bed
FDR = 0.01


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
def compute_real_match_scores(genes, bed_files_dicts):
    print ("start real matrix digitalization...")
    matrix_01_list = signal_digitalisation(genes, bed_files_dicts, areReadsRandomized=False)
    print("digitalization complete...")

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

    print ("start pair comparison...")
    # saves match scores
    match_scores = []
    pair_names_list = []
    for pair in pairs:
        # print ("compare " + pair[0]['file_name'] + " and " + pair[1]['file_name'])
        match_score, pair_names = compare_pair(pair, genes.set_index('GeneID'), gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
        match_scores.append({"match_score": match_score, "pair_name": pair_names})
        pair_names_list.append(pair_names)
        pair_names = pair_names + ".csv"
        match_score.to_csv(os.path.join(match_scores_output_dir, pair_names), index=True, header=True, decimal='.', sep=',', float_format='%.6f')
    print("real comparison complete")

    return gene_list, match_scores, pair_names_list, matrix_01_list


def compare_pair_n_times(parallel_arguments):
    bed_files_pair, genes, gene_list, n = parallel_arguments
    # extract a pair of bed files
    match_scores = []
    for i in range(n):
        print("fake_comp "+str(i))
        # start = time.time()
        matrix_01_pair = []
        for bed_files_dict in bed_files_pair:
            bed_file = bed_files_dict["bed_file"]
            bed_file_name = bed_files_dict["bed_file_name"]
            me = MatricesExtractor(bed_file, genes)
            # extract the matrices
            pd_matrix_coverage, matrix_01 = me.extract_matrices(areReadsRandomized=True)
            matrix_01_pair.append({'matrix': matrix_01, 'file_name': bed_file_name})
        # end = time.time()
        # print('01 completed in ' + str(end - start) + " sec")

        start = time.time()
        match_score, pair_names = compare_pair(matrix_01_pair, genes.set_index('GeneID'), gene_list)
        pair_names = Path(pair_names[0]).stem + ":" + Path(pair_names[1]).stem
        match_scores.append({'pair_name': pair_names, 'match_score': match_score})
        # end = time.time()
        # print('comparisons completed in ' + str(end - start) + " sec")
    return match_scores


def extract_reproducible_sequences(reproducible_genes, matrix_01_list):
    # for each matrix_01 select only the reproducible genes
    reproducible_genes_tables = [matrix_01_struct['matrix'][matrix_01_struct['matrix'].index.isin(reproducible_genes)] for matrix_01_struct in matrix_01_list]
    # select the elements that are one for all the sequences
    sequences_ones_mask = [(f == 1).to_numpy() for f in reproducible_genes_tables]  # FIX ME per fare il confronto qui passo a numpy. Andrebbe verificato che i genei confrontati siano quelli giusti
    sequences_ones_mask = np.stack(sequences_ones_mask)
    all_ones = np.all(sequences_ones_mask, axis=0)
    # select the elements that are minus one for all the sequences
    sequences_minus_ones_mask = [(f == -1).to_numpy() for f in reproducible_genes_tables]  # FIX ME verifica se si può valutare la condizione di uguaglianza con -1 e 1 in un passo solo
    sequences_minus_ones_mask = np.stack(sequences_minus_ones_mask)
    all_minus_ones = np.all(sequences_minus_ones_mask, axis=0)  # FIX ME verificare il comportamento con i nan
    # get a mask with all the elements that are one and minus one for all the sequences
    reproducible_sequence_mask = np.stack([all_ones, all_minus_ones])
    reproducible_sequence_mask = np.any(reproducible_sequence_mask, axis=0)
    return reproducible_sequence_mask, reproducible_genes_tables[0].to_numpy()


def main():
    print("num of core available: " + str(num_cores))

    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes()
    bed_files_dicts = ifm.get_bed_files()
    gene_list, match_scores_real, pair_names_list, matrix_01_list = compute_real_match_scores(genes, bed_files_dicts)

    # create pairs of bed files
    bed_files_pairs = [list(f) for f in combinations(bed_files_dicts, 2)]

    match_scores_list = []
    # # randomize and digitalise each pair
    p = multiprocessing.Pool(processes=num_task)
    # NOTA: processes rappresenta il numero di processi in parallelo che eseguono i calcoli.
    # Idealmente ce ne vorrebbe uno per ogni coppia di file bed (es. con 4 file l'ideale sarebbero 6 processi)
    # La cosa migliore è usare il numero maggiore di processi che possono stare in memoria
    print ("start fake matrix comparisons")
    start = time.time()
    # res = []
    # for bed_files_pair in bed_files_pairs:
    #     res.append(p.apply_async(compare_pair_n_times, [bed_files_pair, genes, gene_list, num_comparison]))
    #
    # p.close()
    # p.join()
    # for i in res:
    #     match_scores_list.append(i.get())
    # end = time.time()

    parallel_arguments = [[bed_files_pair,genes,gene_list,num_comparison] for bed_files_pair in bed_files_pairs]

    res = p.map_async(compare_pair_n_times, parallel_arguments)
    p.close()
    p.join()
    res = res.get()

    for i in res:
        match_scores_list.append(i)
    end = time.time()
    print('fake matrix comparisons completed in ' + str(end - start) + "seconds")

    # compute the match score histograms for the random comparisons
    match_scores_hist = {}
    for fake_match_scores in match_scores_list:
        for fake_match_score in fake_match_scores:
            # fake_match_score contains the scores of one pair
            pair_name = fake_match_score['pair_name']
            match_scores_fake = fake_match_score['match_score']
            gene_hist = {}
            for gene, match_score in match_scores_fake.items():
                gene_hist[gene] = [match_score]

            if pair_name in match_scores_hist:
                for gene, match_score in match_scores_fake.items():
                    match_scores_hist[pair_name][gene].append(match_score)  # = [match_score]
            else:
                match_scores_hist[pair_name] = gene_hist
    p_value_matrix = pd.DataFrame(index=gene_list, columns=pair_names_list)

    # extract pvalues for each gene and dataset pair
    for pair_name in match_scores_hist:
        for gene in match_scores_hist[pair_name]:
            gene_hist = pd.Series(match_scores_hist[pair_name][gene])
            hist_mean = np.mean(gene_hist)
            hist_std = np.std(gene_hist)

            for match_score_real in match_scores_real:
                pair_name_real = match_score_real["pair_name"]
                if pair_name_real == pair_name:
                    real_score = match_score_real["match_score"][gene]
                    z_score = (real_score - hist_mean) / hist_std
                    pvalue = st.norm.sf(abs(z_score))
                    p_value_matrix[pair_name][gene] = pvalue
                    # p_value_matrix[gene][pair_name] = pvalue
            # gene_hist.plot.hist(grid=True, bins=20, rwidth=0.9, color='#607c8e')
            # plt.savefig(histogram_plot_path+"dataset_pair:" + pair_name + " gene:" + gene + ".png")
            # plt.figure(i) #NOTA: fatto solo per evitare di sovrascrivere i plot, credo ci siano modi migliori...
            # i += 1

    reproducible_genes = []
    for gene, pvalue_row in p_value_matrix.iterrows():
        pvalue_row = pvalue_row.to_numpy()
        y = multipletests(pvals=pvalue_row, alpha=FDR, method="fdr_bh")
        number_of_significative_values = len(y[1][np.where(y[1] < FDR)])
        # if all the pvalues are below the threshold for each dataset then the gene can be considered reproducible
        if number_of_significative_values == len(pair_names_list):
            reproducible_genes.append(gene)

    reproducible_sequence_mask, first_matrix_01_with_only_reproducible_genes = extract_reproducible_sequences(reproducible_genes, matrix_01_list)
    # take the first matrix 01 with only reproducible genes and put to zero the non reproducible parts
    first_matrix_01_with_only_reproducible_genes[~reproducible_sequence_mask] = 0
    reproducible_sequence = pd.DataFrame(first_matrix_01_with_only_reproducible_genes, index=reproducible_genes)
    reproducible_sequence.to_csv(os.path.join(reproducible_sequence_output_dir, "reproducible_sequence.csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')

if __name__ == '__main__':
    main()
