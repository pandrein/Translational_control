import multiprocessing
import os
import sys
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from Comparison_digital_profiles_human import compare_pair
from Signal_digitalisation_human import MatricesExtractor
from utils import InputFileManager
from utils import create_dir_if_not_exist

# np.random.seed(10)

np.set_printoptions(threshold=sys.maxsize)

load_intermediate_results_from_csv = False
add_small_random_value_to_real_scores = False


def print_full(x):
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    pd.set_option('display.float_format', '{:20,.2f}'.format)
    pd.set_option('display.max_colwidth', None)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.float_format')
    pd.reset_option('display.max_colwidth')


num_cores = multiprocessing.cpu_count()
num_task = num_cores - 1
plot_data = False
num_bins = 15

# I/O directories
input_dir = os.path.join(os.getcwd(), "./alignToCDS_bed/alignToCDS_bed")  # get the path to the data input directory
match_scores_output_dir = os.path.join(os.getcwd(), "matrix_python/match_scores/")  # Sets the directory where all the saved outputs will be stored
matrix_coverage_real_dir = os.path.join(os.getcwd(), "matrix_python/matrix_coverage_real/")
reproducible_sequence_output_dir = os.path.join(os.getcwd(), "matrix_python/reproducible_sequence/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "./alignToCDS_bed/CDShumanGenesLengths.txt")  # path to upload the file containing each gene's ID and the correspondent gene length
histogram_plot_path = os.path.join(os.getcwd(), "genes_histograms/")  # path to upload the file containing each gene's ID and the correspondent gene length
intermediate_results = os.path.join(os.getcwd(), "intermediate_results/")
plots_folder = os.path.join(os.getcwd(), "plots/")
match_scores_hist_plot_folder = os.path.join(plots_folder, "match_scores_hist/")
match_coverage_hist_plot_folder = os.path.join(plots_folder, "coverage_hist/")
hist_01_plot_folder = os.path.join(plots_folder, "01_hist/")
path_match_score_csv = os.path.join(os.getcwd(), "path_match_score_csv/")
create_dir_if_not_exist([input_dir, match_scores_output_dir, histogram_plot_path, reproducible_sequence_output_dir, intermediate_results, plots_folder, match_scores_hist_plot_folder, path_match_score_csv, matrix_coverage_real_dir, match_coverage_hist_plot_folder])

FDR = 0.01
reproducibility_min_fraction = 0.80

def signal_digitalisation(genes, bed_files_dicts, areReadsRandomized, add_small_random_value):
    matrix_01_list = []
    for bed_files_dict in bed_files_dicts:
        bed_file = bed_files_dict["bed_file"]
        bed_file_name = bed_files_dict["bed_file_name"]
        # extract the matrix coverage and the matrix 01 for each bed file
        me = MatricesExtractor(bed_file, genes, bed_file_name)
        pd_matrix_01 = me.extract_matrices(areReadsRandomized=areReadsRandomized, add_small_random_value=add_small_random_value)
        # if plot_data:
        #     # plot the matrix coverage for some specific genes
        #     for gene, coverage in pd_matrix_coverage.iterrows():
        #         coverage = coverage[~np.isnan(coverage)]
        #
        #         match_scores_hist_pair_plot_folder = os.path.join(match_coverage_hist_plot_folder, bed_file_name)
        #         create_dir_if_not_exist([match_scores_hist_pair_plot_folder])
        #
        #         x = range(0, len(coverage))
        #         if gene == "EG11773" or gene == "EG11983" or gene == "EG14133":
        #             plot = sns.lineplot(x, coverage, color='black')
        #             plot.fill_between(x, coverage, color='black')
        #
        #             plot.set(xticks=((x[0::int(len(coverage) * 0.08)])))
        #             plot.get_figure().savefig(os.path.join(match_scores_hist_pair_plot_folder, "gene:" + gene))
        #             plot.get_figure().clf()
        #     # plot the matrix 01 for some specific genes
        #     for gene, matrix_01 in pd_matrix_01.iterrows():
        #         matrix_01 = matrix_01[~np.isnan(matrix_01)]
        #
        #         hist_01_hist_pair_plot_folder = os.path.join(hist_01_plot_folder, bed_file_name)
        #         create_dir_if_not_exist([hist_01_hist_pair_plot_folder])
        #         # print_full(coverage)
        #         x = range(0, len(matrix_01))
        #         if gene == "EG11773" or gene == "EG11983" or gene == "EG14133":
        #             plot = sns.lineplot(x, matrix_01, color='black')
        #             plot.fill_between(x, matrix_01, color='black')
        #
        #             plot.set(xticks=((x[0::int(len(matrix_01) * 0.08)])))
        #
        #             plot.get_figure().savefig(os.path.join(hist_01_hist_pair_plot_folder, "gene:" + gene))
        #             plot.get_figure().clf()
        #             matrix_01.to_csv(os.path.join(hist_01_hist_pair_plot_folder, "gene" + gene + ".csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')

        matrix_01_list.append({'matrix': pd_matrix_01, 'file_name': bed_file_name})

    return matrix_01_list


# FIX ME al momento esegue la parte modificata comparison_digital_profiles_human, verificare che sia corretta. Per e.coli provare la comparison_digital_profiles.
def compute_real_match_scores(genes, bed_files_dicts, save_results=True):
    print("start real matrix digitalization...")
    matrix_01_list = signal_digitalisation(genes, bed_files_dicts, areReadsRandomized=False, add_small_random_value=add_small_random_value_to_real_scores)
    print("digitalization complete...")
    # first = matrix_01_list[0]["matrix"]
    # second = matrix_01_list[1]["matrix"]

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
        if save_results:
            first_matrix_coverage = pair[0]['matrix_coverage']
            second_matrix_coverage = pair[1]['matrix_coverage']
            first_file_name = pair[0]['file_name']
            second_file_name = pair[1]['file_name']

            first_matrix_coverage.to_csv(os.path.join(matrix_coverage_real_dir, first_file_name), index=True, header=True, decimal='.', sep=',', float_format='%.6f')
            second_matrix_coverage.to_csv(os.path.join(matrix_coverage_real_dir, second_file_name), index=True, header=True, decimal='.', sep=',', float_format='%.6f')

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


def calc_reproducible_sequences(match_scores_list, gene_list, pair_names_list, match_scores_real, matrix_01_list):
    # compute the match score histograms for the random comparisons
    # match_score_list contains all the match scores computed during the random comparisons
    match_scores_hist = {}
    for fake_match_scores in match_scores_list:
        for fake_match_score in fake_match_scores:
            pair_name = fake_match_score['pair_name']
            match_scores_fake = fake_match_score['match_score']
            # match_scores_hist aggregates the match scores indexing by pair name (name of the pair of bed files) and gene
            gene_hist = {}
            for gene, match_score in match_scores_fake.items():
                gene_hist[gene] = [match_score]

            if pair_name in match_scores_hist:
                for gene, match_score in match_scores_fake.items():
                    match_scores_hist[pair_name][gene].append(match_score)  # = [match_score]
            else:
                match_scores_hist[pair_name] = gene_hist

    # the pvale matrix contains the pvalues indexed by the name of pair of bed files and the name of the gene
    p_value_matrix = pd.DataFrame(index=gene_list, columns=pair_names_list)

    # the matrix summary is a pandas dataframe that summarizes the mean, standard deviation, pvalue and zscores of each gene
    matrix_summary_columns = [[c + "_mean", c + "_std", c + "_zscore", c + "_pvalue"] for c in pair_names_list]
    matrix_summary_columns = np.reshape(matrix_summary_columns, (np.shape(matrix_summary_columns)[0] * np.shape(matrix_summary_columns)[1])).T
    matrix_summary = pd.DataFrame(index=gene_list, columns=matrix_summary_columns)
    plot_num = 0
    # extract pvalues for each gene and dataset pair
    for pair_name in match_scores_hist:
        for gene in match_scores_hist[pair_name]:
            # compute mean and standard deviation for each gene
            gene_hist = pd.Series(match_scores_hist[pair_name][gene])
            hist_mean = np.mean(gene_hist)
            hist_std = np.std(gene_hist)

            if plot_data:
                # plot the histogram of match scores for each gene
                match_scores_hist_pair_plot_folder = os.path.join(match_scores_hist_plot_folder, pair_name)
                create_dir_if_not_exist([match_scores_hist_pair_plot_folder])
                sns.set_style('darkgrid')
                plot = sns.distplot(gene_hist, bins=num_bins).set_title("hist_mean: " + str('%.5f' % hist_mean) + "   hist_std: " + str('%.5f' % hist_std))
                plot.get_figure().savefig(os.path.join(match_scores_hist_pair_plot_folder, "gene:" + gene))
                plot.get_figure().clf()

            # compute the pvalues of each gene for each pair of bed files
            for match_score_real in match_scores_real:
                pair_name_real = match_score_real["pair_name"]
                if pair_name_real == pair_name:
                    real_score = match_score_real["match_score"][gene]
                    if hist_std != 0:
                        # compute the zscore and corresponding pvalue
                        z_score = (real_score - hist_mean) / hist_std
                        pvalue = 1 - st.norm.cdf(z_score)
                    else:
                        z_score = np.nan
                        pvalue = np.nan

                    # insert the results in the matrix summary
                    matrix_summary[pair_name + "_mean"][gene] = hist_mean
                    matrix_summary[pair_name + "_std"][gene] = hist_std
                    matrix_summary[pair_name + "_zscore"][gene] = z_score
                    matrix_summary[pair_name + "_pvalue"][gene] = pvalue
                    p_value_matrix[pair_name][gene] = pvalue

            plot_num += 1

    reproducible_genes = []
    matrix_summary.to_csv(os.path.join(reproducible_sequence_output_dir, "matrix_summary_005.csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')
    for gene, pvalue_row in p_value_matrix.iterrows():
        pvalue_row = pvalue_row.to_numpy()

        # test the pvalues with Benjamini/Hochberg

        # python version
        # use multipletest from python stats library, results are the same as in R script: number_of_significative_values_python) ==len(bh_candidates)

        # y = multipletests(pvals=pvalue_row, alpha=FDR, method="fdr_bh")
        # number_of_significative_values_python = len(y[1][np.where(y[1] < FDR)])

        # # if all the pvalues are below the threshold for each dataset then the gene can be considered reproducible
        # if number_of_significative_values == len(pair_names_list):
        #     reproducible_genes.append(gene)

        # R version
        pvalue_row = np.sort(pvalue_row)
        critical_values = ((np.nonzero(pvalue_row >= 0)[0] + 1) / len(pair_names_list)) * FDR

        # select the gene only if the number of selected elements is above the threshold
        reproducibily_threshold = int(round((len(pair_names_list) * reproducibility_min_fraction)))
        #print(reproducibily_threshold)
        if len(critical_values) > 0:
            # remove nan from pvalue row
            pvalue_row = pd.to_numeric(pvalue_row, errors='coerce')
            pvalue_row = pvalue_row[np.logical_not(np.isnan(pvalue_row))]
            bh_candidates = pvalue_row[pvalue_row <= critical_values]
            if len(bh_candidates) > 0:
                idx_of_max_value = np.argwhere(bh_candidates == np.amax(bh_candidates)).flatten().tolist()[-1] + 1
                bh_selected = pvalue_row[np.array(range(0, idx_of_max_value))]
                if len(bh_selected) >= reproducibily_threshold:
                    reproducible_genes.append(gene)

    # extract the reproducible sequences
    reproducible_sequence_mask, first_matrix_01_with_only_reproducible_genes = extract_reproducible_sequences(reproducible_genes, matrix_01_list)
    # take the first matrix 01 with only reproducible genes and put to zero the non reproducible parts
    first_matrix_01_with_only_reproducible_genes[~reproducible_sequence_mask] = 0
    reproducible_sequence = pd.DataFrame(first_matrix_01_with_only_reproducible_genes, index=reproducible_genes)
    reproducible_sequence.to_csv(os.path.join(reproducible_sequence_output_dir, "reproducible_sequence_005.csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')
    p_value_matrix.to_csv(os.path.join(reproducible_sequence_output_dir, "global_matrix_005.csv"), index=True, header=True, decimal='.', sep=',', float_format='%.6f')


def extract_reproducible_sequences(reproducible_genes, matrix_01_list):
    # for each matrix_01 select only the reproducible genes (that are present in reproducible_genes)
    reproducible_genes_tables = [matrix_01_struct['matrix'][matrix_01_struct['matrix'].index.isin(reproducible_genes)] for matrix_01_struct in matrix_01_list]

    # insert padding to have sequences with the same length
    sizes = [f.shape for f in reproducible_genes_tables]
    max_size = max(sizes, key=lambda x: x[1])[1]
    matrix_with_same_col_dim = []
    for m in reproducible_genes_tables:
        size = m.shape[1]
        fill_range = np.arange(size, max_size)
        fill_range = [str(item) for item in fill_range]
        matrix_with_same_col_dim.append(m.reindex(list(m) + fill_range, axis=1))

    # select the elements that are one for all the sequences
    sequences_ones_mask = [(f == 1).to_numpy() for f in matrix_with_same_col_dim]
    sequences_ones_mask = np.stack(sequences_ones_mask)
    all_ones = np.all(sequences_ones_mask, axis=0)

    # select the elements that are minus one for all the sequences
    sequences_minus_ones_mask = [(f == -1).to_numpy() for f in matrix_with_same_col_dim]
    sequences_minus_ones_mask = np.stack(sequences_minus_ones_mask)
    all_minus_ones = np.all(sequences_minus_ones_mask, axis=0)
    # get a mask with all the elements that are one and minus one for all the sequences
    reproducible_sequence_mask = np.stack([all_ones, all_minus_ones])
    reproducible_sequence_mask = np.any(reproducible_sequence_mask, axis=0)
    return reproducible_sequence_mask, matrix_with_same_col_dim[0].to_numpy()


def main():
    print("num of core available: " + str(num_cores) + " used: " + str(num_task))

    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes_human()
    bed_files_dicts = ifm.get_bed_files()

    # the match scores are saved in one or multiple files (i.e. intermediate results)
    # match_scores_list contains a list of all the random match scores
    match_scores_list = []
    if load_intermediate_results_from_csv:
        # load the match scores from a csv it is inefficient and used only for debug
        intermediate_results_path_list = [os.path.abspath(os.path.join(path_match_score_csv, f)) for f in os.listdir(path_match_score_csv) if f.endswith(tuple([".csv"]))]
        for i in intermediate_results_path_list:
            with open(i, 'r') as f:
                match_scores = []
                gene_list = []
                ir = pd.read_csv(f)
                for (columnName, columnData) in ir.iteritems():
                    if columnName == "Unnamed: 0":
                        gene_list.append(columnData.to_numpy())
                    else:
                        match_scores.append(columnData.to_numpy())

                for match_score in match_scores:
                    s = pd.Series(match_score, index=gene_list[0]).rename_axis("GeneID")
                    match_scores_list.append({"pair_name": str(Path(os.path.basename(i)).with_suffix('')), "match_score": s})
        match_scores_list = [match_scores_list]
    else:
        # load the match scores from npy files in to the match_scores_list
        intermediate_results_path_list = [os.path.abspath(os.path.join(intermediate_results, f)) for f in os.listdir(intermediate_results) if f.endswith(tuple([".npy"]))]
        for i in intermediate_results_path_list:
            with open(i, 'rb') as f:
                ir = np.load(f, allow_pickle=True)
                match_scores_list.extend(ir)

    # compute the real match scores
    print("start computing real match scores")
    gene_list, match_scores_real, pair_names_list, matrix_01_list = compute_real_match_scores(genes, bed_files_dicts, save_results=False)
    print("end computing real match scores")

    print("start computing reproducible sequences")
    # compute the reproducible sequences
    calc_reproducible_sequences(match_scores_list, gene_list, pair_names_list, match_scores_real, matrix_01_list)


if __name__ == '__main__':
    main()
