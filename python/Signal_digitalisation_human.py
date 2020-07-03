#!/usr/bin/env python
# coding: utf-8

# This code generates Matrix_coverage (a matrix storing the information about the coverage of each gene in the dataset) and Matrix_01 in which,
# for each position, a 1 if the coverage in that position is higher than the median and -1 if not.

# Uploads the needed modules -----------

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
from os import listdir
from os.path import isfile, join
from pathlib import Path


class MatricesExtractor:
    def __init__(self, bed_table_FP_reduced, genes):
        self.bed_table_FP_reduced = bed_table_FP_reduced
        self.genes = genes

    # defining function to be applied to every bed file
    def extract_matrices(self):

        ######################
        #### FP READS ########
        ######################

        # creates a table joining the list of genes present in bed_table_FP_reduced and the number of reads mapping on each of them

        table_FP = self.bed_table_FP_reduced["Chromosome"].value_counts().sort_index().rename_axis('GeneID').reset_index(name='ReadsCounts')
        table_FP_Geneslengths = pd.merge(table_FP, self.genes, on="GeneID")

        gen_max_lenght = table_FP_Geneslengths["GeneLength"].max()

        bed_table_FP_reduced = self.bed_table_FP_reduced.join(table_FP_Geneslengths.set_index('GeneID'), on='Chromosome').set_index("Chromosome")

        gen_list = [group.to_numpy() for name_of_the_group, group in bed_table_FP_reduced.groupby("Chromosome")]
        gene_name_list = [name_of_the_group for name_of_the_group, group in bed_table_FP_reduced.groupby("Chromosome")]#added
        # extracts the matrix_coverage and matrix_01
        matrix_coverage = self.create_matrix_coverage(gen_list, gen_max_lenght)
        matrix_01 = self.create_matrix_01(matrix_coverage)

        # converts the matrices to pandas dataFrame
        # gene_name_list = bed_table_FP_reduced.index.unique()
        pd_matrix_coverage = pd.DataFrame(matrix_coverage, index=gene_name_list)
        matrix_01 = pd.DataFrame(matrix_01, index=gene_name_list)

        return pd_matrix_coverage, matrix_01

    def create_matrix_coverage(self, gen_list, gen_max_lenght):
        # CREATING MATRIX COVERAGE
        matrix_coverage = []
        for gene in gen_list:
            start, counts_start = np.unique(gene[:, 0], return_counts=True)
            counts_start_dict = dict(zip(start, counts_start))
            end, counts_end = np.unique(gene[:, 1], return_counts=True)
            counts_end_dict = dict(zip(end + 1, counts_end))
            gene_lenght = gene[0][3]
            vector_coverage = np.empty([gen_max_lenght])
            vector_coverage[:] = np.nan
            val = 0
            for i in range(gene_lenght):
                if i in counts_start_dict:
                    val += counts_start_dict[i]
                if i in counts_end_dict:
                    val -= counts_end_dict[i]
                vector_coverage[i] = val

            matrix_coverage.append(vector_coverage)

        return matrix_coverage

    def create_matrix_01(self, matrix_coverage):

        # CREATING MATRIX 01
        # compare the profile heights at each nucleotide position (coverage) with its median value computed along the entire ORF

        rand_add = (np.random.random(np.shape(matrix_coverage)) * 0.001) + 0.001  # add small random values to avoid an artifact that comes out whan we compute the median of a vector having many zeros
        matrix_01 = matrix_coverage + rand_add

        # computes the median of the coverage values at each nucleotide.
        median = np.nanmedian(matrix_01, axis=1) #Computes the median along the entire ORF, while ignoring NaNs.
        median = np.expand_dims(median, axis=1)
        median = median.repeat(np.shape(matrix_01)[1], axis=1)

        matrix_01 = np.subtract(matrix_01, median)  # element by element subtraction between arrays

        # assigns +1 to the positions having a coverage value higher than the median.
        # assigns -1 to the positions having a coverage value lower than the median.
        matrix_01[matrix_01 >= 0] = 1
        matrix_01[matrix_01 < 0] = -1

        return matrix_01


def main():
    # I/O directories
    input_data_dir = "check_reduced"
    output_data_dir = "matrix_python/check_reduced"

    working_dir = os.getcwd()  # Sets the working directory

    input_dir = os.path.join(working_dir, input_data_dir)  # get the path to the data input directory
    output_dir = os.path.join(working_dir, output_data_dir)  # Sets the directory where all the saved outputs will be stored
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)  # create input directory if not already present
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # DATA INPUT =========================================
    # genes_lengths_path = os.path.join(input_dir,"gene_lengths.csv")
    genes_lengths_path = os.path.join(working_dir, "CDShumanGenesLengths.txt")  # path to upload the file containing each gene's ID and the correspondent gene length
    #genes = pd.read_csv(genes_lengths_path)
    genes= pd.read_csv(genes_lengths_path, sep= "\t", header=None)
    # removes a not useful column from the "genes" dataframe
    # renames column from the "genes" dataframe
    #genes = genes.drop(genes.columns[0], axis=1)
    genes = genes.rename(columns={genes.columns[0]: "GeneID", genes.columns[1]: "GeneLength"})

    bed_table_FP_path_list = [os.path.abspath(os.path.join(input_dir, f)) for f in listdir(input_dir) if isfile(os.path.join(input_dir, f)) if f.endswith(".bed")]  # list all the file names in input_dir

    # execute the coverage function for each .bed file
    for bed_table_FP_path in bed_table_FP_path_list:
        bed_file_name = Path(os.path.basename(bed_table_FP_path)).stem
        coverage_matrix_csv_path = os.path.join(output_dir, bed_file_name + "_matrix_coverage.csv")
        matrix_01_csv_path = os.path.join(output_dir, bed_file_name + "_matrix_01.csv")

        ###############################
        #### BED FILE FOOTPRINTS ######
        ###############################

        # uploads the .bed files containing a list of gene's ID that have at least one FP, start and end point of mapping
        bed_table_FP = pr.read_bed(bed_table_FP_path, as_df=True)
        bed_table_FP_reduced = bed_table_FP.iloc[:, [0, 1, 2]]

        me = MatricesExtractor(bed_table_FP_reduced, genes)
        # extract the matrices
        pd_matrix_coverage, matrix_01 = me.extract_matrices()

        # Exports the dataFrames into CSV files
        matrix_01.to_csv(matrix_01_csv_path, index=True)
        pd_matrix_coverage.to_csv(coverage_matrix_csv_path, index=True)


if __name__ == '__main__':
    main()
