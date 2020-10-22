#!/usr/bin/env python
# coding: utf-8

# This code generates Matrix_coverage (a matrix storing the information about the coverage of each gene in the dataset) and Matrix_01 in which,
# for each position, a 1 if the coverage in that position is higher than the median and -1 if not.

# Uploads the needed modules -----------

import os
from os import listdir
from os.path import isfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import sys
import time
from utils import create_dir_if_not_exist

np.random.seed(10)


class MatricesExtractor:
    def __init__(self, bed_table_FP_reduced, genes, bed_file_name=" "):
        self.bed_file_name = bed_file_name

        self.bed_table_FP_reduced = bed_table_FP_reduced
        self.genes = genes
        table_FP = self.bed_table_FP_reduced["Chromosome"].value_counts().sort_index().rename_axis('GeneID').reset_index(name='ReadsCounts')
        table_FP_Geneslengths = pd.merge(table_FP, self.genes, on="GeneID")

        self.gen_max_lenght = table_FP_Geneslengths["GeneLength"].max()

        self.bed_table_FP_reduced = self.bed_table_FP_reduced.join(table_FP_Geneslengths.set_index('GeneID'), on='Chromosome').set_index("Chromosome")

    def randomize_reads(self, bed_table):
        readsLength = bed_table["End"] - bed_table["Start"] + 1
        randomStartInterval = bed_table["GeneLength"] - readsLength
        random_vector = np.random.rand(len(readsLength))
        randomStarts = (random_vector * randomStartInterval).astype(int)
        randomEnds = randomStarts + readsLength
        # insert the random values
        bed_table_rand = bed_table.iloc[:, 2:4]
        bed_table_rand.insert(0, "Start", randomStarts)
        bed_table_rand.insert(1, "End", randomEnds)
        return bed_table_rand

    # defining function to be applied to every bed file
    def extract_matrices(self, areReadsRandomized=False, add_small_random_value=False, rep_num=0):

        ######################
        #### FP READS ########
        ######################
        # creates a table joining the list of genes present in bed_table_FP_reduced and the number of reads mapping on each of them
        bed_table_FP_reduced = self.bed_table_FP_reduced
        if areReadsRandomized:
            bed_table_FP_reduced = self.randomize_reads(bed_table_FP_reduced)  # randomizza tutte le reads

        random_start_reads = bed_table_FP_reduced["Start"]
        group_by_genes = random_start_reads.groupby("Chromosome")

        if areReadsRandomized:
            folder_rand_reads = "./reads_random/"
            create_dir_if_not_exist([folder_rand_reads])

            for name_of_the_group, group in group_by_genes:  # salva le reads random FIX ME andrebbe fatto solo quando areReadsRandomized = True???
                creation_path = os.path.join(folder_rand_reads, name_of_the_group)
                create_dir_if_not_exist([creation_path])
                group.to_csv(os.path.join(creation_path, "bedfile_" + self.bed_file_name + "_rep_" + str(rep_num) + ".csv"))

        gen_list = []
        gene_name_list = []
        for name_of_the_group, group in bed_table_FP_reduced.groupby("Chromosome"):
            gen_list.append(group.to_numpy())
            gene_name_list.append(name_of_the_group)
        matrix_coverage = self.create_matrix_coverage(gen_list, self.gen_max_lenght)  # crea la matrice di coverage
        matrix_01 = self.create_matrix_01(matrix_coverage, add_small_random_value=add_small_random_value) # crea la matrice 01
        matrix_coverage = pd.DataFrame(matrix_coverage, index=gene_name_list)
        matrix_01 = pd.DataFrame(matrix_01, index=gene_name_list)
        return matrix_coverage, matrix_01

    def create_matrix_coverage(self, gen_list, gen_max_lenght):
        # CREATING MATRIX COVERAGE
        matrix_coverage = np.empty([len(gen_list), gen_max_lenght])
        matrix_coverage[:] = np.nan
        # matrix_coverage = cp.array(matrix_coverage)
        # matrix_coverage = gpu.garray(matrix_coverage)
        # gen_list = gpu.garray(gen_list)
        for gene_idx in range(len(gen_list)):
            gene = gen_list[gene_idx]
            start, counts_start = np.unique(gene[:, 0], return_counts=True)
            counts_start_dict = dict(zip(start, counts_start))
            end, counts_end = np.unique(gene[:, 1], return_counts=True)
            counts_end_dict = dict(zip(end + 1, counts_end))
            gene_lenght = gene[0][3]
            # vector_coverage = np.empty([gen_max_lenght])
            # vector_coverage[:] = np.nan
            val = 0
            for i in range(gene_lenght):
                if i in counts_start_dict:
                    val += counts_start_dict[i]
                if i in counts_end_dict:
                    val -= counts_end_dict[i]
                # vector_coverage[i] = val
                matrix_coverage[gene_idx, i] = val
            # matrix_coverage[gene_idx, :] = vector_coverage
            # matrix_coverage.append(vector_coverage)
        # cp.cuda.Stream.null.synchronize()

        return matrix_coverage

    def add_rand_to_matrix(self, matrix_coverage):
        rand_add = (np.random.random(np.shape(matrix_coverage)) * 0.001) + 0.001  # add small random values to avoid an artifact that comes out when we compute the median of a vector having many zeros
        matrix_01 = matrix_coverage + rand_add
        return matrix_01

    def subtract_median_to_matrix(self, matrix_01, median):
        median = median[:, np.newaxis]
        matrix_01 = matrix_01 - median
        return matrix_01

    def create_matrix_01(self, matrix_coverage, add_small_random_value):

        if add_small_random_value:
            # CREATING MATRIX 01
            # compare the profile heights at each nucleotide position (coverage) with its median value computed along the entire ORF
            matrix_coverage = self.add_rand_to_matrix(matrix_coverage)

        # computes the median of the coverage values at each nucleotide.
        median = np.nanmedian(matrix_coverage, axis=1)  # Computes the median along the entire ORF, while ignoring NaNs.
        matrix_01 = self.subtract_median_to_matrix(matrix_coverage, median)

        # assigns +1 to the positions having a coverage value higher than the median.
        # assigns -1 to the positions having a coverage value lower than the median.
        matrix_01[matrix_01 >= 0] = 1
        matrix_01[matrix_01 < 0] = -1

        return matrix_01
