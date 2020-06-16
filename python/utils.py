import pandas as pd
import os
from os import listdir
from os.path import isfile, join
import pyranges as pr
from pathlib import Path


class InputFileManager:
    def __init__(self, gene_lengths_path, input_bed_files_dir):
        self.gene_lengths_path = gene_lengths_path
        self.input_bed_files_dir = input_bed_files_dir

    def get_genes(self):
        genes = pd.read_csv(self.gene_lengths_path)
        # removes a not useful column from the "genes" dataframe
        # renames column from the "genes" dataframe
        genes = genes.drop(genes.columns[0], axis=1)
        genes = genes.rename(columns={genes.columns[0]: "GeneID", genes.columns[1]: "GeneLength"})
        return genes

    def get_bed_files(self):
        bed_table_FP_path_list = [os.path.abspath(os.path.join(self.input_bed_files_dir, f)) for f in listdir(self.input_bed_files_dir) if isfile(os.path.join(self.input_bed_files_dir, f)) if f.endswith(".bed")]  # list all the file names in input_dir
        bed_table_FPs = [{"bed_file": pr.read_bed(f, as_df=True), "bed_file_name": Path(os.path.basename(f)).stem} for f in bed_table_FP_path_list]
        bed_table_FPs_reduced = [{"bed_file": f["bed_file"].iloc[:, [0, 1, 2]], "bed_file_name": f["bed_file_name"]} for f in bed_table_FPs]
        return bed_table_FPs_reduced


def create_dir_if_not_exist(dir_list):
    for directory in dir_list:
        if not os.path.exists(directory):
            os.makedirs(directory)  # create input directory if not already present
