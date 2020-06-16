from utils import InputFileManager
from utils import create_dir_if_not_exist
from Signal_digitalisation import MatricesExtractor
import os

# I/O directories
input_dir = os.path.join(os.getcwd(), "input_data/")  # get the path to the data input directory
output_dir = os.path.join(os.getcwd(), "output_matrix/")  # Sets the directory where all the saved outputs will be stored
genes_lengths_path = os.path.join(os.getcwd(), "gene_lengths.csv")  # path to upload the file containing each gene's ID and the correspondent gene length

addRand = True  # add noise
areReadsRandomized = True  # randomize reads

create_dir_if_not_exist([input_dir, output_dir])


def main():
    ifm = InputFileManager(genes_lengths_path, input_dir)
    genes = ifm.get_genes()
    bed_files_dicts = ifm.get_bed_files()
    for bed_files_dict in bed_files_dicts:
        bed_file = bed_files_dict["bed_file"]
        bed_file_name = bed_files_dict["bed_file_name"]
        me = MatricesExtractor(bed_file, genes)
        # extract the matrices
        pd_matrix_coverage, matrix_01 = me.extract_matrices(addRand=addRand, areReadsRandomized=areReadsRandomized)

        coverage_matrix_csv_path = os.path.join(output_dir, bed_file_name + "_matrix_coverage.csv")
        matrix_01_csv_path = os.path.join(output_dir, bed_file_name + "_matrix_01.csv")

        # Exports the dataFrames into CSV files
        matrix_01.to_csv(matrix_01_csv_path, index=True)
        pd_matrix_coverage.to_csv(coverage_matrix_csv_path, index=True)


if __name__ == '__main__':
    main()
