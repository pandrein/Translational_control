rm(list = ls())

require(Hmisc)
require(data.table)

# SETUP =================================================================================================================================

# Sets the input and output directories
input_dir <- "~/Dropbox (Cambridge University)/Private-work/Ribo-seq/Reproducibility/New_analysis/Outputs/Matrices_coverage_and_01/"
output_dir <- "~/Dropbox (Cambridge University)/Private-work/Ribo-seq/Reproducibility/New_analysis/Outputs/Match_Scores/"

#DATA UPLOAD =================================================================================================================================

setwd(input_dir)

# Inputs a file containing each gene's ID and the correspondent gene length

genes_lengths <- read.csv("~/Dropbox (Cambridge University)/Private-work/Ribo-seq/Datasets/E.coli_gene_lengths.csv", sep=",",stringsAsFactors=FALSE, row.names = 2)
genes_lengths<-genes_lengths[,-1, drop=F] # removes a not useful column from the "genes" dataframe


# Uploading Matrices coverage ----------------------------------------------------------------------------------------------------------

# #Lists the files in the input directory
# files_list <- list.files("./")
# files_list_coverage <- files_list[grep(files_list, pattern = "coverage")]
# 
# # Uploads the files in a list
# matrices_coverage <- lapply(files_list_coverage, read.csv, header = T, stringsAsFactors = FALSE, row.names=1)
# names(matrices_coverage) <- tools::file_path_sans_ext(files_list_coverage)
# 
# # finds the genes (row names) in common between all the datasets
# common_genes = Reduce(intersect, lapply(matrices_coverage, row.names))
# 
# # From each matrix coverage, takes only the rows corresponding to the common genes
# matrices_coverage = lapply(matrices_coverage, function(x) { x[row.names(x) %in% common_genes,] })

# Uploading Matrices 01 ----------------------------------------------------------------------------------------------------------

#Lists the files in the input directory
files_list <- list.files("./")
files_list_01 <- files_list[grep(files_list, pattern = "01")]

# Uploads the files in a list
matrices_01 <- lapply(files_list_01, read.csv, header = T, stringsAsFactors = FALSE, row.names=1)
names(matrices_01) <- tools::file_path_sans_ext(files_list_01)

# finds the genes (row names) in common between all the datasets
common_genes = Reduce(intersect, lapply(matrices_01, row.names))

# From each matrix coverage, takes only the rows corresponding to the common genes
matrices_01 = lapply(matrices_01, function(x) { x[row.names(x) %in% common_genes,] })


# COMPARING PROFILES ##========================================================================================================

setwd(output_dir)

# Compares the matrices 01 pairwise, row by row

pairs <- combn(names(matrices_01), 2)

for(i in 1:ncol(pairs)){
  matrix_comparison <- matrix(nrow = length(common_genes), ncol = ncol(matrices_01[[pairs[1,i]]]))
  rownames(matrix_comparison) <- rownames(matrices_01[[pairs[1,i]]])
  
  for(t in 1:nrow(matrix_comparison)){
    
    string_length <- genes_lengths[rownames(matrix_comparison)[t],]
    
    matrix_comparison[t,][which(matrices_01[[pairs[1,i]]][t, 1:string_length] == matrices_01[[pairs[2,i]]][t, 1:string_length]) ] = 2
    
    matrix_comparison[t,][which(matrices_01[[pairs[1,i]]][t, 1:string_length] != matrices_01[[pairs[2,i]]][t, 1:string_length]) ] = -2
    
    }

    match_scores <- vector(length = nrow(matrix_comparison))
    for(w in 1:length(match_scores)){
      
      match_scores[w] = length(which(matrix_comparison[w,]==2))/ genes_lengths[rownames(matrix_comparison)[w],]
      names(match_scores)[w]<-rownames(matrix_comparison)[w]
      }
   
     assign(paste0("match_scores_", pairs[1,i], "_",  pairs[2,i]), match_scores )
     write.csv(match_scores, file=paste0("match_scores_", pairs[1,i], "_",  pairs[2,i]))
     #print(i)
}



