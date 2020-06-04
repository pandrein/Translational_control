# This code generates Matrix_coverage (a matrix storing the information about the coverage of each gene in the dataset) and Matrix_01 in which,
# for each position, a 1 if the coverage in that position is higher than the median and -1 if not

rm(list = ls())

#### PRELIMINARIES ############################################################################################# 

#*Uploads the needed libraries --------------------------------------------------------------------------------

require(Hmisc)
require(data.table)

#*Sets the number of significant digits for the output --------------------------
sig_dig = 4

#*Sets the working directory ---------------------------------------------------------------------------------

#gets the default wd
default_wd <- getwd()

#Sets the directory where all the saved outputs will be stored

setwd(dir = )# <--- insert here the path to the working (output) directory

new_wd <- getwd()

# DATA INPUT =========================================

#Inputs a file containing each gene's ID and the correspondent gene length
genes <- as.data.frame(read.csv("./gene_lengths.csv", sep=","),stringsAsFactors=FALSE)
genes<-genes[,-1] # removes a not useful column from the "genes" dataframe
colnames(genes) <- c("GeneID", "GeneLength")

#exrtacts as vectors the useful columns of the ''genes'' dataframe
GeneID <- as.vector(genes[,1])
GeneLength <- as.vector(genes[,2])

###############################
#### BED FILE FOOTPRINTS ######
###############################


# Uploads the .bed file containing a list of gene's ID that have at least one FP, start and end point of mapping
bed_table_FP <- as.data.frame(fread(".bed"))
bed_table_FP_reduced <- bed_table_FP[,c(1,2,3)]

######################
#### FP READS ########
######################
#creates a table joining the list of genes present in bed_table_FP_reduced and the number of reads mapping on each of them
table_FP <- as.data.frame(table(bed_table_FP_reduced[,1]))
colnames(table_FP) <- c("GeneID", "ReadsCounts")

table_FP_Geneslengths <- merge(table_FP, genes, by.x="GeneID", by.y="GeneID")

#extracts as vectors the useful columns of the ''table_FP''
#genesID_FP <- as.vector(table_FP[,1])
#readscounts_FP <- as.numeric(table_FP[,2])


#CREATING MATRIX COVERAGE AND MATRIX 01 ====================================================================

#** Matrix coverage ----------------------------------------------------------------
counts <- vector(mode="numeric", length=nrow(table_FP_Geneslengths))
matrix_coverage <- matrix(nrow=nrow(table_FP_Geneslengths), ncol=max(table_FP_Geneslengths$GeneLength))
rownames(matrix_coverage)<- table_FP_Geneslengths[,1]

for (i in 1 : length(table_FP_Geneslengths$GeneID)){
  
  sub_bed_table_FP <- bed_table_FP_reduced[which(bed_table_FP_reduced[,1] %in% table_FP_Geneslengths$GeneID[i]),]
  
  vector_coverage <- vector(mode="numeric", length=table_FP_Geneslengths$GeneLength[i])
  
  for(l in 1:nrow(sub_bed_table_FP)){
    
    vector_coverage[seq(sub_bed_table_FP[l,2]+1,sub_bed_table_FP[l,3]+1)] <- vector_coverage[seq(sub_bed_table_FP[l,2]+1,sub_bed_table_FP[l,3]+1)] + 1
  }
  
  counts[i]<- l
  
  matrix_coverage[i,1:length(vector_coverage)] <- vector_coverage
  
}

write.csv(matrix_coverage, file=" ", quote = FALSE, row.names = TRUE)


#**Matrix 01 -------------------------

matrix_01 <- matrix_coverage 

#add small random values to the matrix 01
randomize=function(x) x+runif(1, min=0.001, max=0.002) #Builds the function for randomizing the matrix
#add random values
matrix_01 <- apply(matrix_01, c(1,2), randomize)

for(m in 1:nrow(matrix_01)){
  
  matrix_01[m,1:table_FP_Geneslengths$GeneLength[m]][which(matrix_01[m, 1: table_FP_Geneslengths$GeneLength[m]]< median(matrix_01[m, 1: table_FP_Geneslengths$GeneLength[m]]))] = -1
  
  matrix_01[m,1:table_FP_Geneslengths$GeneLength[m]][which(matrix_01[m, 1: table_FP_Geneslengths$GeneLength[m]] !=-1)] = 1
  
}

write.csv(matrix_01, file=" ",quote = FALSE, row.names = TRUE)


