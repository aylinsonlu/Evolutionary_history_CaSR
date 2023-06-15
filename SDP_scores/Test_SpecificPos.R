library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
source("./PHACT_ComputeScore.R")
source("./Specifictiy_Code_Generalized.R")
args = commandArgs(trailingOnly=TRUE)

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

###### UPDATE THESE PARTS WRT THE SUBFAMILY CONSIDERED
names <- c("TAS1R3","CaSR_like","GPRC6A","CaSR","TAS1R1","TAS1R2") # The first one is the "MAIN" Subfamily
main <- names[1]
folder_name <- "TAS1R3/"

# Full Fasta, Phylogenetic tree, Probability file
file_fasta <- paste(folder_name, "casr_casr_likes_6a_tastes_nogaptas1r3.fasta", sep="")
file_nwk_main <- paste(folder_name, "ancestral_no_gap_tas1r3.treefile", sep="")
file_rst <- paste(folder_name, "ancestral_no_gap_tas1r3.state", sep="")

############# Nodes Names & Distances are Read from the File

info <- readLines(paste(folder_name, "ancestral_nodes_distances_", main, sep=""))
info_p1 <- info[1:length(names)]
info_p2 <- info[(length(names)+1):length(info)]
nodes <- c()
distances <- c()
for (i in 1:length(names)){
  n <- info[grep(names[i], info_p1)]
  n <- unlist(strsplit(n, ": "))
  n <- n[2]
  nodes <- c(nodes, n)
}
for (i in 2:(length(names))){
  i1 <- paste(main, "-", names[i], sep = "")
  d <- info_p2[grep(i1, info_p2)]
  d <- unlist(strsplit(d, ": "))
  d <- as.numeric(d[2])
  distances <- c(distances, d)
}

####################### PART 1 - PHACT Score Computation #######################

for (i in 1:length(names)){
  name_lower_case <- tolower(names[i])
  file_nwk <- paste(folder_name, names[i], "_from_ancestral.nwk", sep="")
  main_root <- nodes[i]
  print(main_root)
  output_name <- names[i]
  compute_score(file_nwk, file_nwk_main, file_rst, file_fasta, output_name, "Equal", main_root)
}

####################### PART 2 - Specificity #######################

# PHACT Score Files
files <- paste(names, ".csv", sep = "")

# Name of the CSV file
output_name <- paste(names[1], "_Specificity", sep = "")
compute_specificity(files, distances, output_name)
