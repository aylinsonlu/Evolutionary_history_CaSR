library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(caret)
library(data.table)
source("./construct_HMM.R")
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

main_fam <- "CASR"
close_fam <- c("CASRLikes")
rest_fam <- c("GPRC6A", "TAS1R1", "TAS1R2", "TAS1R3")

folder_name <- "/Users/aylin/Documents/2021-10-20-hmm/hmm/no_gap_msa_files/no_gap_casr" # all files are assumed to be in the same folder, the updated profile will be saved to here as well
cfunvals_file <- sprintf("%s.csv", main_fam) # name of the file that contains conservation levels (csv file, names should be CASR CASR_aa GPRC6A GPRC6A_aa etc)

bound1_conslevel <- 0.9 # threshold to say main subfamily is conserved
bound2_conslevel <- 0.9 # threshold to say close/rest are conserved

bound_blosum <- 0.5 # No need to change this

### Call the main function
construct_HMM(folder_name, main_fam, close_fam, rest_fam, cfunvals_file, bound1_conslevel, bound2_conslevel, bound_blosum)






