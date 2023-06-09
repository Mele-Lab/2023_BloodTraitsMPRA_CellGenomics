#This code takes the qtlType_coloc_results_all_summary file generated by summarize_qtl_results.sh and generate a file where each lead SNP-Gene-Tissue result has a single row with PP3, PP4, and PP4/(PP3 + PP4)

library(tidyverse)

#use args to get the summary file name
args = commandArgs(TRUE)
all_summary_coloc_results_file <- args[1]

all_summary_coloc_results_df <- read_tsv(file = all_summary_coloc_results_file, col_names = TRUE)

#Spread the PPID file to make each results 1 row
one_row_results <- all_summary_coloc_results_df %>% spread(PPID,PP)

#calculate PP4 conditioned on two peaks PP4/(PP3+PP4)
one_row_results$cond.PP.H4.abf <-  one_row_results$PP.H4.abf/(one_row_results$PP.H3.abf + one_row_results$PP.H4.abf)

#generate string for output file
outputFilename <- gsub(".txt","spread_condPP4.txt",all_summary_coloc_results_file)

write.table(one_row_results, file= outputFilename, sep="\t", row.names=FALSE, quote=FALSE)

