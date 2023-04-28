#### get sequence logo ######
## misc.
library(BiocParallel)    # parallel computation
library(SummarizedExperiment)  # to manipulate [Ranged]SummarizedExperiment objects
library(stringr)         # manipulating character strings

## plotting
library(ggplot2)         # for plotting

## genomes
library(BSgenome.Hsapiens.UCSC.hg38)   # human genome

## pacakges w/ sequence motif-related functions
library(Biostrings)      # manipualting sequences & motif representations
library(universalmotif)  # manipulating motif representations
library(MotifDb)         # accessing motif DBs
library(TFBSTools)       # accessing motif DBs
library(JASPAR2018)      # JASPAR 2018 database
library(monaLisa)        # motif searching, enrichment, regression

#### select SNPs and TF of interest for PWM motifs #####
### 5:148385938:C:G GATA1, GATA3 REF
### rs10050933 MEF2A ALT
### rs118040162 TBX1 REF
### rs1806920 RARA REF
### rs8057203 TBX19 ALT
## snp in 68

ref <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/FIMO_results/hum_tf/fimo.txt')

res <- query(MotifDb, andStrings=(c("TBX19", "JASPAR2018", "Hsapiens")))
res <- convert_motifs(res)
view_motifs(res, dedup.names = TRUE)
