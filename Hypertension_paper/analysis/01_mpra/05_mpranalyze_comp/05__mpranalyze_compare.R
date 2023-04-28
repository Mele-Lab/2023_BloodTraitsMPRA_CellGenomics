# 05__mpranalyze_compare ####
# in this notebook, i run MPRAnalyze in 'compare' mode to get log2 foldchanges and p-values between  
# sequence orthologs.


  # load the package
library(MPRAnalyze)
library(tidyr)

library(BiocParallel)

## 1. load data ####
## first load data for library depth correction #####

dna_counts_depth <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/dna_counts.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)

# since we only have 1 dna replicate -- add another so code doesn't crash (expects matrix)
dna_counts_depth["dna_2"] <- dna_counts_depth["dna_1"]

row.names(dna_counts_depth) <- dna_counts_depth$element
dna_counts_depth <- dna_counts_depth[ , !(names(dna_counts_depth) %in% c("element")), drop=FALSE]
dna_counts_depth <- as.matrix(dna_counts_depth)

rna_counts_depth <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/rna_counts.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)
row.names(rna_counts_depth) <- rna_counts_depth$element
rna_counts_depth <- rna_counts_depth[ , !(names(rna_counts_depth) %in% c("element")), drop=FALSE]
rna_counts_depth <- as.matrix(rna_counts_depth)

dna_cols_depth <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/dna_col_ann.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)
names(dna_cols_depth) <- c("id", "condition", "sample")

# add second row to dna_cols_depth
row2 <- data.frame(id="dna_2", condition="dna", sample="2")
dna_cols_depth <- rbind(dna_cols_depth, row2)
row.names(dna_cols_depth) <- dna_cols_depth$id

rna_cols_depth <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/rna_col_ann.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)
names(rna_cols_depth) <- c("id", "condition", "sample")
row.names(rna_cols_depth) <- rna_cols_depth$id
rna_cols_depth

# make sure everything is a factor
dna_cols_depth$condition <- as.factor(dna_cols_depth$condition)
rna_cols_depth$condition <- as.factor(rna_cols_depth$condition)
rna_cols_depth$sample <- as.factor(rna_cols_depth$sample)
rna_cols_depth

## then data to model: first, DNA (same for all models) #####

all_comp_dna_counts <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/dna_counts.seq_comp.mpranalyze.CM.5.Mital.new_back.txt", sep="\t", header=TRUE)
row.names(all_comp_dna_counts) <- all_comp_dna_counts$comp_id
all_comp_dna_counts <- all_comp_dna_counts[ , !(names(all_comp_dna_counts) %in% c("comp_id","X"))]
all_comp_dna_counts <- as.matrix(all_comp_dna_counts)

all_comp_dna_cols <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/dna_col_ann.seq_comp.mpranalyze.CM.5.Mital.new_back.txt", sep="\t", header=TRUE)
row.names(all_comp_dna_cols) <- all_comp_dna_cols$X
head(all_comp_dna_cols)

all_comp_dna_cols$barcode <- as.factor(all_comp_dna_cols$barcode)
all_comp_dna_cols$seq <- as.factor(all_comp_dna_cols$seq)
all_comp_dna_cols$condition <- as.factor(all_comp_dna_cols$condition)
all_comp_dna_cols

## then controls (same for all models) ####

all_comp_ctrls <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/ctrl_status.seq_comp.mpranalyze.CM.5.Mital.new_back.txt", sep="\t", header=TRUE)
all_comp_ctrls <- as.logical(all_comp_ctrls$ctrl_status)
head(all_comp_ctrls)

length(all_comp_ctrls)

## then data to model: native effects ####
rna_counts <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/all_rna_counts.seq_comp.mpranalyze.CM.5.Mital.new_back.txt", sep="\t", header=TRUE)
row.names(rna_counts) <- rna_counts$comp_id
rna_counts <- rna_counts[ , !(names(rna_counts) %in% c("comp_id"))]
rna_counts <- as.matrix(rna_counts)

rna_cols <- read.table("../../../../Hypertension/data/01_counts/mpranalyze_files/rna_col_ann.seq_comp.mpranalyze.CM.5.Mital.new_back.txt", sep="\t", header=TRUE)
row.names(rna_cols) <- rna_cols$X
head(rna_cols)

# make sure everything is a factor
rna_cols$barcode <- as.factor(rna_cols$barcode)
rna_cols$seq <- as.factor(rna_cols$seq)
rna_cols$condition <- as.factor(rna_cols$condition)
head(rna_cols)

# 2. estimate library depth for sample/condition pair ####
# create MPRA object
depth_obj <- MpraObject(dnaCounts = dna_counts_depth, rnaCounts = rna_counts_depth, 
                        dnaAnnot = dna_cols_depth, rnaAnnot = rna_cols_depth)

# estimate depth factors using uq -- here, a sample/condition pair == 1 library
depth_obj <- estimateDepthFactors(depth_obj, lib.factor = c("sample", "condition"),  depth.estimator='uq',
                                  which.lib = "dna")
depth_obj <- estimateDepthFactors(depth_obj, lib.factor = c("id"),  
                                  depth.estimator='uq', which.lib = "rna")

rna_depths <- rnaDepth(depth_obj)
rna_depths

rna_cols_depth$depth <- rna_depths
rna_cols_depth

## 3. run model to compare native effects ####
nrow(rna_cols_depth)
nrow(rna_cols)
head(rna_cols)

# first need to set the dnadepths and rnadepths manually
all_comp_dna_cols$depth <- rep(1, nrow(all_comp_dna_cols))

# note 13 will change depending how many barcodes there are per element
rna_cols$depth <- rep(rna_depths, each=100)

# create MPRA object
all_comp_dna_counts[is.na(all_comp_dna_counts)] <- 0
rna_counts[is.na(rna_counts)] <- 0
obj <- MpraObject(dnaCounts = all_comp_dna_counts, rnaCounts = rna_counts, 
                  dnaAnnot = all_comp_dna_cols, rnaAnnot = rna_cols, controls = all_comp_ctrls,
                  BPPARAM = SnowParam(workers=16,type="SOCK"))

obj <- setDepthFactors(obj, dnaDepth = all_comp_dna_cols$depth, rnaDepth = rna_cols$depth)

obj <- analyzeComparative(obj = obj, 
                          dnaDesign = ~ barcode, 
                          rnaDesign = ~ seq, 
                          reducedDesign = ~ 1)

comp_res <- testLrt(obj)
head(comp_res)

hist(comp_res[all_comp_ctrls,]$fdr)

hist(comp_res[!all_comp_ctrls,]$fdr)

write.table(comp_res, file = "../../../data/02_activs/comp_results.CM.new_back.txt", sep = "\t",
            quote = FALSE)
