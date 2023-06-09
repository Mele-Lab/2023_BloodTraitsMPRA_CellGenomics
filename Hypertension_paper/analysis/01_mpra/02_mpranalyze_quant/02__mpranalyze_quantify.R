####
#in this notebook, i run the R package MPRAnalyze on the set of variants in order to quantify transcriptional activities.
#i first have MPRAnalyze estimate the library depth correction factors based on the full set of elements 
#then i have it run in quantify mode, using the RANDOM sequences as negative 
#controls.
####

# # install MPRAnalyze
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("MPRAnalyze", version = "3.12")

install.packages("RCurl")
BiocManager::install("BiocParallel")

# load the mpranalyze package
library(MPRAnalyze)
library(BiocParallel)
options(MulticoreParam=quote(MulticoreParam(workers=8)))

npar.backend <- bpparam()

## 1. Load data ####
#data for library depth correction
## change CM for VSMC to analyze the other cell line

dna_counts_depth <- read.table("../../../data/01_counts/mpranalyze_files/dna_counts.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)

# since we only have 1 dna replicate -- add another so code doesn't crash (expects matrix)
dna_counts_depth["dna_2"] <- dna_counts_depth["dna_1"]

row.names(dna_counts_depth) <- dna_counts_depth$element
dna_counts_depth <- dna_counts_depth[ , !(names(dna_counts_depth) %in% c("element")), drop=FALSE]
dna_counts_depth <- as.matrix(dna_counts_depth)

rna_counts_depth <- read.table("../../../data/01_counts/mpranalyze_files/rna_counts.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)
row.names(rna_counts_depth) <- rna_counts_depth$element
rna_counts_depth <- rna_counts_depth[ , !(names(rna_counts_depth) %in% c("element")), drop=FALSE]
rna_counts_depth <- as.matrix(rna_counts_depth)


dna_cols_depth <- read.table("../../../data/01_counts/mpranalyze_files/dna_col_ann.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)
names(dna_cols_depth) <- c("id", "condition", "sample")

# add second row to dna_cols_depth
row2 <- data.frame(id="dna_2", condition="dna", sample="2")
dna_cols_depth <- rbind(dna_cols_depth, row2)
row.names(dna_cols_depth) <- dna_cols_depth$id

rna_cols_depth <- read.table("../../../data/01_counts/mpranalyze_files/rna_col_ann.for_depth_estimation.mpranalyze.CM.5.txt", sep="\t", header=TRUE)
names(rna_cols_depth) <- c("id", "condition", "sample")
row.names(rna_cols_depth) <- rna_cols_depth$id
dna_cols_depth

rna_cols_depth

# make sure everything is a factor
dna_cols_depth$condition <- as.factor(dna_cols_depth$condition)
rna_cols_depth$condition <- as.factor(rna_cols_depth$condition)
rna_cols_depth$sample <- as.factor(rna_cols_depth$sample)

#Load data to model
dna_counts <- read.table("../../../data/01_counts/mpranalyze_files/dna_counts.mpranalyze.for_quantification.CM.5.txt", sep="\t", header=TRUE)
row.names(dna_counts) <- dna_counts$element
dna_counts <- dna_counts[ , !(names(dna_counts) %in% c("element","X"))]
dna_counts <- as.matrix(dna_counts)

rna_counts <- read.table("../../../data/01_counts/mpranalyze_files/rna_counts.mpranalyze.for_quantification.CM.5.txt", sep="\t", header=TRUE)
row.names(rna_counts) <- rna_counts$element
rna_counts <- rna_counts[ , !(names(rna_counts) %in% c("element"))]
rna_counts <- as.matrix(rna_counts)

dna_cols <- read.table("../../../data/01_counts/mpranalyze_files/dna_col_ann.mpranalyze.for_quantification.CM.5.txt", sep="\t", header=TRUE)
row.names(dna_cols) <- dna_cols$X
rna_cols <- read.table("../../../data/01_counts/mpranalyze_files/rna_col_ann.mpranalyze.for_quantification.CM.5.txt", sep="\t", header=TRUE)
row.names(rna_cols) <- rna_cols$X

# make sure everything is a factor
dna_cols$barcode <- as.factor(dna_cols$barcode)
rna_cols$barcode <- as.factor(rna_cols$barcode)
dna_cols$sample <- as.factor(dna_cols$sample)
dna_cols$condition <- as.factor(dna_cols$condition)
rna_cols$condition <- as.factor(rna_cols$condition)

ctrls <- read.table("../../../data/01_counts/mpranalyze_files/ctrl_status.mpranalyze.for_quantification.CM.5.Mital.txt", sep="\t", header=TRUE)
ctrls <- as.logical(ctrls$ctrl_status)

## 2. estimate library depth for sample/condition pair ####

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

rna_cols_depth

## 3. run MPRAnalyze quantification to get alpha per element ####

# first need to set the dnadepths and rnadepths manually
dna_cols$depth <- rep(1, nrow(dna_cols))

# note 13 will change depending how many barcodes there are per element
rna_cols$depth <- rep(rna_depths, each=25)

# create MPRA object
obj <- MpraObject(dnaCounts = dna_counts, rnaCounts = rna_counts, 
                  dnaAnnot = dna_cols, rnaAnnot = rna_cols, 
                  controls = ctrls, BPPARAM = SnowParam(workers=16,type="SOCK"))

# set depth factors manually
obj <- setDepthFactors(obj, dnaDepth = dna_cols$depth, rnaDepth = rna_cols$depth)

# analyze quantification in unpaired DNA library
obj <- analyzeQuantification(obj = obj, 
                             dnaDesign = ~ barcode,
                             rnaDesign = ~ sample)

alpha <- getAlpha(obj, by.factor = "condition")
head(alpha)

head(alpha[ctrls,])

## 4. find TSSs with significant activities by comparing them to negative controls ####

# test against negative controls
res<- testEmpirical(obj = obj, statistic = alpha$CM)
summary(res)

alpha$CM_pval <- res$pval.mad
head(alpha)

# histogram for negative controls
hist(alpha[ctrls,]$CM_pval)

# histogram for TSSs
hist(alpha[!ctrls,]$CM_pval)

# correct for multiple testing
alpha$CM_padj <- p.adjust(alpha$CM_pval, method = "fdr")
head(alpha)

hist(alpha[ctrls,]$CM_padj)

## 5. Write alphas to file ####
write.table(alpha, file = "../../../data/02_activs/alpha_per_elem.quantification.CM.txt", sep = "\t",
            quote = FALSE)

