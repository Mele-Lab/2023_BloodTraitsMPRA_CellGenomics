##### merge ñiftover annotations with GTEx ------

library(data.table)
#library(R.utils)
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 


### read gtex file 
tissue = args[1]
eQTL_all_header = c("chrom_b37", "chromStart_b37", "chromEnd_b37", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")

all_pairs <- read.table(file = paste0('hg19/',tissue,'.allpairs.tab.gz'), sep = '\t', header = F)
colnames(all_pairs) <- eQTL_all_header
print(nrow(all_pairs))

#### read liftover file 
liftover_hg38 <- read.table(file = paste0('hg19/',tissue,'_temp_hg38.bed'), sep = '\t', header = F)
colnames(liftover_hg38) <- c("chrom_b38", "chromStart_b38", "chromEnd_b38", "chromStart_b37", "chromEnd_b37", "eGeneID", "A1_eqtl", "A2_eqtl")
print(nrow(liftover_hg38))

### merge data 
liftover_hg38 <- as.data.table(liftover_hg38)
all_pairs <- as.data.table(all_pairs)

hg38_all_pairs <- merge(all_pairs, liftover_hg38)
print(nrow(hg38_all_pairs))
hg38_all_pairs <- as.data.frame(hg38_all_pairs)
rm(all_pairs)
rm(liftover_hg38)

### re-arrange columns 
cols = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")
hg38_all_pairs <- hg38_all_pairs[,cols]

#### save file 
print("writing table")
write.table(hg38_all_pairs, paste0(tissue, '.allpairs.tab'), sep = '\t', col.names = F, row.names = F, quote = F)
print("done")
