#### analysis outlier regions clustering + haplotypes ######
#### read data #####

sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

length <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_snps_length_block.txt',
                     sep = '\t', header = T)
head(length)

CM_results_s <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                           header = T)

VSMC_results_s <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                             header = T)


#### repeat density without first and last ######
sentinel_snps$DiffCM <- 'Not Diff'
sentinel_snps$DiffVSMC <- 'Not Diff'
sentinel_snps$DiffCM[sentinel_snps$snp %in% CM_results_s$snp_info] <- 'Diff'
sentinel_snps$DiffVSMC[sentinel_snps$snp %in% VSMC_results_s$snp_info] <- 'Diff'
table(sentinel_snps$DiffCM)

head(length)
sentinel_snps <- sentinel_snps %>%
  separate(change, c("chr", "pos","ref","alt"), ":")
sentinel_snps$pos_sent <- paste0(sentinel_snps$sentinel,':',sentinel_snps$pos)

length$pos_min <- paste0(length$sentinel,':',length$min)
length$pos_max <- paste0(length$sentinel,':',length$max)

head(sentinel_snps)
#first_last_sentinel <- sentinel_snps[!sentinel_snps$pos_sent %in% c(length$pos_min, length$pos_max),]
first_last_sentinel <- sentinel_snps

number_snps <- as.data.frame(table(first_last_sentinel$sentinel[first_last_sentinel$DiffVSMC == 'Diff']))
number_snps_all <- as.data.frame(table(first_last_sentinel$sentinel))
number_snps_all <- merge(number_snps_all, length, by.x='Var1', by.y='sentinel')
number_snps_all$clust <- number_snps_all$Freq/((number_snps_all$length)/1000000)

number_snps <- merge(number_snps, length, by.x='Var1', by.y='sentinel')
number_snps <- merge(number_snps, number_snps_all, by='Var1')

number_snps$ratio <- (number_snps$Freq.x/number_snps$clust)
number_snps$ratio_all <- (number_snps$Freq.y/((number_snps$length.y/1000000)))
number_snps$ratio_reg <- (number_snps$Freq.x/((number_snps$length.y/1000000)))

vector <-number_snps$ratio_reg[number_snps$length.x>300]
hist((vector), col = 'steelblue', breaks = 60, main = 'Nº of reg SNPs / Mb', xlab = 'density of reg variants') + 
  abline(v=median((vector)), col='red')


### get outliers top 10% cluster ####
o <- order(number_snps$ratio,decreasing=TRUE)
ordered <- (number_snps[o,])
top <- ordered %>% top_frac(.1)

#### plot scores for cluster outliers #####
CM_results_s$clustering_group <- 'Bottom density'
CM_results_s$clustering_group[CM_results_s$sentinel %in% top$Var1] <- 'Top density'
write.table(CM_results_s[,c('snp_info','sentinel','clustering_group')], 
            '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_clustering_info.txt',
            quote = F, sep = '\t', row.names = F, col.names = T)

VSMC_results_s$clustering_group <- 'Bottom density'
VSMC_results_s$clustering_group[VSMC_results_s$sentinel %in% top$Var1] <- 'Top density'
write.table(VSMC_results_s[,c('snp_info','sentinel','clustering_group')], 
            '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_clustering_info.txt',
            quote = F, sep = '\t', row.names = F, col.names = T)

library(wesanderson)
library(ggpubr)
group_by(VSMC_results_s, clustering_group) %>%
  summarise(
    count = n(),
    median = median(score_groups_2, na.rm = TRUE),
    IQR = IQR(score_groups_2, na.rm = TRUE)
  )
compare_means(score_groups_2 ~ clustering_group, data = VSMC_results_s)
wilcox.test(score_groups_2 ~ clustering_group, data = VSMC_results_s, 
            exact = FALSE, alternative = "less")
var.test(score_groups_2 ~ clustering_group, data = VSMC_results_s)
t.test(score_groups_2 ~ clustering_group, data = VSMC_results_s,
       var.equal = F, alternative = "less")

ggplot(VSMC_results_s, aes(x = clustering_group, y = score_groups_2, fill=clustering_group)) + 
  geom_boxplot() + 
  theme_classic() + xlab('') + ylab('Score') + 
  scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest1"))+
  annotate("text",
           x = 1:length(table(VSMC_results_s$clustering_group)),
           y = aggregate(score_groups_2 ~ clustering_group, VSMC_results_s, mean)[ , 2]+0.25,
           label = table(VSMC_results_s$clustering_group),
           col = "black",
           vjust =  2) +
  theme(legend.position="none") +
  stat_compare_means(label.x = 1.25, label.y = 7)

###### clustering of reg in genes --- outliers #####

library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#### go_genes #####
reg_blood_pressure <- read.delim('~/marenostrum/Data/Hypertension/GO_reg_of_BP.txt', header = F)

sentinel <- sentinel_snps[,c('sentinel','chr')] %>% distinct()
length_sentinel <- merge(length, sentinel)
snp_coord <- as.list(paste0(length_sentinel$chr,':',length_sentinel$min,':',length_sentinel$max))

names(snp_coord) <- length_sentinel$sentinel
snp_coord_df = lapply(snp_coord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  as.numeric %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  mutate(chrom=paste0('chr', chrom))

rownames(snp_coord_df) <- length_sentinel$sentinel

for (snp in length_sentinel$sentinel) {
  tryCatch({
    mycoords.gr <- makeGRangesFromDataFrame(snp_coord_df[snp,])
    genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), mycoords.gr)
    entrez_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="ENSEMBL", keytype="ENTREZID", multiVals="first")
    symbol_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="SYMBOL", keytype="ENTREZID", multiVals="first")
    
    genes$ensembl <- entrez_DE
    genes$symbol <- symbol_DE
    
    chr <- as.character(unique(seqnames(mycoords.gr)))
    gen <- genome(genes)[1]
    
    options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")
    
    gtrack <- GenomeAxisTrack()
    itrack <- IdeogramTrack(genome = 'hg19', chromosome = chr)
    
    reg_snps <- CM_results_s[CM_results_s$sentinel == snp,c('chr',"start")]
    reg_snps$end <- reg_snps$start+1
    reg_snps$chr <- paste0('chr',reg_snps$chr)
    
    snp_coord_gr <- makeGRangesFromDataFrame(reg_snps)
    atrack <- AnnotationTrack(snp_coord_gr, name = 'reg variants', id = ' ', genome = 'hg19')
    
    ### all snps ####
    all_snps <- sentinel_snps[sentinel_snps$sentinel == snp,c('chr',"pos")]
    all_snps$end <- as.integer(all_snps$pos)+1
    colnames(all_snps) <- c('chr','start','end')
    all_snps$start <- as.numeric(all_snps$start)
    all_snps$chr <- as.numeric(all_snps$chr)
    
    snp_coord_gr_all <- makeGRangesFromDataFrame(all_snps)
    atrack_all <- AnnotationTrack(snp_coord_gr_all, name = 'LD variants', id = ' ')
    
    grtrack <- GeneRegionTrack(genes, genome = gen,
                               chromosome = chr, name = "Genes", 
                               transcriptAnnotation = "symbol",
                               background.panel = "#FFFEDB",
                               background.title = "darkblue")
    
    pdf(file = paste0("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/plots/gen_region/",snp,"gen_track.pdf"),   # The directory you want to save the file in
        width = 6, # The width of the plot in inches
        height = 4) # The height of the plot in inches
    plotTracks(list(itrack,gtrack,atrack, grtrack),
                     featureAnnotation = "id", 
                     fontcolor.feature = "darkblue", 
                     just.feature = 'above',  
                     col = NULL, 
                     shape="arrow") 
    dev.off()

  }, error = function(err){
    
    print("no gene in the region")
    
    
  })
  #plotTracks(list(itrack,atrack,grtrack))
  

}

#### check MPRA bias on haplotype #####
### first get haplotypes ####
library(LDlinkR)
library(ieugwasr)
library(zoo)

for (snp in length_sentinel$sentinel[9:nrow(length_sentinel)]) {
  snps_ld <- sentinel_snps$snp[sentinel_snps$sentinel == snp]
  snps_ld <- snps_ld[grep('rs',snps_ld)]
  
  if (length(snps_ld)<2) {
    snps_ld <- c(snps_ld,snps_ld)
  }
  
  print(snp)
  
  # LD_2 <- LDmatrix(snps = snps_ld,
  #                    pop = c("CEU","TSI","FIN","GBR","IBS"),
  #                    token = "c2ed83367cc2",
  #                    genome_build = "grch38")
  
  if (length(snps_ld)<30) {
    LD_plink<- LDhap(snps = snps_ld, 
                     pop =  c("CEU","TSI","FIN","GBR","IBS"), 
                     token = "c2ed83367cc2",
                     genome_build = "grch38_high_coverage"
    )
    print("done LD")
    
    write.table(as.data.frame(LD_plink), 
                file = paste0('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/haplotype_info/',snp,'_1','_LD_matrix.txt'), 
                quote = F, col.names = T, row.names = F, sep = '\t')
    
  } else {
    
    batches <- rollapply( 1:length(snps_ld), 30, c, by=5, partial=5, align="left" )
    
    for (batch in c(1:nrow(batches))) {
      LD_plink<- LDhap(snps = snps_ld[batches[batch,]], 
                       pop =  c("CEU","TSI","FIN","GBR","IBS"), 
                       token = "c2ed83367cc2",
                       genome_build = "grch38_high_coverage"
      )
      print("done LD")
      
      write.table(as.data.frame(LD_plink), 
                  file = paste0('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/haplotype_info/',snp,'_',batch,'_LD_matrix.txt'), 
                  quote = F, col.names = T, row.names = F, sep = '\t')
      
    }}
  
}


##### read haplotypes and merge them by sentinel ######
sentinel <- sentinel_snps[,c('sentinel','chr')] %>% distinct()
sentinel_hapl <- list()
CM_results_s$sign <- ''
CM_results_s$snp_ref <- paste0(CM_results_s$snp_info,'_',CM_results_s$ref)
CM_results_s$snp_alt <- paste0(CM_results_s$snp_info,'_',CM_results_s$alt)
for (snp in sentinel$sentinel) {
  files <- list.files('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/haplotype_info/', 
                      pattern=paste0(snp,'_*'))
  
  if (length(files)>1) {
    all_hapl <- lapply(c(1:length(files)), function(x) as.data.frame(t(read.table(paste0('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/haplotype_info/',
                                                                           files[x]), header = T, sep = '\t',colClasses = "character"))[,c(1)]))
    
    all_hapl_df <- do.call('rbind.data.frame', all_hapl)
    colnames(all_hapl_df) <- ('all_hapl')
  } else {
  
  all_hapl <- t(read.table(paste0('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/haplotype_info/',
                                                                  files[1]), header = T, sep = '\t',colClasses = "character"))[,c(1)]
  all_hapl_df <- as.data.frame(all_hapl)
  }
  all_hapl_df$snp <- rownames(all_hapl_df)
  
  snps_ld <- sentinel_snps[sentinel_snps$sentinel == snp,c('sentinel','snp','ref','alt')]
  
  all_hapl_df <- all_hapl_df[grep('rs', all_hapl_df$snp),] %>% distinct()
  all_hapl_df <- all_hapl_df[all_hapl_df$snp %in% snps_ld$snp,]
  all_hapl_df <- merge(all_hapl_df, snps_ld, by.x=c('snp','all_hapl'), by.y=c('snp','ref'), all.x=T)
  all_hapl_df <- merge(all_hapl_df, snps_ld, by.x=c('snp','all_hapl'), by.y=c('snp','alt'), all.x=T)
  
  all_hapl_df <- all_hapl_df[!rowSums(is.na(all_hapl_df))>2,]
  tryCatch({
  all_hapl_df$sign_rev <- paste0(all_hapl_df$snp,'_',all_hapl_df$ref)
  all_hapl_df$sign_alt <- paste0(all_hapl_df$snp,'_',all_hapl_df$alt)
  
  CM_results_s$sign[CM_results_s$snp_info %in% all_hapl_df$snp & CM_results_s$snp_ref %in% all_hapl_df$sign_rev] <- '-'
  CM_results_s$sign[CM_results_s$snp_info %in% all_hapl_df$snp & CM_results_s$snp_alt %in% all_hapl_df$sign_alt] <- '+'
  
  }, error = function(err){
    print("not haplotype info")
})
}

#### calculate possible bias ####
sentinel$binomial_pval <- ''
sentinel$n_pos <- ''
sentinel$n_test <- ''
sentinel$total <- ''

sentinel$binomial_pval_2 <- ''
sentinel$n_pos_2 <- ''
sentinel$n_test_2 <- ''
sentinel$total_2 <- ''

for (snp in sentinel$sentinel) {
  snps_ld <- CM_results_s[CM_results_s$sentinel == snp,]
  snps_ld <- snps_ld[snps_ld$sign != '',]
  #snps_ld$logFC_comp[snps_ld$sign == '-'] <- -(snps_ld$logFC_comp[snps_ld$sign == '-'])
  
  n_pos <- length(snps_ld$logFC_comp[snps_ld$logFC_comp>0])
  n_neg <- length(snps_ld$logFC_comp[snps_ld$logFC_comp<0])
  total <- length(snps_ld$logFC_comp)
  tryCatch({
  binomial <- binom.test(n_pos, total, 0.5)$p.value
  } , error = function(err) {
    binomial <- NA
  })
  
  sentinel$binomial_pval_2[sentinel$sentinel == snp] <- binomial
  sentinel$n_pos_2[sentinel$sentinel == snp] <- n_pos
  sentinel$n_test_2[sentinel$sentinel == snp] <- n_neg
  sentinel$total_2[sentinel$sentinel == snp] <- total
  
}

