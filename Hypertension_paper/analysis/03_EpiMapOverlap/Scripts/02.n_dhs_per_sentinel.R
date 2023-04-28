#### read data for numbe rof DHS per sentinel -----
sentinel_snps = read.table('../../../data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

DHS <- read.table('../../../analysis/03_EpiMapOverlap/DHS_overlap_with_info.txt',
sep = '\t', header = T)


### merge data ----
head(DHS)
merged_data <- merge(DHS, sentinel_snps, by.x='V8', by.y='change', all.y=T)

### read length block ----
length <- read.delim('../../../data/design/sentinel_snps_length_block.txt',
                     sep = '\t', header = T)
head(length)

###reg results ####
CM_results_s <- read.table('../../../analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                           header = T)

VSMC_results_s <- read.table('../../../analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                             header = T)

#### summarize sentinel and DHS ----
merged_data <- merged_data[,c('sentinel','V4','enhancers','promoters','dyadic')] %>% distinct()
merged_data <- merged_data[!is.na(merged_data$V4),]
head(merged_data)

number_dhs_per_sentinel <- as.data.frame(table(merged_data$sentinel))

### merge with length of region ------
length <- merge(length, number_dhs_per_sentinel, by.x='sentinel', by.y='Var1', all.x=T)
head(length)
length$Freq[is.na(length$Freq)] <- 0

sentinel_chr <- sentinel_snps[,c('sentinel','Chr')] %>% distinct()
sentinel_chr <- merge(sentinel_chr, length)
sentinel_chr$Chr <- paste0('chr',sentinel_chr$Chr)
write.table(sentinel_chr[,c('Chr','min','max','sentinel')], 
            '../../../data/design/sentinel_coordinates_block.bed',
            quote = F, col.names = F, row.names = F, sep = '\t')

library("ggpubr")
ggscatter(length, x = "length", y = "Freq",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Length LD block", ylab = "Nº DHS elem")


ggplot(length, aes(x=length, y = Freq)) +
  geom_point(color='darkgrey')+
  xlab('Length LD block')+
  ylab('Nº DHS elem')+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_classic()+
  geom_point() +
  geom_text_repel(
    data = length,
    aes(label = sentinel),
    size = 4)


##### desity -----
length$clust <- length$Freq/(length$length/1000000)

ggplot(length[length$length>1,], aes(x=clust)) + 
  geom_density(adjust = 1/3, color='darkgrey')+
  xlab('N of LD SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

length$clust <- as.numeric(length$clust)
length$length <- as.numeric(length$length)
vector <- length$clust[length$length>1]
hist(log10(vector), col = 'steelblue', breaks = 60, main = 'Nº of DHS / Mb', xlab = 'log10(DHS/Mbase)') + 
  abline(v=median(log10(vector)), col='red')

### density regulatory 

#### repeat density without first and last ######
sentinel_snps$DiffCM <- 'Not Diff'
sentinel_snps$DiffVSMC <- 'Not Diff'
sentinel_snps$DiffCM[sentinel_snps$snp %in% CM_results_s$snp_info] <- 'Diff'
sentinel_snps$DiffVSMC[sentinel_snps$snp %in% VSMC_results_s$snp_info] <- 'Diff'
table(sentinel_snps$DiffCM)

head(length)

head(sentinel_snps)

number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$DiffCM == 'Diff']))
number_snps_all <- as.data.frame(table(sentinel_snps$sentinel))
number_snps_all <- merge(number_snps_all, length, by.x='Var1', by.y='sentinel')
number_snps_all$clust <- number_snps_all$Freq/((number_snps_all$length)/1000000)

number_snps <- merge(number_snps, length, by.x='Var1', by.y='sentinel')
number_snps <- merge(number_snps, number_snps_all, by='Var1')

number_snps$ratio <- (number_snps$Freq.x/number_snps$clust)
number_snps$ratio_all <- (number_snps$Freq.y/((number_snps$length.y/1000000)))

merged_data <- merged_data[,c('sentinel','V4','enhancers','promoters','dyadic')] %>% distinct()
merged_data <- merged_data[!is.na(merged_data$V4),]
head(merged_data)

number_dhs_per_sentinel <- as.data.frame(table(merged_data$sentinel))
number_snps <- merge(number_snps, number_dhs_per_sentinel, by='Var1')

number_snps$reg_len <- (number_snps$Freq.x/((number_snps$length.y)/1000000))
number_snps$dhs_len <- (number_snps$Freq/((number_snps$length.y)/1000000))
                        
library("ggpubr")
ggscatter(number_snps[number_snps$length.x>300,], x = "reg_len", y = "dhs_len",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "density regulatory variants (reg var / Kb)", 
          ylab = "density DHS elements")

