#### install and load packages ------
library(coloc)
library(data.table)
library(rjson)
library(tidyverse)
library(ieugwasr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(glue)
library(vroom)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(Homo.sapiens)
library(ggbio)
library(egg)


#Run runsusie on dataset 1, storing the results
#Run runsusie on dataset 2, storing the results
#Run coloc.susie on the two outputs from above

### coloc must be run on a single region

#One common mistake is to use the standard error of beta in place of the variance of beta. 
#If your dataset provides the standard error, simply square it to get the variance.

#Read in the arguments from the config file
source("QTL_config.R")
source(setup_config_R)

#Print all of the config file settings to screen or the stnd out file
print("trait")
print(trait)
print("trait_file_header_info")
print(traitFilePath)
print(trait_A1col)
print(trait_A2col)
print(trait_SNPcol)
print(trait_CHRcol)
print(trait_BPcol)
print(trait_Pcol)
print(trait_Ncol)
print(trait_MAFcol)
print("traitType:")
print(traitType)
print(traitProp)
print("Locus Info:")
print(chrom)
print(colocStart)
print(colocStop)
print(lead_SNP)
print("QTL type:")
print(qtlType)

## testing region rs6026739	57701709	57767727	66019	20

# Set up QTL variables for either eQTL or sQTL
qtlType == "eqtl"
  sig_qtl_tabix_dir = eQTL_sig_qtl_tabix_dir
  sig_geneID_col = eQTL_sig_geneID_col
  all_qtl_tabix_dir = eQTL_all_qtl_tabix_dir
  all_header = eQTL_all_header
  all_geneID = eQTL_all_geneID
  all_chrom = eQTL_all_chrom
  all_chromEnd = eQTL_all_chromEnd
  all_pvalue = eQTL_all_pvalue
  QTL_tissue_table = eQTL_tissue_table
  
#### First get the locus from GWAS data
#minimum_data=D1[c("beta","varbeta","snp","position","type","sdY")]
#nosdY_data=D1[c("beta","varbeta","snp","position","type","N","MAF")]
  #type=quant
#read in the tissue specific eQLT summary file with the file names added
trait_region = fread(file=traitFilePath, sep="\t", header=TRUE)
print("trait input file successfully loaded")

trait_region <- trait_region[trait_region$CHR == chrom & trait_region$POS >= colocStart & trait_region$POS <= colocStop,]
head(trait_region)

trait_A1col_str = paste(trait_A1col,"_trait",sep="")
trait_A2col_str = paste(trait_A2col,"_trait",sep="")
trait_SNPcol_str = paste(trait_SNPcol,"_trait",sep="")
trait_MAFcol_str = paste(trait_MAFcol,"_trait",sep="")

colnames(trait_region)[colnames(trait_region)== trait_A1col] <- trait_A1col_str
colnames(trait_region)[colnames(trait_region)== trait_A2col] <- trait_A2col_str
colnames(trait_region)[colnames(trait_region)== trait_SNPcol] <- trait_SNPcol_str
colnames(trait_region)[colnames(trait_region)== trait_MAFcol] <- trait_MAFcol_str

#grab rs numbers if present
trait_region_rs = trait_region[[trait_SNPcol_str]]
head(trait_region_rs, 3)

#use hash tables to find chromosome positions
hash_table_file = paste0(hash_table_dir, "chr_", chrom, "_snp151_hash_table.json")
print(hash_table_file)

hash_table <- fromJSON(file = hash_table_file) 
head(hash_table, 3)
print("hash table loaded")

trait_chrom_pos = hash_table[trait_region_rs]

#read in tissue table from setup_config.R
tissueTable = read.table(file=QTL_tissue_table, sep=",", header=TRUE)
print(QTL_tissue_table)
print(head(tissueTable))
colnames(tissueTable) <- c('Tissue', 'NumberRNASeqandGTSamples', 'NumberRNASeqSamples',
                           'Number_of_eGenes','allPairsTabixFilename','sigPairsTabixFilename')

#Convert to strings from factors
tissueTable$Tissue = as.character(tissueTable$Tissue)
tissueTable$allPairsTabixFilename = as.character(tissueTable$allPairsTabixFilename)
tissueTable$sigPairsTabixFilename = as.character(tissueTable$sigPairsTabixFilename)

#Format Tissue names
tissueTable_tissue_noSpace = gsub("\\(","",tissueTable$Tissue)
tissueTable_tissue_noSpace = gsub("\\)","",tissueTable_tissue_noSpace)
tissueTable_tissue_noSpace = gsub("[[:space:]]","_",tissueTable_tissue_noSpace)
tissueTable_tissue_noSpace = gsub("-_", "", tissueTable_tissue_noSpace)
tissueTable$Tissue = tissueTable_tissue_noSpace

#create csv file from significant pair files
lead_SNP_pos = hash_table[[lead_SNP]]
print(lead_SNP_pos)
#convert format for tabix
lead_SNP_pos_tabix_with_chr = paste0(gsub("_",":",lead_SNP_pos), "-", gsub("^.*?_","",lead_SNP_pos))
lead_SNP_pos_tabix_without_chr = paste0(gsub("chr", "", gsub("_",":",lead_SNP_pos)), "-", gsub("^.*?_","",lead_SNP_pos))

for (i in 1:nrow(tissueTable)) {
  sigpair_filename = tissueTable$sigPairsTabixFilename[i]
  if (is.na(sigpair_filename)) {
    print(paste(tissueTable$Tissue[i], "is not available in the significant pairs files"))
    next
  }
  
  file = paste0(sig_qtl_tabix_dir, "/", sigpair_filename) 
  
  system(paste("tabix", file, lead_SNP_pos_tabix_with_chr, ">", paste0(lead_SNP, "_temp.csv"))) 
  system(paste("tabix", file, lead_SNP_pos_tabix_without_chr, ">>", paste0(lead_SNP, "_temp.csv"))) 
  
  system(paste0("sed -i \"s/$/\t", tissueTable$Tissue[i], "/\" ", lead_SNP, "_temp.csv"))
  system(paste0("cat ", lead_SNP, "_temp.csv >> ", lead_SNP, ".csv"))
}

#read in the csv file of eGene-Tissue pairs 
eGenes = tryCatch({
  read.table(file=paste0(lead_SNP, ".csv"), sep="\t", stringsAsFactors=FALSE)
}, error = function(err) {
  print(paste(lead_SNP, "is not a significant", qtlType, "in any tissue in the", qtlType, "dataset."))
  quit(status=0)
})

#loop through the eGene-Tissue pairs in eGenes and prep running COLOC
library(clusterProfiler)
for(i in 1:nrow(eGenes)){
  geneID <- eGenes[i, sig_geneID_col]
  
  geneID_noDOT <-  gsub("\\..*","", geneID)
  
  geneSymbol <- tryCatch({
    
    bitr(geneID_noDOT, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
    
  }, error = function(err){
    
    print("ENSEMBL ID could not be converted to HGNC Symbol")
    print(paste("geneSymbol will be set to the ENSEMBL ID",geneID_noDOT))
    
    return(geneID_noDOT)
    
  })
  
  tissue <- eGenes[i, length(eGenes)]
  
  print(geneID)
  print(geneSymbol)
  
  if (length(geneSymbol) > 1) {
    geneSymbol <- geneSymbol[1]
    print(paste("Using first gene symbol: ", geneSymbol))
  }
  
  print(tissue)
  
  # find the all pair file that contains the tissue of interest
  
  tissueLine <- tissueTable[tissueTable$Tissue == tissue,] 
  allpair_filename <- tissueLine$allPairsTabixFilename
  
  #if the Filename field is NA then skip this eGene-Tissue pair
  if (is.na(allpair_filename)){
    print(tissue)
    print("This tissue is not available in the all pairs files currently")
    next
  }
  
  tabix_allpair_path = paste0(all_qtl_tabix_dir, allpair_filename)
  qtl_N <- tissueLine$NumberRNASeqandGTSamples
  
  #parentheses are causing issues too
  tissue_noSpace = gsub("\\(","",tissue)
  tissue_noSpace = gsub("\\)","",tissue_noSpace)
  # The tissue names have any whitespace in them and we want to use these in the output file names so replace " " with "_"
  tissue_noSpace = gsub("[[:space:]]","_",tissue_noSpace)
  
  # make a eGene-Tissue and trait prefix for file names
  if (qtlType == "eqtl") {
    out_prefix = paste(geneSymbol,geneID,tissue_noSpace,trait,sep="_")
  }	
  #run liftover on colocStart and colocStop if in HG19
  if (build == "hg19") {	
    repeat {
      print("running liftOver")
      bed_liftover = data.frame("chr" = c(paste0("chr", chrom), paste0("chr", chrom)), "bp1" = c(colocStart - 1, colocStop - 1), "bp2" = c(colocStart, colocStop)) 
      write.table(bed_liftover,file=paste0("temp_hg19.bed"),sep="\t",quote = FALSE,row.names=FALSE,col.names=FALSE)
      
      #generate liftOver command
      liftOver_command = paste( "liftOver temp_hg19.bed", liftOver_chain, "temp_hg38.bed temp_hg19.unmapped -bedPlus=3 -tab", sep=" ")
      system(liftOver_command)
      
      hg38_positions = as.data.frame(read.table("temp_hg38.bed", header=FALSE, sep="\t")) 
      
      if (!is.na(hg38_positions[2,3])) {
        break
      } else {
        unmapped = fread(file="temp_hg19.unmapped", sep='\t', header=FALSE)
        if (unmapped[1,3] == colocStart) {
          print(paste("Could not map to hg38:", colocStart))
          colocStart = colocStart + 5000
        } else if (unmapped[1,3] == colocStop) {
          print(paste("Could not map to hg38:", colocStop))
          colocStop = colocStop - 5000
        }
        print(paste("Trying new region:", colocStart, "-", colocStop)) 
        
        if (colocStart >= colocStop) {
          print("Could not perform liftover")
          quit(status=0)
        }
      }
    }
  }
  
  print("Grabbing the all pairs data")    
  #Use tabix to grab data, try both with and without "chr"
  eGeneTissueInputFile = paste(geneSymbol,tissue_noSpace,chrom,colocStart,colocStop,".txt", sep="_")
  if (build == "hg38") {
    system(paste0("tabix ", tabix_allpair_path, " ", chrom, ":", colocStart, "-", colocStop, " >> ", eGeneTissueInputFile))
    system(paste0("tabix ", tabix_allpair_path, " chr", chrom, ":", colocStart, "-", colocStop, " >> ", eGeneTissueInputFile))
  } else if (build == "hg19") {
    system(paste0("tabix ", tabix_allpair_path, " ", chrom, ":", hg38_positions[1,3], "-", hg38_positions[2,3], " >> ", eGeneTissueInputFile))    
    system(paste0("tabix ", tabix_allpair_path, " chr", chrom, ":", hg38_positions[1,3], "-", hg38_positions[2,3], " >> ", eGeneTissueInputFile))    
  } else {
    print("ERROR: Please specify build: \"hg19\" or \"hg38\"")
    quit()
  }
  
  print("reading the all pairs data into R")
  #read the file we just generated from the grep command into R
  eGeneTissueInput = fread(file = eGeneTissueInputFile, header=FALSE)
  
  #check if tabix result is empty
  print(dim(eGeneTissueInput))
  if (dim(eGeneTissueInput)[1] == 0) {
    print("Did not find gene in all pairs data")
    if (qtlType == "eqtl") {
      quit(status=0)
    } else if (qtlType == "sqtl") {
      next
    }
  }	
  
  #remove the eGeneTissueInputFile after has been read into R to save disk space
  system(paste0("rm ", eGeneTissueInputFile, sep=""))
  
  colnames(eGeneTissueInput) = all_header
  eGeneTissueInput$eGeneID <- gsub('\\..*','', eGeneTissueInput$eGeneID)
  
  print("Filtering on the geneID")
  #filter for the geneID of interest
  eGeneTissue_region = eGeneTissueInput[eGeneTissueInput[[all_geneID]] == geneID_noDOT,]
  
  if (nrow(eGeneTissue_region) == 0) {
    print("Warning: There was not an exact match on Ensembl ID. Likely this is due to a GTEX version  update.")
    #make a string that removes the everything after in the geneID
    noDecimalGeneID = gsub("\\..*","",geneID)
    
    #recreated eGeneTissue_region
    eGeneTissue_region = eGeneTissueInput
    
    #grep the simplified Ensembl ID
    possible_Ensembl_gene_lines <- eGeneTissue_region[grepl(noDecimalGeneID, eGeneTissue_region[[all_geneID]]),]
    
    #check to make sure there is just one other Ensembl ID
    possible_Ensembl_genes <- unique(possible_Ensembl_gene_lines[[all_geneID]])
    
    if(length(possible_Ensembl_genes) == 1){
      print("Found a unique Ensembl ID so this analysis will continue using the Ensembl ID:")
      print(possible_Ensembl_genes)
      
      #this will be a single Ensembl ID string
      geneID = possible_Ensembl_genes
      eGeneTissue_region = eGeneTissue_region[eGeneTissue_region[[all_geneID]] == geneID,]
    } else {
      print("The Ensembl ID from your GTEx csv was not able to be reliably mapped to an Enseml ID in the GTEx database, so this gene will be skipped:")
      print(geneID)
      next
    }    
  }
  
  if (qtlType == "eqtl") {
    #Keep only the columns that are needed
    eGeneTissue_region <- eGeneTissue_region %>% dplyr::select(all_of(eQTL_all_chrom), all_of(eQTL_all_chromEnd), all_of(eQTL_all_geneID), all_of(eQTL_all_pvalue), all_of('maf'), all_of('ma_samples'), all_of('slope'),all_of('slope_se'))
  } else if (qtlType == "sqtl") {
    #Keep only the columns that are needed
    eGeneTissue_region <- eGeneTissue_region %>% dplyr::select(all_of(sQTL_all_chrom), all_of(sQTL_all_chromEnd), all_of(sQTL_all_geneID), all_of(sQTL_all_pvalue), all_of(sQTL_all_intron_chr), all_of(sQTL_all_intron_bp_first), all_of(sQTL_all_intron_bp_end), all_of(sQTL_all_intron_clu))
  }
  
  print("adding rs numbers to the QTL data")
  #add rs genegene,,numbers to the eGeneTissue_region DF
  
  #Add "chr" prefix to chromosome number if not present
  eGeneTissue_region[[all_chrom]] = paste0("chr", sub("chr", "", eGeneTissue_region[[all_chrom]]))
  
  #make chromosome_position column for merging
  eGeneTissue_region$chromosome_position <- paste(eGeneTissue_region[[all_chrom]],eGeneTissue_region[[all_chromEnd]]+1,sep="_")
  
  #create data frame with rs numbers associated with chromosome_position and add to eGeneTissue region
  uniqID_DF = as.data.frame(as.data.frame(unlist(trait_chrom_pos)))
  uniqID_DF$SNP <- rownames(uniqID_DF)
  colnames(uniqID_DF) <- c("chromosome_position", "SNP")
  
  eGeneTissue_region = merge(eGeneTissue_region, uniqID_DF, by = "chromosome_position")  
  
  ################################ eQTL colocalization and RA Plots ################################ 
  
  if (qtlType == "eqtl") {
    print("merging the trait and eqtl data on unique ID")
    #merge the trait and eGeneTissue region DFs on rs numbers
    colocInputFile = merge(eGeneTissue_region, trait_region, by.x="SNP", by.y=trait_SNPcol_str)
    
    #remove any NAs
    colocInputFile = colocInputFile[complete.cases(colocInputFile), ]
    
    #check for 0s in the trait_Pcol
    if (0 %in% colocInputFile[[trait_Pcol]]){
      
      print("WARNING: THERE ARE SNPS WITH P-VALUES OF 0 AT THIS LOCUS. These SNPs have been removed for the Colocalization anlysis and may lead to unusual regional association plots")
      
      #remove SNPs who's trait P-value is 0 
      colocInputFile = colocInputFile[colocInputFile[[trait_Pcol]] != 0,]
      
    }
    
    #write colocInputFile to file for making locus zoom plots
    colocInputFile_outputStr = paste(out_prefix,"coloc_input_data.txt",sep="_")
    write.table(colocInputFile, file= colocInputFile_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
    
    print("Running coloc")
    #run coloc
    if (traitType == "cc"){
      coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], N=colocInputFile[[trait_Ncol]], type=traitType, s=traitProp), dataset2=list(pvalues=colocInputFile[[eQTL_all_pvalue]], N=qtl_N, type="quant"),MAF=colocInputFile[[trait_MAFcol_str]])
    } else {
      coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile[[trait_Pcol]], 
                                               N=colocInputFile[[trait_Ncol]], type=traitType), 
                                 dataset2=list(pvalues=colocInputFile[[eQTL_all_pvalue]], N=qtl_N, type="quant"),
                                 MAF=colocInputFile[['maf']])
    }
    
    #prepare useful outputs
    coloc_results_summary = coloc_results$summary
    coloc_results_full = coloc_results$results
    
    #calculate pp4 / pp3 + pp4
    PP3andPP4 = coloc_results_summary[5] + coloc_results_summary[6]
    
    pp4_conditional = coloc_results_summary[6] / PP3andPP4
    
    #prep coloc output strings
    coloc_results_summary_outputStr = paste(out_prefix,"coloc_results_summary.txt",sep="_")
    coloc_results_full_outputStr = paste(out_prefix,"coloc_results_full.txt",sep="_")
    coloc_results_pp4_cond_outputStr = paste(out_prefix,"coloc_results_pp4_cond.txt",sep="_")
    
    #write to file
    write.table(coloc_results_summary, file=coloc_results_summary_outputStr, sep="\t", row.names=TRUE, quote=FALSE)
    write.table(coloc_results_full, file=coloc_results_full_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
    write.table(pp4_conditional, file=coloc_results_pp4_cond_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
    
    ### now run sussie
    ## read LD matrix
    tryCatch({
      
      LD_matrix <- read.delim(paste0(path_ld,lead_SNP,'_LD_matrix.txt'), row.names = 1)
      rownames(LD_matrix) <- gsub('_.*','', rownames(LD_matrix))
      colnames(LD_matrix) <- gsub('_.*','', colnames(LD_matrix))
      colocInputFile <- colocInputFile[colocInputFile$SNP %in% colnames(LD_matrix),]
      LD_matrix <- LD_matrix[colocInputFile$SNP, colocInputFile$SNP]
  
      d1 = list( N=min(colocInputFile[[trait_Ncol]]), 
                type=traitType, beta=setNames(colocInputFile[['BETA']], colocInputFile$SNP), 
                varbeta=setNames(colocInputFile[['SE']]^2, colocInputFile$SNP),
                MAF=setNames(colocInputFile$maf, colocInputFile$SNP),LD=as.matrix(LD_matrix),
                snp=colocInputFile$SNP, position=colocInputFile$chromEnd_b38)
      
      S1=runsusie(d1)
      print(summary(S1))
      
      d2 = list( N=max(colocInputFile$ma_samples), 
                 type=traitType, beta=setNames(colocInputFile$slope, colocInputFile$SNP), 
                 varbeta=setNames(colocInputFile[['slope_se']]^2, colocInputFile$SNP),
                 sdY=1,LD=as.matrix(LD_matrix),MAF=setNames(colocInputFile$maf, colocInputFile$SNP),
                 snp=colocInputFile$SNP, position=colocInputFile$chromEnd_b38)
      
      S2=runsusie(d2)
      print(summary(S2))
      
      if(requireNamespace("susieR",quietly=TRUE)) {
        susie.res=coloc.susie(S1,S2)
        print(susie.res$summary)
      }
      
      ## save susie results 
      #prepare useful outputs
      susie_results_summary = susie.res$summary
      susie_results_full = susie.res$results
      
      #calculate pp4 / pp3 + pp4
      PP3andPP4 = susie_results_summary$PP.H3.abf + susie_results_summary$PP.H4.abf
      
      pp4_conditional = susie_results_summary$PP.H4.abf / PP3andPP4
      
      #prep coloc output strings
      susie_results_summary_outputStr = paste(out_prefix,"susie_results_summary.txt",sep="_")
      susie_results_full_outputStr = paste(out_prefix,"susie_results_full.txt",sep="_")
      susie_results_pp4_cond_outputStr = paste(out_prefix,"susie_results_pp4_cond.txt",sep="_")
      
      #write to file
      write.table(susie_results_summary, file=susie_results_summary_outputStr, sep="\t", row.names=TRUE, quote=FALSE)
      write.table(susie_results_full, file=susie_results_full_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
      write.table(pp4_conditional, file=susie_results_pp4_cond_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
        
    }, error = function(err){
      
      print("susie can't run in these pair")
      

    })
    

    
    #generate regional association plot
    
    # leadSNP_DF = colocInputFile#[colocInputFile$SNP == lead_SNP,]
    # leadSNP_DF[[eQTL_all_chrom]] = as.integer(gsub('[a-zA-Z]', '', leadSNP_DF[[eQTL_all_chrom]])) 
    # leadSNP_DF = leadSNP_DF %>% dplyr::select(SNP, all_of(eQTL_all_chrom), all_of(trait_BPcol), all_of(eQTL_all_pvalue), all_of(trait_Pcol))
    # 
    # library(locuscomparer)
    # 
    # gwas <- as.data.frame(leadSNP_DF[,c('SNP','P')])
    # colnames(gwas) <- c('rsid','pval')
    # eqtl <- as.data.frame(leadSNP_DF[,c('SNP','pvalue_eQTL')])
    # colnames(eqtl) <- c('rsid','pval')
    # 
    # locuscomp_gwas_str = paste(out_prefix,"gwas_locuscomp.txt",sep="_")
    # locuscomp_eqtl_str = paste(out_prefix,"eqtl_locuscomp.txt",sep="_")
    # 
    # write.table(gwas, file=locuscomp_gwas_str, sep="\t", row.names=F, quote=F, col.names = T)
    # write.table(eqtl, file=locuscomp_eqtl_str, sep="\t", row.names=F, quote=F, col.names = T)
    # 
    # RA_plot <- locuscompare(in_fn1 = locuscomp_gwas_str, 
    #              in_fn2 = locuscomp_eqtl_str, 
    #              title = 'GWAS', title2 = 'eQTL')
    # 
    # pdf(file = paste0(out_prefix, "_", geneSymbol, "_", tissue,".pdf"), paper = 'USr', width = 15, height = 20)  
    # print(RA_plot)
    # dev.off()
    
    
  } 
  }


  