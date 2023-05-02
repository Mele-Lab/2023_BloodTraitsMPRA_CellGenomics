 #!/bin/sh

#BSUB -cwd /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/
#BSUB -J bedtools
#BSUB -e bedtools.err
#BSUB -o bedtools.out
#BSUB -q bsc_ls
#BSUB -W 24:00
#BSUB -n 1
#BSUB -M 7000
#tot 378


module purge && module load java/1.8.0u66 intel/2017.4 impi/2017.4 MKL/2017.4 gcc/7.2.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip SAMTOOLS/1.9 subread GATK/4.1.2.0 BEDTOOLS/2.25.0 R/3.6.3

data='/gpfs/projects/bsc83/Data/gene_annotation/gencode/release_34/'
fantom='/gpfs/projects/bsc83/Data/FANTOM/'

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.chr.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_34/genes_sorted.hg19.bed \
#     -sorted > SNPs_master_list_genes_overlap.bed

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_34/exon_sorted.hg19.bed \
#     -sorted > SNPs_master_list_exons_overlap.bed

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.chr.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_34/intergenic_sorted.hg19.bed \
#     -sorted > SNPs_master_list_intergenic_overlap.bed

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/gene_annotation/gencode/release_34/gencode.v34.introns.hg19.sorted.bed \
#     -sorted > SNPs_master_list_introns_overlap.bed

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.chr.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/FANTOM/human_permissive_enhancers_phase_1_and_2.sorted.bed \
#     -sorted > SNPs_master_list_enhancers_overlap.bed

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.chr.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/FANTOM/TSS_human.sorted.bed \
#     -sorted > SNPs_master_list_TSS_overlap.bed

# bedtools intersect -wa -wb \
#     -a snp_coord_master.hg38.simple.chr.sorted.bed \
#     -b /gpfs/projects/bsc83/Data/FANTOM/hg19.cage_peak_phase1and2combined_coord.sorted.bed \
#     -sorted > SNPs_master_list_peaks_overlap.bed

bedtools intersect -wa -wb \
    -a snp_coord_master.hg19.simple.chr.sorted.bed \
    -b ${data}genes_sorted.hg19.bed ${data}exon_sorted.hg19.chr.bed ${data}intergenic_sorted.hg19.bed ${data}gencode.v34.introns.hg19.chr.sorted.bed ${fantom}human_permissive_enhancers_phase_1_and_2.sorted.bed ${fantom}TSS_human.sorted.bed ${fantom}hg19.cage_peak_phase1and2combined_coord.sorted.bed\
    -names Gene Exon Intergenic Intronic Enhancer TSS CAGE_peak\
    -sorted > all_results_SNPS_overlap_genome.bed

