#!/bin/sh
#SBATCH --job-name=bwtool
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/Breast/
#SBATCH --cpus-per-task=24  
#SBATCH --ntasks=1
#SBATCH --qos=bsc_ls
#SBATCH --time=24:00:00
#SBATCH --error=ANALYSIS/Hypertension/analysis/04_PhyloP/bwtool_submit.err
#SBATCH --output=ANALYSIS/Hypertension/analysis/04_PhyloP/bwtool_submit.out
# normal queue bsc_ls

##/gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/bed_oligos_testtiles.sorted.bed

/gpfs/projects/bsc83/utils/bwtool/bwtool summary /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/window_10_bp_before_snp.sorted.bed \
/gpfs/projects/bsc83/Data/PhyloP/hg19/placentalMammals.phyloP46way.bw ANALYSIS/Hypertension/analysis/04_PhyloP/hg19_summary_window_b.phyloP46placentaway.bed \
-header -with-sum -fill=0

#/gpfs/projects/bsc83/utils/bwtool/bwtool summary /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/data/design/window_10_bp_after_snp.sorted.bed \
#/gpfs/projects/bsc83/Data/PhyloP/hg19/placentalMammals.phyloP46way.bw ANALYSIS/Hypertension/analysis/04_PhyloP/hg19_summary_window_a.phyloP46placentaway.bed \
#-header -with-sum -fill=0

#placentalMammals.phyloP46way.bw
#hg19.100way.phyloP100way.bw

#bed_oligos_testtiles.sorted.bed
#snp_coord_master.hg19.simple.chr.sortedplus1.bed
#window_10_bp_before_snp.bed
#window_10_bp_after_snp.bed