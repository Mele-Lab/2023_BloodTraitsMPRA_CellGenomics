#!/bin/bash
#Code to collect the QTL colocalization results after they have all been generated.

#set variables for date and time
today=`date +%F`
now=`date +%T`

pwd

#load config file
source /gpfs/projects/bsc83/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/colocQuial_res/qtl_config.sh
source $setup_config_sh

#write header to the summary results file
#nsnps   hit1    hit2    PP.H0.abf       PP.H1.abf       PP.H2.abf       PP.H3.abf       PP.H4.abf       idx1    idx2
printf "SNP\tGene\tGeneID-Tissue\tTrait\tnsnps\thit1\thit2\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tidx1\tidx2\n" > $trait"_"$qtlType"_susie_results_all_summary_"$today"_"$now".txt" 

#go through each lead SNP directory
cat $leadSNPsFilePath | while read line
do 
	echo $line

	#parse the SNP
	dir=`echo $line | cut -f 1 -d " "`
	echo $dir

    cd  $dir

    #for each GTEx coloc output in the directory
    for colocOut in *susie_results_summary.txt
    do
        echo $colocOut
        colocfilename=`basename $colocOut`
        #echo $colocfilename

        #save a "trait" string for formatting
        trait_str="_"$trait"_"

        #grab the information we need from the locus coloc results file
        sed "s/^/$colocfilename /" $colocOut | sed "s/^/$dir /" | sed "s,$trait_str, $trait," | sed "s/susie_results_summary.txt//" | sed "s/_ENSG/ ENSG/"| tr " " "\t" | tail -n+2 >> ../$trait"_"$qtlType"_susie_results_all_summary_"$today"_"$now".txt"

    done

    cd ..

done