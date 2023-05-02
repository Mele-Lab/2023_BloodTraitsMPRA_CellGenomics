cription in all 4 raQTL files:
1) 	"chr"	chromosome on which SNP lies
2) 	"SNP_ID"	SNP identifyer
3) 	"SNPabspos"	genomic position of SNP (hg19)
4) 	"ref.element.count"	number of DNA fragments that contain the REF allele for this SNP
5) 	"alt.element.count"	number of DNA fragments that contain the VAR allele for this SNP

In k562.matched.control.LP190708.txt.gz/k562.sign.id.LP190708.txt.gz:
6) 	"k562.ref.mean"	mean expression in K562 of all DNA fragments that contain the REF allele for this SNP
7) 	"k562.alt.mean"	mean expression in K562 of all DNA fragments that contain the VAR allele for this SNP
8) 	"k562.wilcox.p.value"	P-value by a Wilcoxon rank sum test in K562, REF DNA fragments versus VAR DNA fragments
9) 	"k562.wilcox.p.value.random"	same as "k562.wilcox.p.value", but after shuffeling the expression values of the DNA fragments

In hepg2.matched.control.LP190708.txt.gz/hepg2.sign.id.LP190708.txt.gz:
6) 	"hepg2.ref.mean"	same as "k562.ref.mean", but for HepG2
7) 	"hepg2.alt.mean"	same as "k562.alt.mean", but for HepG2
8) 	"hepg2.wilcox.p.value"	same as "k562.wilcox.p.value", but for HepG2
9) 	"hepg2.wilcox.p.value.random"	same as "k562.wilcox.p.value.random", but for HepG2

In all 4 raQTL files:
10)     "ref"	reference allele for this SNP
11)     "alt"	alternative allele for this SNP
