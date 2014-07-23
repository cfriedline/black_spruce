library(topGO)
countfile = "seqclean/all_ests.fa.clean_output/contig_member.counts"
counts = read.table(countfile, header=T, row.names=1, sep="\t")
gofile = "all_annot_GOs_20140722_2254.txt"
go = read.table(gofile, header=F, sep="\t")
library(ALL)
data(ALL)
data(geneList)