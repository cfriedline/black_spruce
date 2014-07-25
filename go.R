library(topGO)
rm(list=ls())
count_file = "seqclean/all_ests.fa.clean_output/contig_member.counts"
counts = read.table(count_file, header=T, row.names=1, sep="\t")

gene_names = rownames(counts)

cambium_sums = rowSums(counts[,c("P32C","P40C")])
needle_sums = rowSums(counts[,c("P32N","P40N")])

cambium_interesting = cambium_sums[cambium_sums>1]
needle_interesting = needle_sums[needle_sums>1]

go_file = "all_annot_GOs_20140722_2254.txt"
go = read.table(go_file, header=F, sep="\t")
gene_id_2go <- readMappings(file="all_annot_GOs_20140722_2254.txt_topGO.txt")

interesting = list()
interesting$cambium = cambium_interesting
interesting$needle = needle_interesting
godata = list()
onts = c("BP","CC", "MF")
sigs = list()
for (i in 1:length(onts)) {
    for (j in 1:length(interesting)) {
        interest = interesting[[j]]
        gene_list <- factor(as.integer(gene_names %in% names(interest)))
        names(gene_list) <- gene_names
        
        GOdata = new("topGOdata",
                     description=paste(names(interesting)[j], onts[i], sep="-"),
                     ontology = onts[i], 
                     allGenes = gene_list, 
                     annot = annFUN.gene2GO, 
                     gene2GO = gene_id_2go,
                     nodeSize=5)
        print(GOdata)
        godata = append(godata, GOdata)
        test.stat = runTest(GOdata, algorithm = "classic", statistic = "fisher")
        sigs = append(sigs, test.stat)
        printGraph(GOdata, 
                   test.stat, 
                   firstSigNodes = 5, 
                   fn.prefix = paste("tGO", "for", description(GOdata)), 
                   useInfo = "all", 
                   pdfSW = TRUE)
        
        gt = GenTable(GOdata, classic=test.stat, topNodes=20)
        write.table(gt, file=paste(description(GOdata), ".txt", sep=""), row.names=F)
    }    
}

