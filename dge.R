f = "~/projects/black_spruce/seqclean/all_ests.fa.clean_output/contig_member.counts"
counts = read.table(f, header=T, row.names=1)
env = data.frame(row.names=colnames(counts), 
                    tissue=c("cambium", "cambium", "needle", "needle"),
                    parent=c("mom","dad","mom","dad"))
design = model.matrix(~env$parent + env$tissue)
y = DGEList(counts=counts)
y = estimateGLMCommonDisp(y, design, verbose=T)
y = estimateGLMTrendedDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)
fit = glmFit(y, design)
lrt = glmLRT(fit)
print(topTags(lrt))