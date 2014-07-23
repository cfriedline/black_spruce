f = "~/projects/black_spruce/seqclean/all_ests.fa.clean_output/contig_member.counts"
counts = read.table(f, header=T, row.names=1)
tissue = c("cambium", "cambium", "needle", "needle")
parent = c("mom","dad","mom","dad")
env = data.frame(row.names=colnames(counts), 
                    tissue=tissue,
                    parent=parent)
design = model.matrix(~parent + tissue)
rownames(design)=colnames(counts)
print(design)
y = DGEList(counts=counts)
y = calcNormFactors(y)
y = estimateGLMCommonDisp(y, design, verbose=T)
y = estimateGLMTrendedDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)
fit = glmFit(y, design)
lrt = glmLRT(fit)
tt = topTags(lrt, n=nrow(counts))
write.csv(tt, file="~/projects/black_spruce/toptags.csv")