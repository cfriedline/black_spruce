library(edgeR)
f = "~/projects/black_spruce/seqclean/all_ests.fa.clean_output/contig_member.counts"
counts = read.table(f, header=T, row.names=1)
#counts = counts[sort(colnames(counts))]
# counts = counts[c("P32N", "P40N", "P32C", "P40C")]
tissue = factor()
parent = factor()
for (k in 1:ncol(counts)) {
    c = colnames(counts)[k]
    if (length(i <- grep("32", c))) {
        parent = c(parent, "mom")
    } else {
        parent = c(parent, "dad")
    }
    
    if (length(i <- grep("N", c))) {
        tissue = c(tissue, "needle")
    } else {
        tissue = c(tissue, "cambium")
    }
}
tissue = as.factor(tissue)
parent = as.factor(parent)
env = data.frame(row.names=colnames(counts), 
                    tissue=tissue,
                    parent=parent)
design = model.matrix(~parent + tissue, env)
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
print(tt[1:10,])
write.csv(tt, file="~/projects/black_spruce/toptags.csv")