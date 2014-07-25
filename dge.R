library(edgeR)
rm(list=ls())
f = "~/projects/black_spruce/seqclean/all_ests.fa.clean_output/contig_member.counts"
counts = read.table(f, header=T, row.names=1)
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
DENameList = list()
for (i in 1:length(levels(tissue))) {
    ref = levels(tissue)[i]
    tissue = relevel(tissue, ref=ref)
    design = model.matrix(~tissue+parent)
    rownames(design)=colnames(counts)
    print(design)

    y = DGEList(counts=counts,group=tissue)
    y = calcNormFactors(y)
    plotMDS(y)
    y = estimateGLMCommonDisp(y, design)
    y = estimateGLMTrendedDisp(y,design)
    y = estimateGLMTagwiseDisp(y,design)
    plotBCV(y)
    fit = glmFit(y, design)
    lrt = glmLRT(fit,coef=2)
    tt = topTags(lrt, n=150)
    print(tt[1:10,])
    print(summary(dt <- decideTestsDGE(lrt)))
    isDE <- as.logical(dt)
    DEnames <- rownames(y)[isDE]
    DENameList[[i]] = DEnames
    write.csv(tt[DEnames,], file=sprintf("~/projects/black_spruce/toptags_%s.csv",ref))
}
DEunion = union(DENameList[[1]], DENameList[[2]])
DEintersect = intersect(DENameList[[1]], DENameList[[2]])
DEdiff = setdiff(DENameList[[1]], DENameList[[2]])
