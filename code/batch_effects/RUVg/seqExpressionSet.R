plot.pca <- function(dataToPca, speciesCol, namesPch, sampleNames){
    pca <- prcomp(t(assay(se)), scale=T, center=T)
    pca.var <- pca$sdev^2/sum(pca$sdev^2)
    for (i in 1:9){
        pca_axis1=i
        pca_axis2=i+1
        plot(pca$x[,pca_axis1], pca$x[,pca_axis2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
        points(pca$x[,pca_axis1][which(colnames(lcpm) %in% allreps)], pca$x[,pca_axis2][which(colnames(lcpm) %in% allreps)], col="black", pch=8, cex=2)
        text(pca$x[,pca_axis1][which(samplenames %in% allreps)], pca$x[,pca_axis2][which(samplenames %in% allreps)], labels=samplenames[which(samplenames %in% allreps)], pos=3)
        #legend(legend=unique(sampleNames), col=unique(speciesCol), pch=unique(namesPch), x="bottomright", cex=0.6)
        legend(legend=unique(sampleNames), pch=16, x="bottomright", col=unique(speciesCol), cex=0.6, title=name, border=F, bty="n")
        legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 15, title="Batch", cex=0.6, border=F, bty="n")
    }

    return(pca)
}

speciesCol=as.numeric(get(name)), namesPch=y$samples$batch + 15
rv = rowVars(exprs(x)) 
select = order(rv, decreasing=TRUE)[seq_len(ntop)] 
pca = prcomp(t(exprs(x)[select,]))



set2=makeSummarizedExperimentFromExpressionSet(set1)
function (object, intgroup = "Island", ntop = 500, returnData = FALSE) 
{
    rv <- rowVars(assay(se))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
    pca2 <- prcomp(t(assay(se)[select, ]))
    percentVar <- pca2$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(set2)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(set2)[, "Island", drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        colData(set2)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, name = colnames(set2))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }

    plot(pca$x[,1], pca$x[,2], col=speciesCol, pch=namesPch, cex=2, xlab=paste0("PC", pca_axis1, " (", round(pca.var[pca_axis1]*100, digits=2), "% of variance)"), ylab=paste0("PC", pca_axis2, " (", round(pca.var[pca_axis2]*100, digits=2), "% of variance)", sep=""), main=name)
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed() +
        annotate("text", x=pca$x[,pca_axis1][which(samplenames %in% allreps)], y=pca$x[,pca_axis2][which(samplenames %in% allreps)], label= samplenames[which(samplenames %in% allreps)])
}



which(apply(t(assay(se)), 2, var)==0)


#######

y <- gsePathway(geneList, nPerm=10000,pvalueCutoff=0.2,pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(y)
head(res)
