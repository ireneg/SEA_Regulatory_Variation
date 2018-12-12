# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

# load dependencies: libraries, human count data, plasmodium data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/RUVs")

# load extra libraries
library(NineteenEightyR)
library(ReactomePA)

# DE analysis using RUVs  ----------------------------------------------------------------------------------------------------

# set up colors
colors <- electronic_night(n=5)
lightcolors=c("thistle", "darkblue")

# RUVSeq works with a genes-by-samples numeric matrix or a SeqExpressionSet object containing read counts. Let's set up a SeqExpressionSet with our counts matrix
set <- newSeqExpressionSet(as.matrix(y$counts), phenoData = data.frame(Island, row.names=colnames(y)))

# exploratory analyis of results before normalisation
pdf("no_Normalisation.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[batch])
plotPCA(set, col=colors[batch], cex=1.2, pch=as.numeric(batch) + 14, labels=F)
legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 14, title="Batch", cex=0.6, border=F, bty="n", col=unique(colors[batch]))
dev.off()

# Normalisation can be performed using "median","upper", or "full", however when passing to edgeR's normalisation method (below), the only options are "TMM","RLE", and "upperquartile". In order to keep consistency, we'll go ahead and choose the "upper" method, since it's in both.
set <- betweenLaneNormalization(set, which="upper")

# identify which samples are replicated
allreplicated=as.factor(samplenames %in% allreps)

# exploratory analysis after normalisation
pdf("upperQuartileNOrmalisation.pdf")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[batch])
plotPCA(set, col=colors[batch], cex=1.2, pch=as.numeric(batch) + 14, labels=F, main="Batch")
# plot PCA to see how closely replicates sit
plotPCA(set, col=lightcolors[allreplicated], cex=1.2, main="replicates")
legend(legend=unique(y$samples$batch), "topright", pch=unique(y$samples$batch) + 14, title="Batch", cex=0.6, border=F, bty="n", col=unique(colors[batch]))
legend(legend=unique(allreplicated), pch=16, x="bottomright", col=unique(colors[allreplicated]), cex=0.6, title="replicate", border=F, bty="n")
dev.off()

# First, construct a matrix specifying the replicates. 
# create matrix filled with -1s 
replicates=matrix(-1, nrow=length(allreps), ncol=3)
rownames(replicates)=unique(samplenames[samplenames %in% allreps])
for (i in 1:nrow(replicates)){
    replicates[i,1:length(grep(rownames(replicates)[i], samplenames))] = grep(rownames(replicates)[i], samplenames)
}

# set up genes used for RUVs
genes <- rownames(y)

# set RUVs. We'll first set this to 5 and see how it looks
pdf("RUVsNormalisation_choosingK.pdf")
for (i in c(1,2,3,4,5,10)){
    set1 <- RUVs(set, genes, k=i, replicates)
    plotRLE(set1, main=paste("Housekeeping Controls",i, sep="\n"), col=as.numeric(batch), outline=FALSE)
    # plot PCA to see how closely replicates sit
    plotPCA(set1, col=lightcolors[allreplicated], cex=1.2, main=paste(i, "replicates", sep="\n"))
    # plot pca to see if batch effect is eliminated/minimized
    plotPCA(set1, col=as.numeric(batch), cex=1.2, main="Batch")
}
dev.off()


# Set up list of housekeeping genes as controls (from Eisenberg and Levanon, 2003) for volcano plots
housekeeping=read.table("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/BatchEffects/Housekeeping_ControlGenes.txt", as.is=T, header=F)
# if this is broken, use host = "uswest.ensembl.org"
ensembl.mart.90 <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'www.ensembl.org', ensemblRedirect = FALSE)
biomart.results.table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), mart = ensembl.mart.90,values=housekeeping, filters="hgnc_symbol")
hkGenes=as.vector(biomart.results.table[,1])
hkControls=hkGenes[which(hkGenes %in% rownames(y$counts))]

# look at the p-value distributions and choose the best k
deGenes=vector()
k=vector()
for (j in 1:3){
    counter=0
    pdf(paste0("pvalueDist_choosingK_RUVs_",j,".pdf"))
    for (i in c(1:6)){
        counter=counter+1
        set1 <- RUVs(set, genes, k=i, replicates)
        design <- cbind(model.matrix(~0 + Island), pData(set1)[2:(i+1)])
        colnames(design)=gsub("Island", "", colnames(design))
        colnames(design)[3]="Mappi"
        z <- DGEList(counts=counts(set1), group=Island)
        # here, we use uppserquarile normalisation since "upper" is used for between lane normalisation of our 'set' object
        z <- calcNormFactors(z, method="upperquartile")
        contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))
        v <- voom(z, design, plot=FALSE, normalize.method = "quantile")
        vfit <- lmFit(v, design)
        vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
        efit <- eBayes(vfit)
        topTable <- topTable(efit, coef=j, p.value=0.01, lfc=1, n=Inf)
        deGenes[counter]=nrow(topTable)
        k[counter]=i
        # get p-value distribution from raw pvalues
        hist(efit$p.value[,j], main=paste(colnames(efit)[j],i,sep="\n"), ylim=c(0,max(table(round(efit$p.value[,j], 1)))+1000), xlab="p-value")
        # get volcano plots
        plot(efit$coef[,j], -log10(as.matrix(efit$p.value)[,j]), pch=20, main=paste(colnames(efit)[j],i,sep="\n"), xlab="log2FoldChange", ylab="-log10(pvalue)")
        points(efit$coef[,j][which(names(efit$coef[,j]) %in% hkControls)], -log10(as.matrix(efit$p.value)[,j][which(names(efit$coef[,j]) %in% hkControls)]) , pch=20, col=4, xlab="log2FoldChange", ylab="-log10(pvalue)")
        legend("topleft", c("genes", "hk genes"),fill=c("black",4))
        abline(v=c(-1,1), lty=2)
    }
    dev.off()
    pdf(paste("numberDeGenes_choosingK_RUVs_",colnames(efit)[j],".pdf"))
    plot(k, deGenes, main=colnames(efit)[j])
    dev.off()
}

# we'll set k to 5 for now
set1 <- RUVs(set, genes, k=5, replicates)

# make a matrix for labels of replicates
xcord=c(-0.14,-0.03, 0.02, 0.02, 0.05)
ycord=c(0.1,-0.09,0.0,0.087,0.03)
textcords=cbind(xcord,ycord)

# now plot all of the covariates and see how the batch correction worked
for (name in c(covariate.names[c(1:10,18)], "allreplicated")){
    pdf(paste0("batchCorrectedPCA_RUVs_",name,".pdf"))
    plotPCA(set1, labels=F, pch=as.numeric(batch) + 14, col=as.numeric(get(name)), main=name)
    # text(textcords, c("SMB-WNG-021", "MPI-381", "SMB-ANK-016", "SMB-ANK-027", "MTW-013"), cex=0.8)
    legend(legend=unique(get(name)), pch=16, x="bottomright", col=unique(as.numeric(get(name))), cex=0.6, title=name, border=F, bty="n")
    legend(legend=unique(y$samples$batch), "topright", pch=unique(as.numeric(batch)) + 14, title="Batch", cex=0.6, border=F, bty="n")
    dev.off()
}

# perform DE using a Limma pipeline --------------------------------------------------------------------------------------------------

design <- model.matrix(~0 + Island + W_1 + W_2 + W_3 + W_4 + W_5, data=pData(set1))
colnames(design)=gsub("Island", "", colnames(design))
colnames(design)[3]="Mappi"
z <- DGEList(counts=counts(set1), group=Island)
# add in genes to the DGEList for downstream use
geneid <- rownames(z)
genes <- select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),]
z$genes <- genes
z <- calcNormFactors(z, method="upperquartile")

# set up contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))

# now perform voom and see how it performs
pdf("Limma_voom_upperquartilenormalisation.pdf", height=8, width=15)
par(mfrow = c(1,2))
v <- voom(z, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final Model: Mean-variance trend")
dev.off()

# get number of DE genes that pass our theshold and write to file
dt <- decideTests(efit, p.value=0.01, lfc=1)
write.table(summary(dt), file="numberSigDEgenes_RUVs_voom.txt")

# for all contrasts, get names of significant genes.
pdf("MDPlot_housekeeping_TopGenes_FDRpval01_LFC01_.pdf", height=15,width=10)
par(mfrow = c(3,1))
for (i in 1:ncol(efit)){
    o <- which(names(efit$Amean) %in% names(which(abs(dt[,i])==1)))
    x <- efit$Amean
    m <- efit$coefficients[,i]
    t=which(names(efit$coefficients[,i]) %in% names(which(abs(dt[,i])==1)))
    G <- y$genes[names(which(dt[,i]==1)),]$SYMBOL
    plotMD(efit, column=i, status=dt[,i], main=colnames(efit)[i], hl.col=c("blue","red"), values=c(-1,1))
    abline(h=c(1,-1), lty=2)
    legend(legend=paste(names(summary(dt)[,i]), summary(dt)[,i], sep="="), x="bottomright", border=F, bty="n")
    text(x[o], m[t], labels=G)
}
dev.off()

# plot log2 fold change between islands
pdf("log2FC_IslandComparisons_pval01.pdf", height=15, width=10)
# note 'p.value' is the cutoff value for adjusted p-values
topTable <- topTable(efit, coef=1, n=Inf, p.value=0.01)
plot(density(topTable$logFC), col=1, xlim=c(-2,2), main="LogFC Density", xlab="LogFC", ylab="Density")
abline(v=c(-1,-0.5,0.5,1), lty=3)
for (i in 2:ncol(efit)){
    topTable <- topTable(efit, coef=i, n=Inf, p.value=0.01)
    lines(density(topTable$logFC), col=i, xlim=c(-2,2))
}
legend(x="topright", bty="n", col=1:ncol(efit), legend=colnames(efit), lty=1, lwd=2)
dev.off()


# We can also look at the top DE genes with a heatmap of logCPM values for the top 100 genes. Each gene (or row) is scaled so that mean expression is zero and the standard deviation is one (we're using 'E' from the voom object which is a numeric matrix of normalized expression values on the log2 scale). Samples with relatively high expression of a given gene are marked in red and samples with relatively low expression are marked in blue. Lighter shades and white represent genes with intermediate expression levels. Samples and genes are reordered by the method of hierarchical clustering
# set up colors
mycol <- colorpanel(1000,"blue","white","red")
colors=c("steelblue","goldenrod1","red2")
cc=colors[Island]
pdf("topGenes_Heatmap.pdf", height=15, width=15)
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf)
    if (nrow(topTable) > 0){
        index <- which(v$genes$ENSEMBL %in% topTable$ENSEMBL[1:100])
        heatmap.2(v$E[index,], scale="row",labRow=v$genes$SYMBOL[index], labCol=Island, col=mycol, trace="none", density.info="none", margin=c(8,6), lhei=c(2,10), dendrogram="column", ColSideColors=cc, keysize=1, cexRow=1.5, main=colnames(efit)[i])
    }
    # write out topTable genes for the whole gene list
    write.table(topTable, file=paste0("topTable_",colnames(efit)[i],".txt"))
    # get what these genes are doing from Biomart and save to file
    sig.Genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=rownames(topTable), filters="ensembl_gene_id")
    write.table(sig.Genes, file=paste0("significantGenes_RUVs_",colnames(efit)[i],".txt")) 
}
dev.off()

# show the number of DE genes between all islands
pdf("vennDiagram_allSigDEGenes_05FDR.pdf", height=15, width=15)
vennDiagram(dt[,1:3], circle.col=c("red", "blue", "green"))
dev.off()

# get DE genes in common with all islands
de.common <- which(dt[,1]!=0 & dt[,2]!=0 & dt[,3]!=0)

# find out what these genes are doing
commonGenes.results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=names(de.common), filters="ensembl_gene_id")
write.table(commonGenes.results, file="allCommonGenes.txt")

# get which genes are in common between populations with mappi (i.e., SMBvsMPI and MTWvsMPI)
de.common.MPI <- which(dt[,2]!=0 & dt[,3]!=0)

# find out what these genes are doing
commonGenes.MPI <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "interpro","interpro_description"), mart = ensembl.mart.90,values=names(de.common.MPI), filters="ensembl_gene_id")
write.table(commonGenes.MPI, file="allCommonGenes_MPI.txt")


# Let's see how the expression levels of each of our significant genes are distributed within each island. First, assign our top genes and ensembl IDs to variables
topGenes=c("MARCO", "SIGLEC6", "SIGLEC7","SIGLEC14", "KLRF1", "IFI27", "RNF182", "MDGA1", "UTS2", "MYOM2", "LOC102724159")
topEnsembl=c("ENSG00000019169", "ENSG00000105492", "ENSG00000168995", "ENSG00000254415","ENSG00000150045", "ENSG00000165949", "ENSG00000180537", "ENSG00000112139", "ENSG00000049247", "ENSG00000036448", "ENSG00000275464")

# To visualise distributions, we'll be making violin plots using ggpubr which needs p-value labels. Let's go ahead and make a matrix to input this into ggpubr
# first set up matrix
topGenes.pvalue=matrix(nrow=length(topEnsembl), ncol=ncol(efit))
rownames(topGenes.pvalue)=topEnsembl
colnames(topGenes.pvalue)=colnames(efit)
for (i in 1:ncol(efit)){
    # get significant genes over a logFC of 1 for all Island comparisons
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf)
    for(j in topEnsembl){
        # input the adjusted p.value for each gene
        topGenes.pvalue[j,i]=topTable[j,"adj.P.Val"]
    }
}

# make pvalues into scientific notation with max 3 digits
topGenes.pvalue=formatC(topGenes.pvalue, format="e", digits=2, drop0trailing=T)
# convert e notation to base 10 notation
topGenes.pvalue=sub("e", "x10^", topGenes.pvalue)

# finally, let's make the violin plots using fancy ggpubr
pdf("TopGenes_ggboxplot_Island.pdf", height=8, width=10)
counter=0
for(ensembl in topEnsembl){
    counter=counter+1
    gene.df <- data.frame(cpm(y$counts[which(y$genes$ENSEMBL==ensembl),],log=T),Island)
    colnames(gene.df)=c("CPM", "Island")
    annotation_df <- data.frame(start=c("Sumba","Sumba", "Mentawai"), end=c("Mentawai","West Papua","West Papua"), y=c(max(gene.df[,1]+4),max(gene.df[,1]+5),max(gene.df[,1]+6)), label=paste("limma p-value =",topGenes.pvalue[ensembl,],sep=" "))
    print(ggviolin(gene.df, x = "Island", y = "CPM", fill="Island", add=c("jitter","boxplot"), main=topGenes[counter], palette=1:3, add.params = c(list(fill = "white"), list(width=0.05))) + geom_signif(data=annotation_df,aes(xmin=start, xmax=end, annotations=label, y_position=y),textsize = 5, vjust = -0.2,manual=TRUE) + ylim(NA, max(gene.df[,1])+7))
}
dev.off()


# Enrichment analysis for Gene Ontology ----------------------------------------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(z)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
z$entrezID=entrez[match(rownames(z), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=z$entrezID) 
for (i in 1:ncol(contr.matrix)){
    camera.matrix=camera(v,idx,design,contrast=contr.matrix[,i])
    write.table(camera.matrix, file=paste0("cameraMatrix_",colnames(contr.matrix)[i],".txt"))
}

# gene set testing with goSeq
for(pop in 1:ncol(efit)){
    for(pval in c(0.05, 0.01)){
        topTable <- topTable(efit, coef=pop, n=Inf, p.value=pval, lfc=1)
        gene.vector=as.integer(rownames(z) %in% rownames(topTable))
        names(gene.vector)=rownames(z)

        # set the probability weighting fcuntion, i.e., implement a weight for each gene dependent on its length
        pwf=nullp(gene.vector,"hg19","ensGene")
        # use  the  default  method  to  calculate  the  over  and  under  expressed  GO categories among DE genes
        GO.wall=goseq(pwf,"hg19","ensGene",use_genes_without_cat=TRUE)

        # now let's interpret the results. First we need to apply a multiple hypothesis testing correction set at 5% (BH method)
        enriched.GO=GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05,]
        write.table(enriched.GO, file=paste("enrichedGOterms",colnames(efit)[pop],pval,".txt", sep="_"), quote=F, row.names=F)
        if(nrow(enriched.GO)>0){
            zz=file(paste("topTen_enrichedGOterms",pop,pval,".txt", sep="_"), open="wt")
            sink(zz)
            sink(zz, type = "message")
            # get GO terms for top ten enriched GO terms - write output with the sink() function
            for(go in 1:length(enriched.GO$category)){
                print(GOTERM[[enriched.GO$category[go]]])
            }
        }
        sink(type = "message")
        sink()

        # KEGG pathway analysis
        en2eg=as.list(org.Hs.egENSEMBL2EG)
        # Get the mapping from Entrez 2 KEGG
        eg2kegg=as.list(org.Hs.egPATH)
        # Define a function which gets all unique KEGG IDs
        # associated with a set of Entrez IDs
        grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
        # Apply this function to every entry in the mapping from
        # ENSEMBL 2 Entrez to combine the two maps
        kegg=lapply(en2eg,grepKEGG,eg2kegg)

        # produce PWF as before
        pwf=nullp(gene.vector,"hg19","ensGene")
        KEGG=goseq(pwf,gene2cat=kegg, use_genes_without_cat=TRUE)
        enriched.GO.kegg=KEGG[p.adjust(KEGG$over_represented_pvalue, method="BH")<.05,]
        write.table(enriched.GO.kegg, file=paste("enrichedGOkegg",pop,pval,".txt", sep="_"))
    }
}

# Enrichment Visualisation ----------------------------------------------------------------------------------------------------


# Having played around with many visualisation options, it seems as though using Revigo is the best and most user-friendly. Using the output from the GoSeq enriched Go results (FDR <0.05), we can plug our GO IDs into Revigo (http://revigo.irb.hr/) and then output this as an R script.
# First save  output (remember that SMBvsMTW has no significantly enriched GO pathways so we'll only do this for SMBvsMPI and MTWvsMPI)  

# Revigo for MTW vs MPI (pval of 0.05 FDR)- Treemap
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/REVIGO_MTWvsMPI_Treemap.r")
# scatterplot
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/REVIGO_MTWvsMPI_Scatterplot.r")

# Revigo for SMB vs MPI (pval of 0.05 FDR)- Treemap
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/REVIGO_SMBvsMPI_Treemap.r")
# scatterplot
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_RegulatoryVariation/code/Differential_Expression/123_combined/REVIGO_SMBvsMPI_Scatterplot.r")

# Reactome --------------------------------------------------
# get what the top genes are doing from Reactome
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf)
    entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene"),values= rownames(topTable),mart=ensembl.mart.90)
    de <- entrez[,2]
    x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
    if(nrow(x) > 0){
        pdf(paste0("enriched_pathways_reactome_",colnames(efit)[i],"PopComparisons.pdf"), height=5, width=10)
        print(barplot(x, showCategory=nrow(x), font.size = 8, title = colnames(efit)[i]))
        dev.off()
    }
}

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
