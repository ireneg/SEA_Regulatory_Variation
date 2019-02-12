# script created by KSB, 08.08.18
# Perform DE analysing relationship between islands

# load dependencies: libraries, human count data, and data preprocessing
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/countData_123_combined.R")
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/dataPreprocessing_123_combined.R")

library(EGSEA)

# set working directory
setwd("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Output/DE_Analysis/123_combined/DE_Island/LM_allCovarPlusBlood")

# Removing heteroscedascity with voom and fitting linear models -----------------------------------------------------------------------

# First, set up design matrix
# We don't know what the age is for SMB-PTB028 (#116) so we will just add in the median age of Sumba (44.5)
y$samples$Age[which(is.na(y$samples$Age) == T)]=45
design <- model.matrix(~0 + Island + y$samples$Age + batch + y$samples$RIN + y$sample$CD8T + y$sample$CD4T + y$sample$NK + y$sample$Bcell + y$sample$Mono + y$sample$Gran)
colnames(design)=gsub("Island", "", colnames(design))
#rename columns to exclude spaces and unrecognised characters
colnames(design)[c(3,4,7:13)]=c("Mappi","Age","RIN", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")

# set up contrast matrix
contr.matrix <- makeContrasts(SMBvsMTW=Sumba - Mentawai, SMBvsMPI=Sumba - Mappi, MTWvsMPI=Mentawai - Mappi, levels=colnames(design))

v <- voom(y, design, plot=F)
# fit linear models
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust=T)
plotSA(efit, main="Mean-variance trend elimination")
dev.off()

dt <- decideTests(efit, p.value=0.01, lfc=1)

# Enrichment analysis for Gene Ontology ----------------------------------------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
y$entrezID=entrez[match(rownames(y), entrez[,1]), 2]

# gene set testing with Camera
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata")) 
idx <- ids2indices(Hs.c2,id=y$entrezID) 
for (i in 1:ncol(contr.matrix)){
    camera.matrix=camera(v,idx,design,contrast=contr.matrix[,i])
    write.table(camera.matrix, file=paste0("cameraMatrix_",colnames(contr.matrix)[i],".txt"))
}

# gene set testing with goSeq
for(pop in 1:ncol(efit)){
    for(pval in c(0.05, 0.01)){
        topTable <- topTable(efit, coef=pop, n=Inf, p.value=pval, lfc=1)
        gene.vector=as.integer(rownames(y) %in% rownames(topTable))
        names(gene.vector)=rownames(y)

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

# Revigo for MTW vs MPI (pval of 0.05 FDR)
# scatterplot
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/REVIGO/Limma/REVIGO-MTWvsMPI.r")
# treemap
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/REVIGO/Limma/REVIGO_treemap-MTWvsMPI.r")

# Revigo for SMB vs MPI (pval of 0.05 FDR)
# scatterplot
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/REVIGO/Limma/REVIGO-SMBvsMPI.r")
# treemap
source("/Users/katalinabobowik/Documents/UniMelb_PhD/Analysis/UniMelb_Sumba/Scripts/GIT/SEA_Regulatory_Variation/code/Differential_Expression/123_combined/REVIGO/Limma/REVIGO_treemap-SMBvsMPI.r")

# Reactome ---------------------------------------------------------------------

# get what the top genes are doing from Reactome
for (i in 1:ncol(efit)){
    topTable <- topTable(efit, coef=i, p.value=0.01, lfc=1, n=Inf)
    entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene"),values= rownames(y),mart=ensembl.mart.90)
    de <- entrez[,2]
    x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
    if(nrow(x) > 0){
        assign(colnames(efit)[i], print(barplot(x, showCategory=nrow(x), font.size = 8, title = colnames(efit)[i])))
    }
}

pdf("enriched_pathways_reactome_PopComparisons.pdf", height=5, width=15)
ggarrange(SMBvsMPI, MTWvsMPI)
dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# EGSEA -----------------------------------------------------------------------

# transform ensembl IDs to entrez IDs to be compatible with human c2 dataset (below)
ensembl_genes=rownames(y)
entrez=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "description"),values= ensembl_genes,mart=ensembl.mart.90)
v$genes$entrezID=entrez[match(rownames(y), entrez[,1]), 2]


#subset voom object to contain only genes with an Entrez Gene ID
v <- v[which(!is.na(v$genes$entrezID)),]
v <- v[-c(which(duplicated(v$genes$entrezID))),]
rownames(v) <- v$genes$entrezID
rownames(v$genes) <- v$genes$entrezID

#build gene-set collections to compute enrichment for
gs.annots = buildIdx(entrezIDs=v$genes$entrezID, species="human",msigdb.gsets=c("c2", "c5"), go.part = TRUE)

#generate a symbolsMap
symbolsMap = v$genes[, c(4, 2)]
colnames(symbolsMap) = c("FeatureID", "Symbols")
symbolsMap[, "Symbols"] = as.character(symbolsMap[, "Symbols"])

#establish the list of independent gene-set enrichment methods for EGSEA to use
baseMethods = egsea.base()[-2]

# Ensemble testing with EGSEA
gsa = egsea(voom.results=v, contrasts=contr.matrix, gs.annots=gs.annots, symbolsMap=symbolsMap, baseGSEAs=baseMethods, sort.by="p.adj", num.threads = 8, report = FALSE)

# Visualisation of results
# first make function for plotting
all.gsa.vis=function(db,pop,sorter){
    plotSummaryHeatmap(gsa, gs.label=db,hm.vals = "avg.logfc.dir",file.name=paste("plotSummaryHeatmap",db,".pdf",sep="_"), sort.by="p.adj")
    plotSummary(gsa, gs.label = db, contrast = pop, file.name = paste(pop,db,sorter,"plotSummary.pdf",sep="_"), x.cutoff=100, sort.by=sorter)
    plotMethods(gsa, gs.label = db, contrast = pop,file.name = paste(pop,db,"mds.pdf",sep="_"))
    plotBars(gsa, gs.label = db, contrast = pop, file.name=paste(pop,db,sorter,"plotBars.pdf",sep="_"), sort.by=sorter)
    write.table(topSets(gsa, contrast = pop,names.only=FALSE, number = Inf, verbose = FALSE, gs.label=db), file=paste(pop,db,"topSets.txt",sep="_"))
}

# now plot for all pops
for (pop in colnames(efit)){
    for (db in c("c2","c5BP","kegg")){
        all.gsa.vis(db,pop,sorter="p.adj")
    }
}

# we can also plot GO graphs
for (pop in colnames(efit)){
    for (label in c("c5BP", "c5CC")){
        plotGOGraph(gsa, gs.label=label, contrast = pop, file.name=paste(pop,label,"GOGraph.pdf"sep="_"))
    }
}

# interstingly, we see malaria as a signature. Let's plot the pathway and a heatmap of this

# plot results at the gene level
for (pop in colnames(efit)){
    plotPathway(gsa, gene.set = "Malaria",contrast = pop, gs.label = "kegg",file.name = paste(pop,"pathwayMalaria.pdf",sep="_"))
    plotHeatmap(gsa, gene.set="Malaria", gs.label="kegg",contrast = pop, file.name = paste(pop,"kegg_Heatmap_Malaria.pdf",sep="_"))   
}
