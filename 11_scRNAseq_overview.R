## ----load_416b------------------------------------------------------------------------
library("scRNAseq")
library("SingleCellExperiment")
library("AnnotationHub")
library("scater")

## Load the data set
sce.416b <- LunSpikeInData(which = "416b")

## We convert the blocking factor to a factor so that downstream steps do not treat it as an integer.
sce.416b$block <- factor(sce.416b$block)

## rename the rows with the symbols, reverting to Ensembl identifiers
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID", column = "SYMBOL"
)
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID", column = "SEQNAME"
)

rownames(sce.416b) <- uniquifyFeatureNames(
    rowData(sce.416b)$ENSEMBL,
    rowData(sce.416b)$SYMBOL
)


## ----sce_basics-----------------------------------------------------------------------
## Look at your SCE
sce.416b

## Get in the slot "assay", in the count matrix
## [genes, cells]
assay(sce.416b, "counts")[110:113, 1:2] # gene, cell

## We can do it like this too
counts(sce.416b)[110:113, 1:2]

## We could add more assays to our SCE
sce.416b <- logNormCounts(sce.416b)
sce.416b

## Acces to the column names (cell identifyers)
head(colnames(sce.416b))

## Acces to the column data (cell information)
head(colData(sce.416b))

## Acces to the row names (gene names)
head(rownames(sce.416b))

## Acces to the row data (gene information)
head(rowData(sce.416b))

## We can create another SCE subsetitng the first one
sce_2 <- sce.416b[110:130, 1:2]
sce_2


## ----$_operator-----------------------------------------------------------------------
head(sce.416b$`cell type`)


## ----reducedDimNames------------------------------------------------------------------
## This is empty
reducedDimNames(sce_2)
## Compute PCA
sce_2 <- runPCA(sce_2)
## Check again
reducedDimNames(sce_2)


## ----identify_mito_percellQC----------------------------------------------------------
library("scuttle")

## Identify mitochondrial genes (those with SEQNAME equal to "MT") in the row data
mito <- which(rowData(sce.416b)$SEQNAME == "MT")

## Compute per-cell QC metrics, including a subset for mitochondrial genes
stats <- perCellQCMetrics(sce.416b, subsets = list(Mt = mito))

summary(stats$sum) # total library sizes for all cells
summary(stats$detected) # detected features (genes)
summary(stats$subsets_Mt_percent) # percentage of reads mapping to mitochondrial genes
summary(stats$altexps_ERCC_percent) # percentage of reads mapping to spike-in controls


## ----addpercellQC_subMito-------------------------------------------------------------
## Compute addPerCellQCMetrics, including a subset for mitochondrial genes
sce.416b <- addPerCellQCMetrics(sce.416b, subsets = list(Mito = mito))
colnames(colData(sce.416b))


## ----QC_fixedThresholds---------------------------------------------------------------
## Using our previous perCellQCMetrics data:

## Identify cells with a total library size (sum of counts) less than 100,000
c.lib <- stats$sum < 1e5
## Identify cells with fewer than 5,000 detected features (genes)
qc.nexprs <- stats$detected < 5e3
## Identify cells with more than 10% of reads mapping to spike-in controls (e.g., ERCC)
qc.spike <- stats$altexps_ERCC_percent > 10
## Identify cells with more than 10% of reads mapping to mitochondrial genes
qc.mito <- stats$subsets_Mt_percent > 10

## Create a combined logical vector that marks cells to discard if they meet any of the above criteria
discard <- c.lib | qc.nexprs | qc.spike | qc.mito

## Summarize the number of cells removed for each reason.
DataFrame(
    LibSize = sum(c.lib), # Number of cells removed due to low library size
    NExprs = sum(qc.nexprs), # Number of cells removed due to low number of detected features
    SpikeProp = sum(qc.spike), # Number of cells removed due to high spike-in proportion
    MitoProp = sum(qc.mito), # Number of cells removed due to high mitochondrial proportion
    Total = sum(discard) # Total number of cells removed
)


## ----QC_adaptiveThreshold-------------------------------------------------------------
## Identify cells that are outlier
reasons <- perCellQCFilters(stats,
    sub.fields = c("subsets_Mt_percent", "altexps_ERCC_percent")
) # No transformation

colSums(as.matrix(reasons))

## Extract the exact filter thresholds
attr(reasons$low_lib_size, "thresholds")
attr(reasons$low_n_features, "thresholds")


## ----QC_sce416b_plots-----------------------------------------------------------------
library("scater")

## Add the information to the SCE columns
colData(sce.416b) <- cbind(colData(sce.416b), stats)
sce.416b$block <- factor(sce.416b$block)
sce.416b$phenotype <- ifelse(grepl("induced", sce.416b$phenotype), "induced", "wild type")
sce.416b$discard <- reasons$discard

## Plot
gridExtra::grid.arrange(
    ## Diccard low total counts
    plotColData(sce.416b,
        x = "block", y = "sum", colour_by = "discard",
        other_fields = "phenotype"
    ) + facet_wrap(~phenotype) +
        scale_y_log10() + ggtitle("Total count"),
    ## Discard low detected genes
    plotColData(sce.416b,
        x = "block", y = "detected", colour_by = "discard",
        other_fields = "phenotype"
    ) + facet_wrap(~phenotype) +
        scale_y_log10() + ggtitle("Detected features"),
    ## Discard high mitocondrial percentage
    plotColData(sce.416b,
        x = "block", y = "subsets_Mito_percent",
        colour_by = "discard", other_fields = "phenotype"
    ) +
        facet_wrap(~phenotype) + ggtitle("Mito percent"),
    ## Discard high
    plotColData(sce.416b,
        x = "block", y = "altexps_ERCC_percent",
        colour_by = "discard", other_fields = "phenotype"
    ) +
        facet_wrap(~phenotype) + ggtitle("ERCC percent"),
    ncol = 1
)


## ----eval=interactive()---------------------------------------------------------------
library("iSEE")
iSEE(sce.416b)


## ----discard_sce416b------------------------------------------------------------------
## Keep the columns we DON'T want to discard.
filtered <- sce.416b[, !reasons$discard]


## ----load_sceZeisel-------------------------------------------------------------------
library("scRNAseq")
library("scater")

## Load dataset
sce.zeisel <- ZeiselBrainData()
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel,
    ids = sub("_loc[0-9]+$", "", rownames(sce.zeisel))
)
## Compute perCellQCMetrics
stats <- perCellQCMetrics(sce.zeisel, subsets = list(
    Mt = rowData(sce.zeisel)$featureType == "mito"
))
## Compute quickPerCellQC
qc <- quickPerCellQC(stats, percent_subsets = c(
    "altexps_ERCC_percent",
    "subsets_Mt_percent"
))
## Discard low quality cells
sce.zeisel <- sce.zeisel[, !qc$discard]


## ----ibrarySizeFactors_zeisel---------------------------------------------------------
library("scater")

## Compute librarySizeFactors
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)
summary(lib.sf.zeisel)


## ----hist_libSizeFactors--------------------------------------------------------------
## Plot the library size factors differences
hist(log10(lib.sf.zeisel), xlab = "Log10[Size factor]", col = "grey80")


## ----deconvolution_norm---------------------------------------------------------------
library("scran")

## Compute quickCluster + calculateSumFactor for deconvolution normalization
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)
table(clust.zeisel)

deconv.sf.zeisel <- calculateSumFactors(sce.zeisel, clusters = clust.zeisel)
summary(deconv.sf.zeisel)


## ----spike-ins_norm-------------------------------------------------------------------
library("scRNAseq")
sce.richard <- RichardTCellData()
sce.richard <- sce.richard[, sce.richard$`single cell quality` == "OK"]
sce.richard


## ----computeSpikeFactors--------------------------------------------------------------
## computeSpikeFactors() to estimate spike-in size factors
sce.richard <- computeSpikeFactors(sce.richard, "ERCC")
summary(sizeFactors(sce.richard))


## ----logNormCounts_zeisel-------------------------------------------------------------
## Compute normalized expression values and log-transformation
sce.zeisel <- logNormCounts(sce.zeisel)

assayNames(sce.zeisel)


## ----modelGeneVar_zeisel--------------------------------------------------------------
library("scran")

## Model the mean-variance relationship
dec.zeisel <- modelGeneVar(sce.zeisel)

## Plot the fit
fit.zeisel <- metadata(dec.zeisel)
plot(fit.zeisel$mean, fit.zeisel$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
curve(fit.zeisel$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


## ----HVGs_zeisel----------------------------------------------------------------------
## Order by most interesting genes for inspection
dec.zeisel[order(dec.zeisel$bio, decreasing = TRUE), ]


## ----modelGeneVarWithSpikes_416b------------------------------------------------------
## Fit a mean-dependent trend to the variance of the spike-in transcripts
dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")
## Order by most interesting genes for inspection
dec.spike.416b[order(dec.spike.416b$bio, decreasing = TRUE), ]

## Plot the fit
plot(dec.spike.416b$mean, dec.spike.416b$total,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
fit.spike.416b <- metadata(dec.spike.416b)
points(fit.spike.416b$mean, fit.spike.416b$var, col = "red", pch = 16)
curve(fit.spike.416b$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


## ----modelGeneVarByPoisson_zeisel-----------------------------------------------------
## construct a mean-variance trend in the log-counts
set.seed(0010101)
dec.pois.zeisel <- modelGeneVarByPoisson(sce.zeisel)

## Order by most interesting genes for inspection
dec.pois.zeisel <- dec.pois.zeisel[order(dec.pois.zeisel$bio, decreasing = TRUE), ]
head(dec.pois.zeisel)

## Plot the fit
plot(dec.pois.zeisel$mean, dec.pois.zeisel$total,
    pch = 16, xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
curve(metadata(dec.pois.zeisel)$trend(x), col = "dodgerblue", add = TRUE)


## ----modelGeneVar_batch---------------------------------------------------------------
## Fit a mean-dependent trend to the variance of the spike-in transcripts
## Independently for each batch (block)
dec.block.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", block = sce.416b$block) # block=sce.416b$block
head(dec.block.416b[order(dec.block.416b$bio, decreasing = TRUE), 1:6])

## Plot the fit by batch (block)
par(mfrow = c(1, 2))
blocked.stats <- dec.block.416b$per.block
for (i in colnames(blocked.stats)) {
    current <- blocked.stats[[i]]
    plot(current$mean, current$total,
        main = i, pch = 16, cex = 0.5,
        xlab = "Mean of log-expression", ylab = "Variance of log-expression"
    )
    curfit <- metadata(current)
    points(curfit$mean, curfit$var, col = "red", pch = 16)
    curve(curfit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
}


## ----Select_TopHVGs-------------------------------------------------------------------
## Top 1000 genes
hvg.zeisel.var <- getTopHVGs(dec.zeisel, n = 1000)
str(hvg.zeisel.var)


## ----getHVGs_zeisel-------------------------------------------------------------------
library("scran")

## Top 2000 HVGs
top.zeisel <- getTopHVGs(dec.zeisel, n = 2000)

## Principal component analysis using top 2000 HVGs, 50 PCs
set.seed(100)
sce.zeisel <- fixedPCA(sce.zeisel, subset.row = top.zeisel)

reducedDimNames(sce.zeisel)


## ----VarExplained_PCs-----------------------------------------------------------------
## Variance explained by PCs
percent.var <- attr(reducedDim(sce.zeisel), "percentVar")
plot(percent.var, log = "y", xlab = "PC", ylab = "Variance explained (%)")


## ----PCs_zeisel-----------------------------------------------------------------------
library("scater")
## Plot PCA (Top 2 PCs for 2 dimentional visualization)
plotReducedDim(sce.zeisel, dimred = "PCA", colour_by = "level1class")


## ----Plot_multiplePCA_PCs-------------------------------------------------------------
## plot top 4 PCs against each other in pairwise plots
plotReducedDim(sce.zeisel, dimred = "PCA", ncomponents = 4, colour_by = "level1class")


## ----runTSNE_zeisel-------------------------------------------------------------------
## TSNE using runTSNE() stores the t-SNE coordinates in the reducedDims
set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA")
## Plot TSNE
plotReducedDim(sce.zeisel, dimred = "TSNE", colour_by = "level1class")


## ----TSNE_perplexity_plots------------------------------------------------------------
## run TSNE using diferent perplexity numbers and plot

## TSNE using perplexity = 5
set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA", perplexity = 5)
out5 <- plotReducedDim(sce.zeisel,
    dimred = "TSNE",
    colour_by = "level1class"
) + ggtitle("perplexity = 5")
## TSNE using perplexity = 20
set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA", perplexity = 20)
out20 <- plotReducedDim(sce.zeisel,
    dimred = "TSNE",
    colour_by = "level1class"
) + ggtitle("perplexity = 20")
## TSNE using perplexity = 80
set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred = "PCA", perplexity = 80)
out80 <- plotReducedDim(sce.zeisel,
    dimred = "TSNE",
    colour_by = "level1class"
) + ggtitle("perplexity = 80")
## Combine plots
gridExtra::grid.arrange(out5, out20, out80, ncol = 3)


## ----Umap_zeisel----------------------------------------------------------------------
## UMAP using runUMAP() stores the coordinates in the reducedDims
set.seed(100)
sce.zeisel <- runUMAP(sce.zeisel, dimred = "PCA")
## Plot UMAP
plotReducedDim(sce.zeisel, dimred = "UMAP", colour_by = "level1class")


## ----Clustering_clusterCells----------------------------------------------------------
library("scran")
## Cluster using "scran::clusterCells"
nn.clusters <- clusterCells(sce.zeisel, use.dimred = "PCA")
## Cluster assignments
table(nn.clusters)


## ----plot_clusters_zeisel-------------------------------------------------------------
## Save the cluster assignments
colLabels(sce.zeisel) <- nn.clusters
## Plot TSNE coloured by cluster assignments
plotReducedDim(sce.zeisel, "TSNE", colour_by = "label")


## ----clustering_specify_K-------------------------------------------------------------
library(bluster)

## Clustering using k=10
nn.clusters2 <- clusterCells(sce.zeisel,
    use.dimred = "PCA",
    BLUSPARAM = SNNGraphParam(k = 10, type = "rank", cluster.fun = "walktrap")
)

table(nn.clusters2)


## ----ClusterCells_graph---------------------------------------------------------------
## Obtain the graph
nn.clust.info <- clusterCells(sce.zeisel, use.dimred = "PCA", full = TRUE)
head(nn.clust.info$objects$graph)


## ----moreRes_clustering---------------------------------------------------------------
## More resolved clustering using a smaller k (k=5)
clust.5 <- clusterCells(sce.zeisel, use.dimred = "PCA", BLUSPARAM = NNGraphParam(k = 5))
table(clust.5)


## ----lessRes_clustering---------------------------------------------------------------
## Less resolved clustering using a larger k (k=50)
clust.50 <- clusterCells(sce.zeisel, use.dimred = "PCA", BLUSPARAM = NNGraphParam(k = 50))
table(clust.50)


## ----weighting_scheme-----------------------------------------------------------------
## Cluster using the number of shared nearest neighbors (type="number")
clust.num <- clusterCells(sce.zeisel,
    use.dimred = "PCA",
    BLUSPARAM = NNGraphParam(type = "number")
)
table(clust.num)

## Cluster using the Jaccard index (similarity between sample sets)
clust.jaccard <- clusterCells(sce.zeisel,
    use.dimred = "PCA",
    BLUSPARAM = NNGraphParam(type = "jaccard")
)
table(clust.jaccard)

## Cluster without specifying a graph type (default method-KNNGraphParam)
clust.none <- clusterCells(sce.zeisel,
    use.dimred = "PCA",
    BLUSPARAM = KNNGraphParam()
)
table(clust.none)


## ----community_detection, eval=FALSE--------------------------------------------------
## clust.walktrap <- clusterCells(sce.zeisel,
##     use.dimred = "PCA",
##     BLUSPARAM = NNGraphParam(cluster.fun = "walktrap")
## )
## 
## clust.louvain <- clusterCells(sce.zeisel,
##     use.dimred = "PCA",
##     BLUSPARAM = NNGraphParam(cluster.fun = "louvain")
## )
## 
## clust.infomap <- clusterCells(sce.zeisel,
##     use.dimred = "PCA",
##     BLUSPARAM = NNGraphParam(cluster.fun = "infomap")
## )
## 
## clust.fast <- clusterCells(sce.zeisel,
##     use.dimred = "PCA",
##     BLUSPARAM = NNGraphParam(cluster.fun = "fast_greedy")
## )
## 
## clust.labprop <- clusterCells(sce.zeisel,
##     use.dimred = "PCA",
##     BLUSPARAM = NNGraphParam(cluster.fun = "label_prop")
## )
## 
## clust.eigen <- clusterCells(sce.zeisel,
##     use.dimred = "PCA",
##     BLUSPARAM = NNGraphParam(cluster.fun = "leading_eigen")
## )


## ----Hierarchical_clust---------------------------------------------------------------
library("scran")
## Top 2000 HVGs
top.416b <- getTopHVGs(sce.416b, n = 2000)
## Principal component analysis using top 2000 HVGs, 50 PCs
set.seed(100)
sce.416b <- fixedPCA(sce.416b, subset.row = top.416b)
## TSNE
sce.416b <- runTSNE(sce.416b, dimred = "PCA")


## ----plot_dendogram-------------------------------------------------------------------
library("dendextend")

## Perform hierarchical clustering on the PCA-reduced data from sce.416b
## The BLUSPARAM argument specifies the clustering method (here "ward.D2").
## The full=TRUE argument ensures that additional objects related to clustering are returned.
hclust.416b <- clusterCells(sce.416b,
    use.dimred = "PCA",
    BLUSPARAM = HclustParam(method = "ward.D2"), full = TRUE
)

## Extract the hierarchical clustering tree from the clustering result
tree.416b <- hclust.416b$objects$hclust

## Customize the dendrogram for better visualization
tree.416b$labels <- seq_along(tree.416b$labels)
## Convert the hierarchical clustering tree to a dendrogram object
dend <- as.dendrogram(tree.416b, hang = 0.1)
combined.fac <- paste0(
    sce.416b$block, ".",
    sub(" .*", "", sce.416b$phenotype)
)

labels_colors(dend) <- c(
    "20160113.wild" = "blue",
    "20160113.induced" = "red",
    "20160325.wild" = "dodgerblue",
    "20160325.induced" = "salmon"
)[combined.fac][order.dendrogram(dend)]

## Plot the dendrogram
plot(dend)


## ----cut_dendogram--------------------------------------------------------------------
library("dynamicTreeCut")

## Perform hierarchical clustering with dynamic tree cut on the PCA
## The BLUSPARAM argument specifies the clustering method (here "ward.D2"),
## and enables dynamic tree cut (cut.dynamic=TRUE) with specific parameters.
hclust.dyn <- clusterCells(sce.416b,
    use.dimred = "PCA",
    BLUSPARAM = HclustParam(
        method = "ward.D2", cut.dynamic = TRUE,
        cut.params = list(minClusterSize = 10, deepSplit = 1)
    )
)
table(hclust.dyn)

## Plot dendogram
labels_colors(dend) <- as.integer(hclust.dyn)[order.dendrogram(dend)]
plot(dend)

## Obtain assignations and plot TSNE
colLabels(sce.416b) <- factor(hclust.dyn)
plotReducedDim(sce.416b, "TSNE", colour_by = "label")


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------


## -------------------------------------------------------------------------------------

