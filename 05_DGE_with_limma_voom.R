## ----download_data_DGE, warning=FALSE, message=FALSE----------------------------------
## Load the container package for RSE
library("SummarizedExperiment")

## Connect to ExperimentHub
library("ExperimentHub")
eh <- ExperimentHub::ExperimentHub()

## Load package datasets
myfiles <- query(eh, "smokingMouse")

## Download the mouse gene data
rse_gene <- myfiles[["EH8313"]]

## Samples from the nicotine experiment and from pups only
rse_gene_nic <- rse_gene[, which(rse_gene$Expt == "Nicotine" & rse_gene$Age == "Pup")]

## Retain only expressed genes (passed the filtering step)
rse_gene_filt <- rse_gene_nic[
    rowData(rse_gene_nic)$retained_after_feature_filtering == TRUE,
]


## ----explore_data---------------------------------------------------------------------
## Data dimensions: number of genes and samples
dim(rse_gene_filt)

## Raw counts for first 3 genes in the first 5 samples
assays(rse_gene_filt)$counts[1:3, 1:5]

## Log-normalized counts for first 3 genes in the first 5 samples
assays(rse_gene_filt)$logcounts[1:3, 1:5]

## Data for the first 2 samples
head(colData(rse_gene_filt), 2)


## ----model.matrix()-------------------------------------------------------------------
## Define formula
formula <- ~ Group + Sex + flowcell + mitoRate + overallMapRate + totalAssignedGene + detected + ERCCsumLogErr

## Model matrix
model <- model.matrix(formula, data = colData(rse_gene_filt))
head(model)


## ----coeff----------------------------------------------------------------------------
## Comparison of interest: Group
coef <- "GroupExperimental"


## ----voom-----------------------------------------------------------------------------
library("limma")

## voom():
#   1. Transform counts to log2(cpm)
#      ----------------------------------------------------------------------------
# .     |   Note we passed voom() raw counts as input, not the lognorm counts!!!   |
#      ----------------------------------------------------------------------------
#   2. Estimate mean-variance relationship for each gene
#   3. Compute observation weights for limma (next step)

vGene <- voom(assay(rse_gene_filt), design = model, plot = TRUE)


## ----voom_outs------------------------------------------------------------------------
## Returned data
names(vGene)

## E: contains the computed log(cpm)
dim(vGene$E)
vGene$E[1:5, 1:5]

## weights: contains the computed variance weight for each observation
dim(vGene$weights)
vGene$weights[1:5, 1:5]

## design: is the provided design matrix
head(vGene$design)

## targets: the sample library sizes used to compute log(cpm) in the first step
dim(vGene$targets)
head(vGene$targets)

identical(vGene$targets$lib.size, colSums(assay(rse_gene_filt)))


## ----lmFit----------------------------------------------------------------------------
## lmFit():
#   1. Fit linear model for each gene to estimate logFCs
fitGene <- lmFit(vGene)

## Corroborate "ls" method was applied
fitGene$method

## Explore outputs: estimated coefficients (logFCs)
head(fitGene$coefficients)


## ----lm-------------------------------------------------------------------------------
## Lognorm expression of first gene
rse_gene_one_gene <- rse_gene_filt[1, ]
colData(rse_gene_one_gene) <- cbind(colData(rse_gene_one_gene),
    "lognorm_expr" = assays(rse_gene_one_gene)$logcounts[1, ]
)


## Fit simple linear model
formula <- lognorm_expr ~ Group
lm <- lm(formula, data = colData(rse_gene_one_gene))
summary(lm)

## Two sample t-test
t.test(formula, data = colData(rse_gene_one_gene), var.equal = TRUE)


## ----eBayes---------------------------------------------------------------------------
## eBayes()
## 1. Compute the empirical Bayes statistics for DE
eBGene <- eBayes(fitGene)

## Outputs of interest:
## s2.prior -> prior residual variance (prior mean 1/s0^2)
##             in prior distribution of residual variances
eBGene$s2.prior

## df.prior -> degrees of freedom d0 in prior distribution
##             of residual variances
eBGene$df.prior

## s2.post -> posterior residual sample variances of the genes (~sg^2)
length(eBGene$s2.post)
head(eBGene$s2.post)

##  t -> moderated t-stats of the genes for each contrast
dim(eBGene$t)
eBGene$t[1:5, 1:5]

## p.value: corresponding unadjusted p-values of moderated t-stats
dim(eBGene$p.value)
eBGene$p.value[1:5, 1:5]


## ----topTable-------------------------------------------------------------------------
## topTable()
## 1. Obtain gene-wise DE stats for Group (Nicotine vs Ctrl)
top_genes <- topTable(eBGene, coef = coef, p.value = 1, number = nrow(rse_gene_filt), sort.by = "none")

## Outputs for each gene and for the coeff selected (Group):
##  logFC: log2-fold-changes
head(top_genes$logFC)


## ----betas_logFCs---------------------------------------------------------------------
setdiff(top_genes$logFC, eBGene$coefficients[, "GroupExperimental"])


## ----topTable2------------------------------------------------------------------------
##  t: moderated t-stats
head(top_genes$t)

## . P.value: unadjusted p-values of t-stats
head(top_genes$P.Value)

##  adj.P.Val: p-values adjusted to control the FDR
head(top_genes$adj.P.Val)


## ----hist_p---------------------------------------------------------------------------
## Histogram of unadjusted p-values
hist(top_genes$P.Value, xlab = "p-values", main = "")


## ----quantify_DEGs--------------------------------------------------------------------
## DEGs for FDR<0.05
de_genes <- top_genes[which(top_genes$adj.P.Val < 0.05), ]

## Number of DEGs
dim(de_genes)


## ----volcano plot, warning=FALSE, message=FALSE---------------------------------------
library("ggplot2")

## Define up- and down-regulated DEGs, and non-DEGs
FDR <- 0.05
DE <- vector()
for (i in 1:dim(top_genes)[1]) {
    if (top_genes$adj.P.Val[i] > FDR) {
        DE <- append(DE, "n.s.")
    } else {
        if (top_genes$logFC[i] > 0) {
            DE <- append(DE, "Up")
        } else {
            DE <- append(DE, "Down")
        }
    }
}
top_genes$DE <- DE

## Colors, sizes and transparencies for up & down DEGs and non-DEGs
cols <- c("Up" = "indianred2", "Down" = "steelblue2", "n.s." = "grey")
sizes <- c("Up" = 1.3, "Down" = 1.3, "n.s." = 0.8)
alphas <- c("Up" = 0.4, "Down" = 0.6, "n.s." = 0.5)

## Plot volcano plot
ggplot(
    data = top_genes,
    aes(
        x = logFC, y = -log10(adj.P.Val),
        color = DE,
        fill = DE,
        size = DE,
        alpha = DE
    )
) +
    geom_point(shape = 21) +
    geom_hline(
        yintercept = -log10(FDR),
        linetype = "dashed", color = "gray35", linewidth = 0.5
    ) +
    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed", color = "gray35", linewidth = 0.5
    ) +
    labs(y = "-log10(FDR)", x = "logFC(Nicotine vs Control)") +
    theme_bw() +
    scale_color_manual(values = cols, name = "Differential expression") +
    scale_fill_manual(values = cols, name = "Differential expression") +
    scale_size_manual(values = sizes, name = "Differential expression") +
    scale_alpha_manual(values = alphas, name = "Differential expression") +
    theme(
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.key.height = unit(0.15, "cm"),
        axis.title = element_text(size = (13)),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12)
    )


## ----complexHeatmap, warning=FALSE,message=FALSE--------------------------------------
library("ComplexHeatmap")

## We plot lognorm counts
lognorm_data <- assays(rse_gene_filt)$logcounts

## Subset to DEGs only
lognorm_data <- lognorm_data[rownames(de_genes), ]

## Define column (sample) names to display
colnames(lognorm_data) <- paste0("Pup_", 1:dim(lognorm_data)[2])


## ----scale_data-----------------------------------------------------------------------
## Center and scale the data to make differences more evident
lognorm_data <- (lognorm_data - rowMeans(lognorm_data)) / rowSds(lognorm_data)

## Sample annotation: Sex, Group, and library size
col_anno <- HeatmapAnnotation(
    df = as.data.frame(colData(rse_gene_filt)[, c("Sex", "Group")]),
    library_size = anno_barplot(colData(rse_gene_filt)$sum, gp = gpar(fill = "lightyellow2")),
    col = list(
        "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
        "Group" = c("Control" = "gray68", "Experimental" = "gold2")
    )
)

## Gene annotation: logFC and biotype
de_genes$logFC_binary <- sapply(de_genes$logFC, function(x) {
    if (x > 0) {
        ">0"
    } else {
        "<0"
    }
})
de_genes$protein_coding_gene <- sapply(rowData(rse_gene_filt[rownames(de_genes), ])$gene_type, function(x) {
    if (x == "protein_coding") {
        "TRUE"
    } else {
        "FALSE"
    }
})

gene_anno <- rowAnnotation(
    df = as.data.frame(cbind(
        "logFC" = de_genes$logFC_binary,
        "protein_coding_gene" = de_genes$protein_coding_gene
    )),
    col = list(
        "logFC" = c("<0" = "deepskyblue3", ">0" = "brown2"),
        "protein_coding_gene" = c("TRUE" = "darkseagreen3", "FALSE" = "magenta")
    )
)


## ----heat map, warning=FALSE, message=FALSE, fig.width=10,fig.height=6.9--------------
library("circlize")

## Plot
Heatmap(lognorm_data,
    name = "lognorm counts",
    show_row_names = FALSE,
    top_annotation = col_anno,
    left_annotation = gene_anno,
    row_km = 2,
    column_km = 2,
    col = colorRamp2(c(-4, -0.0001, 00001, 4), c("darkblue", "lightblue", "lightsalmon", "darkred")),
    row_title = "DEGs",
    column_title = "Samples",
    column_names_gp = gpar(fontsize = 7),
    heatmap_width = unit(12.5, "cm"),
    heatmap_height = unit(12.5, "cm")
)

