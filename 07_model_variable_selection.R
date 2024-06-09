## ----download_data, warning=FALSE, message=FALSE--------------------------------------
## Load the container package for this type of data
library("SummarizedExperiment")

## Connect to ExperimentHub
library("ExperimentHub")
eh <- ExperimentHub::ExperimentHub()

## Load the datasets of the package
myfiles <- query(eh, "smokingMouse")

## Download the mouse gene data
rse_gene <- myfiles[["EH8313"]]

## Keep samples from nicotine experiment and pups only
rse_gene_nic <- rse_gene[
    ,
    which(rse_gene$Expt == "Nicotine" & rse_gene$Age == "Pup")
]

## Use expressed genes only (i.e. that passed the filtering step)
rse_gene_filt <- rse_gene_nic[
    rowData(rse_gene_nic)$retained_after_feature_filtering,
]

## Keep samples that passed QC and manual sample filtering steps (all passed)
rse_gene_filt <- rse_gene_filt[
    ,
    rse_gene_filt$retained_after_QC_sample_filtering &
        rse_gene_filt$retained_after_manual_sample_filtering
]


## ----CCA, message=FALSE, warning=FALSE, out.height=7, out.width=7---------------------
library("variancePartition")
library("pheatmap")

## Plot heatmap of correlations
## Define all variables to examine; remove those with single values
formula <- ~ Group + Sex + plate + flowcell + mitoRate + overallMapRate + totalAssignedGene + rRNA_rate + sum + detected + ERCCsumLogErr

## Measure correlations
CCA <- canCorPairs(formula, colData(rse_gene_filt))
## Heatmap
pheatmap(
    CCA, ## data
    color = hcl.colors(50, "YlOrRd", rev = TRUE), ## color scale
    fontsize = 8, ## text size
    border_color = "black", ## border color for heatmap cells
    cellwidth = unit(0.4, "cm"), ## height of cells
    cellheight = unit(0.4, "cm") ## width of cells
)


## ----message=FALSE, warning=FALSE-----------------------------------------------------
library("ggplot2")
library("cowplot")

## Boxplots/Scatterplots/Barplots for each pair of correlated variables

corr_plots <- function(sample_var1, sample_var2, sample_color) {
    ## Define sample colors by variable
    colors <- list(
        "Group" = c("Control" = "brown2", "Experimental" = "deepskyblue3"),
        "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
        "plate" = c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1"),
        "flowcell" = c(
            "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta",
            "HKCNKDSXX" = "turquoise3", "HKCTMDSXX" = "tomato"
        )
    )

    data <- colData(rse_gene_filt)

    ## a) Barplots for categorical variable vs categorical variable
    if (class(data[, sample_var1]) == "character" & class(data[, sample_var2]) == "character") {
        ## y-axis label
        y_label <- paste("Number of samples from each ", sample_var2, sep = "")

        ## Stacked barplot with counts for 2nd variable
        plot <- ggplot(data = as.data.frame(data), aes(
            x = !!rlang::sym(sample_var1),
            fill = !!rlang::sym(sample_var2)
        )) +
            geom_bar(position = "stack") +
            ## Colors by 2nd variable
            scale_fill_manual(values = colors[[sample_var2]]) +
            ## Show sample counts on stacked bars
            geom_text(aes(label = after_stat(count)),
                stat = "count",
                position = position_stack(vjust = 0.5), colour = "gray20", size = 3
            ) +
            theme_bw() +
            labs(
                subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits = 3)),
                y = y_label
            ) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }


    ## b) Boxplots for categorical variable vs continuous variable
    else if (class(data[, sample_var1]) == "character" & class(data[, sample_var2]) == "numeric") {
        plot <- ggplot(data = as.data.frame(data), mapping = aes(
            x = !!rlang::sym(sample_var1),
            y = !!rlang::sym(sample_var2),
            color = !!rlang::sym(sample_var1)
        )) +
            geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = NA) +
            geom_jitter(width = 0.15, alpha = 1, size = 1.5) +
            stat_smooth(method = "lm", geom = "line", alpha = 0.6, size = 0.4, span = 0.3, aes(group = 1), color = "orangered3") +
            scale_color_manual(values = colors[[sample_var1]]) +
            theme_bw() +
            guides(color = "none") +
            labs(
                subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits = 3)), y = gsub("_", " ", sample_var2),
                x = sample_var1
            ) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }


    ## c) Scatterplots for continuous variable vs continuous variable
    else if (class(data[, sample_var1]) == "numeric" & class(data[, sample_var2]) == "numeric") {
        plot <- ggplot(as.data.frame(data), aes(
            x = !!rlang::sym(sample_var1),
            y = !!rlang::sym(sample_var2),
            color = !!rlang::sym(sample_color)
        )) +
            geom_point(size = 2) +
            stat_smooth(method = "lm", geom = "line", alpha = 0.6, size = 0.6, span = 0.25, color = "orangered3") +
            ## Color by sample_color variable
            scale_color_manual(name = sample_color, values = colors[[sample_color]]) +
            theme_bw() +
            labs(
                subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits = 3)),
                y = gsub("_", " ", sample_var2),
                x = gsub("_", " ", sample_var1)
            ) +
            theme(
                axis.title = element_text(size = (7)),
                axis.text = element_text(size = (6)),
                plot.subtitle = element_text(size = 7, color = "gray40"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 7)
            )
    }

    return(plot)
}


## ----message=FALSE, warning=FALSE-----------------------------------------------------
## Correlation plot for Group and plate
p <- corr_plots("Group", "plate", NULL)
p + theme(plot.margin = unit(c(1, 5.5, 1, 5.5), "cm"))


## ----message=FALSE, warning=FALSE-----------------------------------------------------
## Correlation plot for overallMapRate and rRNA_rate
p <- corr_plots("overallMapRate", "rRNA_rate", "Group")
p + theme(plot.margin = unit(c(2, 3.5, 2, 3.5), "cm"))


## ----message=FALSE, warning=FALSE-----------------------------------------------------
## Correlation plot for overallMapRate and plate
p <- corr_plots("plate", "overallMapRate", NULL)
p + theme(plot.margin = unit(c(2, 5, 2, 5), "cm"))


## ----message=FALSE, warning=FALSE-----------------------------------------------------
## Correlation plot for overallMapRate and flowcell
p <- corr_plots("flowcell", "overallMapRate", NULL)
p + theme(plot.margin = unit(c(2, 5, 2, 5), "cm"))


## ----message=FALSE, warning=FALSE-----------------------------------------------------
## Correlation plots for sum and detected
p <- corr_plots("sum", "detected", "Group")
p + theme(plot.margin = unit(c(2, 3.5, 2, 3.5), "cm"))


## ----message=FALSE, warning=FALSE-----------------------------------------------------
p <- corr_plots("Group", "flowcell", NULL)
plots <- plot_grid(p)
plots + theme(plot.margin = unit(c(0.5, 5, 0.5, 5), "cm"))


## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------
## ## Fit a linear mixed model (LMM) that takes continuous variables as fixed effects and categorical variables as random effects
## 
## varPartAnalysis <- function(formula) {
##     ## Ignore genes with variance 0
##     genes_var_zero <- which(apply(assays(rse_gene_filt)$logcounts, 1, var) == 0)
##     if (length(genes_var_zero) > 0) {
##         rse_gene_filt <- rse_gene_filt[-genes_var_zero, ]
##     }
## 
##     ## Loop over each gene to fit the model and extract variance explained by each variable
##     varPart <- fitExtractVarPartModel(assays(rse_gene_filt)$logcounts, formula, colData(rse_gene_filt))
## 
##     # Sort variables by median fraction of variance explained (FVE)
##     vp <- sortCols(varPart)
##     p <- plotVarPart(vp)
## 
##     return(list(p, vp))
## }


## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------
## #####  Fit model with all variables  #####
## 
## # sum, detected, and ERCCsumLogErr are not included as they are in very different scales!
## formula <- ~ (1 | Group) + (1 | Sex) + (1 | plate) + (1 | flowcell) + mitoRate + overallMapRate +
##     totalAssignedGene + rRNA_rate
## plot <- varPartAnalysis(formula)[[1]]
## plot + theme(
##     plot.margin = unit(c(1, 1, 1, 1), "cm"),
##     axis.text.x = element_text(size = (7)),
##     axis.text.y = element_text(size = (7.5))
## )


## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------
## #####  Fit model without correlated variables  #####
## 
## ## Pup plots without overallMapRate and plate
## formula <- ~ (1 | Group) + (1 | Sex) + (1 | flowcell) + mitoRate + overallMapRate + totalAssignedGene
## varPart <- varPartAnalysis(formula)
## varPart_data <- varPart[[2]]
## plot <- varPart[[1]]
## plot + theme(
##     plot.margin = unit(c(1, 1, 1, 1), "cm"),
##     axis.text.x = element_text(size = (7)),
##     axis.text.y = element_text(size = (7.5))
## )


## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------
## library("rlang")
## 
## ## Plot of gene expression lognorm counts vs. sample variable
## plot_gene_expr <- function(sample_var, gene_id) {
##     colors <- list(
##         "Group" = c("Control" = "brown2", "Experimental" = "deepskyblue3"),
##         "Age" = c("Adult" = "slateblue3", "Pup" = "yellow3"),
##         "Sex" = c("F" = "hotpink1", "M" = "dodgerblue"),
##         "Pregnancy" = c("Yes" = "darkorchid3", "No" = "darkolivegreen4"),
##         "plate" = c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1"),
##         "flowcell" = c(
##             "HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta", "HKCNKDSXX" = "turquoise3",
##             "HKCTMDSXX" = "tomato"
##         )
##     )
## 
##     ## Lognorm counts of the gene across samples
##     data <- colData(rse_gene_filt)
##     data$gene_expr <- assays(rse_gene_filt)$logcounts[gene_id, ]
## 
##     ## Percentage of variance explained by the variable
##     percentage <- 100 * signif(varPart_data[gene_id, sample_var], digits = 3)
## 
##     ## Boxplots for categorical variables
##     if (class(data[, sample_var]) == "character") {
##         plot <- ggplot(data = as.data.frame(data), mapping = aes(
##             x = !!rlang::sym(sample_var),
##             y = gene_expr, color = !!rlang::sym(sample_var)
##         )) +
##             geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
##             geom_jitter(width = 0.15, alpha = 1, size = 1) +
##             stat_smooth(geom = "line", alpha = 0.6, size = 0.4, span = 0.3, method = "lm", aes(group = 1), color = "orangered3") +
##             scale_color_manual(values = colors[[sample_var]]) +
##             theme_bw() +
##             guides(color = "none") +
##             labs(
##                 title = gene_id,
##                 subtitle = paste0("Variance explained: ", percentage, "%"),
##                 y = "lognorm counts", x = sample_var
##             ) +
##             theme(
##                 axis.title = element_text(size = (7)),
##                 axis.text = element_text(size = (6)),
##                 plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
##                 plot.subtitle = element_text(size = 7, color = "gray40"),
##                 legend.text = element_text(size = 6),
##                 legend.title = element_text(size = 7)
##             )
##     }
## 
##     ## Scatterplots for continuous variables
##     else {
##         colors <- c(
##             "mitoRate" = "khaki3", "overallMapRate" = "turquoise", "totalAssignedGene" = "plum2", "rRNA_rate" = "orange3",
##             "sum" = "palegreen3", "detected" = "skyblue2", "ERCCsumLogErr" = "slateblue1"
##         )
## 
##         plot <- ggplot(as.data.frame(data), aes(x = eval(parse_expr(sample_var)), y = gene_expr)) +
##             geom_point(color = colors[[sample_var]], size = 2) +
##             stat_smooth(geom = "line", alpha = 0.4, size = 0.4, span = 0.25, method = "lm", color = "orangered3") +
##             theme_bw() +
##             guides(color = "none") +
##             labs(
##                 title = gene_id,
##                 subtitle = paste0("Variance explained: ", percentage, "%"),
##                 y = "lognorm counts", x = gsub("_", " ", sample_var)
##             ) +
##             theme(
##                 plot.margin = unit(c(0.4, 0.1, 0.4, 0.1), "cm"),
##                 axis.title = element_text(size = (7)),
##                 axis.text = element_text(size = (6)),
##                 plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
##                 plot.subtitle = element_text(size = 7, color = "gray40"),
##                 legend.text = element_text(size = 6),
##                 legend.title = element_text(size = 7)
##             )
##     }
## 
##     return(plot)
## }


## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------
## ## Function to plot gene expression vs sample variable data for top 3 most affected genes
## 
## plot_gene_expr_sample <- function(sample_var) {
##     ## Top 3 genes most affected by sample variable
##     affected_genes <- rownames(varPart_data[order(varPart_data[, sample_var], decreasing = TRUE), ][1:3, ])
## 
##     ## Plots
##     plots <- list()
##     for (i in 1:length(affected_genes)) {
##         plots[[i]] <- plot_gene_expr(sample_var, affected_genes[i])
##     }
##     plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
## }


## ----message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------
## ## Plots for top affected genes by 'overallMapRate'
## plots <- plot_gene_expr_sample("overallMapRate")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'totalAssignedGene'
## plots <- plot_gene_expr_sample("totalAssignedGene")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'Group'
## plots <- plot_gene_expr_sample("Group")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))
## 
## ## Plots for top affected genes by 'Sex' (genes in sexual chrs)
## plots <- plot_gene_expr_sample("Sex")
## plots + theme(plot.margin = unit(c(3, 1, 2, 3), "cm"))


## ----exercise1_varPart, message=FALSE, warning=FALSE, echo=FALSE, eval=FALSE----------
## ## Solution
## 
## ## Gene ID
## gene_id <- "ENSMUSG00000042348.10"
## ## % of variance explained by Group
## percentage <- 100 * signif(varPart_data[gene_id, "Group"], digits = 3)
## ## Sample colors
## colors <- c("Control" = "brown2", "Experimental" = "deepskyblue3")
## ## Gene expression logcounts
## rse_gene_filt$gene_expr <- assays(rse_gene_filt)$logcounts[gene_id, ]
## 
## ## Plot
## plot <- ggplot(
##     data = as.data.frame(colData(rse_gene_filt)),
##     mapping = aes(x = Group, y = gene_expr, color = Group)
## ) +
##     geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
##     geom_jitter(width = 0.15, alpha = 1, size = 1) +
##     scale_color_manual(values = colors) +
##     theme_bw() +
##     guides(color = "none") +
##     labs(
##         title = gene_id,
##         subtitle = paste0("Variance explained: ", percentage, "%"),
##         y = "lognorm counts"
##     ) +
##     theme(
##         plot.margin = unit(c(2, 6, 2, 6), "cm"),
##         axis.title = element_text(size = (7)),
##         axis.text = element_text(size = (6)),
##         plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
##         plot.subtitle = element_text(size = 7, color = "gray40"),
##         legend.text = element_text(size = 6),
##         legend.title = element_text(size = 7)
##     )
## 
## plot


## ----exercise2_varPart, message=FALSE, warning=FALSE, echo=FALSE, eval=FALSE----------
## ## Solution
## 
## ## Gene ID
## gene_id <- "ENSMUSG00000064372.1"
## ## % of variance explained by Group
## percentage <- 100 * signif(varPart_data[gene_id, "Group"], digits = 3)
## ## Sample colors
## colors <- c("Control" = "brown2", "Experimental" = "deepskyblue3")
## ## Gene expression logcounts
## rse_gene_filt$gene_expr <- assays(rse_gene_filt)$logcounts[gene_id, ]
## 
## ## Plot
## plot <- ggplot(
##     data = as.data.frame(colData(rse_gene_filt)),
##     mapping = aes(x = Group, y = gene_expr, color = Group)
## ) +
##     geom_boxplot(size = 0.25, width = 0.32, color = "black", outlier.color = "#FFFFFFFF") +
##     geom_jitter(width = 0.15, alpha = 1, size = 1) +
##     scale_color_manual(values = colors) +
##     theme_bw() +
##     guides(color = "none") +
##     labs(
##         title = gene_id,
##         subtitle = paste0("Variance explained: ", percentage, "%"),
##         y = "lognorm counts"
##     ) +
##     theme(
##         plot.margin = unit(c(2, 6, 2, 6), "cm"),
##         axis.title = element_text(size = (7)),
##         axis.text = element_text(size = (6)),
##         plot.title = element_text(hjust = 0.5, size = 7.5, face = "bold"),
##         plot.subtitle = element_text(size = 7, color = "gray40"),
##         legend.text = element_text(size = 6),
##         legend.title = element_text(size = 7)
##     )
## 
## plot

