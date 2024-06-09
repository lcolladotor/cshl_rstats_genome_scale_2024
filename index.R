## ----install, eval = FALSE------------------------------------------------------
## ## For installing Bioconductor packages
## if (!requireNamespace("BiocManager", quietly = TRUE)) {
##     install.packages("BiocManager")
## }
## 
## ## Install required packages
## BiocManager::install(
##     c(
##         "usethis", ## Utilities
##         "BiocFileCache",
##         "RefManageR",
##         "gitcreds",
##         "gert",
##         "gh",
##         "here",
##         "Hmisc",
##         "biocthis",
##         "lobstr",
##         "postcards",
##         "scater",
##         "sessioninfo",
##         "stringr",
##         "SummarizedExperiment", ## Main containers / vis
##         "iSEE",
##         "edgeR", ## RNA-seq
##         "ExploreModelMatrix",
##         "limma",
##         "smokingMouse",
##         "recount3",
##         "rlang",
##         "scRNAseq",
##         "pheatmap", ## Visualization
##         "ggplot2",
##         "ggrepel",
##         "patchwork",
##         "RColorBrewer",
##         "ComplexHeatmap",
##         "cowplot",
##         "Polychrome",
##         "spatialLIBD", ## Advanced
##         "variancePartition"
##     )
## )


## ----session_packages, eval = TRUE, message = FALSE-----------------------------
## Load the package at the top of your script
library("sessioninfo")

## Utilities
library("BiocFileCache")
library("BiocStyle")
library("biocthis")
library("gitcreds")
library("gert")
library("gh")
library("here")
library("lobstr")
library("postcards")
library("usethis")
library("sessioninfo")

## Data
library("smokingMouse")
library("scRNAseq")

## Main containers / vis
library("SummarizedExperiment")
library("iSEE")

## RNA-seq
library("edgeR")
library("ExploreModelMatrix")
library("limma")
library("recount3")

## QCA
library("scater")

## Variance Partition
library("variancePartition")

## Visualization: plots & text
library("ComplexHeatmap")
library("ggplot2")
library("patchwork")
library("pheatmap")
library("RColorBrewer")
library("Hmisc")
library("stringr")
library("cowplot")
library("rlang")
library("ggrepel")
library("Polychrome")

## Spatial transcriptomics
library("spatialLIBD")


## ----session_info---------------------------------------------------------------
## Reproducibility information
options(width = 120)
session_info()
proc.time()


## ----"check_curl"---------------------------------------------------------------
curl::curl_version()

