## ----load_iSEE-------------------------------------------
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("iSEE")

packageVersion("iSEE")
library("iSEE")


## ----vignettes_iSEE, eval=FALSE--------------------------
## browseVignettes("iSEE")


## ----quick_launch, eval=interactive()--------------------
## Launch iSEE for the se ("SummarizedExperiment" object)
iSEE(se)

## Launch iSEE for the sce ("SingleCellExperiment" object)
iSEE(sce)


## ----download_sce_layer----------------------------------
## Lets get some data using spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

## We can check how big the object is with lobstr
lobstr::obj_size(sce_layer)


## --------------------------------------------------------
curl::curl_version()$version


## ## Install homebrew from https://brew.sh/
## brew install curl

## ----eval = FALSE----------------------------------------
## Sys.setenv(PKG_CONFIG_PATH = "/opt/homebrew/opt/curl/lib/pkgconfig")
## install.packages("curl", type = "source")


## ----"sce_layer_manual_workaround"-----------------------
tmp_sce_layer <- tempfile("sce_layer.RData")
download.file(
    "https://www.dropbox.com/s/bg8xwysh2vnjwvg/Human_DLPFC_Visium_processedData_sce_scran_sce_layer_spatialLIBD.Rdata?dl=1",
    tmp_sce_layer,
    mode = "wb"
)
load(tmp_sce_layer, verbose = TRUE)
sce_layer


## ----explore_iSEE, eval = FALSE--------------------------
## ## Load library
## library("iSEE")
## 
## ## Deploy
## iSEE(sce_layer)


## ----load_scRNAseq_data----------------------------------
library("scRNAseq")
library("scater")
library("iSEE")

# Load the dataset
sce <- ReprocessedAllenData(assays = "tophat_counts")

# Normalize counts and perform PCA
sce <- logNormCounts(sce, exprs_values = "tophat_counts")
sce <- runPCA(sce, ncomponents = 4)


## ----single_geneExpr-------------------------------------
## Initial settings for a single gene expression
initial_single <- list(
    FeatureAssayPlot(Assay = "logcounts", YAxisFeatureName = "Serpine2"),
    ReducedDimensionPlot(Type = "PCA", ColorBy = "Column selection", ColumnSelectionSource = "FeatureAssayPlot1")
)

## Launch iSEE with the initial settings
if (interactive()) {
    iSEE(sce, initial = initial_single)
}


## ----2_genesExpr-----------------------------------------
## Initial settings for 2 genes expression on the same "FeatureAssayPlot"
initial_combined <- list(
    FeatureAssayPlot(Assay = "logcounts", XAxis = "Feature name", XAxisFeatureName = "Serpine2", YAxisFeatureName = "Bcl6"),
    ReducedDimensionPlot(Type = "PCA", ColorBy = "Column selection", ColumnSelectionSource = "FeatureAssayPlot1")
)

## Launch iSEE with the initial settings
if (interactive()) {
    iSEE(sce, initial = initial_combined)
}


## ----chain_FeatureAssayPlots-----------------------------
## Initial settings chainning multiple "FeatureAssayPlot"
initial_double <- list(
    FeatureAssayPlot(Assay = "logcounts", YAxisFeatureName = "Serpine2"),
    FeatureAssayPlot(Assay = "logcounts", YAxisFeatureName = "Bcl6", ColumnSelectionSource = "FeatureAssayPlot1", ColumnSelectionRestrict = TRUE),
    ReducedDimensionPlot(Type = "PCA", ColorBy = "Column selection", ColumnSelectionSource = "FeatureAssayPlot2")
)

## Launch iSEE with the initial settings
if (interactive()) {
    iSEE(sce, initial = initial_double)
}

