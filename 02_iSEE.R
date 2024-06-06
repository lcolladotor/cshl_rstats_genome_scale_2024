## ----load_iSEE, eval=FALSE, warning=FALSE, message=FALSE------------------------------------
## # if (!require("BiocManager", quietly = TRUE))
## #     install.packages("BiocManager")
## #
## # BiocManager::install("iSEE")
## 
## library(iSEE)


## ----vignettes_iSEE, eval=FALSE, warning=FALSE, message=FALSE-------------------------------
## browseVignettes("iSEE")


## ----quick_launch, eval=FALSE, warning=FALSE, message=FALSE---------------------------------
## ## Launch iSEE for the se ("SummarizedExperiment" object)
## iSEE(se)
## 
## ## Launch iSEE for the sce ("SingleCellExperiment" object)
## iSEE(sce)


## ----download_iSEEdata, eval=FALSE, warning=FALSE, message=FALSE----------------------------
## ## Lets get some data using spatialLIBD
## sce_layer <- spatialLIBD::fetch_data("sce_layer")
## sce_layer
## # class: SingleCellExperiment
## # dim: 22331 76
## # metadata(0):
## # assays(2): counts logcounts
## # rownames(22331): ENSG00000243485 ENSG00000238009 ... ENSG00000278384 ENSG00000271254
## # rowData names(10): source type ... is_top_hvg is_top_hvg_sce_layer
## # colnames(76): 151507_Layer1 151507_Layer2 ... 151676_Layer6 151676_WM
## # colData names(13): sample_name layer_guess ... layer_guess_reordered_short spatialLIBD
## # reducedDimNames(6): PCA TSNE_perplexity5 ... UMAP_neighbors15 PCAsub
## # mainExpName: NULL
## # altExpNames(0):
## # 33.99 MB
## 
## ## We can check how big the object is with lobstr
## lobstr::obj_size(sce_layer)
## # 33.99 MB


## -------------------------------------------------------------------------------------------
curl::curl_version()$version


## ## Install homebrew from https://brew.sh/
## brew install curl

## -------------------------------------------------------------------------------------------
Sys.setenv(PKG_CONFIG_PATH = "/opt/homebrew/opt/curl/lib/pkgconfig")
install.packages("curl", type = "source")


## ----explore_iSEE, eval = FALSE-------------------------------------------------------------
## ## Load library
## library("iSEE")
## 
## ## Deploy
## iSEE(sce_layer)

