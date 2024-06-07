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


## ----download_iSEEdata----------------------------------------------------------------------
## Lets get some data using spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

## We can check how big the object is with lobstr
lobstr::obj_size(sce_layer)


## -------------------------------------------------------------------------------------------
curl::curl_version()$version


## ## Install homebrew from https://brew.sh/
## brew install curl

## ----eval = FALSE---------------------------------------------------------------------------
## Sys.setenv(PKG_CONFIG_PATH = "/opt/homebrew/opt/curl/lib/pkgconfig")
## install.packages("curl", type = "source")


## ----explore_iSEE, eval = FALSE-------------------------------------------------------------
## ## Load library
## library("iSEE")
## 
## ## Deploy
## iSEE(sce_layer)

