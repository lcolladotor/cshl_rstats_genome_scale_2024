## ----"download_modeling_results"----------------------------------------------------
## get reference layer enrichment statistics
layer_modeling_results <- spatialLIBD::fetch_data(type = "modeling_results")


## ----"layer_modeling_results_workaround"--------------------------------------------
tmp_modeling_results <- tempfile("modeling_results.RData")
download.file(
    "https://www.dropbox.com/s/se6rrgb9yhm5gfh/Human_DLPFC_Visium_modeling_results.Rdata?dl=1",
    tmp_modeling_results,
    mode = "wb"
)
load(tmp_modeling_results, verbose = TRUE)

## Let's rename the object into the name used in the
## spatial registration vignette (from spatialLIBD)
layer_modeling_results <- modeling_results

