rmds <- dir(pattern = "^[[:digit:]]+.*\\.Rmd$")
rmd_num <- as.integer(substr(rmds, 0, 2))
start_i <- 9
i <- seq(start_i, length(rmds))
rmd_num[i] <- rmd_num[i] - 1

rmds_new <-
    paste0(formatC(
        rmd_num,
        width = 2,
        format = "d",
        flag = "0"
    ),
        gsub("^[[:digit:]]+", "", rmds))

data.frame(rmds, rmds_new)

mapply(file.rename, rmds, rmds_new)
