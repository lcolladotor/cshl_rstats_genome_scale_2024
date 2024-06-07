## ----eval=FALSE, warning=FALSE, message=FALSE-----------------------------------------------
## ## Install the package manually
## # install.packages("here")
## 
## ## Load "here" (previously installed)
## library("here")


## ----here::here, eval=FALSE, warning=FALSE, message=FALSE-----------------------------------
## here::here()


## ----check_paths, eval=FALSE, warning=FALSE, message=FALSE----------------------------------
## getwd() # returns the current path
## setwd("desired/directory") # changes to the specified path


## ----bestpractice, eval=FALSE, warning=FALSE, message=FALSE---------------------------------
## ## Instead of "C:/Users/user/Desktop/data/myfile.csv"
## 
## ## Use here to construct file paths
## file_path <- here("Users", "user", "Desktop", "data", "myfile.csv")
## # file_path <- here:here("Users", "user", "Desktop","data", "myfile.csv")
## data <- read.csv(file_path)


## ----here_examples, eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## ## Example: save data to a file and load it
## a <- 1
## c <- 23
## 
## save(a, c, file = here("test-data.RData"))
## # save(a, c, file = here:here("test-data.RData"))
## load(here("test-data.RData"))
## # load(here:here("test-data.RData"))
## 
## ## Create a directory
## dir.create(here("subdirectory"), showWarnings = FALSE)
## # dir.create(here:here("subdirectory"), showWarnings = FALSE)
## 
## ## Create a file, indicating the subdirectory (the first argument in this case)
## file.create(here("subdirectory", "filename"))
## # file.create(here:here("subdirectory", "filename"))
## 
## ## Open the new created file
## file.show(here("subdirectory", "filename"))
## # file.show(here:here("subdirectory", "filename"))
## 
## ## For example, if we want to see our files in the directory
## list.files(here(), recursive = TRUE)
## # list.files(here:here(), recursive = TRUE)


## ----load_usethis, eval=FALSE, warning=FALSE, message=FALSE---------------------------------
## ## Install the package manually
## # install.packages("usethis")
## 
## ## Load "usethis (previously installed)
## library("usethis")


## ----use_*, eval=FALSE, warning=FALSE, message=FALSE----------------------------------------
## ## usethis::use_*()
## usethis::use_r()
## usethis::use_git()
## usethis::use_readme_md()


## ----usethis_example, eval=FALSE, warning=FALSE, message=FALSE------------------------------
## ## For example, create a README file
## usethis::use_readme_md()
## #> ✔ Writing 'README.md'
## #> ● Edit 'README.md'


## ----install_gitreq, eval=FALSE, warning=FALSE, message=FALSE-------------------------------
## ## Packages we will need
## install.packages(c("gitcreds", "gert", "gh"))
## 
## ## Load them separately
## library("gitcreds")
## library("gert")
## library("gh")


## ----token, eval=FALSE, warning=FALSE, message=FALSE----------------------------------------
## ## Initiate connection with GitHub
## usethis::create_github_token() # redirects to GitHub where you'll choose a specific name for the token


## ----gitcreds, eval=FALSE, warning=FALSE, message=FALSE-------------------------------------
## gitcreds::gitcreds_set() # here you place the token (NOT your GitHub password!!!)


## ----git_config, eval=FALSE, warning=FALSE, message=FALSE-----------------------------------
## ## Configure GitHub user
## usethis::edit_git_config() # opens the .gitconfig file
## 
## ## Place the name and email of your GitHub account.
## ## JUST remove the "#" and respect the other spaces
## 
## # [user]
## #   name = N A M E
## #   email = github_email


## ----git_repo, eval=FALSE, warning=FALSE, message=FALSE-------------------------------------
## ## Initialize the Git repository
## usethis::use_git()
## 
## ## Connect your local Git repository with GitHub servers
## usethis::use_github()


## ----gh_whoami, eval=FALSE, warning=FALSE, message=FALSE------------------------------------
## gh::gh_whoami()


## ----git_commands, eval=FALSE, warning=FALSE, message=FALSE---------------------------------
## ## Write a new file, using here::here to specify the path
## writeLines("hello", here::here("R", "test-here.R"))
## 
## ## Another way is to use use_r
## usethis::use_r("test-file-github.R") # adds file to the project's R directory
## 
## ## For example, we might try adding something new
## gert::git_add("R/test-file-github.R")
## 
## ## Add commit of what was done
## gert::git_commit("uploaded test file")
## 
## ## Gives info about the commits
## gert::git_log()
## 
## ## Upload your changes from the local repo to GitHub
## gert::git_push() # IMPORTANT COMMAND


## ----render_site, eval=FALSE, warning=FALSE, message=FALSE----------------------------------
## rmarkdown::render_site()


## ----install_postcards, eval=FALSE, warning=FALSE, message=FALSE----------------------------
## ## You can install Postcards with the following command:
## # install.packages("postcards")
## 
## ## Or you can install the latest development version (not recommended):
## # remotes::install_github("seankross/postcards@main")


## ----postcards_project, eval=FALSE, warning=FALSE, message=FALSE----------------------------
## ## Create a new project
## usethis::create_project("Your_Username.github.io")


## ----setup_gitpcds, eval=FALSE, warning=FALSE, message=FALSE--------------------------------
## ## Set up Git and GitHub
## usethis::use_git() # Restart the session
## usethis::use_github()


## ----choose_template, eval=FALSE, warning=FALSE, message=FALSE------------------------------
## ## Choose only one template (the one you like the most)
## postcards::create_postcard(template = "jolla")
## postcards::create_postcard(template = "jolla-blue")
## postcards::create_postcard(template = "trestles")
## postcards::create_postcard(template = "onofre")
## postcards::create_postcard(template = "solana")


## ----deploy_postcards, eval=FALSE, warning=FALSE, message=FALSE-----------------------------
## ## Deploy the GitHub page
## rmarkdown::render("index.Rmd")

