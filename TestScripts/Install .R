# Create package
#devtools::create('SPATAwrappers')

library(devtools)
setwd("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/SPATAwrappers")

# Updata namespace
document() #<- Namespace usw
load_all() # <- library
library(SPATAwrappers)


devtools::install_github("heilandd/SPATAwrappers")

SPATAwrappers::

#### GIT Tocken

usethis::use_git_config(user.name = "heilandd", user.email = "dieter.henrik.heiland@uniklinik-freiburg.de")
## create a personal access token for authentication:
usethis::create_github_token() 

## set personal access token:
credentials::set_github_pat("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")

## or store it manually in '.Renviron':
usethis::edit_r_environ()
## store your personal access token with: GITHUB_PAT=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
## and make sure '.Renviron' ends with a newline
#Check 
usethis::git_sitrep()

#New Tocken
#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
#R:GITHUB_PAT







