# Create package
#devtools::create('SPATAwrappers')

library(devtools)
setwd("~/Desktop/Projekt_Metabolom/Bioinformatic Tests/SPATAwrappers")

# Updata namespace
document() #<- Namespace usw
load_all() # <- library
library(SPATAwrappers)


devtools::install_github("heilandd/SPATAwrappers")

object@spatial$`334_T`$Cell_coords
