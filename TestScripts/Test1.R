

# CellProfiler Functions --------------------------------------------------
library(tidyverse)

#Import Feature data
data1 <- CellProfilerImport("MyExpt_FilterObjects.csv" )

#Import Object data
data <- rbind(
  CellProfilerImportObjectPositions("MyExpt_IdentifyPrimaryObjects_T.csv", type="T"),
  CellProfilerImportObjectPositions("MyExpt_IdentifyPrimaryObjects_M.csv", type="M"),
  CellProfilerImportObjectPositions("MyExpt_IdentifyPrimaryObjects_Tumor.csv", type="TT")
)

#Quandify distance
data.dist <- CellProfilerCellDistance(data)

#Plot Distance plot
p=CellProfilerPlotDistance(data.dist, image = 1, dist.threshold = 300, pt.size = 6, lt.size = 4, lt.alpha = 0.3)
CellProfilerPlotDistanceAddLayer(p,image = 1, data = data1, color.to = "HMOX")


#Get cell type
# Example HMOX+/- stained cells

mean <- data1[[1]] %>% pull(HMOX) %>% mean()
data1[[1]] <- 
data1[[1]] %>% mutate(type=ifelse(HMOX>mean, "H+", "H-"))
CellProfilerPlotDistanceAddLayer(p,image = 1, data = data1, color.to = "type")

# Distance plots (density)
plot <- CellProfilerDistancePlot(data=data1, prefix="M", data.dist, exclude="M")
library(patchwork)
plot$T+plot$TT



# Test SPATA juxtaposition functions --------------------------------------

samplr.read <- "275_T"
object.spata <- readRDS(paste0("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions/", samplr.read, "_SPATA_CNV_Pred.RDS"))
#spata.ind.plot <-   runAutoencoderDenoising(spata.ind.plot, activation="relu", bottleneck=c(32), layers = c(128, 64, 32), dropout = 0.1, epochs = 20)
object.spata@used_genesets <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/GS_new.RDS")


object <- object.spata

#Position connection

parameter= "HM_G2M_CHECKPOINT"
parameter = "HM_HYPOXIA"
df <- 
  SPATA2::joinWithGeneSets(object, gene_sets = parameter) %>% 
  as.data.frame() %>% 
  dplyr::select(barcodes,x,y,{{parameter}})


VF <- SPATAwrappers::inferVectorFields(object, df, parameter, dist.spot=60)


plotJuxtapositionSPATA(object, feature.source="MES", feature.target = "HM_HYPOXIA",thr=0.9, dist=100, lt.alpha = 0.1)
p1 <- plotJuxtapositionSPATA(object, feature.source="MES", feature.target = "HM_HYPOXIA",thr=0.9, dist=100, add=T, color="black")
p2 <- plotJuxtapositionSPATA(object, feature.source="MES", feature.target = "AC",thr=0.85, dist=200, add=T, color="red")
stream <- SPATAwrappers::plotStreamlines(VF, parameter, res=2,skip = 2,pt.alpha = 0.5, color.extern = object %>% SPATA2::getFeatureDf() %>% pull(seurat_clusters))

stream+p2+p1

data <- plotJuxtapositionSPATA(object, feature.source="MES", feature.target = "HM_HYPOXIA",thr=0.9, dist=100, data=T)


a=matrix(1,3,3)
b=matrix(2,3,3)
c=matrix(3,3,3)

e <- array(c(a,b,c), dim=c(3,3,3))
dim(e)




# get VF
parameter= "HM_G2M_CHECKPOINT"
parameter = "HM_HYPOXIA"
df <- 
  SPATA2::joinWithGeneSets(object, gene_sets = parameter) %>% 
  as.data.frame() %>% 
  dplyr::select(barcodes,x,y,{{parameter}})


VF <- SPATAwrappers::inferVectorFields(object, df, parameter, dist.spot=60)
SPATAwrappers::plotStreamlines(VF, parameter, res=2)
SPATAwrappers::plotVectorFields(VF, parameter, skip=1)
  
  
  




inferSpatial.plot(object, feature = c("AC", "MES", "NPC", "OPC"))

object.spata %>% SPATA2::plotCnvResults(., across="seurat_clusters")


object@cnv$`275_T`$regions_df

saveRDS(object,"UCSF.RDS")

plotCNV(object, across="CAF")

# IMC data Functions ------------------------------------------------------



  















