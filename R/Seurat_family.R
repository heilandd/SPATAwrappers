#' @title  Juxtaposition in stRNA-seq
#' @author Dieter Henrik Heiland
#' @description Juxtaposition in stRNA-seq
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
run.SNN.stability <- function(object, 
                              assay="RNA", 
                              reduction="pca", 
                              dims = 1:30,
                              resolution = seq(from=0.1, to=1.5, by=0.025), 
                              cluster_id="new", 
                              algorithm=1,
                              plot=T){
  
  #SNN-Analysis:
  # Find optimal resolution
  Seurat::DefaultAssay(object) <- assay
  
  object <- Seurat::FindNeighbors(object, assay=assay, dims = dims, reduction=reduction)
  if(is.null(object@graphs$SCT_snn)){
    object@graphs$SCT_snn <- object@graphs$integrated_snn
    object@graphs$SCT_nn <- object@graphs$integrated_nn
  }
  object <- Seurat::FindClusters(object, resolution = resolution, algorithm=algorithm)
  
  object@meta.data
  # Validate resolution by cluster stability and Recall
  
  prefix=paste0(assay,"_snn_res.")
  
  library(tidyverse)
  library(ggraph)
  clt <- clustree::clustree(object, prefix = prefix, node_colour = "sc3_stability")
  clt.plot <- clt$data %>% 
    dplyr::select(!!sym(prefix), sc3_stability) %>% 
    dplyr::group_by(!!sym(prefix)) %>% 
    dplyr::summarise(mean(sc3_stability)) %>% 
    as.data.frame()
  names(clt.plot)[2] <- "sc3_stability"
  clt.plot[,prefix]<-as.numeric(as.character(clt.plot[,prefix]))
  
  if(plot==T){
  val_plot=ggplot(clt.plot, mapping=aes(x=!!sym(prefix), y=sc3_stability))+
    geom_point()+
    theme_classic()+ 
    geom_line()+
    geom_vline(xintercept = clt.plot %>% arrange(desc(sc3_stability)) %>% head(1) %>% pull(!!sym(prefix)),
               color="red", size=2)
  
  print(val_plot)
  }
  
  #Optimal Cluster and Run cluster Analysis
  object <- Seurat::FindClusters(object, resolution = clt.plot %>% dplyr::arrange(desc(sc3_stability)) %>% head(1) %>% dplyr::pull(!!sym(prefix)))
  
  print(paste0("Optimal Cluster resolution: ", clt.plot %>% dplyr::arrange(desc(sc3_stability)) %>% head(1) %>% dplyr::pull(!!sym(prefix)) ))
  
  cr.op <- 
    paste0(prefix, 
           clt.plot %>% 
             dplyr::arrange(desc(sc3_stability)) %>% 
             head(1) %>% 
             dplyr::pull(!!sym(prefix))
  )
  
  object@meta.data$new <- object@meta.data[,cr.op] 
  names(object@meta.data)[names(object@meta.data)=="new"] <- cluster_id
  
  return(object)
}



