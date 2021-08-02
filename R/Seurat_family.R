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



#' @title  MergeInferCNVSeurat
#' @author Dieter Henrik Heiland
#' @description MergeInferCNVSeurat
#' @param object Seurat Object
#' @param results The InferCNV results folder
#' @inherit 
#' @return 
#' @examples 
#' @export

MergeInferCNVSeurat <- function(object, results=getwd()){
  
  library(infercnv)
  library(SPATA2)
  library(tidyverse)
  
  infer.obj <- readRDS(stringr::str_c(results, "/infercnv-obj.RDS"))
  gene_pos_df <-
    infer.obj@gene_order %>% 
    dplyr::select(chr, start,stop) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("hgnc_symbol") %>% 
    dplyr::rename(chromosome_name:=chr)
  
  result_dir <-stringr::str_c(results, "/infercnv.observations.txt")
  results <- utils::read.table(result_dir)
  barcodes <- base::colnames(results)
  confuns::give_feedback(msg = "Summarizing cnv-results by chromosome.")
  
  
  
  
  # join cnv results (per gene) with chromosome positions and summarize by chromosome
  cnv_prefix="Chr"
  
  ordered_cnv_df <-
    base::as.data.frame(results) %>%
    tibble::rownames_to_column("hgnc_symbol") %>%
    dplyr::left_join(., gene_pos_df, by = "hgnc_symbol") %>%
    dplyr::group_by(chromosome_name) %>%
    dplyr::select(chromosome_name, dplyr::any_of(barcodes)) %>%
    dplyr::summarise_all(base::mean) %>%
    dplyr::mutate(Chr = stringr::str_c(cnv_prefix, chromosome_name)) %>%
    dplyr::select(Chr, dplyr::everything()) %>%
    dplyr::ungroup() %>%
    dplyr::select(-chromosome_name)
  
  cnames <- c("barcodes", ordered_cnv_df$Chr)
  
  ordered_cnv_df2 <-
    dplyr::select(ordered_cnv_df, -Chr) %>%
    base::t() %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(var = "barcodes") %>%
    magrittr::set_colnames(value = cnames) %>%
    dplyr::mutate(barcodes = stringr::str_replace_all(string = barcodes, pattern = "\\.", replacement = "-")) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(match = cnv_prefix), .fns = base::as.numeric)) %>%
    tibble::as_tibble()
  
  # add results to spata object
  confuns::give_feedback(msg = "Adding results to the Seurat object's feature data.")
  
  # feature variables
  
  join <- left_join(object@meta.data %>% rownames_to_column("barcodes") , ordered_cnv_df2 , by="barcodes")
  object@meta.data <- cbind(object@meta.data, join[,ordered_cnv_df$Chr])
  
  return(object)
  
}


