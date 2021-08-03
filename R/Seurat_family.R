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

MergeInferCNVSeurat <- function(object, results=getwd(), remove.prefix=NULL){
  
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
  if(!is.null(remove.prefix)){barcodes <- stringr::str_remove_all(barcodes, "X");colnames(results) <-barcodes }
  
  
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
  
  object@meta.data <- dplyr::left_join(object@meta.data %>% tibble::rownames_to_column("barcodes") , ordered_cnv_df2 , by="barcodes")
  
  # cnv matrix
  base::colnames(results) <- stringr::str_replace_all(string = base::colnames(results), pattern = "\\.", replacement = "-")
  cnv_mtr <- base::as.matrix(results)
  
  cnv_list <- list(prefix = cnv_prefix,
                   gene_pos_df = gene_pos_df,
                   cnv_df = ordered_cnv_df2,
                   cnv_mtr = cnv_mtr,
                   regions_df = SPATA2::cnv_regions_df
  )
  
  
  object@assays$CNV <- cnv_list
  
  return(object)
  
}

#' @title  plotCNV
#' @author Dieter Henrik Heiland
#' @description plotCNV
#' @param object Seurat Object
#' @param results The InferCNV results folder
#' @inherit 
#' @return 
#' @examples 
#' @export
#' 
plotCNV <- function(object,
                    across = NULL,
                    across_subset = NULL,
                    relevel = NULL,
                    smooth_span = 0.08,
                    nrow = NULL,
                    ncol = NULL,
                    clr = "blue",
                    of_sample = NA,
                    verbose = NULL
){
  
  if(is.null(object@assays$CNV)) stop("No CNV assay stored in Seurat object, call ::MergeInferCNVSeurat ")
  
  # cnv results
  cnv_results <- object@assays$CNV
  cnv_data <- cnv_results$cnv_mtr
  
  if(base::is.null(across)){
    
    confuns::give_feedback(msg = "Plotting cnv-results for whole sample.", verbose = verbose)
    
    plot_df <-
      base::data.frame(
        mean = base::apply(cnv_data, MARGIN = 1, FUN = stats::median),
        sd = base::apply(cnv_data, MARGIN = 1, FUN = stats::sd)
      ) %>%
      tibble::rownames_to_column(var = "hgnc_symbol") %>%
      dplyr::left_join(x = ., y = cnv_results$gene_pos_df, by = "hgnc_symbol") %>%
      dplyr::mutate(
        chromosome_name = base::factor(chromosome_name, levels = base::as.character(0:23))
      ) %>%
      tibble::as_tibble()
    
    line_df <-
      dplyr::count(x = plot_df, chromosome_name) %>%
      dplyr::mutate(
        line_pos = base::cumsum(x = n),
        line_lag = dplyr::lag(x = line_pos, default = 0) ,
        label_breaks = (line_lag + line_pos) / 2
      ) %>%
      tidyr::drop_na()
    
    final_plot <-
      ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = 1:base::nrow(plot_df), y = mean)) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, span = 0.08, se = FALSE, color = clr) +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = mean-sd, ymax = mean + sd),
                           alpha = 0.2) +
      ggplot2::geom_vline(data = line_df, mapping = ggplot2::aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(breaks = line_df$label_breaks, labels = line_df$chromosome_name) +
      ggplot2::labs(x = "Chromosomes", y = "CNV-Results")
    
  } else if(base::is.character(across)){
    
    confuns::give_feedback(
      msg = glue::glue("Plotting cnv-results across '{across}'. This might take a few moments."),
      verbose = verbose
    )
    
    gene_names <- base::rownames(cnv_data)
    
    prel_df <-
      base::as.data.frame(cnv_data) %>%
      base::t() %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column(var = "barcodes") %>% 
      dplyr::left_join(., object@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("barcodes") %>% select(barcodes, {{across}}), by="barcodes") %>% 
      confuns::check_across_subset(df = ., across = across, across.subset = across_subset, relevel = relevel) %>% 
      tidyr::pivot_longer(
        cols = dplyr::all_of(gene_names),
        names_to = "hgnc_symbol",
        values_to = "cnv_values"
      ) %>%
      dplyr::left_join(x = ., y = cnv_results$gene_pos_df, by = "hgnc_symbol") %>%
      dplyr::mutate(
        chromosome_name = base::factor(chromosome_name, levels = base::as.character(0:23))
      ) %>%
      tibble::as_tibble()
    
    summarized_df <-
      dplyr::group_by(prel_df, !!rlang::sym(x = across), chromosome_name, hgnc_symbol) %>%
      dplyr::summarise(
        cnv_mean = stats::median(x = cnv_values, na.rm = TRUE),
        cnv_sd = stats::sd(x = cnv_values, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(!!rlang::sym(x = across)) %>%
      dplyr::mutate(x_axis = dplyr::row_number(),across = !!rlang::sym(across))
    
    line_df <-
      dplyr::count(x = summarized_df, chromosome_name) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(!!rlang::sym(across)) %>%
      dplyr::mutate(
        line_pos = base::cumsum(x = n),
        line_lag = dplyr::lag(x = line_pos, default = 0) ,
        label_breaks = (line_lag + line_pos) / 2
      )  %>%
      tidyr::drop_na()
    

    final_plot <- 
      ggplot2::ggplot(data = summarized_df, mapping = ggplot2::aes(x = x_axis, y = cnv_mean)) + 
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, span = smooth_span, se = FALSE, color = clr) + 
      ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = cnv_mean - cnv_sd, ymax = cnv_mean + cnv_sd), alpha = 0.2) + 
      ggplot2::geom_vline(data = line_df, mapping = ggplot2::aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5) + 
      ggplot2::facet_wrap( ~ across, nrow = nrow, ncol = ncol) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(strip.background = ggplot2::element_blank()) + 
      ggplot2::scale_x_continuous(breaks = line_df$label_breaks, labels = line_df$chromosome_name) + ggplot2::labs(x = "Chromosomes", y = "CNV-Results")
    
  }
  
  confuns::give_feedback(msg = "Done.", verbose = verbose)
  
  base::return(final_plot)
  
}

