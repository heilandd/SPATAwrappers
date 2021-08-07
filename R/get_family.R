#' @title  getDistance
#' @author Dieter Henrik Heiland
#' @description getDistance
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

getDistance <- function(x1, y1, x2, y2) {sqrt((x2 - x1)^2 + (y2 - y1)^2)}



#' @title  check_color_to
#' @author Dieter Henrik Heiland
#' @description check_color_to
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 
check_color_to <- function(object,
                           color_to,
                           color_by,
                           max_length = NULL){
  
  require(SPATA2)
  all_genes = getGenes(object.spata)
  all_gene_sets = getGeneSets(object.spata)
  all_features = getFeatureNames(object.spata)
  
  if(!base::is.null(max_length)){
    base::warning("max_length is deprecated. ")
  }
  
  confuns::is_vec(color_to, "character", "color_to")
  
  return_list <- list()
  
  if(base::any(color_to %in% all_genes)){
    
    return_list[["genes"]] <-
      confuns::check_vector(
        input = color_to,
        against = all_genes,
        verbose = TRUE,
        ref.input = "input for argumet 'color to'",
        ref.against = "all known genes"
      )
    
  } else if(base::any(color_to %in% all_features)){
    
    if(base::length(color_to) != 1){
      
      base::stop("Features have to be specified as a single character value.")
      
    }
    
    return_list[["features"]] <- all_features[all_features %in% color_to]
    
  } else if(base::any(color_to %in% all_gene_sets)){
    
    if(base::length(color_to) != 1){
      
      base::stop("Gene-sets have to be specified as a single character value.")
    }
    
    return_list["gene_sets"] <- all_gene_sets[all_gene_sets %in% color_to]
    
  } else {
    
    base::stop(glue::glue("Could not find '{color_to}' among all genes, gene-sets and features."))
    
  }
  
  base::return(return_list)
  
}


#' @title  getSurroundedSpots
#' @author Dieter Henrik Heiland
#' @description getSurroundedSpots
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 
getSurroundedSpots <- function(object){
  of_sample <- SPATA2::getSampleNames(object)
  coords <- SPATA2::getCoordsDf(object)
  bc_origin <- coords$barcodes
  bc_destination <- coords$barcodes
  
  # get grouped data.frame with all barcode combinations
  
  dfgr <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
    dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>%
    dplyr::group_by(bc_origin)
  
  # filter barcodes that are surrounded by a total of six spots
  
  sufficient_bc <-
    dplyr::slice_min(.data = dfgr, order_by = distance, n = 7) %>%
    dplyr::group_by(bc_origin, distance) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(distance != 0 & count == 6) %>%
    dplyr::pull(bc_origin)
  
  # filter for barcodes 
  
  final_df <-
    dplyr::filter(dfgr, bc_origin %in% sufficient_bc) %>%
    dplyr::slice_min(order_by = distance, n = 7) %>%
    dplyr::mutate(sample = {{of_sample}}) %>%
    dplyr::select(sample, dplyr::everything()) 
  
  base::return(final_df)
  
}