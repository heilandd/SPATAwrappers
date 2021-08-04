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