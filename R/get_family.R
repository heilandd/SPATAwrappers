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


#' @title  getSurroundedSpots
#' @author Dieter Henrik Heiland
#' @description getSurroundedSpots
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 
getSurroundedSpots.sc <- function(object){
  of_sample <- SPATA2::getSampleNames(object)
  coords <- SPATA2::getCoordsDf(object)
  bc_origin <- coords$barcodes
  bc_destination <- coords$barcodes
  
  # get grouped data.frame with all barcode combinations
  dfgr <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>%
    dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 10)) %>%
    dplyr::group_by(bc_origin)
  
  # filter barcodes that are surrounded by a total of six spots
  
  sufficient_bc <-
    dplyr::slice_min(.data = dfgr, order_by = distance, n = 7) %>%
    dplyr::group_by(bc_origin, distance) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(distance != 0 & count == 1) %>%
    dplyr::pull(bc_origin)
  
  # filter for barcodes 
  
  final_df <-
    dplyr::filter(dfgr, bc_origin %in% sufficient_bc) %>%
    dplyr::slice_min(order_by = distance, n = 7) %>%
    dplyr::mutate(sample = {{of_sample}}) %>%
    dplyr::select(sample, dplyr::everything()) 
  
  base::return(final_df)
  
}


#' @title  getMoransI
#' @author Dieter Henrik Heiland
#' @description getMoransI
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 

getMoransI <- function(object, features, smooth=F, normalize=F){
  
  color_var <- SPATA2::hlpr_join_with_aes(object, 
                                          df = SPATA2::getCoordsDf(object), 
                                          color_by = features, 
                                          normalize = normalize, 
                                          smooth = smooth)
  coords <- color_var[, c("x", "y")] %>% as.data.frame()
  data <- color_var[, features] %>% as.data.frame()
  rownames(coords) = rownames(data) <- color_var$barcodes
  
  adj <- as.matrix(dist(coords))
  adj <- 1/adj
  diag(adj)=0
  data <- map_dfr(.x=features, function(feat){
    moran <- Rfast2::moranI(data[,feat], adj)
    moran <- data.frame(MoranI=moran[1], p_value=moran[2])
    rownames(moran) <- feat
    return(moran)
  })
  
  return(data)
}



#' @title  getAdjacentSpots
#' @author Dieter Henrik Heiland
#' @description getAdjacentSpots
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#'
getAdjacentSpots <- function(object, export.object=T, workers=16, ram=40){
  
  color_var <- SPATA2::getCoordsDf(object)
  coords <- color_var[, c("x", "y")] %>% as.data.frame()
  rownames(coords) <- color_var$barcodes
  
  adj <- as.matrix(dist(coords)) %>% reshape2::melt() 

  
  #calculate the mean min dist
  min_adj <- adj %>% dplyr::filter(value!=0) %>% dplyr::pull(value) %>% min()
  min_adj <- min_adj+(min_adj/2)
  
  
  base::options(future.fork.enable = TRUE)
  future::plan("multisession", workers = workers)
  future::supportsMulticore()
  base::options(future.globals.maxSize = ram* 10* 1024^2)
  message("... Run multicore ... ")
  
  
  message("---- Get neighbors ---- ")
  map.adj <- furrr::future_map(.x=SPATA2::getBarcodes(object), .f=function(bc){
    
    adj %>% 
      filter(Var1==bc) %>% 
      filter(value<min_adj & value>0) %>% 
      pull(Var2) %>% 
      as.character()
    
  }, .progress = T)
  
  names(map.adj) <- SPATA2::getBarcodes(object)
  
  
  if(export.object==T){
    object@spatial[[SPATA2::getSampleNames(object)]]$AdjacentSpots <- map.adj
    
  }else{
    object <- map.adj
  }
  
  return(object)
}



#' @title  Spatial Adjacent Regression
#' @author Dieter Henrik Heiland
#' @description Spatial Adjacent Regression
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#'

getSpatialAdjacentRegression <- function(object,
                                         feature1,
                                         feature2,
                                         smooth=F,
                                         normalize=F,
                                         workers=16, 
                                         ram=40){
  
  #Get data
  
  features <- c(feature1, feature2)
  
  color_var <- SPATA2::hlpr_join_with_aes(object, 
                                          df = SPATA2::getCoordsDf(object), 
                                          color_by = features, 
                                          normalize = normalize, 
                                          smooth = smooth)
  coords <- color_var[, c("x", "y")] %>% as.data.frame()
  data <- color_var[, features] %>% as.data.frame()
  rownames(coords) = rownames(data) <- color_var$barcodes
  
  if(is.null(object@spatial[[SPATA2::getSampleNames(object)]]$AdjacentSpots)) stop(
    "Adjacent Spots are not defined, please run getAdjacentSpots() before"
    )

  map.adj <- object@spatial[[SPATA2::getSampleNames(object)]]$AdjacentSpots
  
  base::options(future.fork.enable = TRUE)
  future::plan("multisession", workers = workers)
  future::supportsMulticore()
  base::options(future.globals.maxSize = ram* 10* 1024^2)
  message("... Run multicore ... ")
  
  
  message("---- Get Neighbor Expression ---- ")
  adj_exp <- furrr::future_map_dfr(.x=color_var$barcodes, .f=function(bc){
    
      color_var %>% 
      dplyr::filter(barcodes %in% as.character(map.adj[[bc]])) %>% 
      dplyr::select({{features}}) %>% 
      colMeans() %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>% 
      mutate(barcodes=bc)
      
      
    
  }, .progress = T)
  rownames(adj_exp)=adj_exp$barcodes
  adj_exp$barcodes=NULL
  names(adj_exp) <- paste0("adj_",names(adj_exp))
  
  data <- cbind(data[rownames(adj_exp), ], adj_exp)
  
  #get adjacent mean expression
  
  model_1 <- lm(data[,feature1]~data[,feature2]) %>% summary()
  model_2 <- lm(data[,feature2]~data[,feature1]) %>% summary()


  res_1 <- model_1$residuals %>% as.numeric()
  res_2 <- model_2$residuals %>% as.numeric()
  adj.r.squared <- model_1$adj.r.squared
  cor.out <- cor(data[,feature1] , data[,feature2])
  
  
  model_res_1 <- lm(data[,paste0("adj_",feature2)]~res_1) %>% summary()
  model_res_2 <- lm(data[,paste0("adj_",feature1)]~res_2) %>% summary()
  
  p_res_1 <- model_res_1$coefficients[2,4]
  p_res_2 <- model_res_2$coefficients[2,4]
  
  moran <- getMoransI(object, features)
  
  out <- data.frame(
    feature1 =feature1,
    feature2 = feature2,
    PearsonCor=cor.out,
    p=model_1$coefficients[2,4],
    adj.r.squared=adj.r.squared,
    Adj_feature1=model_res_1$adj.r.squared,
    Adj_feature2=model_res_2$adj.r.squared,
    p_Adj_feature1=p_res_1,
    p_Adj_feature2=p_res_2,
    moranI_feature1=moran$MoranI[1],
    moranI_feature1=moran$MoranI[2]
    
  )
  
  return(out)
  
  
}














