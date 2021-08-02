#' @title  Juxtaposition in stRNA-seq
#' @author Dieter Henrik Heiland
#' @description Juxtaposition in stRNA-seq
#' @param source Variable define basline of state from which the distance is computet
#' @param target Variable "to" for distance computation
#' @param source.type and target.type can be "feature", "gene" or "geneset"
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

inferJuxtaposition <- function(object, 
                               source,
                               source.type="feature",
                               target,
                               target.type="feature"){
  
  ###From
  
  if(source.type=="feature"){from <- 
    object %>% 
    SPATA2::joinWithFeatures(features = source) %>% 
    dplyr::mutate(source=SPATA2::hlpr_normalize_vctr(!!sym(source)) ) %>% 
    dplyr::arrange(dplyr::desc(source))}
  
  if(source.type=="gene"){from <- 
    object %>% 
    SPATA2::joinWithGenes(genes = source, average_genes = T) %>% 
    dplyr::mutate(source=SPATA2::hlpr_normalize_vctr(mean_genes) ) %>% 
    dplyr::arrange(dplyr::desc(source))}
  
  if(source.type=="geneset"){from <- 
    object %>% 
    SPATA2::joinWithGeneSets(gene_sets = source) %>% 
    dplyr::mutate(source=SPATA2::hlpr_normalize_vctr(!!sym(source)) ) %>% 
    dplyr::arrange(dplyr::desc(source))}
  
  #### to
  
  
  if(target.type=="feature"){to <- 
    object %>% 
    SPATA2::joinWithFeatures(features = target) %>% 
    dplyr::mutate(target=SPATA2::hlpr_normalize_vctr(!!sym(target)) ) %>% 
    dplyr::arrange(dplyr::desc(target))}
  
  if(target.type=="gene"){to <- 
    object %>% 
    SPATA2::joinWithGenes(genes = target, average_genes = T) %>% 
    dplyr::mutate(target=SPATA2::hlpr_normalize_vctr(mean_genes) ) %>% 
    dplyr::arrange(dplyr::desc(target))}
  
  if(target.type=="geneset"){to <- 
    object %>% 
    SPATA2::joinWithGeneSets(gene_sets = target) %>% 
    dplyr::mutate(target=SPATA2::hlpr_normalize_vctr(!!sym(target)) ) %>% 
    dplyr::arrange(dplyr::desc(target))}
  
  
  to <- 
    to %>% 
    dplyr::select(barcodes, target, x, y) %>% 
    dplyr::rename("x.to":=x) %>% 
    dplyr::rename("y.to":=y) %>% 
    dplyr::mutate(aligned.target=target)
  
  from <- 
    from %>% 
    dplyr::select(barcodes, source, x, y) %>% 
    dplyr::rename("x.from":=x) %>% 
    dplyr::rename("y.from":=y)

  dist <- 
    cbind(from, to %>% dplyr::select(-aligned.target)) %>% 
    dplyr::select(-5) %>% 
    dplyr::mutate(dist=SPATAwrappers::getDistance(.$x.from, .$y.from,.$x.to, .$y.to )) %>% 
    dplyr::left_join(., to %>% dplyr::select(barcodes, aligned.target), by="barcodes")
  
  
  
  return(dist)
  
  
  
}

#' @title  inferSpatial.mc
#' @author Dieter Henrik Heiland
#' @description inferSpatial.mc
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
inferSpatial.mc <- function(P1,P2, n = 599){
  # Monte-Carlo Simulation of spatial correlation
  if(length(P2)!=length(P1)) stop("Unequal Inputs")
  message(paste0(Sys.time(), "Start Model"))
  M <- glm(P1 ~ P2)
  #coef(M)[2]
  message(paste0(Sys.time(), "Start MC Simulation"))
  #MC
  I.r <- purrr::map_dbl(.x=1:n, .f=function(i){
    x <- sample(P1, replace=FALSE)
    y <- P2
    # Compute new set of lagged values
    #x.lag <- lag.listw(lw, x)
    # Compute the regression slope and store its value
    M.r    <- glm(y ~ x)
    I.r <- coef(M.r)[2]
    return(I.r)
  })
  
  #hist(I.r, main=NULL, xlab="Spatial-Cor-MC", las=1, xlim=c(-0.5,0.5))
  #abline(v=coef(M)[2], col="red")
  
  #p-value
  message(paste0(Sys.time(), "P-Value Extraction"))
  N.greater <- sum(coef(M)[2] > I.r)
  p <- min(N.greater + 1, n + 1 - N.greater) / (n + 1)
  p
  out <- list(coef(M)[2], p)
  names(out) <- c("Cor", "p")
  return(out)
}

#' @title  inferSpatial.mc
#' @author Dieter Henrik Heiland
#' @description inferSpatial.mc
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
inferSpatial.ac <- function(object,feature){
  
  message(paste0(Sys.time(), "  Start Autocorrelation"))
  
  tissue.pos <- SPATA2::getCoordsDf(object)
  plot.new()
  message(paste0(Sys.time(), "  Transforme in spatial dataset"))
  segments <- pbmcapply::pbmclapply(1:nrow(tissue.pos), function(xx){
    
    segment_c <-plotrix::draw.circle(x=tissue.pos$x[xx], 
                                     y=tissue.pos$y[xx],
                                     radius=5*2) 
    
    segment <- as.data.frame(do.call(cbind, segment_c))
    
    segment <- sp::Polygon(cbind(segment$x, segment$y))
    segment <- sp::Polygons(list(segment), tissue.pos$barcodes[xx])
    return(segment)
  }, mc.cores = 8)
  
  SpP = sp::SpatialPolygons(segments, 1:length(segments))
  
  attr <- SPATA2::getFeatureDf(object) %>% as.data.frame()
  rownames(attr) <- attr$barcodes
  attr <- attr[,feature]
  
  SrDf = sp::SpatialPolygonsDataFrame(SpP, attr)
  
  message(paste0(Sys.time(), "  Define neighboring Spots"))
  
  nb <- spdep::poly2nb(SrDf, queen=F, snap=25)
  lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
  
  message(paste0(Sys.time(), "  Computing the Moran’s I statistic"))
  
  stat <- purrr::map(.x=feature, .f=function(x){
    print(x)
    model <- spdep::moran.test(SrDf@data[,x],lw, zero.policy=T)
    model$estimate[[1]]
  })
  
  names(stat) <- feature
  
  return(stat)
  
  
}

#' @title  inferSpatial.mc
#' @author Dieter Henrik Heiland
#' @description inferSpatial.mc
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
inferSpatial.plot <- function(object,feature, plot=T){
  
  
  #Run analysis 
  message(paste0(Sys.time(), "  Run MC Simulation for Spatial  Correlations analysis"))
  
  cor.mat <- matrix(NA, length(feature),length(feature));colnames(cor.mat) = rownames(cor.mat) <- feature
  cor.mat <- reshape2::melt(cor.mat)
  
  # fill data 
  cor.mat$value <- pbmcapply::pbmclapply(1:nrow(cor.mat), function(i){
    cor.out <- inferSpatial.mc(P1=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,1]),
                               P2=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,2]),
                               n=200)
    return(cor.out$Cor)
    
  }, mc.cores = 8) %>% unlist()
  
  message(paste0(Sys.time(), "  Run Auto Corelation analysis"))
  
  cor.Auto <- SPATAwrappers::inferSpatial.ac(object, feature)
  cor.Auto <- as.data.frame(do.call(rbind,cor.Auto)) %>% 
    tibble::rownames_to_column("Var2") %>% 
    mutate(Var1="Autocorelation",
           value=V1) %>%
    select(Var1,Var2,value)
  
  cor.mat <- rbind(cor.mat, cor.Auto)
  
  
  message(paste0(Sys.time(), " Plotting "))
  
  if(plot==T){
    corrplot::corrplot.mixed(t(reshape2::acast(cor.mat, Var1~Var2, value.var="value")))
  }
  
  
  return(cor.mat)
  
}

#' @title  jaccard
#' @author Dieter Henrik Heiland
#' @description jaccard
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#' @title  inferSpotlight
#' @author Dieter Henrik Heiland
#' @description inferSpotlight
#' @inherit 
#' @return 
#' @examples 
#' @export
inferSpotlight <- function(spata.obj, seurat.obj, feature, marker.genes=NULL){
  
  #Set up Seurat features
  seurat.obj <- Seurat::SetIdent(seurat.obj, value=feature)
  if(is.null(marker.genes)){
    marker.spotlight <- Seurat::FindAllMarkers(seurat.obj)
  }else{marker.spotlight <- marker.genes}
  
  
  #Run Spotlight
  spot <- 
    SPOTlight::spotlight_deconvolution(se_sc= seurat.obj,
                                       counts_spatial=SPATA2::getCountMatrix(spata.obj),
                                       clust_vr=feature,
                                       cluster_markers=marker.spotlight)
  
  
  
  # Transfer back
  spata.obj <- 
    SPATA2::addFeatures(spata.obj, 
                        feature_df = spot[[2]] %>% as.data.frame() %>% dplyr::mutate(barcodes=SPATA2::getBarcodes(spata.obj)) %>% dplyr::select(barcodes, 1,2),,
                        overwrite = T
                        )
  
  
  return(spata.obj)
  
  
}










