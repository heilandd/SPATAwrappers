#' @title  Juxtaposition in stRNA-seq
#' @author Dieter Henrik Heiland
#' @description Juxtaposition in stRNA-seq
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

inferJuxtaposition <- function(object, genset=NULL, feature=NULL){
  
  
  #Check input
  if(is.null(genset) & is.null(feature)) stop ("No feature selected")
  if(is.null(genset) & length(feature)!=2) stop ("Not enough feature selected: First source second target")
  if(is.null(feature) & length(genset)!=2) stop ("Not enough feature selected: First source second target")
  
  #Select Modus
  if(length(genset)==1 & length(feature)==1){
    
    #from GS to Feature
    
  }
  if(length(genset)==2 ){
    
    #from GS to GS
    
  }
  if(length(feature)==2){
    
    #from Feature to Feature
    
  }
  
  
  
  
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
                                     radius=5*0.4) 
    
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
    model <- moran.test(SrDf@data[,x],lw)
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
    cor.out <- spatial.mc(P1=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,1]),
                          P2=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,2]),
                          n=200)
    return(cor.out$Cor)
    
  }, mc.cores = 8) %>% unlist()
  
  message(paste0(Sys.time(), "  Run Auto Corelation analysis"))
  
  cor.Auto <- spatial.ac(object, feature)
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

#' @title  jaccard
#' @author Dieter Henrik Heiland
#' @description jaccard
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
inferSpotlight <- function(spata.obj, seurat.obj, feature){
  
  #Set up Seurat features
  seurat.obj <- Seurat::SetIdent(caf.object.spotlight, value=feature)
  marker.spotlight <- Seurat::FindAllMarkers(caf.object.spotlight)
  
  #Run Spotlight
  spot <- 
    SPOTlight::spotlight_deconvolution(se_sc= seurat.obj,
                                       counts_spatial=SPATA2::getCountMatrix(spata.obj),
                                       clust_vr=feature,
                                       cluster_markers=marker.spotlight)
  
  
  # Transfer back
  spata.obj <- 
    SPATA2::addFeatures(spata.obj, 
                        feature_df = spot[[2]] %>% as.data.frame() %>% dplyr::mutate(barcodes=SPATA2::getBarcodes(spata.obj)) %>% dplyr::select(barcodes, 1,2)
                        )
  
  
  return(spata.obj)
  
  
}










