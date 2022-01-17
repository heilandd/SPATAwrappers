
#' @title  BaySpace Cluster
#' @author Dieter Henrik Heiland
#' @description BaySpace Cluster
#' @inherit 
#' @return 
#' @examples 
#' @export

runBaySpace <- function(object, pathToOuts, Spatial.enhancer=T, max.cluster=13, return.model=T, empty.remove=F){
  

  sce <- BayesSpace::readVisium(pathToOuts)
  sample <- SPATA2::getSampleNames(object)
  object@data[[sample]]$SCE <- list(colData=SingleCellExperiment::colData(sce),
                                    rowData=SingleCellExperiment::rowData(sce))

  set.seed(102)

  if(empty.remove==T){
    sce <- sce[, colSums(counts(sce)) > 0]
  }else{
    spots.no.read= colData(sce[, colSums(counts(sce)) == 0]) %>% rownames()
    #add 2 random reads
    
    sce@assays@data$counts[runif(n=length(spots.no.read), 
                               min = 1, 
                               max=nrow(sce@assays@data$counts)) %>% round(), spots.no.read] <-  1
  }
  
  

  
  
  
  space <- BayesSpace::spatialPreprocess(sce , platform="Visium", n.PCs=30, n.HVGs=2000, log.normalize=T)
  space <- BayesSpace::qTune(space, qs=seq(2, max.cluster), platform="Visium")
  #calculate the ellbow point

  logliks <- attr(space, "q.logliks")
  optimal.cluster <- akmedoids::elbow_point(x=logliks$q, y=logliks$loglik)$x %>% round()
  space <- BayesSpace::spatialCluster(space, 
                                      q=optimal.cluster, 
                                      platform="Visium", 
                                      d=7,
                                      init.method="mclust", 
                                      model="t", 
                                      gamma=2,
                                      nrep=1000, 
                                      burn.in=100,
                                      save.chain=TRUE)
  
  cluster.df <- 
  space@colData %>% 
  as.data.frame() %>% 
  dplyr::select(spot, spatial.cluster) %>% 
  dplyr::rename("barcodes":=spot) %>% 
  dplyr::rename("BayesSpace":=spatial.cluster)
  object <- SPATA2::addFeatures(object, cluster.df)
  if(Spatial.enhancer==T){

    space.enhanced <- BayesSpace::spatialEnhance(space, q=optimal.cluster, platform="Visium", d=7,
                                                 model="t", gamma=2,
                                                 jitter_prior=0.3, jitter_scale=3.5,
                                                 nrep=1000, burn.in=100,
                                                 save.chain=TRUE)

    space.data <- as.data.frame(space.enhanced@colData)
    scCoords <- SPATAwrappers::getNucleusPosition(object)
    x <- space.data %>% pull(col)
    y <- space.data %>% pull(row)
    z <- space.data %>% pull(spatial.cluster)
    s1 =  akima::interp(x = x, 
                        y = y, 
                        z = z, 
                        nx = nrow(space.data), 
                        ny = nrow(space.data), 
                        xo = seq(min(x), max(x), length = nrow(space.data)), 
                        yo = seq(min(y), max(y), length = nrow(space.data)))
    message(paste0(Sys.time(), " ---- ", "Predict single-cell expression levels ", " ----"))
    r.pred <- raster::raster(t(s1$z[,ncol(s1$z):1]), 
                             xmn = min(scCoords$x), xmx = max(scCoords$x),
                             ymn = min(scCoords$y), ymx = max(scCoords$y))
    pts <- sp::SpatialPointsDataFrame(scCoords[,c('x','y')], scCoords)
    scCoords$pred <- raster::extract(r.pred, pts)
    scCoords <- scCoords %>% filter(!is.na(pred))
    scCoords$pred <- round(scCoords$pred, digits = 0) %>% as.character()
    object@spatial[[sample]]$BayesSpace <- scCoords
    if(return.model==T){
      object@spatial[[sample]]$SpaceEnhancer <-space.enhanced
      object@spatial[[sample]]$Space <-space

    }
  }
  return(object)
}



#' @title  Enhanced Nucleus Feature
#' @author Dieter Henrik Heiland
#' @description Enhanced Nucleus Feature
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotEnhancedNucleusFeature <- function(object, 
                                       genes,
                                       pt.size=0.5, 
                                       pt.alpha=1,
                                       alpha2pred=F,
                                       addImage=F,
                                       Palette=NULL,
                                       pt_clrsp="Reds"){

  message(paste0(Sys.time(), " ---- ", "Get BayerSpace Enhanced Features ", " ----"))
  

  #Get Space enhanced Feature
  

  sample <- SPATA2::getSampleNames(object)
  

  space.enhanced <- object@spatial[[sample]]$SpaceEnhancer
  

  space <- object@spatial[[sample]]$Space
  

  space.enhanced <- BayesSpace::enhanceFeatures(space.enhanced, 

                                                space,

                                                feature_names=genes,

                                                nrounds=0)

  space.data <- as.data.frame(space.enhanced@colData)
  
  space.data$z <- as.numeric(space.enhanced@assays@data$logcounts[genes, ])

  scCoords <- SPATAwrappers::getNucleusPosition(object)

  x <- space.data %>% pull(col)
  

  y <- space.data %>% pull(row)
  

  z <- space.data %>% pull(z)

  s1 =  akima::interp(x = x, 

                      y = y, 

                      z = z, 

                      nx = nrow(space.data), 

                      ny = nrow(space.data), 

                      xo = seq(min(x), max(x), length = nrow(space.data)), 

                      yo = seq(min(y), max(y), length = nrow(space.data)))

  message(paste0(Sys.time(), " ---- ", "Predict single-cell expression levels ", " ----"))

  r.pred <- raster::raster(t(s1$z[,ncol(s1$z):1]), 
                           xmn = min(scCoords$x), xmx = max(scCoords$x),
                           ymn = min(scCoords$y), ymx = max(scCoords$y))

  pts <- sp::SpatialPointsDataFrame(scCoords[,c('x','y')], scCoords)
  
  scCoords$pred <- raster::extract(r.pred, pts)
  scCoords <- scCoords %>% filter(!is.na(pred))

  if(addImage==T){p=SPATA2::plotSurface(object, display_image=T, pt_alpha = 0)}else{p=ggplot()+theme_void()}
  
  
  if(alpha2pred==T){pt.alpha <- scCoords$pred}

  p=p+geom_point(data=scCoords, aes(x=x, y=y, color=pred), size=pt.size, alpha=pt.alpha)

  if(is.null(Palette)){p=p+SPATA2::scale_color_add_on(aes = "color", 
                                                      clrsp = pt_clrsp)}else{p=p+scale_colour_gradientn(colours = Palette(50), oob=scales::squish,...)}
  

  p=p+ggplot2::coord_equal()
  
  return(p)
}


#' @title  plotEnhancedCluster
#' @author Dieter Henrik Heiland
#' @description plotEnhancedCluster
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
plotEnhancedCluster <- function(object, 
                                pt.size=0.5, 
                                pt.alpha=1,
                                addImage=F){
  
  sample <- SPATA2::getSampleNames(object)
  scCoords <- object@spatial[[sample]]$BayesSpace
  if(addImage==T){p=SPATA2::plotSurface(object, display_image=T, pt_alpha = 0)}else{p=ggplot()+theme_void()}
  p=p+geom_point(data=scCoords, 
                 aes(x=x, y=y, color=pred), 
                 size=pt.size, 
                 alpha=pt.alpha)

  p=p+ggplot2::coord_equal()

  return(p)

}

