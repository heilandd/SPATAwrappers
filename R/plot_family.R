#' @title  plotJuxtaposition
#' @author Dieter Henrik Heiland
#' @description plotJuxtaposition
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export


plotJuxtaposition <- function(dist, points=T, span=0.2,rescale=F, quantile.cor=T, q = 0.1){
  
  if(quantile.cor==T){
    dist <- NFCN2::getCleaned(dist, feat = "aligned.target",q = q)
  }
  
  if(rescale==T){
    dist$aligned.target <- scales::rescale(dist$aligned.target, c(0,1))
  }
  
  p=ggplot2::ggplot()+ ggplot2::theme_classic()
  if(points==T){
    p=p+ggplot2::geom_point(data=dist, mapping=aes(x=dist, y=aligned.target))
  }
  p=p+geom_smooth(data=dist, mapping=aes(x=dist, y=aligned.target), se=F, span=span)
  return(p)
  
}


#' @title  plotJuxtapositionSPATA
#' @author Dieter Henrik Heiland
#' @description plotJuxtapositionSPATA
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotJuxtapositionSPATA <- function(object,
                                   feature.source="NPC",
                                   feature.target="OPC",
                                   color_by="seurat_clusters",
                                   pt_alpha=0.3,
                                   lt.alpha=0.1,
                                   pt_size=6,
                                   lt.size=0.2,
                                   thr=0.8,
                                   dist=200){
  
  
  
  
  
  get.source <- SPATAwrappers::check_color_to(object,color_to = feature.source)
  get.target <- SPATAwrappers::check_color_to(object,color_to = feature.target)
  source <- object %>% SPATA2::joinWithVariables(variables = get.source, verbose = F) %>% dplyr::rename("bc_origin":=barcodes)
  target <- object %>% SPATA2::joinWithVariables(variables =  get.target ,verbose = F)%>% dplyr::rename("bc_destination":=barcodes)
  
  
  coords <- SPATA2::getCoordsDf(object)
  bc_origin <- coords$barcodes
  bc_destination <- coords$barcodes
  
  df_distance <-
    tidyr::expand_grid(bc_origin, bc_destination) %>%
    dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), key = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), key = "bc_destination") %>%
    dplyr::left_join(x = ., y = dplyr::select(source, bc_origin, xo = x, yo = y, {{feature.source}}), key = "bc_origin") %>%
    dplyr::left_join(x = ., y = dplyr::select(target, bc_destination, xd = x, yd = y, {{feature.target}}), key = "bc_destination") %>%
    dplyr::mutate(distance = base::round(sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0))
  
  
  df_distance.processed <- 
    df_distance %>% 
    dplyr::mutate(source=SPATA2::hlpr_normalize_vctr(df_distance %>% dplyr::pull(!!sym(feature.source))), 
                  target=SPATA2::hlpr_normalize_vctr(df_distance %>% dplyr::pull(!!sym(feature.target)))) %>%  
    dplyr::filter(source>thr) %>% 
    dplyr::filter(target>thr) %>% 
    dplyr::filter(distance<dist) %>% 
    dplyr::mutate(size=SPATA2::hlpr_normalize_vctr(max(distance)-(distance)))
  
  
  p <- SPATA2::plotSurface(object, color_by = color_by, pt_size = pt_size, verbose = F, pt_alpha =pt_alpha )
  p <- p+geom_segment(data=df_distance.processed, mapping = aes(x=xo, y=yo, xend=xd, yend=yd,size=size), alpha=lt.alpha)
  p <- p+Seurat::NoLegend()
  
  return(p)
  
}
