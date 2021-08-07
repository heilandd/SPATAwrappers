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
                                   color="black",
                                   pt_alpha=0.3,
                                   lt.alpha=0.1,
                                   pt_size=6,
                                   lt.size=0.2,
                                   thr=0.8,
                                   dist=200,
                                   data=F,
                                   add=F){
  
  
  
  
  
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
  
  if(add==T){p <- geom_segment(data=df_distance.processed, mapping = aes(x=xo, y=yo, xend=xd, yend=yd,size=size), alpha=lt.alpha, color=color)}else{
  p <- SPATA2::plotSurface(object, color_by = color_by, pt_size = pt_size, verbose = F, pt_alpha =pt_alpha )
  p <- p+geom_segment(data=df_distance.processed, mapping = aes(x=xo, y=yo, xend=xd, yend=yd,size=size), alpha=lt.alpha,color=color)
  p <- p+Seurat::NoLegend()
  }
  
  if(data==F){out <- p}else{
  out=list(df_distance.processed,p)
  names(out)=c("data", "plot")}
  
  return(out)
  
}



#' @title  plotStreamlines
#' @author Dieter Henrik Heiland
#' @description plotStreamlines
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
plotStreamlines <- function(VF,
                            parameter,
                            size.arrow=1,
                            grid=c(50,50),
                            gamma.u=0.2,
                            gamma.v=0.2,
                            arrow.angle=25,
                            arrow.length=0.6,
                            L=100,
                            color.extern=NULL,
                            res=1,
                            pt.size=6,
                            pt.alpha=0.8,
                            skip=2.5,
                            xwrap=NULL,
                            ywrap=NULL){
  VF <- 
    VF %>% 
    dplyr::select(x,y,{{parameter}}, t.x, t.y) %>% 
    rename("parameter":=!!sym(parameter))
  
  drifter.split.sf <-  
    VF %>% 
    sf::st_as_sf(coords = c("x", "y"))
  drifter.grid <-  drifter.split.sf %>% 
    sf::st_make_grid(n = grid)%>%
    sf::st_sf()
  drifter.split.sf.se <- 
    drifter.split.sf %>% 
    dplyr::filter(parameter!=0)
  drifter.gridded <-  
    drifter.grid %>% 
    mutate(id = 1:n(), contained = lapply(sf::st_contains(sf::st_sf(geometry),drifter.split.sf.se),identity),
           obs = sapply(contained, length),
           u = sapply(contained, function(x) {median(drifter.split.sf.se[x,]$t.x, na.rm = TRUE)}),
           v = sapply(contained, function(x) {median(drifter.split.sf.se[x,]$t.y, na.rm = TRUE)})) 
  drifter.gridded = drifter.gridded %>% dplyr::select(obs, u, v) %>% na.omit()
  
  coordinates <-  
    drifter.gridded %>% 
    sf::st_centroid() %>% 
    sf::st_coordinates() %>% 
    tibble::as_tibble() %>% 
    dplyr::rename(x = X, y = Y)
  sf::st_geometry(drifter.gridded) = NULL
  
  current.gridded.se <-  
    coordinates %>% 
    dplyr::bind_cols(drifter.gridded) %>% 
    dplyr::mutate(season = "SE")
  
  drifter.current.gridded <-  
    current.gridded.se %>% 
    dplyr::bind_rows(current.gridded.se)
  
  drf.se <-  
    drifter.current.gridded %>% 
    dplyr::filter(season == "SE")
  
  u.se <-  oce::interpBarnes(x = drf.se$x, y = drf.se$y, z = drf.se$u, gamma=gamma.u)
  
  dimension = data.frame(lon = u.se$xg, u.se$zg) %>% dim()
  
  u.tb <-  
    data.frame(lon = u.se$xg, u.se$zg) %>% 
    tidyr::gather(key = "lata", value = "u", 2:dimension[2]) %>% 
    dplyr::mutate(lat = rep(u.se$yg, each = dimension[1])) %>% 
    dplyr::select(lon,lat, u) %>% 
    tibble::as_tibble()
  
  v.se = oce::interpBarnes(x = drf.se$x, y = drf.se$y, z = drf.se$v,gamma=gamma.v)
  
  v.tb = data.frame(lon = v.se$xg, v.se$zg) %>% 
    tidyr::gather(key = "lata", value = "v", 2:dimension[2]) %>% 
    dplyr::mutate(lat = rep(v.se$yg, each = dimension[1])) %>% 
    dplyr::select(lon,lat, v) %>% 
    tibble::as_tibble()
  
  uv.se <-  
    u.tb %>% 
    dplyr::bind_cols(v.tb %>% dplyr::select(v)) %>% 
    dplyr::mutate(vel = sqrt(u^2+v^2))
  
  
  color.points <- VF$parameter
  if(!is.null(color.extern)){color.points <- color.extern}
  
  if(color.points %>% class()=="factor"){
    p= ggplot2::ggplot()+ggplot2::theme_void()
    p=p+geom_point(data=VF, mapping=aes(x,y, color=color.points), size=pt.size, alpha=pt.alpha)
    p=p+metR::geom_streamline(data = uv.se, aes(x = lon, y = lat, dx = u, dy = v),
                              size=size.arrow,
                              arrow.length = arrow.length,
                              arrow.angle = arrow.angle,
                              arrow.type = "closed",
                              L = L, res =res,skip=skip,xwrap =xwrap,ywrap = ywrap,
                              lineend = "round")
    
  }else{
    
    p= ggplot2::ggplot()+ggplot2::theme_void()
    p=p+geom_point(data=VF, mapping=aes(x,y, color=color.points), size=pt.size, alpha=pt.alpha)
    p=p+ggplot2::scale_color_viridis_c(guide = "none")
    p=p+metR::geom_streamline(data = uv.se, aes(x = lon, y = lat, dx = u, dy = v),
                              size=size.arrow,
                              arrow.length = arrow.length,
                              arrow.angle = arrow.angle,
                              arrow.type = "closed",
                              L = L, res =res,skip=skip,xwrap =xwrap,ywrap = ywrap,
                              lineend = "round")
  }
  
  return(p)
  
  
}

#' @title  plotVectorFields
#' @author Dieter Henrik Heiland
#' @description plotVectorFields
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
plotVectorFields <- function(VF, parameter, pt.size=6,pt.alpha=0.8,color.extern=NULL,skip=1){
  
  VF <- 
    VF %>% 
    dplyr::select(x,y,{{parameter}}, t.x, t.y) %>% 
    rename("parameter":=!!sym(parameter))
  
  
  color.points <- VF$parameter
  if(!is.null(color.extern)){color.points <- color.extern}
  
  if(color.points %>% class()=="factor"){
  p <- 
    ggplot2::ggplot(data=VF, aes(x,y))+
    ggplot2::geom_point(data=VF , mapping=aes(x,y, color=color.points), size=pt.size, alpha=pt.alpha)+
    metR::geom_vector(aes(dx = t.x, dy = t.y),skip=skip) +
    metR::scale_mag()+
    ggplot2::theme_void()+
    Seurat::NoLegend()
  }else{
    p <- 
      ggplot2::ggplot(data=VF, aes(x,y))+
      ggplot2::geom_point(data=VF , mapping=aes(x,y, color=color.points), size=pt.size, alpha=pt.alpha)+
      ggplot2::scale_color_viridis_c(guide = "none")+
      metR::geom_vector(aes(dx = t.x, dy = t.y),skip=skip) +
      ggplot2::theme_void()+
      Seurat::NoLegend()+
      metR::scale_mag()
  }
  
  return(p)
}








