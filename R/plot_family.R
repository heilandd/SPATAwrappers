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
                            surface=T,
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
  
  color.points <- parameter
  if(!is.null(color.extern)){color.points <- color.extern }
  
  VF <- 
    VF %>% 
    dplyr::select(x,y,{{parameter}},{{color.points}}, t.x, t.y) %>% 
    dplyr::rename("parameter":=!!sym(parameter))
  
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
  
  
  
  if(color.points %>% class()=="factor"){
    p= ggplot2::ggplot()+ggplot2::theme_void()
    if(surface==T){
    p=p+geom_point(data=VF, mapping=aes(x,y, color=.data[[color.points]]), size=pt.size, alpha=pt.alpha)
    }
    p=p+metR::geom_streamline(data = uv.se, aes(x = lon, y = lat, dx = u, dy = v),
                              size=size.arrow,
                              arrow.length = arrow.length,
                              arrow.angle = arrow.angle,
                              arrow.type = "closed",
                              L = L, res =res,skip=skip,xwrap =xwrap,ywrap = ywrap,
                              lineend = "round")
    
  }else{
    
    p= ggplot2::ggplot()+ggplot2::theme_void()
    if(surface==T){
      if(color.points==parameter){ p=p+geom_point(data=VF, mapping=aes(x,y, color=parameter), size=pt.size, alpha=pt.alpha) }else{
        p=p+geom_point(data=VF, mapping=aes(x,y, color=.data[[color.points]]), size=pt.size, alpha=pt.alpha)
      }
      
      p=p+ggplot2::scale_color_viridis_c(guide = "none")
    }
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
    dplyr::rename("parameter":=!!sym(parameter))
  
  
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


#' @title  plotHeatmap
#' @author Dieter Henrik Heiland
#' @description plotHeatmap
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plotHeatmap <-  function(object, 
                         across=NULL,
                         gs, 
                         thr=c(0.8,0.7,0.7,0.6,0.5),
                         nr.sample=NULL,
                         plot.type=c("heatmap", "bar.comp", "point", "line", "stat"),
                         size=2,
                         sample.type="sample.type",
                         color="black",
                         alpha=1,
                         feature=NULL){
    
    
    sc.spata <- object
    spata2_obj <-   object
    
    # Prepare data ------------------------------------------------------------
    
    # Check if GS is in the GS list of spata.obj
    
    if (length(c(intersect(gs, SPATA2::getGeneSets(sc.spata)) == gs ) %>% unique()) ==2 ) stop("Input Gene Set do not match gene Sets in spata object")
    if (c(intersect(gs, SPATA2::getGeneSets(sc.spata)) == gs ) %>% unique() != T ) stop("Input Gene Set do not match gene Sets in spata object")
    if(length(gs)!=length(thr)) stop("thr and gs do not match. Different length")
    
    # Run analysis of the basline gs
    message(paste0(Sys.time(),  "  Run analysis and extract spots  "))
    df <- SPATA2::joinWithGeneSets(spata2_obj, gene_sets = gs, smooth=F, ignore = T) 
    
    #Extract Freatures
    if(!is.null(feature)){
      # If not Element of the feature DF
      if(!is.element(feature, c(SPATA2::getFeatureNames(spata2_obj) %>% as.character()))){
        
        # Test if feature is element of GS
        if(is.element(feature, SPATA2::getGeneSets(spata2_obj) %>% as.character())){
          gene_df=c(SPATA2::joinWithGeneSets(spata2_obj, gene_sets = feature) %>% dplyr::select(barcodes, {{feature}})) %>% as.data.frame() 
          spata2_obj <- SPATA2::addFeatures(spata2_obj,gene_df, overwrite = T)}
        
        if(is.element(feature, SPATA2::getGenes(spata2_obj) %>% as.character())){
          gene_df=c(SPATA2::joinWithGenes(spata2_obj, genes = feature) %>% dplyr::select(barcodes, {{feature}})) %>% as.data.frame() 
          spata2_obj <- SPATA2::addFeatures(spata2_obj,gene_df, overwrite = T)}
        
        element <- c(is.element(feature, SPATA2::getGeneSets(spata2_obj) %>% as.character()),is.element(feature, SPATA2::getGenes(spata2_obj) %>% as.character()))
        if(length(unique(element))==1) stop("feature not found")
        
      }}
    
    
    # Run analysis ------------------------------------------------------------
    
    if(plot.type=="heatmap"){
      
      if(!is.null(across)){
        split.f <- SPATA2::joinWithFeatures(spata2_obj, features = across) %>% dplyr::select(barcodes, !!sym(across))
        df <- df %>% left_join(., split.f, by="barcodes")
        
        message(paste0(df %>% pull(!!sym(across)) %>% unique()))
        
        p <- purrr::map(.x=unique(df %>% pull(!!sym(across))), .f=function(sep.i){
          
          df <- df %>% filter(!!sym(across)==sep.i)
          
          # get thresholded Spots 
          gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
            if(!is.null(nr.sample)){
              df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
            }else{
              df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
            }
          })
          
          df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
          rownames(df.1) <- df$barcodes
          
          type <- purrr::map(.x=1:length(gs), .f=function(i){
            rep(gs[i], length(gl.ls[[i]]))
          }) %>% unlist()
          
          bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type) #%>% filter(dup==F)
          bc <- bc.1$inp
          df.2 <- df.1[bc, gs] %>% as.matrix() %>% t()
          
          col= colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(200)
          myBreaks <- c(seq(0, 0.59, length.out=25), 
                        seq(0.6, 1, length.out=25))
          
          p=pheatmap::pheatmap(df.2 , 
                               cluster_rows = F,
                               cluster_cols = F, 
                               show_colnames = F,
                               color=col,
                               silent=T)
          
          return(p)
          
          
        })
        
        
        
      }
      else{
        
        # get thresholded Spots 
        gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
          if(!is.null(nr.sample)){
            df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
          }else{
            df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
          }
        })
        
        
        
        
        df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
        rownames(df.1) <- df$barcodes
        
        type <- purrr::map(.x=1:length(gs), .f=function(i){
          rep(gs[i], length(gl.ls[[i]]))
        }) %>% unlist()
        
        bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type) #%>% filter(dup==F)
        bc <- bc.1$inp
        df.2 <- df.1[bc, gs] %>% as.matrix() %>% t()
        
        col= colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(200)
        myBreaks <- c(seq(0, 0.59, length.out=25), 
                      seq(0.6, 1, length.out=25))
        
        p=pheatmap::pheatmap(df.2 , 
                             cluster_rows = F,
                             cluster_cols = F, 
                             show_colnames = F,
                             color=col,
                             silent=T)
        p=list(p, df.2)
        return(p) 
      }
      return(p) 
      
    }
    if(plot.type=="stat"){
      
      if(!is.null(across)){
        split.f <- SPATA2::joinWithFeatures(spata2_obj, features = across) %>% dplyr::select(barcodes, !!sym(across))
        df <- df %>% left_join(., split.f, by="barcodes")
        
        
        p <- purrr::map(.x=unique(df %>% pull(!!sym(across))), .f=function(sep.i){
          
          df <- df %>% filter(!!sym(across)==sep.i)
          
          # get thresholded Spots 
          # get thresholded Spots 
          gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
            if(!is.null(nr.sample)){
              df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
            }else{
              df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
            }
          })
          df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
          rownames(df.1) <- df$barcodes
          
          type <- purrr::map(.x=1:length(gs), .f=function(i){
            rep(gs[i], length(gl.ls[[i]]))
          }) %>% unlist()
          
          bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
          
          rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
          bc.1$module_score <- 0
          scores <- df.1[bc.1$inp, gs]
          
          for(i in names(scores)[1:5]){
            bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
              scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
          }
          bc <- bc.1$inp
          TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
          
          TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
          
          p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=1, fill=!!sym(feature)))+geom_col(size=1)+ theme_void()
          
          if( c(TC.cont %>% pull(!!sym(feature)) %>% class()) == "numeric"){
            message("feature is numeric")
            p=p+scale_fill_viridis_c()
          }
          
          df.to.test <- TC.cont %>% dplyr::select(module,!!sym(feature))
          names(df.to.test) <- c("variable", "sample")
          t=chisq.test(table(df.to.test$variable, df.to.test$sample),simulate.p.value = T)
          
          p=list(p,t)
          
          return(p)
          
          
        })
        
        
        
      }else{
        
        # get thresholded Spots 
        gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
          if(!is.null(nr.sample)){
            df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
          }else{
            df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
          }
        })
        df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
        rownames(df.1) <- df$barcodes
        
        type <- purrr::map(.x=1:length(gs), .f=function(i){
          rep(gs[i], length(gl.ls[[i]]))
        }) %>% unlist()
        
        bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
        
        rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
        bc.1$module_score <- 0
        scores <- df.1[bc.1$inp, gs]
        
        for(i in names(scores)[1:5]){
          bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
            scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
        }
        bc <- bc.1$inp
        TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
        
        TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
        
        p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=1, fill=!!sym(feature)))+geom_col(size=size)+ theme_void()
        
        df.to.test <- TC.cont %>% dplyr::select(module,!!sym(feature))
        names(df.to.test) <- c("variable", "sample")
        t=chisq.test(table(df.to.test$variable, df.to.test$sample),simulate.p.value = T)
        
        p=list(p,t)
        
      }
      
      
      
      
      
    }
    if(plot.type=="line"){
      
      if(!is.null(across)){
        split.f <- SPATA2::joinWithFeatures(spata2_obj, features = across) %>% dplyr::select(barcodes, !!sym(across))
        df <- df %>% left_join(., split.f, by="barcodes")
        
        
        p <- purrr::map(.x=unique(df %>% pull(!!sym(across))), .f=function(sep.i){
          
          df <- df %>% filter(!!sym(across)==sep.i)
          
          # get thresholded Spots 
          gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
            if(!is.null(nr.sample)){
              df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
            }else{
              df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
            }
          })
          df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
          rownames(df.1) <- df$barcodes
          
          type <- purrr::map(.x=1:length(gs), .f=function(i){
            rep(gs[i], length(gl.ls[[i]]))
          }) %>% unlist()
          
          bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
          
          rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
          bc.1$module_score <- 0
          scores <- df.1[bc.1$inp, gs]
          
          for(i in names(scores)[1:5]){
            bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
              scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
          }
          bc <- bc.1$inp
          TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
          
          TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
          
          p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=!!sym(feature),col=bc.1$module))+geom_line(size=size, alpha=alpha)+ theme_classic()
          
          df.to.test <- TC.cont %>% dplyr::select(module,!!sym(sample.type), !!sym(feature))
          names(df.to.test) <- c("variable", "sample", "value")
          anova <- aov( value~ variable * sample, data = df.to.test)
          t <- TukeyHSD(anova, which = "variable")$variable
          p=list(p,t)
          
          
        })
        
        
        
      }else{
        
        # get thresholded Spots 
        gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
          if(!is.null(nr.sample)){
            df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
          }else{
            df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
          }
        })
        df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
        rownames(df.1) <- df$barcodes
        
        type <- purrr::map(.x=1:length(gs), .f=function(i){
          rep(gs[i], length(gl.ls[[i]]))
        }) %>% unlist()
        
        bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
        
        rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
        bc.1$module_score <- 0
        scores <- df.1[bc.1$inp, gs]
        
        for(i in names(scores)[1:5]){
          bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
            scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
        }
        bc <- bc.1$inp
        TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
        
        TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
        
        p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=!!sym(feature),col=bc.1$module))+geom_line(size=size, alpha=alpha)+ theme_classic()
        
        df.to.test <- TC.cont %>% dplyr::select(module,!!sym(sample.type), !!sym(feature))
        names(df.to.test) <- c("variable", "sample", "value")
        anova <- aov( value~ variable * sample, data = df.to.test)
        t <- TukeyHSD(anova, which = "variable")$variable
        p=list(p,t)
        
      }
      
      
      
      
      
    }
    if(plot.type=="point"){
      
      if(!is.null(across)){
        split.f <- SPATA2::joinWithFeatures(spata2_obj, features = across) %>% dplyr::select(barcodes, !!sym(across))
        df <- df %>% left_join(., split.f, by="barcodes")
        
        
        p <- purrr::map(.x=unique(df %>% pull(!!sym(across))), .f=function(sep.i){
          
          df <- df %>% filter(!!sym(across)==sep.i)
          
          # get thresholded Spots 
          gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
            if(!is.null(nr.sample)){
              df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
            }else{
              df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
            }
          })
          df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
          rownames(df.1) <- df$barcodes
          
          type <- purrr::map(.x=1:length(gs), .f=function(i){
            rep(gs[i], length(gl.ls[[i]]))
          }) %>% unlist()
          
          bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
          
          rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
          bc.1$module_score <- 0
          scores <- df.1[bc.1$inp, gs]
          
          for(i in names(scores)[1:5]){
            bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
              scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
          }
          bc <- bc.1$inp
          TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
          
          TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
          
          p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=!!sym(feature),col=bc.1$module))+geom_point(size=size, alpha=alpha)+ theme_classic()
          
          df.to.test <- TC.cont %>% dplyr::select(module,!!sym(sample.type), !!sym(feature))
          names(df.to.test) <- c("variable", "sample", "value")
          anova <- aov( value~ variable * sample, data = df.to.test)
          t <- TukeyHSD(anova, which = "variable")$variable
          p=list(p,t)
          
          
        })
        
        
        
      }
      else{
        
        # get thresholded Spots 
        gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
          if(!is.null(nr.sample)){
            df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
          }else{
            df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
          }
        })
        df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
        rownames(df.1) <- df$barcodes
        
        type <- purrr::map(.x=1:length(gs), .f=function(i){
          rep(gs[i], length(gl.ls[[i]]))
        }) %>% unlist()
        
        bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
        
        rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
        bc.1$module_score <- 0
        scores <- df.1[bc.1$inp, gs]
        
        for(i in names(scores)[1:5]){
          bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
            scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
        }
        bc <- bc.1$inp
        TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
        
        TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
        
        #Problem of duplicated barcodes add .1 to duplicated bcs
        TC.cont$barcodes <- TC.cont$barcodes %>% make.unique(., ".")
        TC.cont$barcodes <- factor(TC.cont$barcodes, levels=TC.cont$barcodes)
        
        if(is.null(size)){
          p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=!!sym(feature),col=bc.1$module, size=!!sym(feature)))+geom_point(alpha=alpha)+ theme_classic()
        }else{
          p=ggplot(TC.cont, mapping=aes(x=1:nrow(TC.cont), y=!!sym(feature), col=bc.1$module))+geom_point(size=size, alpha=alpha)+ theme_classic()
        }
        
        
        df.to.test <- TC.cont %>% dplyr::select(module,!!sym(sample.type), !!sym(feature))
        names(df.to.test) <- c("variable", "sample", "value")
        anova <- aov( value~ variable * sample, data = df.to.test)
        t <- TukeyHSD(anova, which = "variable")$variable
        p=list(p,t, TC.cont)
        
        return(p)
      }
      
      return(p)
      
    
    }
    if(plot.type=="bar.comp"){
      
      if(!is.null(across)){
        split.f <- SPATA2::joinWithFeatures(spata2_obj, features = across) %>% dplyr::select(barcodes, !!sym(across))
        df <- df %>% left_join(., split.f, by="barcodes")
        
        
        p <- purrr::map(.x=unique(df %>% pull(!!sym(across))), .f=function(sep.i){
          
          df <- df %>% filter(!!sym(across)==sep.i)
          
          # get thresholded Spots 
          # get thresholded Spots 
          gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
            if(!is.null(nr.sample)){
              df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
            }else{
              df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
            }
          })
          df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
          rownames(df.1) <- df$barcodes
          
          type <- purrr::map(.x=1:length(gs), .f=function(i){
            rep(gs[i], length(gl.ls[[i]]))
          }) %>% unlist()
          
          bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
          
          rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
          bc.1$module_score <- 0
          scores <- df.1[bc.1$inp, gs]
          
          for(i in names(scores)[1:5]){
            bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
              scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
          }
          bc <- bc.1$inp
          TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
          
          TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(module, desc(module_score))
          
          p=ggplot(data=TC.cont, aes(x=module, y=1, fill=!!sym(feature)))+
            geom_bar(position="fill", stat="identity")+
            theme_classic()+
            #scale_fill_brewer(Feat,type='qual')+
            ylab("Percentage")+
            theme(
              plot.margin = margin(t = 50, r = 100, b = 50, l = 100, unit = "pt"),
              axis.text.y = element_text(color="black"),
              axis.text.x = element_text(color="black", angle = 75, vjust = .5))
          
          df.to.test <- TC.cont %>% dplyr::select(module,!!sym(feature))
          names(df.to.test) <- c("variable", "sample")
          t=chisq.test(table(df.to.test$variable, df.to.test$sample),simulate.p.value = T)
          
          p=list(p,t)
          
          return(p)
          
          
        })
        
        
        
      }else{
        
        # get thresholded Spots 
        gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
          if(!is.null(nr.sample)){
            df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
          }else{
            df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
          }
        })
        df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
        rownames(df.1) <- df$barcodes
        
        type <- purrr::map(.x=1:length(gs), .f=function(i){
          rep(gs[i], length(gl.ls[[i]]))
        }) %>% unlist()
        
        bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
        
        rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
        bc.1$module_score <- 0
        scores <- df.1[bc.1$inp, gs]
        
        for(i in names(scores)[1:5]){
          bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
            scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
        }
        bc <- bc.1$inp
        TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
        
        TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(module, desc(module_score))
        
        p=ggplot(data=TC.cont, aes(x=module, y=1, fill=!!sym(feature)))+
          geom_bar(position="fill", stat="identity")+
          theme_classic()+
          #scale_fill_brewer(Feat,type='qual')+
          ylab("Percentage")+
          theme(
            plot.margin = margin(t = 50, r = 100, b = 50, l = 100, unit = "pt"),
            axis.text.y = element_text(color="black"),
            axis.text.x = element_text(color="black", angle = 75, vjust = .5))
        
        df.to.test <- TC.cont %>% dplyr::select(module,!!sym(feature))
        names(df.to.test) <- c("variable", "sample")
        t=chisq.test(table(df.to.test$variable, df.to.test$sample),simulate.p.value = T)
        
        p=list(p,t)
        
      }
      
      
      
      
      
    }
    if(plot.type=="circular"){
      
      if(!is.null(across)){
        split.f <- SPATA2::joinWithFeatures(spata2_obj, features = across) %>% dplyr::select(barcodes, !!sym(across))
        df <- df %>% left_join(., split.f, by="barcodes")
        
        
        p <- purrr::map(.x=unique(df %>% pull(!!sym(across))), .f=function(sep.i){
          
          df <- df %>% filter(!!sym(across)==sep.i)
          
          # get thresholded Spots 
          # get thresholded Spots 
          gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
            if(!is.null(nr.sample)){
              df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
            }else{
              df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
            }
          })
          df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
          rownames(df.1) <- df$barcodes
          
          type <- purrr::map(.x=1:length(gs), .f=function(i){
            rep(gs[i], length(gl.ls[[i]]))
          }) %>% unlist()
          
          bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
          
          rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
          bc.1$module_score <- 0
          scores <- df.1[bc.1$inp, gs]
          
          for(i in names(scores)[1:5]){
            bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
              scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
          }
          bc <- bc.1$inp
          TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
          
          TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(module, desc(module_score))
          
          p=ggplot(data=TC.cont, aes(x=module, y=1, fill=!!sym(feature)))+
            geom_bar(position="fill", stat="identity")+
            theme_classic()+
            #scale_fill_brewer(Feat,type='qual')+
            ylab("Percentage")+
            theme(
              plot.margin = margin(t = 50, r = 100, b = 50, l = 100, unit = "pt"),
              axis.text.y = element_text(color="black"),
              axis.text.x = element_text(color="black", angle = 75, vjust = .5))
          
          df.to.test <- TC.cont %>% dplyr::select(module,!!sym(feature))
          names(df.to.test) <- c("variable", "sample")
          t=chisq.test(table(df.to.test$variable, df.to.test$sample),simulate.p.value = T)
          
          p=list(p,t)
          
          
          
        })
        
        
        
      }else{
        
        # get thresholded Spots 
        gl.ls <- purrr::map(.x=1:length(gs), .f=function(i){
          if(!is.null(nr.sample)){
            df %>%  arrange(desc(!!sym(gs[i]))) %>% head(nr.sample)  %>% pull(barcodes)
          }else{
            df %>% filter(!!sym(gs[i])>thr[i]) %>% arrange(desc(!!sym(gs[i])))  %>% pull(barcodes)
          }
        })
        df.1 <- df %>% dplyr::select({{gs}}) %>% as.data.frame()
        rownames(df.1) <- df$barcodes
        
        type <- purrr::map(.x=1:length(gs), .f=function(i){
          rep(gs[i], length(gl.ls[[i]]))
        }) %>% unlist()
        
        bc.1 <- data.frame(inp=gl.ls %>% unlist(), dup=duplicated(gl.ls %>% unlist()), module=type)
        
        rownames(bc.1) <- paste0(bc.1$inp,"_",1:nrow(bc.1))
        bc.1$module_score <- 0
        scores <- df.1[bc.1$inp, gs]
        
        for(i in names(scores)[1:5]){
          bc.1[bc.1 %>% filter(module==i) %>% rownames(), ]$module_score <- 
            scales::rescale(scores[bc.1 %>% filter(module==i) %>% pull(inp), i], c(0,1))
        }
        bc <- bc.1$inp
        TC.cont <- SPATA2::getFeatureDf(spata2_obj) %>% mutate(inp=barcodes) %>% as.data.frame() %>% filter(inp %in% bc.1$inp)
        
        TC.cont <- TC.cont %>% dplyr::left_join(.,bc.1, by="inp") %>% arrange(match(.$module, gs), desc(module_score))
        
        cir <-  
          TC.cont %>% 
          select( c(1:ncol(TC.cont))[str_detect(names(TC.cont), pattern = "Chr")]) %>% 
          reshape2::melt() %>% 
          rename("y":=value) %>% 
          mutate(x=1:nrow(.)) %>% 
          arrange(match(.$variable, paste0("Chr", 0:24)))
        
        colors <- scales::hue_pal()(5)
        colors <- scales::hue_pal()(8) %>% sample()
        cir$variable <- as.character(cir$variable)
        cir <- cir %>% 
          filter(!variable %in% c("Chr0","Chr24")) %>% 
          mutate(module=rep(TC.cont$module, length(unique(.$variable))))
        var <- TC.cont$module
        var2 <- TC.cont$Phase
        color=lapply(1:length(unique(var)), function(i){rep(scales::hue_pal()(length(unique(var)))[i], length(var[var==unique(var)[i]])) }) %>% unlist()
        color2=lapply(1:length(var2), function(i){var2[i] <- colors[match(var2[i], unique(var2))]}) %>% unlist()
        
        cir$y[is.na(cir$y)] <- 1
        cir$variable <- as.factor(cir$variable)
        cir$Prolif <- rep(TC.cont$G2M.Score, length(unique(cir$variable)))
        
        correspondance
        
        
        #Reorder levels
        cir$variable <- factor(cir$variable,paste0("Chr", 1:23))
        
        #Layer1
        circlize::circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0))
        circlize::circos.initialize(cir$variable, x=cir$x)
        #Layer1
        circlize::circos.track(cir$variable, x=cir$x,y=cir$y,ylim=c(0.9,1.1), panel.fun = function(x,y) {
          circlize::circos.points(x, y, col = color, pch = 16, cex = 0.1)
          circlize::circos.lines(x, y=rep(1,length(x)), col = "black")
          circlize::circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + circlize::mm_y(2), CELL_META$sector.index, niceFacing = TRUE)})
        
        #Layer2
        zoom <- cir %>% filter(variable=="Chr1"); zoom$variable <- as.character(zoom$variable) 
        circlize::circos.par(gap.degree = 2)
        circlize::circos.initialize(zoom$module, x=zoom$x)
        #Layer1
        circlize::circos.track(zoom$module, x=zoom$x,y=zoom$Prolif, panel.fun = function(x,y) {
          circlize::circos.lines(x, y=rep(0,length(x)), col = "black")
          circlize::circos.points(x, y, col = color2, cex=0.5, pch=16)})
        
        
        
        df.2 <- df.1[bc, gs] %>% as.matrix()
        col_fun1 = colorRamp2(c(0, 0.5, 1), rev(RColorBrewer::brewer.pal(3, "RdBu")))
        col= colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(3)
        circos.clear()
        circos.heatmap(df.2,col = col_fun1, cluster = F, split = zoom$module)
        
        
        p=ggplot(data=TC.cont, aes(x=module, y=1, fill=!!sym(feature)))+
          geom_bar(position="fill", stat="identity")+
          theme_classic()+
          #scale_fill_brewer(Feat,type='qual')+
          ylab("Percentage")+
          theme(
            plot.margin = margin(t = 50, r = 100, b = 50, l = 100, unit = "pt"),
            axis.text.y = element_text(color="black"),
            axis.text.x = element_text(color="black", angle = 75, vjust = .5))
        
        df.to.test <- TC.cont %>% dplyr::select(module,!!sym(feature))
        names(df.to.test) <- c("variable", "sample")
        t=chisq.test(table(df.to.test$variable, df.to.test$sample),simulate.p.value = T)
        
        p=list(p,t)
        
      }
      
      
      
      
      
    }
    
    
    
    return(p)
    
  }

#' @title  plotCNVPoints
#' @author Dieter Henrik Heiland
#' @description plotCNVPoints
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
plotCNVPoints <- 
  function(object,
           across = NULL,
           across_subset = NULL,
           relevel = NULL,
           clr = "blue",
           ...,
           of_sample = NA,
           verbose = NULL){
    
    # 1. Control --------------------------------------------------------------
    
    #hlpr_assign_arguments(object)
    
    of_sample <- SPATA2::check_sample(object, of_sample = of_sample, of.length = 1)
    
    # -----
    
    
    # 2. Data preparation -----------------------------------------------------
    
    # cnv results
    cnv_results <- getCnvResults(object, of_sample = of_sample)
    
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
        tibble::as_tibble() %>% 
        dplyr::mutate(mean=runif(n=nrow(.), 
                                 min=c(mean-sd), 
                                 max=c(mean+sd)) )
      
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
        ggplot2::geom_point(color=clr)+
        ggplot2::geom_vline(data = line_df, mapping = ggplot2::aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5) +
        ggplot2::theme_classic() +
        ggplot2::scale_x_continuous(breaks = line_df$label_breaks, labels = line_df$chromosome_name) +
        ggplot2::labs(x = "Chromosomes", y = "CNV-Results")
      
    } else if(base::is.character(across)){
      
      confuns::give_feedback(msg = glue::glue("Plotting cnv-results across '{across}'. This might take a few moments."),verbose = verbose)
      
      gene_names <- base::rownames(cnv_data)
      
      prel_df <-
        base::as.data.frame(cnv_data) %>%
        base::t() %>%
        base::as.data.frame() %>%
        tibble::rownames_to_column(var = "barcodes") %>%
        joinWith(object = object, spata_df = ., features = across) %>%
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
          cnv_mean = mean(x = cnv_values, na.rm = TRUE),
          cnv_sd = stats::sd(x = cnv_values, na.rm = TRUE)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(cnv_mean=runif(n=nrow(summarized_df), 
                                     min=c(cnv_mean-cnv_sd), 
                                     max=c(cnv_mean+cnv_sd)) ) %>% 
        dplyr::group_by(!!rlang::sym(x = across)) %>%
        dplyr::mutate(x_axis = dplyr::row_number())
      
      
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
      
      names(summarized_df)[1] <- "across"
      
      final_plot <-
        ggplot2::ggplot(data = summarized_df, mapping = ggplot2::aes(x = x_axis, y = cnv_mean)) +
        ggplot2::geom_point()+
        ggplot2::geom_vline(data = line_df,
                            mapping = ggplot2::aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5
        ) +
        ggplot2::facet_wrap(facets = ~ across) +
        ggplot2::theme_classic() +
        ggplot2::theme(strip.background = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(breaks = line_df$label_breaks, labels = line_df$chromosome_name) +
        ggplot2::labs(x = "Chromosomes", y = "CNV-Results")
      
    }
    
    confuns::give_feedback(msg = "Done.", verbose = verbose)
    
    base::return(final_plot)
    
  }



#' @title  addNucleusPosition
#' @author Dieter Henrik Heiland
#' @description addNucleusPosition
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

addNucleusPosition <- function(object, segment.df, scale_file, tissue_positions.file){
  
  #get position and rescale
  message(paste0("Process and scale:  ",sample.sp))
  sample <- SPATA2::getSampleNames(object)
  tissue.pos <- read.csv(tissue_positions.file, header=F)
  scale <- as.numeric(str_remove(read.table(scale_file,header=F)$V5, pattern=","))
  pos.data <- 
    tissue.pos %>% 
    dplyr::filter(V2==1) %>% 
    dplyr::mutate(x.pos=as.numeric(V6)*scale, 
                  y.pos=as.numeric(V5)*scale) %>% 
    dplyr::rename("barcodes":=V1) %>% 
    dplyr::select(barcodes,x.pos,y.pos)
  
  #rescale position.df
  coords <- SPATA2::getCoordsDf(object)
  rescale.df <- pos.data %>% left_join(.,coords,by="barcodes")
  rescale.factor <- rescale.df$x[1]/rescale.df$x.pos[1]
  segment.df.scaled <- 
    segment.df %>% 
    dplyr::mutate(x=x*rescale.factor, 
                  y=y*rescale.factor,
                  Cell=paste0("Cell_", ObjectNumber)) %>% 
    dplyr::select(Cell,x,y)
  
  
  object@spatial[[sample]]$Cell_coords <- segment.df.scaled
  
  return(object)
  
}

#' @title  getNucleusPosition
#' @author Dieter Henrik Heiland
#' @description getNucleusPosition
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

getNucleusPosition <- function(object){
  sample <- SPATA2::getSampleNames(object)
  if(!is.null(object@spatial[[sample]]$Cell_coords)){return(object@spatial[[sample]]$Cell_coords)} else{message("Cell_coords do not exist")}
}

#' @title  plot2DInterpolation
#' @author Dieter Henrik Heiland
#' @description plot2DInterpolation
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

plot2DInterpolation <- function(object, 
                                color_by,  
                                pt.size=0.5, 
                                pt.alpha=1,
                                normalize=T,
                                alpha2pred=T,
                                smooth=F, 
                                smooth_span=NULL, 
                                addImage=F,
                                Palette=NULL,
                                pt_clrsp="Reds",
                                ...){
  
  
  
  # Get Data ----------------------------------------------------------------
  
  scCoords <- getNucleusPosition(object)
  
  
  coords_df <- 
    getCoordsDf(object) %>% 
    SPATA2::hlpr_join_with_color_by(object = object, 
                            df = ., 
                            color_by = color_by, 
                            normalize = normalize, 
                            smooth = smooth, 
                            smooth_span = smooth_span)
  
  if(is.numeric(coords_df %>% pull(!!sym(color_by)))==T){
    
    message(paste0(Sys.time(), " ---- ", "Start 2D interpolation ", " ----"))
    
    x <- coords_df %>% pull(x)
    y <- coords_df %>% pull(y)
    z <- coords_df %>% pull(!!sym(color_by))
    
    
    
    
    s1 =  akima::interp(x = x, 
                        y = y, 
                        z = z, 
                        nx = nrow(coords_df), 
                        ny = nrow(coords_df), 
                        xo = seq(min(x), max(x), length = nrow(coords_df)), 
                        yo = seq(min(y), max(y), length = nrow(coords_df)))
    
    message(paste0(Sys.time(), " ---- ", "Predict single-cell expression levels ", " ----"))
    
    r.pred <- raster::raster(t(s1$z[,ncol(s1$z):1]), 
                             xmn = min(scCoords$x), xmx = max(scCoords$x),
                             ymn = min(scCoords$y), ymx = max(scCoords$y))
    
    pts <- sp::SpatialPointsDataFrame(scCoords[,c('x','y')], scCoords)
    scCoords$pred <- raster::extract(r.pred, pts)
    scCoords <- scCoords %>% filter(!is.na(pred))
    
    
    if(addImage==T){p=SPATA2::plotSurface(object, display_image=T, pt_alpha = 0)}else{p=ggplot()+theme_void()}
    if(alpha2pred==T){pt.alpha <- scCoords$pred}
    
    p=p+geom_point(data=scCoords, 
                   aes(x=x, y=y, color=pred), 
                   size=pt.size, 
                   alpha=pt.alpha)
    
    if(is.null(Palette)){p=p+SPATA2::scale_color_add_on(aes = "color", 
                                                        clrsp = pt_clrsp)}else{
                                                          p=p+scale_colour_gradientn(colours = Palette(50), oob=scales::squish,...)
                                                        }
    p=p+ggplot2::coord_equal()
    
    
    
    
    
    
  }else{
    
    message("missing")
    
    x <- coords_df %>% pull(x)
    y <- coords_df %>% pull(y)
    z <- coords_df %>% pull(!!sym(color_by))
    
    z <- as.factor(z) %>% as.numeric()
    
    s1 =  akima::interp(x = x, 
                        y = y, 
                        z = z, 
                        nx = nrow(coords_df), 
                        ny = nrow(coords_df), 
                        xo = seq(min(x), max(x), length = nrow(coords_df)), 
                        yo = seq(min(y), max(y), length = nrow(coords_df)))
    
    message(paste0(Sys.time(), " ---- ", "Predict single-cell expression levels ", " ----"))
    
    r.pred <- raster::raster(t(s1$z[,ncol(s1$z):1]), 
                             xmn = min(scCoords$x), xmx = max(scCoords$x),
                             ymn = min(scCoords$y), ymx = max(scCoords$y))
    
    pts <- sp::SpatialPointsDataFrame(scCoords[,c('x','y')], scCoords)
    scCoords$pred <- raster::extract(r.pred, pts)
    scCoords <- scCoords %>% filter(!is.na(pred))
    
    
    if(addImage==T){p=SPATA2::plotSurface(object, display_image=T, pt_alpha = 0)}else{p=ggplot()+theme_void()}
    scCoords$type <- round(scCoords$pred, digits = 0) 
    scCoords$pred <- abs(1-c(c(scCoords$pred-scCoords$type)))*pt.alpha
    
    levels <- data.frame(type=unique(z), real=coords_df[,color_by] %>% unique() %>% pull(!!sym(color_by)))
    levels <- levels %>% filter(type %in% unique(scCoords$type))
    scCoords$type <- as.factor(scCoords$type)
    scCoords$type <- factor(scCoords$type,  labels = levels$real)
    

    
    p=p+geom_point(data=scCoords, 
                   aes(x=x, y=y, color=type), 
                   size=pt.size,
                   alpha=scCoords$pred)
    p=p+ggplot2::coord_equal()
    
    
  }
  
  
  
  
  return(p)
  
}














