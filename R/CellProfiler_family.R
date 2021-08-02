#' @title ImportCellProfiler data
#' @author Dieter Henrik Heiland
#' @description ImportCellProfiler data
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

CellProfilerImport <- function(csv, 
                               type=NULL, 
                               var="Mean",
                               rename=T,
                               rescale=T){
  
  data <- utils::read.csv(csv) 
  
  #Quality check
  if(!any(names(data) %in% c("ImageNumber","ObjectNumber"))) stop("Variable: ImageNumber and / or ObjectNumber are not exist in the csv file")
  
  
  #check position
  if(any(names(data) %in% c("Location_Center_X","Location_Center_Y"))){
    position <- data %>% dplyr::select(ImageNumber,ObjectNumber,Location_Center_X, Location_Center_Y )
    names(position)[3:4] <- c("x", "y")
      
  }else{position <- NA}
  
  # check for variable name
  var.get <- names(data)[stringr::str_detect(names(data), pattern = var)]
  
  message(paste0(length(var.get), " matches of variables are detected: ", var.get))
  
  data <- data %>% dplyr::select(ImageNumber,ObjectNumber,{{var.get}})
  
  if(!is.null(type)){
    data <- data %>% dplyr::mutate(type={{type}})
  }else{
    data <- data %>% dplyr::mutate(type=NA)
    }
  
  if(rename==T){
    
    var.remove <- paste0("Intensity_",var,"Intensity_")
    names.new <- var.get %>% stringr::str_remove_all(., pattern = var.remove)
    names(data)[stringr::str_detect(names(data), pattern = var)] <- names.new
    
  }
  
  if(rescale==T){
    
    var.rescale <- data %>% dplyr::select(-ImageNumber,-ObjectNumber,-type) %>% names()
    
    for(i in var.rescale){data[, i] <- scales::rescale(data[, i], c(0,1))}
    
    }
  
  

  out <- list(data, position)
  names(out) <- c("data", "coords")
  return(out)
  
  
}



#' @title ImportCellProfiler data
#' @author Dieter Henrik Heiland
#' @description ImportCellProfiler data
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

CellProfilerImportObjectPositions <- function(csv, 
                                              type=NULL){
  
  data <- utils::read.csv(csv) 
  
  #Quality check
  if(!any(names(data) %in% c("ImageNumber","ObjectNumber"))) stop("Variable: ImageNumber and / or ObjectNumber are not exist in the csv file")
  
  
  #check position
  if(any(names(data) %in% c("Location_Center_X","Location_Center_Y"))){
    position <- 
      data %>% 
      dplyr::select(ImageNumber,ObjectNumber,Location_Center_X, Location_Center_Y )
    names(position)[3:4] <- c("x", "y") 
    } else{stop("No x,y variavles found")}
  
  if(!is.null(type)){
    position <- position %>% dplyr::mutate(type={{type}})
  }else{
    position <- position %>% dplyr::mutate(type=NA)
  }
  
  return(position)
  
  
}


#' @title ImportCellProfiler data
#' @author Dieter Henrik Heiland
#' @description ImportCellProfiler data
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
CellProfilerCellDistance <- function(data){
  
  
  #check number of images
  images <- unique(data$ImageNumber)
  cell.types <- unique(data$type)
  
  data.selected <- 
  purrr::map(.x=images, function(img){
    
    img.selected <- purrr::map(.x=cell.types, function(c.t){
      
      selected <- 
        data %>% 
        dplyr::filter(type==c.t) %>% 
        dplyr::filter(ImageNumber==img) %>% 
        mutate(ID=paste0(c.t,"_", 1:nrow(.)))
      
      return(selected)
    })
    
    permutations <- 
      permute::allPerms(length(img.selected)) %>% 
      reshape2::melt() %>% select(Var2,value) %>% 
      mutate(sum=Var2+value, 
             duplicated=duplicated(sum), 
             dup2=c(Var2==value) ) %>% 
      filter(duplicated==F) %>% 
      filter(dup2==F) %>% 
      select(Var2, value)
    
    names(permutations) <- c("From", "To")
    
    distance.list <- purrr::map(.x=1:nrow(permutations), function(i){
      
      S1 <- permutations$From[i]
      S2 <- permutations$To[i]
      #matrix 
      mat=matrix(NA, nrow=nrow(img.selected[[S1]]), ncol=nrow(img.selected[[S2]]))
      rownames(mat) <- img.selected[[S1]]$ID
      colnames(mat) <- img.selected[[S2]]$ID
      
      #loop distance
      for(ix in 1:nrow(img.selected[[S1]])){
        S1.df=img.selected[[S1]][ix, ]
        for(ii in 1:nrow(img.selected[[S2]])){
          S2.df <- img.selected[[S2]][ii, ]
          
          distance <- function(x1, y1, x2, y2) {
            sqrt((x2 - x1)^2 + (y2 - y1)^2)}
          mat[ix,ii] <- distance(S1.df$x,S1.df$y, S2.df$x, S2.df$y)
          
          
        }
      }
      
      all.con <- reshape2::melt(mat)
      names(all.con) <- c("From", "To", "dist")
      
      #Add coordinaes
      all.cor <- as.data.frame(do.call(rbind, img.selected))
      
      cor_from <- all.cor %>% dplyr::left_join(data.frame(ID=all.con$From), ., by="ID") %>% dplyr::select(x,y)
      names(cor_from) <- c("x1", "y1")
      cor_to <- all.cor %>% dplyr::left_join(data.frame(ID=all.con$To), ., by="ID") %>% dplyr::select(x,y)
      names(cor_to) <- c("x2", "y2")
      
      all.con <- cbind(all.con, cor_from, cor_to)
      
      return(all.con)
      
    })
    
    
    out <- list(img.selected, permutations, distance.list)
    names(out) <- c("img.selected", "permutations", "distance.list")
    
    return(out)
    
    
  })
  
  
  
}



#' @title ImportCellProfiler data
#' @author Dieter Henrik Heiland
#' @description ImportCellProfiler data
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
CellProfilerPlotDistance <- function(data.dist, 
                                     image=1, 
                                     color=RColorBrewer::brewer.pal(8, "Dark2"), 
                                     pt.size=4,lt.size=2, 
                                     pt.alpha=0.8,
                                     lt.alpha=0.8, 
                                     dist.threshold=200 ){
  
  data <- data.dist[[image]]
  permutations <- data$permutations
  
  p <- ggplot2::ggplot() + theme_classic()
  
  #Add layer connections
  for(i in 1:length(data$distance.list)){
    print(i)
    relationship <- data$distance.list[i] %>% as.data.frame() %>% dplyr::filter(dist<dist.threshold)
    p=p+geom_segment(data=relationship, mapping = aes(x=x1, y=y1, xend=x2, yend=y2),color=color[i], alpha=lt.alpha, size=lt.size)
  }
  
  #Add layer of cells
  p=p+ggplot2::geom_point(data=as.data.frame(do.call(rbind,data$img.selected)), mapping=aes(x,y,color=type), size=pt.size, alpha=pt.alpha)
  
  return(p)
}

#' @title ImportCellProfiler data
#' @author Dieter Henrik Heiland
#' @description ImportCellProfiler data
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
CellProfilerPlotDistanceAddLayer <- function(p, image=1, data, color.to,pt.size=4,pt.alpha=0.8,limit=c(0,1)){
  
  plot.df=cbind(data$data, data$coords[,c("x", "y")]) %>% dplyr::filter(ImageNumber==image)
  
  if(is.numeric(plot.df[,color.to])){
  p=p+ggplot2::geom_point(data=plot.df, mapping=aes(x,y, fill=!!sym(color.to)),colour = "black", size=pt.size,alpha=pt.alpha, shape=21)+
    scale_fill_viridis_c(limit=limit, option="C")
  }else{
    p=p+ggplot2::geom_point(data=plot.df, mapping=aes(x,y, fill=!!sym(color.to)),colour = "black", size=pt.size,alpha=pt.alpha, shape=21)
  }
  
  
  
  return(p)
}


#' @title ImportCellProfiler data
#' @author Dieter Henrik Heiland
#' @description ImportCellProfiler data
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
CellProfilerDistancePlot <- function(data,prefix, data.dist, exclude=NULL){
  
  message("The distance between the target cells in 'data' will be compared to cell types in the 'data.dist' file ")
  
  images <- 1:length(data.dist)
  merged.cell.pos <- as.data.frame(do.call(rbind,data.dist[[1]]$img.selected))
  cell.types <-  merged.cell.pos %>% dplyr::pull(type) %>% unique()
  
  if(!is.null(exclude)){cell.types <- cell.types[!cell.types %in% exclude ] }
  
  distance <- 
  purrr::map(.x=images, .f=function(img){
    
    img.dist <- purrr::map(.x=1:length(data.dist[[img]]$distance.list), .f=function(type.cell){
      
      distance <- data.dist[[img]]$distance.list[[type.cell]]
      data.join <- 
        data$data %>% 
        dplyr::filter(ImageNumber==img) %>% 
        dplyr::mutate(ID = paste0(prefix,"_",ObjectNumber))
      
      if(!stringr::str_detect(distance$From[1], pattern = prefix)){
        a <- distance$To; distance$To <- distance$From;distance$From <- a
      }
      
      distance <- 
      distance %>% 
        dplyr::mutate(ID=From) %>% 
        dplyr::select(ID, To, dist) %>% 
        dplyr::left_join(.,  data.join, by="ID" )
    
    return(distance)
    })
    vv <- purrr::map(.x=img.dist, .f=function(i){is.na(i$ImageNumber[1])}) %>% unlist()
    img.dist <- img.dist[!vv]
  
  
    return(img.dist)
  })
  
  
  #Sum up images
  sum.dist <- purrr::map(.x=1:length(cell.types), .f=function(i){
    
    plot.df <- as.data.frame(do.call(rbind,lapply(images, function(ii){distance[[ii]][[i]] })))
    p=ggplot2::ggplot()+theme_classic()
    p=p+ggplot2::geom_density(data=plot.df, mapping = aes(x=dist, fill=type),alpha=0.3)
    return(p)
  })
  
  names(sum.dist) <- cell.types
  
  return(sum.dist)
  
  
   
}





