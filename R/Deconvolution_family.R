


#' @title  runRCTD cell-type deconvolution
#' @author Dieter Henrik Heiland
#' @description Run the a cell type deconvolution using the SpaceXR/RCTD algorithm from SPATA objects
#' @inherit 
#' @param ref The reference seurat object
#' @param cell_type_var The variale in the Seurat meta.data file containing the cell type group
#' @param object The variale in the Seurat meta.data file containing the cell type group
#' @return 
#' @examples 
#' @export
#' 
runRCTD <- function(object,ref, cell_type_var, overwrite=T){
  
  #Some check up
  if(!any("cell_type_var" %in% names(ref@meta.data))) stop(paste0("The variable: ",cell_type_var, " was not found in the seurat object"))
  SPATA2::check_object(object)
  
  #Load some packages
  library(spacexr)
  library(Matrix)
  
  #Get reference data
  counts <- ref@assays$RNA@counts %>% as.matrix()
  meta_data <- 
    ref@meta.data %>% 
    as.data.frame() %>% 
    dplyr::select({{cluster}}, nFeature_RNA)
  
  
  #Creat RCTD object
  cell_types <- meta_data[, cluster]
  names(cell_types) <- rownames(meta_data)
  cell_types <- as.factor(cell_types)
  nUMI <- colSums(counts)
  reference <- Reference(counts, cell_types, nUMI)
  counts <- SPATA2::getCountMatrix(object) %>% as.matrix()
  coords <- SPATA2::getCoordsDf(object) %>% dplyr::select(barcodes,x, y) %>% as.data.frame()
  rownames(coords) <- coords$barcodes
  coords <- coords[, c(2:3)]
  names(coords) <- c("xcoord", "ycoord")
  nUMI <- colSums(counts)
  
  #Run Deconv
  puck <- spacexr::SpatialRNA(coords, counts, nUMI)
  myRCTD <- spacexr::create.RCTD(puck, reference, max_cores = 5)
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = "doublet")
  
  
  
  # Get Results
  results <- myRCTD@results
  norm_weights = spacexr::normalize_weights(results$weights)
  cell_type_names <- myRCTD@cell_type_info$info[[2]]
  spatialRNA <- myRCTD@spatialRNA
  
  #Create output
  
  norm_weights <- 
    norm_weights %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("barcodes") %>% 
    dplyr::left_join(object %>% 
                       SPATA2::getCoordsDf() %>% 
                       dplyr::select(barcodes),., by = "barcodes")
  
  norm_weights[is.na(norm_weights)] <- 0
  
  object <- object %>% SPATA2::addFeatures(., norm_weights, overwrite = overwrite)
  return(object)
}


#' @title  runRCTD cell-type deconvolution
#' @author Dieter Henrik Heiland
#' @description Run the a cell type deconvolution using the SpaceXR/RCTD algorithm from SPATA objects
#' @inherit 
#' @object SPATA2 object
#' @deconv_cell_types Character value of all celltypes (from the feature data.frame) should be used
#' @spot_extension Numeric value between 0 and 0.7. Specifies extension of the spot radius to include more cells. For example, a spot extension of 0.2 will increase the spot radius of 20%) 
#' @multicore Logical. If TRUE will be run on all cores (recommended)
#' @flip.y Logical. If TRUE y axis will be flipped
#' @scale_factor Numeric value, if adjustment of the SPATAwrappers::getNucleusPosition() is required
#' @return 
#' @examples 
#' @export
#'


getSingleCellDeconv <- function(object,
                                deconv_cell_types,
                                spot_extension=NULL,
                                multicore=T,
                                flip.y=T,
                                workers=16,
                                scale_factor=NULL){
  
  SPATA2::check_method(object)
  
  if(spot_extension>0.7) stop("The spot extension will cause overlap in segmentation ")
  
  
  
  # First quantify the numbers of cell types in each spot and quantify the celltypes
  sc_dat <- SPATAwrappers::getNucleusPosition(object)
  
  #If required flipping and scaling of single cell coords
  if(flip.y==T){
    message("Y-axis will be fliped")
    yrange <- range(sc_dat$y)
    sc_dat$y <- yrange[2] - sc_dat$y + yrange[1]
    
  }
  if(!is.null(scale_factor)){
    sc_dat[,c("x","y")]=sc_dat[,c("x", "y")]*scale_factor
    
    
  }
  
  #Get cells per spot
  
  #-define the spotsize radius
  d <- SPATAwrappers::getSurroundedSpots(object) %>% filter(distance!=0) %>%  pull(distance) %>% min()
  r=(d*c(55/100))/2
  
  if(!is.null(spot_extension)){r=r+(r*spot_extension)}
  
  grid.plot <- SPATA2::getCoordsDf(object)
  
  
  if(multicore==T){
    
    # Run multicore
    base::options(future.fork.enable = TRUE)
    future::plan("multiprocess", workers = workers)
    future::supportsMulticore()
    base::options(future.globals.maxSize = 600 * 1024^2)
    message("... Run multicore ... ")
    
    plot(sc_dat$x, sc_dat$y, pch=".")
    segments <- furrr::future_map_dfr(.x=1:nrow(grid.plot), 
                                      .f=function(x){
                                        
                                        
                                        
                                        #print(x)
                                        #Create Segment
                                        segment <- swfscMisc::circle.polygon(x=grid.plot$x[x],
                                                                             y=grid.plot$y[x],
                                                                             radius=r,
                                                                             poly.type= "cartesian") %>% as.data.frame()
                                        #polygon(segment, border="red")
                                        
                                        #Get Objects in Segment
                                        nuc <- sp::point.in.polygon(pol.x=segment$x, 
                                                                    pol.y=segment$y, 
                                                                    point.x=sc_dat$x, 
                                                                    point.y=sc_dat$y)
                                        
                                        cells <- sc_dat[nuc==1, ]$Cell
                                        cells_in_spot<-sum(nuc)
                                        return <- 
                                          data.frame(barcodes=grid.plot$barcodes[x], 
                                                     Nr_of_cells=cells_in_spot,
                                                     cells=cells)
                                        
                                        
                                        
                                        return(return)
                                        
                                        
                                      },
                                      .progress = T)
    
    
    
  }else{
    
    plot(sc_dat$x, sc_dat$y, pch=".")
    segments <- map_dfr(.x=nrow(grid.plot), .f=function(x){
      
      
      
      #print(x)
      #Create Segment
      segment <- swfscMisc::circle.polygon(x=grid.plot$x[x],
                                           y=grid.plot$y[x],
                                           radius=r,
                                           poly.type= "cartesian") %>% as.data.frame()
      polygon(segment, border="red")
      
      #Get Objects in Segment
      nuc <- sp::point.in.polygon(pol.x=segment$x, 
                                  pol.y=segment$y, 
                                  point.x=sc_dat$x, 
                                  point.y=sc_dat$y)
      
      cells <- sc_dat[nuc==1, ]$Cell
      cells_in_spot<-sum(nuc)
      return <- 
        data.frame(barcodes=grid.plot$barcodes[x], 
                   Nr_of_cells=cells_in_spot,
                   cells=cells)
      
      
      
      return(return)
      
      
    })
    
  }
  
  
  deconv_df <- SPATA2::getFeatureDf(object) %>% dplyr::select(barcodes,{{deconv_cell_types}})
  
  # Run a cell type quantification per spot
  Cell_types <- furrr::future_map_dfr(.x=1:nrow(grid.plot), 
                                      .f=function(x){
                                        
                                        #print(x)
                                        
                                        out <- 
                                          segments %>% 
                                          dplyr::filter(barcodes==grid.plot$barcodes[x])
                                        
                                        #Get the best cells
                                        nr.cells <- nrow(out)
                                        
                                        scores <- 
                                          deconv_df %>% 
                                          dplyr::filter(barcodes==grid.plot$barcodes[x]) %>% 
                                          dplyr::select(-barcodes) %>% 
                                          t() %>% as.data.frame() %>% 
                                          dplyr::arrange(desc(V1)) %>% 
                                          dplyr::mutate(score=scales::rescale(V1, c(0,1)) %>% round(.,digits = 3)) %>% 
                                          dplyr::mutate(cells=round(score*c(nr.cells/sum(score))+0.5 ,digits = 0)) %>%
                                          dplyr::select(-V1)
                                        
                                        i <- 1
                                        cells_select=0
                                        while (cells_select < nr.cells) {
                                          cells_select <- cells_select+scores$cells[i]
                                          #print(cells_select)
                                          i=i+1
                                        }
                                        scores <- scores[1:i-1, ]
                                        
                                        if(sum(scores$cells)>nr.cells){
                                          scores$cells[nrow(scores)]= scores$cells[nrow(scores)]-1
                                        }
                                        
                                        cells_add <- map(.x=1:nrow(scores),.f=function(i){
                                          rep(rownames(scores)[i], scores$cells[i])}) %>% unlist()
                                        
                                        out$celltypes <- cells_add
                                        
                                        return(out)
                                      },
                                      .progress = T)
  
  
  #Align Coords
  
  Cell_types <- Cell_types %>% dplyr::left_join(., sc_dat %>% rename("cells":=Cell), by="cells")
  
  return(Cell_types)
  
}
















