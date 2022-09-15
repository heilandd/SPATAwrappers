


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
runRCTD <- function(object,ref, cell_type_var, overwrite=T, return.RCTD=F){
  
  #Some check up
  if(!any(cell_type_var %in% names(ref@meta.data))) stop(paste0("The variable: ",cell_type_var, " was not found in the seurat object"))
  SPATA2::check_object(object)
  cluster <- cell_type_var
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
  
  if(return.RCTD==T){
    object <- myRCTD
  }else{
    
    
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

    
  }
  

  
  return(object)
}


#' @title  runRCTD cell-type deconvolution
#' @author Dieter Henrik Heiland
#' @description Run the a cell type deconvolution using the SpaceXR/RCTD algorithm from SPATA objects
#' @inherit 
#' @param object SPATA2 object
#' @param deconv_cell_types Character value of all celltypes (from the feature data.frame) should be used
#' @param spot_extension Numeric value between 0 and 0.7. Specifies extension of the spot radius to include more cells. For example, a spot extension of 0.2 will increase the spot radius of 20%) 
#' @param multicore Logical. If TRUE will be run on all cores (recommended)
#' @param flip.y Logical. If TRUE y axis will be flipped
#' @param scale_factor Numeric value, if adjustment of the SPATAwrappers::getNucleusPosition() is required
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

#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @inherit 
#' @param object SPATA2 object
#' @param scDF Data.frame; Output of the getSingleCellDeconv()
#' @param ref.new Seurat Object; Reference dataset used for runRCTD()
#' @param cell_type_var Character value; The col of the Seurat meta data indicating the cell type annotations
#' @param model Character value; Model "BF" randomly select cell compositions and select the best match. Model "oGA" will perform a conditional genetic algorithm to find the best match.
#' @param subset_ref Logical, if TRUE, Seurat object will we downsized to increase speed.
#' @param only_var Logical, if TRUE only consider variable genes for mapping
#' @param n_feature Integer value: Number of variable genes
#' @param max_cells Integer value: Number of cells for downsampling
#' @param AE_norm  Logical, if TRUE use an autoencoder for data integration and normalization (recommended)
#' @param multicore Logical, if TRUE use multicore
#' @param workers Integer, Number of cors 
#' @param iter Integer: For model BF: Number of random spot compositions; For model oGA size of initial population
#' @param iter_GA Integer: Number of iterations of the oGA model
#' @param nr_mut Integer: Number of mutations
#' @param nr_offsprings Integer: Number of offsprings
#' @param cross_over_point Numeric value: Percentage of cross-over cutt-off
#' @param ram Integer: GB of ram can be used for multicore session
#' @return 
#' @examples 
#' @export
#'


runMapping <- function(object,
                       scDF, 
                       ref.new, 
                       cell_type_var="annotation_level_4",
                       model=c("oGA", "BF"),
                       subset_ref=T,
                       only_var=T,
                       n_feature=3000,
                       max_cells=10000,
                       AE_norm=F,
                       dropout=0,
                       bottleneck=16,
                       activation="relu",
                       layers=c(128, 64, 32),
                       epochs=20,
                       multicore=T,
                       workers=16,
                       iter=200,
                       iter_GA=20,
                       nr_mut=2,
                       nr_offsprings=7,
                       cross_over_point=0.5,
                       ram=60
){
  
  
  if(subset_ref==T){
    
    message("--- Downscale Ref Dataset ----")
    ref.new <- DownScaleSeurat(ref.new, 
                               maintain_var=cell_type_var,
                               n_feature=n_feature,
                               only_var=only_var,
                               max=max_cells)
  }
  
  if(!is.null(multicore)){
    
    base::options(future.fork.enable = TRUE)
    future::plan("multisession", workers = workers)
    future::supportsMulticore()
    base::options(future.globals.maxSize = ram* 10* 1024^2)
    message("... Run multicore ... ")
  }
  
  spots <- unique(scDF$barcodes)
  
  #Load data that are required
  message("--- Load data ----")
  mat.ref <- ref.new %>% Seurat::GetAssayData()
  genes.ref <- rownames(mat.ref)
  
  object <- SPATA2::setActiveExpressionMatrix(object, "scaled")
  mat.spata <- SPATA2::getExpressionMatrix(object)
  genes.spata <- rownames(mat.spata)
  
  message("--- Merge Data ----")
  genes <- intersect(genes.ref, genes.spata)
  mat.ref <- mat.ref[genes, ]
  mat.spata <- mat.spata[genes, ]
  
  
  if (AE_norm==T){
    
    message("--- Normalisation and data integration using AE ----")
    
    
    x_train <- cbind(mat.spata, mat.ref)
    x_train <- scales::rescale(x_train, c(0,1))
    
    input_layer <- keras::layer_input(shape = c(ncol(x_train)))
    
    encoder <- 
      input_layer %>% 
      keras::layer_dense(units = layers[1], activation = activation) %>% 
      #keras::layer_batch_normalization() %>% 
      keras::layer_dropout(rate = dropout) %>% 
      keras::layer_dense(units = layers[2], activation = activation) %>% 
      keras::layer_dropout(rate = dropout) %>% 
      keras::layer_dense(units = layers[3], activation = activation) %>% 
      keras::layer_dense(units = bottleneck)
    
    
    decoder <- 
      encoder %>% 
      keras::layer_dense(units = layers[3], activation = activation) %>% 
      keras::layer_dropout(rate = dropout) %>% 
      keras::layer_dense(units = layers[2], activation = activation) %>% 
      keras::layer_dropout(rate = dropout) %>% 
      keras::layer_dense(units = layers[1], activation = "sigmoid") %>% 
      keras::layer_dense(units = c(ncol(x_train)))
    
    
    autoencoder_model <- keras::keras_model(inputs = input_layer, outputs = decoder)
    
    autoencoder_model %>% keras::compile(loss = "mean_squared_error", optimizer = "adam", metrics = c("accuracy"))
    
    history <- autoencoder_model %>% keras::fit(x_train, x_train, 
                                                epochs = epochs, shuffle = T, 
                                                verbose = T)
    
    reconstructed_points <- autoencoder_model %>% keras::predict_on_batch(x = x_train)
    
    mat.spata_ref <- reconstructed_points[,1:ncol(mat.spata)]
    mat.ref_ref <- reconstructed_points[,1:ncol(mat.ref)]
    rownames(mat.ref_ref) <- genes
    rownames(mat.spata_ref) <- genes
    
    
  }
  
  
  
  
  
  message("--- Run Mapping ----")
  
  nr_of_random_spots=iter
  
  ref_meta <- 
    ref.new@meta.data[colnames(mat.ref),c("nCount_RNA",cell_type_var)] %>% 
    as.data.frame()
  
  if(model=="BF"){
    data.new <- 
      furrr::future_map_dfr(.x=1:length(spots), 
                            .f=function(i){
                              
                              #print(i)
                              #Create 1000 random spot compositions
                              bc_run <- spots[i]
                              data <- scDF %>% dplyr::filter(barcodes==bc_run) %>% arrange(celltypes)
                              nr_cells <- nrow(data)
                              n_select <- data %>% count(celltypes)
                              
                              random_spots <- list()
                              cor.list <- list()
                              
                              for( i in 1:nr_of_random_spots){
                                random_spots[[i]] <- 
                                  map(.x=1:nrow(n_select), 
                                      .f=function(x){
                                        sample(rownames(ref_meta[ref_meta[,cell_type_var]==n_select$celltypes[x], ]), n_select$n[x]) }) %>% unlist()
                                if(nr_cells==1){
                                  cor.list[[i]] <- cor(as.numeric(mat.ref[,random_spots[[i]]]),
                                                       as.numeric(mat.spata[,bc_run]))
                                }else{
                                  cor.list[[i]] <- cor(as.numeric(mat.ref[,random_spots[[i]]] %>% rowMeans()),
                                                       as.numeric(mat.spata[,bc_run]))
                                }
                                
                                
                              }
                              
                              
                              validate_randoms_select <- 
                                cor.list %>% 
                                unlist() %>% 
                                as.data.frame() %>% 
                                rownames_to_column("order") %>% 
                                rename("cor":=.) %>% 
                                arrange(desc(cor))
                              
                              select_cells <- random_spots[[as.numeric(validate_randoms_select$order[1])]]
                              
                              data$best_match <- select_cells
                              
                              return(data)
                              
                              
                            }, .progress = T, .options = furrr::furrr_options(seed = TRUE))
  }else{
    
    if(model=="oGA"){
      
      #nested groups
      nested_ref_meta <- ref_meta %>% rownames_to_column("cells") %>% group_by(!!sym(cell_type_var)) %>% nest()
      
      
      ## Functions
      
      fitness <- function(x,nr_cells){
        
        if(nr_cells==1){
          y <- cor(as.numeric(mat.ref[,x==1]),
                   as.numeric(mat.spata[,bc_run]))
        }else{
          y <- cor(as.numeric(mat.ref[,x==1] %>% rowMeans()),
                   as.numeric(mat.spata[,bc_run]))
        }
        return(y)
      }
      
      initiate_Population <- function(nr_of_random_spots, n_select,nested_ref_meta){
        
        mat <- matrix(0, nrow=nr_of_random_spots, ncol=ncol(mat.ref))
        colnames(mat)=colnames(mat.ref)
        select <- nested_ref_meta %>% filter(!!sym(cell_type_var) %in%  n_select$celltypes)
        
        for(x in 1:nr_of_random_spots){
          selected_cells <- map(.x=1:nrow(select), .f=function(i){
            sample(select$data[i] %>% as.data.frame() %>% pull(cells), n_select$n[i])
          }) %>% unlist()
          mat[x,selected_cells] <- 1
        }
        
        return(mat)
        
      }
      
      cross_over <- function(parents,cross_over_point=0.5){
        
        parents_out <- 
          parents %>%
          apply(1, function(x) which(x == 1)) %>% t()
        #as.data.frame() %>% 
        #filter(!V1 %in% intersect(V1,V2)) %>% 
        #filter(!V2 %in% intersect(V1,V2)) %>% 
        #t()
        
        #message(length(intersect(parents_out[1,], parents_out[2,])))
        
        cross_over_select <- c(ncol(parents_out)*cross_over_point) %>% round()
        parents_out <- parents_out[, 1:cross_over_select]
        
        
        # cross values
        parents[1, parents_out[1,]]=0
        parents[1, parents_out[2,]]=1
        
        parents[2, parents_out[2,]]=0
        parents[2, parents_out[1,]]=1
        
        #message(length(which(parents[1, ]==1)))
        #message(length(which(parents[2, ]==1)))
        
        
        
        return(parents)
      }
      
      mutation_GA <- function(offspring_2,nr_mut, nr_offsprings){
        
        #Create random selection
        selector <- runif(nr_offsprings, 1, 2) %>% round(digits = 0)
        
        # Create a mutated and non-mutated output
        offsprings_mut <- matrix(0, ncol = ncol(offspring_2), nrow=nr_offsprings)
        colnames(offsprings_mut) <- colnames(offspring_2)
        
        
        #Run loop for offsprings
        
        for(u in 1:nr_offsprings){
          
          
          #Select random cell and mutate from same cell type
          cells_off <- which(offspring_2[selector[u], ]==1) %>% names()
          cells_mut <- sample(cells_off, nr_mut)
          # Get cell type
          
          ref_sub <- ref_meta[!rownames(ref_meta) %in% cells_off, ]
          cell_type_select <- ref_meta[cells_mut, cell_type_var]
          
          new <- list()
          for(z in 1:nr_mut){
            new[[z]] <- sample(ref_sub %>% 
                                 filter(!!sym(cell_type_var)==cell_type_select[z]) %>% 
                                 rownames(),1)
            
          }
          
          offsprings_inter <- offspring_2
          
          #message(length(cells_mut)==length(unlist(new)))
          
          offsprings_inter[selector[u],cells_mut]=0
          offsprings_inter[selector[u],unlist(new)]=1
          
          #message(length(which(offsprings_inter[selector[u], ]==1)))
          #message(any((cells_off %in% unlist(new)) == T))
          
          offsprings_mut[u, ] <- offsprings_inter[selector[u], ]
          
          
        }
        
        return(offsprings_mut)
      }
      
      
      data.new <- 
        furrr::future_map_dfr(.x=1:length(spots), 
                              .f=function(i){
                                message("--- Model: oGA will be applied")
                                #print(i)
                                
                                
                                
                                bc_run <- spots[i]
                                data <- 
                                  scDF %>% 
                                  dplyr::filter(barcodes==bc_run) %>% 
                                  dplyr::mutate(celltypes=as.character(celltypes)) %>% 
                                  dplyr::arrange(celltypes) 
                                
                                nr_cells <- nrow(data)
                                n_select <- data %>% count(celltypes)
                                
                                #Initiate population
                                pop <- initiate_Population(nr_of_random_spots, n_select,nested_ref_meta)
                                
                                qc <- list()
                                
                                for(zz in 1:iter_GA){
                                  
                                  # Validate the initial Pop
                                  validate_randoms_select <- 
                                    map(.x=1:nr_of_random_spots, function(j){fitness(pop[j,], nr_cells)}) %>% 
                                    unlist() %>% 
                                    as.data.frame() %>% 
                                    rownames_to_column("order") %>% 
                                    rename("cor":=.) %>% 
                                    arrange(desc(cor))
                                  
                                  #select parents
                                  parents <- pop[as.numeric(validate_randoms_select$order[1:2]),]
                                  
                                  
                                  
                                  qc[[zz]] <- mean(validate_randoms_select$cor[1:2])
                                  print(qc[[zz]])
                                  
                                  #Create Children
                                  offspring_2 <- cross_over(parents,cross_over_point)
                                  #message(length(which(offspring_2[1, ]==1)))
                                  #message(length(which(offspring_2[2, ]==1)))
                                  
                                  offspring <- mutation_GA(offspring_2,nr_mut, nr_offsprings)
                                  
                                  
                                  # remove old parents
                                  remove <- as.numeric(validate_randoms_select %>% tail(dim(offspring)[1]) %>% pull(order))
                                  pop.new <- rbind(pop[-remove, ], offspring)
                                  
                                  #update pop
                                  pop <- pop.new
                                  
                                }
                                
                                validate_randoms_select <- 
                                  map(.x=1:nr_of_random_spots, function(j){fitness(pop[j,], nr_cells)}) %>% 
                                  unlist() %>% 
                                  as.data.frame() %>% 
                                  rownames_to_column("order") %>% 
                                  rename("cor":=.) %>% 
                                  arrange(desc(cor))
                                
                                pop_select <- pop[as.numeric(validate_randoms_select$order[1]), ]
                                
                                
                                select_cells <- names(which(pop_select==1))
                                
                                data$best_match <- select_cells
                                
                                #check
                                #data$anno <- ref.new@meta.data[data$best_match, ]$annotation_level_4
                                
                                
                                return(data)
                                
                                
                              }, .progress = T, .options = furrr::furrr_options(seed = TRUE))
    } else{
      
      stop("Model unknown")
    }
    
  }
  
  
  
  
  return(data.new)
  
  
}




#' @title  Run single-cell mapping
#' @author Dieter Henrik Heiland
#' @description Mapping the output of "getSingleCellDeconv()"to the best matching cell from the reference dataset
#' @inherit 
#' @param tab Input data of the class data.frame
#' @param class Character value; The col containing the main class
#' @param subclass Character value; The col containing the subclass
#' @param pal Character value; Color pal from brewer.pal()
#' @param random Logical. If TRUE random mixing colors from pal
#' @param seed Integer value. set.seed() for constant color alignment
#' @param into Character value; The color for non-classified samples
#' @param add_n Integer value. Adopting the contrast of the subclass colors
#' @return 
#' @examples 
#' @export
#'
getSubColors <- function(tab,
                         class="annotation_level_2", 
                         subclass="annotation_level_4", 
                         pal="Set3", 
                         max_pal=12,
                         random=T,
                         seed=200,
                         into="#EDEDED",
                         add_n=1){
  
  
  out <- 
    tab %>% 
    as.data.frame() %>% 
    dplyr::group_by(!!sym(class), !!sym(subclass)) %>% 
    dplyr::summarise(n=length(!!sym(class)))
  
  class_unique <- unique(tab %>% pull(!!sym(class))) %>% as.character()
  
  if(length(class_unique)>max_pal){
    color <- RColorBrewer::brewer.pal(max_pal,pal)
    color_L1 <- colorRampPalette(color)(length(class_unique))
  }else{
    color_L1 <- RColorBrewer::brewer.pal(length(class_unique),pal)
  }
  
  if(random==T){
    set.seed(seed)
    color_L1 <- sample(color_L1)}
  
  
  out2 <- map_dfr(.x=1:length(class_unique), .f=function(i){
    x <- 
      out %>% 
      dplyr::filter(!!sym(class)==class_unique[i])
    a <- nrow(x)
    x <- 
      x %>% 
      dplyr::mutate(colors=c(colorRampPalette(color=c(into,color_L1[i]))(a+add_n)[c(1+add_n):c(a+add_n)] %>% rev()) )
  })
  
  all <- out %>% dplyr::pull(!!sym(subclass)) %>% unique()
  withcolor <- out2 %>% dplyr::pull(!!sym(subclass)) %>% unique()
  inter <- intersect(withcolor, all)
  
  if(length(all[!all %in% inter])>0){
    out3 <-data.frame(a="NaN", b=all[!all %in% inter], n=1, colors=into)
    names(out3) <- names(out2)
    out_4 <- rbind(out2,out3) %>% dplyr::ungroup()
  }else{
    out_4 <- out2
  }
  
  
  return(out_4)
  
}


