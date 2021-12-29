

# Helper Functions --------------------------------------------------------



#' @title  hlpr_merge.pos
#' @author Dieter Henrik Heiland
#' @description hlpr_merge.pos
#' @inherit 
#' @param Chr Chromosome numeric
#' @param start.pos Start position Chromosom
#' @param end.pos End position Chromosom
#' @param ref.bin Output of the Create.ref.bins() functions
#' @return 
#' @examples 
#' 
#' @export
#' 

hlpr_merge.pos <- function(chr, start.pos, end.pos, ref.bin){
  
  bins <- purrr::map_chr(.x=1:length(chr), .f=function(i){
    out <- 
      ref.bin %>% 
      dplyr::filter(Chr == {chr[i]}) %>% 
      dplyr::filter(start<=start.pos[i]) %>% 
      tail(1) %>% 
      dplyr::pull(bin)
    if(is.null(out)){out <- "NA"}
    
    return(out)
  })
  
  return(bins)
}


#' @title  hlpr_bins
#' @author Dieter Henrik Heiland
#' @description hlpr_bins
#' @inherit 
#' @param chr Chromosome numeric
#' @param start Start position Chromosom
#' @param end End position Chromosom
#' @param bin.size Size of Chromosomal bins
#' @return 
#' @examples 
#' 
#' @export
#' 

hlpr_bins <- function(Chr, start, end, bin.size){
  #calculate the size
  nr.bins <- c(c(end/bin.size) %>% round(digits = 0)+1)
  seq <- seq(from=start, to=c(bin.size*(nr.bins-1)), length.out = nr.bins)
  return(data.frame(bin=paste0(Chr, "_Bin_", 1:length(seq)), start=seq, end=c(seq[2:length(seq)], end) ))
  
}

#' @title  Create.ref.bins
#' @author Dieter Henrik Heiland
#' @description Create.ref.bins is a function to setup a reference data.frame for chromosomal bins
#' @inherit 
#' @param regions  Input from the SPATA2 data.frame SPATA2::cnv_regions_df
#' @param bin.size Size of Chromosomal bins
#' @return 
#' @examples 
#' 
#' @export
#'
Create.ref.bins <- function(object, bin.size){
  
  #Input
  regions <- SPATA2::cnv_regions_df
  sample <- SPATA2::getSampleNames(object)
  
  #Create Bins per chromosome
  ref.bins <- map_df(.x=1:nrow(regions), .f=function(i){
    
    dat <-SPATAwrappers::hlpr_bins(
      Chr=rownames(regions)[i], 
      start=regions$Start[i], 
      end=regions$End[i], 
      bin.size=bin.size)
    
    dat$Chr=regions$Chrom[i]
    dat$Chr.arm=rownames(regions)[i]
    
    return(dat)
    
  })
  
  object@cnv[[sample]]$ref.bins <- ref.bins
  
  return(object)
  
}

#' @title  Create.ref.bins
#' @author Dieter Henrik Heiland
#' @description Create.ref.bins is a function to setup a reference data.frame for chromosomal bins
#' @inherit 
#' @param regions  Input from the SPATA2 data.frame SPATA2::cnv_regions_df
#' @param bin.size Size of Chromosomal bins
#' @return 
#' @examples 
#' 
#' @export
#'


# Map the input InferCNV output to referens (Bins) ------------------------

#' @title  runCNV.Coverage
#' @author Dieter Henrik Heiland
#' @description runCNV.Coverage
#' @inherit 
#' @param object SPATA2 object
#' @return 
#' @examples 
#' 
#' @export
#' 
runCNV.Coverage <- function(object){
  
  #Input
  sample <- SPATA2::getSampleNames(object)
  ref.bin <- SPATA2::getCnvResults(object)[["ref.bins"]]
  
  
  #Get gene_pos_df
  
  cnv.genes <- 
    SPATA2::getCnvResults(object)[["gene_pos_df"]] %>% 
    dplyr::filter(chromosome_name %in% unique(ref.bin$Chr)) %>% 
    dplyr::mutate(bins= SPATAwrappers::hlpr_merge.pos(chr = .$chromosome_name, 
                                                      start.pos = .$start_position,
                                                      end.pos = .$end_position,
                                                      ref.bin=ref.bin))
  
  
  object@cnv[[sample]]$gene_pos_df <- cnv.genes
  
  
  #Merge Bins and count Coverage
  coverage <- 
    cnv.genes %>% 
    dplyr::count(bins) %>% 
    dplyr::rename("bin":=bins)
  
  ref.bin.cov <- 
    ref.bin %>% 
    dplyr::left_join(., coverage, by="bin")
  
  ref.bin.cov$n[is.na(ref.bin.cov$n)]=0
  ref.bin.cov$xaxsis <- 1:nrow(ref.bin.cov)
  
  ref.bin.cov$Arm <- "q"
  ref.bin.cov[ref.bin.cov$Chr.arm %>% str_detect(., patter="p"), ]$Arm <- "p"
  
  object@cnv[[sample]]$ref.bin.cov <- ref.bin.cov
  
  return(object)
  
}


#' @title  plotCoverage
#' @author Dieter Henrik Heiland
#' @description plotCoverage
#' @inherit 
#' @param object SPATA2 object
#' @return 
#' @examples 
#' 
#' @export
#'
plotCoverage <- function(object){
  
  cnv.genes <- SPATA2::getCnvResults(object)[["ref.bin.cov"]]

  #Plot Coverage
  line_df <-
    dplyr::count(x = ref.bin.cov, Chr) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      line_pos = base::cumsum(x = n),
      line_lag = dplyr::lag(x = line_pos, default = 0) ,
      label_breaks = (line_lag + line_pos) / 2
    )  %>%
    tidyr::drop_na()
  
  
  p=ggplot2::ggplot(ref.bin.cov, mapping=aes(x=xaxsis, y=n, color=Arm))+
    ggplot2::geom_point(size=3)+
    ggplot2::theme_classic()+
    ggplot2::geom_vline(data = line_df, mapping=aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5)
  
  return(p)
}


#' @title  runCNVNormalization
#' @author Dieter Henrik Heiland
#' @description runCNVNormalization
#' @inherit 
#' @param object SPATA2 object
#' @return 
#' @examples 
#' 
#' @export
#'
runCNV.Normalization <- function(object,
                                window.k=10,
                                noise=0.035,
                                coverage.model=NULL){
  
  
  
  # Get Data ----------------------------------------------------------------
  
  mat <- SPATA2::getCnvResults(object)[["cnv_mtr"]]  
  cnv.genes <- SPATA2::getCnvResults(object)[["gene_pos_df"]]  
  sample <- SPATA2::getSampleNames(object)
  ref.bin <- SPATA2::getCnvResults(object)[["ref.bins"]]  
  
  # Merge sata to referens bins ---------------------------------------------
  
  message(paste0(Sys.time(), " --- ", "Merge data to the reference bins", " ---"))
  
  data.cnv <- 
    reshape2::melt(mat) %>% 
    dplyr::rename("hgnc_symbol":=Var1) %>% 
    dplyr::rename("barcodes":=Var2) %>% 
    dplyr::left_join(., cnv.genes, by="hgnc_symbol") %>% 
    dplyr::filter(!is.na(bins))
  
  
  message(paste0(Sys.time(), " --- ", "Summarize bins", " ---"))
  
  data.cnv <- 
    data.cnv %>% 
    dplyr::group_by(barcodes,bins) %>% 
    dplyr::summarise(CNV=mean(value), CNV.var=var(value))
  
  data.cnv$CNV.var[is.na(data.cnv$CNV.var)]=0
  
  if(!is.null(coverage.model)){
    
    data.cnv <- 
      data.cnv %>% 
      dplyr::mutate(CNV=predict(coverage.model, data.frame(InferCNV.val = CNV ))
      )
    
  }
  
  
  
  message(paste0(Sys.time(), " --- ", "Noise injection", " ---"))
  
  
  CNV.mat <- matrix(NA, nrow=nrow(ref.bin), ncol=ncol(mat))
  colnames(CNV.mat) <- colnames(mat)
  rownames(CNV.mat) <- ref.bin$bin
  
  CNV.mat.melt <- 
    reshape2::melt(CNV.mat) %>% 
    dplyr::rename("bins":=Var1) %>% 
    dplyr::rename("barcodes":=Var2) %>% 
    dplyr::left_join(.,data.cnv, by=c("bins", "barcodes") )
  
  
  CNV.mat.melt.2 <- 
    CNV.mat.melt %>% 
    dplyr::mutate(impute=imputeTS::na_ma(CNV, k = window.k, weighting = "simple"))
  
  
  
  CNV.mat.melt.3 <- 
    CNV.mat.melt.2 %>% 
    dplyr::mutate(impute.var=
                    imputeTS::na_ma(CNV.var, k = window.k, weighting = "simple")
    ) %>% 
    dplyr::mutate(noise=
                    seq(from=1-noise, to=1+noise, length.out=nrow(.))
                  %>% sample()
    ) %>% 
    dplyr::mutate(CNV.out.noise=
                    c(impute+noise+impute.var) %>% 
                    scales::rescale(c(min(data.cnv$CNV),max(data.cnv$CNV)))
    )
  

  
  
  CNV.mat.out <- reshape2::acast(data=CNV.mat.melt.3, bins~barcodes, value.var = "CNV.out.noise")
  
  
  object@cnv[[sample]]$Normalized.bin.matrix <- CNV.mat.out
  
  return(object)
  
}


#' @title  plotCNV.RefMode
#' @author Dieter Henrik Heiland
#' @description plotCNV.RefMode
#' @inherit 
#' @param object SPATA2 object
#' @return 
#' @examples 
#' 
#' @export
#'
plotCNV.RefMode <- function(object, 
                            across=NULL,
                            sub.across=NULL,
                            pt.size=1.5,
                            sample.order=T){
  
  #get data
  
  ref.bin.cov <- SPATA2::getCnvResults(object)[["ref.bin.cov"]]
  mat <- SPATA2::getCnvResults(object)[["Normalized.bin.matrix"]]
  
  line_df <-
    dplyr::count(x = ref.bin.cov, Chr) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      line_pos = base::cumsum(x = n),
      line_lag = dplyr::lag(x = line_pos, default = 0) ,
      label_breaks = (line_lag + line_pos) / 2
    )  %>%
    tidyr::drop_na()
  
  
  
  if(is.null(across)){
    
    #Prepare plot
    #CNV
    ref.bin.cov.plot <- 
      ref.bin.cov %>% 
      dplyr::left_join(.,mat %>% 
                         rowMeans() %>% 
                         as.data.frame() %>% 
                         tibble::rownames_to_column("bin") %>% 
                         dplyr::rename("cnv_mean":=.),
                       by="bin") %>% 
      dplyr::left_join(., apply(mat,1,var) %>% 
                         as.data.frame() %>% 
                         tibble::rownames_to_column("bin") %>% 
                         dplyr::rename("cnv_sd":=.),
                       by="bin") %>% 
      dplyr::mutate(cnv_mean=runif(n=nrow(.), 
                                   min=c(cnv_mean-cnv_sd), 
                                   max=c(cnv_mean+cnv_sd)) )
    
    p=ggplot2::ggplot(ref.bin.cov.plot, mapping=aes(x=xaxsis, y=cnv_mean, color=Arm))+
      ggplot2::geom_point(size=pt.size)+
      ggplot2::theme_classic()+
      ggplot2::geom_vline(data = line_df, mapping=aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5)
    
  }else{
    
    message(paste0(Sys.time(), " --- ", "Plot across: ",across , " ---"))
    
    across.f <- SPATA2::joinWithFeatures(object, features = across)
    
    if(!is.null(sub.across)){across.f <- across.f %>% dplyr::filter(!!sym(across) %in% sub.across)}
    
    across.list <- 
      purrr::map(.x=unique(across.f %>% pull(!!sym(across))),
                 .f=function(f){
                   
                   mat.select <- mat[,across.f %>% filter(!!sym(across)==f) %>% pull(barcodes)]
                   
                   ref.bin.cov.plot <- 
                     ref.bin.cov %>% 
                     dplyr::left_join(.,mat.select %>% 
                                        rowMeans() %>% 
                                        as.data.frame() %>% 
                                        tibble::rownames_to_column("bin") %>% 
                                        dplyr::rename("cnv_mean":=.),
                                      by="bin") %>% 
                     dplyr::left_join(., apply(mat.select,1,var) %>% 
                                        as.data.frame() %>% 
                                        tibble::rownames_to_column("bin") %>% 
                                        dplyr::rename("cnv_sd":=.),
                                      by="bin") %>% 
                     dplyr::mutate(cnv_mean=runif(n=nrow(.), 
                                                  min=c(cnv_mean-cnv_sd), 
                                                  max=c(cnv_mean+cnv_sd)) ) %>% 
                     dplyr::mutate(across=f)
                   
                   return(ref.bin.cov.plot) 
                   
                 })
    
    ref.bin.cov.plot <- as.data.frame(do.call(rbind, across.list))
    
    
    if(sample.order==T){ref.bin.cov.plot <- ref.bin.cov.plot[sample(1:nrow(ref.bin.cov.plot)), ]}
    
    p=ggplot2::ggplot(ref.bin.cov.plot, mapping=aes(x=xaxsis, y=cnv_mean, color=across))+
      ggplot2::geom_point(size=pt.size)+
      ggplot2::theme_classic()+
      ggplot2::geom_vline(data = line_df, mapping=aes(xintercept = line_pos), linetype = "dashed", alpha = 0.5)
    
    
    
    
    
    
  }
  
  
  
  
  
  
  return(p)
  
  
  
}


#' @title  ggLayerGenSet
#' @author Dieter Henrik Heiland
#' @description ggLayerGenSet
#' @inherit 
#' @param object SPATA2 object
#' @return 
#' @examples 
#' 
#' @export
#'
ggLayerGenSet <- function(object, gene_set, y.pram=1, pt.size=1, pt.shape=15, color="black"){
  
  #Data
  genes.df <- SPATA2::getCnvResults(object)[["gene_pos_df"]]
  coords.df <- SPATA2::getCnvResults(object)[["ref.bin.cov"]]
  
  
  bin.select <- 
    genes.df %>% 
    dplyr::filter(hgnc_symbol %in% 
             SPATA2::getGenes(object, of_gene_sets = gene_set) ) %>% 
    dplyr::pull(bins)
  
  ggplot2::geom_point(data=
                        coords.df %>% 
                        dplyr::filter(bin %in% bin.select), 
                      mapping=aes(x=xaxsis, y=y.pram), 
                      color=color, 
                      shape=pt.shape, 
                      size=pt.size)
  
}

#' @title  getCNVStatistics
#' @author Dieter Henrik Heiland
#' @description getCNVStatistics
#' @inherit 
#' @param object SPATA2 object
#' @return 
#' @examples 
#' 
#' @export
#'
getCNVStatistics <- function(object, 
                             groups,
                             across=NULL,
                             sub.across=NULL){
  ref.bin.cov <- SPATA2::getCnvResults(object)[["ref.bin.cov"]]
  mat <- SPATA2::getCnvResults(object)[["Normalized.bin.matrix"]]
  
    message(paste0(Sys.time(), " --- ", "Plot across: ",across , " ---"))
    
    across.f <- SPATA2::joinWithFeatures(object, features = across)
    
    if(!is.null(sub.across)){across.f <- across.f %>% dplyr::filter(!!sym(across) %in% sub.across)}
    
    across.list <- 
      purrr::map(.x=unique(across.f %>% pull(!!sym(across))),
                 .f=function(f){
                   
                   mat.select <- mat[,across.f %>% filter(!!sym(across)==f) %>% pull(barcodes)]
                   
                   ref.bin.cov.plot <- 
                     ref.bin.cov %>% 
                     dplyr::left_join(.,mat.select %>% 
                                        rowMeans() %>% 
                                        as.data.frame() %>% 
                                        tibble::rownames_to_column("bin") %>% 
                                        dplyr::rename("cnv_mean":=.),
                                      by="bin") %>% 
                     dplyr::left_join(., apply(mat.select,1,var) %>% 
                                        as.data.frame() %>% 
                                        tibble::rownames_to_column("bin") %>% 
                                        dplyr::rename("cnv_sd":=.),
                                      by="bin") %>% 
                     dplyr::mutate(cnv_mean=runif(n=nrow(.), 
                                                  min=c(cnv_mean-cnv_sd), 
                                                  max=c(cnv_mean+cnv_sd)) ) %>% 
                     dplyr::mutate(across=f)
                   
                   return(ref.bin.cov.plot) 
                   
                 })
    
    ref.bin.cov.plot <- 
      as.data.frame(do.call(rbind, across.list)) %>% 
      dplyr::filter(across %in% c(groups)) %>% 
      dplyr::mutate(across=as.character(across)) %>% 
      dplyr::group_by(Chr.arm) %>% 
      dplyr::do(broom::tidy(stats::t.test(cnv_mean~across, data = .))) %>% 
      dplyr::mutate(p.adj=stats::p.adjust(p.value, n=nrow(.))) %>% 
      dplyr::arrange(p.adj)

  return(ref.bin.cov.plot)
    
}


