#' @title  importSNVmatrix
#' @author Dieter Henrik Heiland
#' @description importSNVmatrix
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 

importSNVmatrix <- function(spata.obj,path.vcf, path.snv, count.filter=1){
  
  variants <- vcfR::read.vcfR(path.vcf, verbose = FALSE)@fix %>% as.data.frame()
  variants.string <- data.frame(SNP =paste0("Chr_",variants$CHROM,"_", variants$POS,":", variants$REF, ">", variants$ALT),
                               Pos= paste0(variants$CHROM,":",variants$POS),
                               Mut=paste0(variants$REF, "<", variants$ALT))
  

  
  
  
  
  message(paste0(Sys.time(), "----", " Import csv ----"))
  counts.alt <- read.csv(path.snv)
  message(paste0(Sys.time(), "----", " Reshape and Filter ----"))
  
  filter.alt <- rowSums(counts.alt[,-1]) != 0
  counts.alt <- counts.alt[filter.alt, ]

  re.counts.alt <- 
    reshape2::melt(counts.alt) %>% 
    dplyr::rename("Count":=value) %>% 
    dplyr::rename("Pos":=SNV) %>% 
    dplyr::rename("barcodes":=variable) %>% 
    dplyr::filter(Count>=count.filter) %>% 
    dplyr::left_join(.,variants.string, by="Pos" ) %>% 
    dplyr::mutate(barcodes=str_replace_all(as.character(barcodes), pattern="[.]", "-"))
  
  
  
  
  message(paste0(Sys.time(), "----", " Add to SPATA object ----"))
  sample <- SPATA2::getSampleNames(spata.obj)
  spata.obj@data[[sample]]$SNV <- list()
  spata.obj@data[[sample]]$SNV$SNV.filtered <- re.counts.alt
  
  return(spata.obj)
}


#' @title  getSNV
#' @author Dieter Henrik Heiland
#' @description getSNV
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 

getSNV <- function(spata.obj){
  sample <- SPATA2::getSampleNames(spata.obj)
  if(is.list(spata.obj@data[[sample]]$SNV) ){return(spata.obj@data[[sample]]$SNV$SNV.filtered)}else(message("No SNV list in object"))
}

#' @title  getSNV
#' @author Dieter Henrik Heiland
#' @description getSNV
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 

getSNV.spotwise <- function(spata.obj){
  sample <- SPATA2::getSampleNames(spata.obj)
  if(is.list(spata.obj@data[[sample]]$SNV) ){return(spata.obj@data[[sample]]$SNV$SNV.spotwise)}else(message("No SNV list in object"))
}

#' @title  runSpotwiseSNV
#' @author Dieter Henrik Heiland
#' @description getSNV
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#'
runSpotwiseSNV <- function(spata.obj, as.features=T){
  #Get Mutations per spot
  
  re.counts.alt <- SPATAwrappers::getSNV(spata.obj)
  
  message(paste0(Sys.time(), "----", " Run and Filter Spotwise SNV Profiles  ----"))
  
  mutations <- re.counts.alt$Mut %>% unique()
  
  #save the mutations in SPATA object
  sample <- SPATA2::getSampleNames(spata.obj)
  spata.obj@data[[sample]]$SNV$Conversion.type <- mutations
  
  mut.spot <- 
    re.counts.alt %>% 
    tidyr::pivot_wider(names_from = Mut, values_from = Count, values_fill = 0) %>% 
    dplyr::group_by(barcodes) %>% 
    dplyr::summarise_at(.vars={{mutations}},.funs=sum) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(snp.var= map_dbl(.x=barcodes, .f=function(bc){
        {re.counts.alt} %>% 
        dplyr::filter(barcodes %in% bc) %>% 
        dplyr::pull(SNP) %>% 
        unique() %>% 
        length() })) %>% 
    dplyr::left_join(.,SPATA2::getFeatureDf(spata.obj) %>% 
                       dplyr::select(barcodes, nCount_Spatial), by="barcodes") %>% 
    as.data.frame() %>% 
    dplyr::mutate(Total.snp=base::rowSums(.[,mutations])) %>% 
    dplyr::mutate(var.score=c((Total.snp)/(snp.var+1)))
  
  
  
  message(paste0(Sys.time(), "----", " Export Results  ----"))
  sample <- SPATA2::getSampleNames(spata.obj)
  spata.obj@data[[sample]]$SNV$SNV.spotwise <- mut.spot
  
  
  if(as.features==T){
    spata.obj <- spata.obj %>% SPATA2::addFeatures(., mut.spot, feature_names = c(mutations, "Total.snp", "snp.var", "var.score"), overwrite = T)
    spata.obj@fdata[[sample]][is.na(spata.obj@fdata[[sample]])]=0
  }
  
  
  return(spata.obj)
  
  
}


#' @title  plotSNVtype
#' @author Dieter Henrik Heiland
#' @param type Can be "fill" or "stack
#' @description plotSNVtype
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#'

plotSNVtype <- function(spata.obj, feature, type="stack", order=NULL){
  sample <- SPATA2::getSampleNames(spata.obj)
  conversions <- spata.obj@data[[sample]]$SNV$Conversion.type
  SNV <- SPATAwrappers::getSNV(spata.obj)
  #SNV[,conversions] <- SNV[,conversions]*SNV$var.score
  factor <- SPATAwrappers::getSNV.spotwise(spata.obj) %>% dplyr::select(barcodes, var.score)
  
  SNV <- 
    SNV %>% 
    dplyr::left_join(., factor %>% dplyr::select(barcodes, var.score), by="barcodes") %>% 
    dplyr::mutate(Count=Count*var.score)
  
  
  if(is.null(order)){
  #ggplot
  feature.df <- joinWith(spata.obj, feature=feature)
  SNV <- 
      SNV %>% 
      dplyr::left_join(., 
                       feature.df %>% dplyr::select(barcodes,{feature}), 
                       by="barcodes")
    
  p=ggplot(SNV, aes(fill=Mut, y=Count, x=!!sym(feature))) + 
    geom_bar(position=type, stat="identity")+theme_classic()+
    scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(13))
  }else{
    snv.order <- 
      data.frame(barcodes=order) %>% 
      tidyr::separate(. , barcodes, into=c("barcodes", "X"), sep="[.]") %>% 
      dplyr::left_join(.,SNV %>% dplyr::select(barcodes, Mut), by="barcodes")
    snv.order[which(snv.order$X!="<NA>"), "barcodes"] <- paste0(snv.order[which(snv.order$X!="<NA>"), "barcodes"], ".1")
    snv.order <- 
      snv.order %>% 
      dplyr::count(barcodes,Mut)
    
    #Bring it to right order
    snv.order$barcodes %>% duplicated()
    snv.order$barcodes <- factor(snv.order$barcodes, levels = order )
    
    
    
    p=ggplot(snv.order, aes(fill=Mut, y=n, x=barcodes)) + 
      geom_bar(position=type, stat="identity")+theme_classic()+
      scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(13))
    
    
    
  }
  return(p)
  }












