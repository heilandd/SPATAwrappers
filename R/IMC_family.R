#' @title  createInternalReferenceDataset
#' @author Dieter Henrik Heiland
#' @description createInternalReferenceDataset
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

createInternalReferenceDataset <- function(object, marker.list, q=.95, split_by="ID.int"){
  

# Select the top cells for the reference dataset --------------------------

  cells.list <- purrr::map(.x=marker.list, function(list){
    
    
    select.cells.pos <- 
      SPATA2::joinWithGenes(object, genes=list[[1]], average_genes = T, verbose = F) %>% 
      dplyr::rename("mean_genes":=5) %>% 
      dplyr::arrange(desc(mean_genes)) %>% 
      dplyr::select(barcodes, mean_genes)
    
    if(!is.null(list[[2]])){
      select.cells.pos <- 
        SPATA2::joinWithGenes(object, genes=list[[2]], average_genes = T, verbose = F) %>% 
        dplyr::rename("neg.select":=5) %>% 
        dplyr::left_join(select.cells.pos, . , by="barcodes") %>% 
        dplyr::mutate(mean_genes=mean_genes-neg.select) %>% 
        dplyr::arrange(desc(mean_genes)) %>% 
        dplyr::select(barcodes, mean_genes)
    }
    out <- 
      select.cells.pos %>% 
      dplyr::filter(mean_genes > stats::quantile(mean_genes, q))
    
    return(out)
    
  })
  names(cells.list) <- names(marker.list)

# Remove double annotation ------------------------------------------------

  
  cells.list <- purrr::map_df(.x=1:length(cells.list), function(i){
    target <- cells.list[[i]] %>% dplyr::mutate(cell.type=names(cells.list)[i])
    
    #vv <- (c(1:length(cells.list))[-c(i)])
    #other <- purrr::map_df(.x=vv, .f=function(ii){cells.list[[ii]] %>% dplyr::mutate(cell.type=names(cells.list)[ii]) })
    
    #merged <- 
    #dplyr::left_join(target, other, by="barcodes") %>% 
    #dplyr::arrange((mean_genes.y)) %>% 
    #dplyr::mutate(sum=mean_genes.x-mean_genes.y) %>% 
    #dplyr::filter(sum>0.1) %>% 
    #dplyr::pull(barcodes)
    
    return(target)

    })
  cleaned <- purrr::map_df(.x=unique(cells.list$barcodes), .f=function(bc){
    out <- 
      cells.list %>% 
      dplyr::filter(barcodes==bc) %>% 
      dplyr::arrange(desc(mean_genes)) %>% 
      head(1)
    
    return(out)
  })
  

# Subset object -----------------------------------------------------------
  seurat <- SPATA2::transformSpataToSeurat(object,assay_name = "RNA", RunPCA=F, FindNeighbors=F, FindClusters=F, RunTSNE=F, RunUMAP=F)
  seurat <- subset(seurat, cells=cleaned$barcodes)

  #Cell type
  seurat@meta.data$cell.type <- cleaned$cell.type

  seurat <- Seurat::SetIdent(seurat, value=split_by) 
  object.list=Seurat::SplitObject(seurat)
  seurat <- SeuratWrappers::RunFastMNN(object.list, features=30)
  
  #seurat <- Seurat::FindVariableFeatures(seurat)
  #seurat <- Seurat::ScaleData(seurat)
  #seurat <- Seurat::RunPCA(seurat)
  #seurat <- Seurat::RunUMAP(seurat, dim=1:20, reduction="pca")
  #Seurat::DimPlot(seurat, group.by = "cell.type")
  seurat <- Seurat::RunUMAP(seurat, dim=1:29, reduction="mnn")
  #Seurat::DimPlot(seurat, group.by = "cell.type")
  #Seurat::DimPlot(seurat, group.by = "cell.type")
  #Seurat::FeaturePlot(seurat, features = "EGFR")

  return(seurat)

}

#' @title  createReferenceDIM
#' @author Dieter Henrik Heiland
#' @description createReferenceDIM
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export

createReferenceDIM <- function(object, marker.list, n=500, booster=0.5, batch=50000){
  
  
  # Select the top cells for the reference dataset --------------------------
  
  cells.list <- purrr::map(.x=1:length(marker.list), function(i){
    
    list <- marker.list[[i]]
    
    select.cells.pos <- 
      SPATA2::joinWithGenes(object, genes=list[[1]], average_genes = T, verbose = F) %>% 
      dplyr::rename("mean_genes":=5) %>% 
      dplyr::arrange(desc(mean_genes)) %>% 
      dplyr::select(barcodes, mean_genes)
    
    if(!is.null(list[[2]])){
      select.cells.pos <- 
        SPATA2::joinWithGenes(object, genes=list[[2]], average_genes = T, verbose = F) %>% 
        dplyr::rename("neg.select":=5) %>% 
        dplyr::left_join(select.cells.pos, . , by="barcodes") %>% 
        dplyr::mutate(mean_genes=mean_genes-neg.select) %>% 
        dplyr::arrange(desc(mean_genes)) %>% 
        dplyr::select(barcodes, mean_genes)
    }
    out <- 
      select.cells.pos
    names(out)[2] <- names(marker.list)[i]
    
    
    
    return(out)
    
  })
  names(cells.list) <- names(marker.list)
  
  df <- cells.list[[1]]
  for(i in 2:length(cells.list)){
    df <- df %>% 
      dplyr::left_join(., cells.list[[i]], by="barcodes")
  }
  
  #clean
  vv <- names(marker.list)
  bc <- purrr::map_df(.x=1:length(vv), .f=function(i){
    df$clean <- rowMeans(df[, c(vv[-c(i)])])
    bc <- 
      df %>% 
      dplyr::mutate(Comp=!!sym(vv[i])-clean) %>% 
      dplyr::top_n(n=n, wt=Comp) %>% 
      dplyr::mutate(cell.type=vv[i], cell.score=!!sym(vv[i])) %>% 
      dplyr::select(barcodes, cell.type,cell.score, Comp)
    return(bc)
  }) 
  df <- df %>% dplyr::left_join(data.frame(barcodes=bc$barcodes), ., by="barcodes")
  #pheatmap::pheatmap(df[2:6] %>% as.matrix() %>% t(), cluster_rows = F, cluster_cols = F)
  
  #booster
  for(i in 2:ncol(df)-1){
    #print(i)
    start <- ((n*i)-n)+1; end <- (n*(i+1))-n
    #print(c(start, end))
    df[start:end, i+1] <- df[start:end, i+1]+booster
  }
  
  umap <- Seurat::RunUMAP(df[2:ncol(df)] %>% as.matrix())
  df.umap=umap@cell.embeddings %>% as.data.frame()
  names(df.umap) <- c("UMAP_1", "UMAP_2")
  df.umap$barcodes <- bc$barcodes
  df.umap$cell.types <- bc$cell.type
  #Clean
  df.umap <- df.umap[!duplicated(df.umap$barcodes), ]
  #df.umap <- NFCN2::getCleaned(df.umap, feat = "UMAP_1", q=.005)
  #df.umap <- NFCN2::getCleaned(df.umap, feat = "UMAP_2", q=.005)
  
  
  #ggplot(data=df.umap, mapping=aes(x=UMAP_1, y=UMAP_2, color=cell.types))+
  #  geom_point()+
  #  theme_classic()
  
  
  seurat <- SPATA2::transformSpataToSeurat(object,assay_name = "RNA", RunPCA=F, FindNeighbors=F, FindClusters=F, RunTSNE=F, RunUMAP=F)
  seurat <- subset(seurat, cells=df.umap$barcodes)
  
  
  #Cell type
  seurat@meta.data$cell.type <- df.umap$cell.type
  seurat <- Seurat::FindVariableFeatures(seurat)
  seurat <- Seurat::ScaleData(seurat)
  seurat <- Seurat::RunPCA(seurat)
  
  # Add Assay DimRed Ref
  Ref <- seurat@reductions$pca
  df.new <- df %>% dplyr::left_join(data.frame(barcodes=df.umap$barcodes), ., by="barcodes")
  df.new <- df.new[!duplicated(df.new$barcodes), ]
  rownames(df.new) <- df.new$barcodes
  df.new$barcodes <- NULL
  names(df.new) <- paste0("REF_", 1:ncol(df.new))
  Ref@cell.embeddings <- df.new %>% as.matrix()
  Ref@key <- "REF_"
  seurat@reductions$ref <- Ref
  
  #UMAP
  seurat <- Seurat::RunUMAP(seurat, dim=1:c(seurat@reductions$ref %>% length()), reduction="ref", return.model=T, reduction.name = "ref.umap")
  #Seurat::DimPlot(seurat, group.by = "cell.type", reduction = "umap")
  
  
  
  #Integrate to reference
  
  all.data <- SPATA2::transformSpataToSeurat(object,assay_name = "RNA", RunPCA=F, FindNeighbors=F, FindClusters=F, RunTSNE=F, RunUMAP=F)
  all.data <- subset(all.data, cells=sample(rownames(all.data@meta.data), size=batch))
  all.data <- Seurat::FindVariableFeatures(all.data)
  all.data <- Seurat::ScaleData(all.data)
  all.data <- Seurat::RunPCA(all.data) 
  all.data <- Seurat::RunUMAP(all.data, dims=1:10) 
  
  anchors <- Seurat::FindTransferAnchors(
    reference = seurat,
    query = all.data ,
    reduction = "pcaproject",
    reference.assay="RNA",
    normalization.method = "LogNormalize",
    reference.reduction = "ref",
    dims = 1:c(seurat@reductions$ref %>% length())
  )
  
  all.data <- Seurat::MapQuery(
    anchorset = anchors,
    query = all.data,
    refdata = list(
      cell.type = "cell.type"
    ),
    reference = seurat,
    reference.reduction = "ref", 
    reduction.model = "ref.umap"
  )
  
  
  #Seurat::DimPlot(all.data, reduction = "ref.umap", group.by = "predicted.cell.type")
  
  return(all.data)
  
}
