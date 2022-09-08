#' @title  run.fscore
#' @author Dieter Henrik Heiland
#' @description run.fscore
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
runFscore <- function(df){
  
  res.cluster <- names(df)
  
  mat <- matrix(0,length(res.cluster), length(res.cluster));colnames(mat)=rownames(mat) <- res.cluster
  mat <- mat %>% reshape2::melt()
  purrr::map(.x=1:nrow(mat), .f=function(x){
    mat[x,3]<<-FlowSOM::FMeasure(as.numeric(df[,res.cluster[mat[x,1]]]), as.numeric(df[,res.cluster[mat[x,2]]]))
  })
  mat.m <- mat %>% reshape2::dcast(Var1~Var2, value = "value")
  p <- corrplot::corrplot(mat.m[,-c(1)] %>% as.matrix(), method="circle")
  
  names(mat) <- c("Cluster.M_1","Cluster.M_2", "F-score" )
  
  return(p)
  
}
#' @title  run.fscore.mat
#' @author Dieter Henrik Heiland
#' @description run.fscore.mat
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
runFscore.mat <- function(df){
  
  res.cluster <- names(df)
  
  mat <- matrix(0,length(res.cluster), length(res.cluster));colnames(mat)=rownames(mat) <- res.cluster
  mat <- mat %>% reshape2::melt()
  purrr::map(.x=1:nrow(mat), .f=function(x){
    mat[x,3]<<-FlowSOM::FMeasure(as.numeric(df[,res.cluster[mat[x,1]]]), as.numeric(df[,res.cluster[mat[x,2]]]))
  })
  mat.m <- mat %>% reshape2::dcast(Var1~Var2, value = "value")
  names(mat) <- c("Cluster.M_1","Cluster.M_2", "F-score" )
  
  return(mat)
  
}
#' @title  run.full.cluster.validation (for Spatial Paper Beta Version!!!!!)
#' @author Dieter Henrik Heiland
#' @description run.full.cluster.validation
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
runFullclustervalidation <- function(seuratObj, plot.folder){
  
  require(Seurat)
  
  #Subset Tumor
  tumor <- seuratObj@meta.data %>% filter(cell.type=="Tumor") %>% rownames()
  seuratObj <- subset(seuratObj, cells=tumor)
  
  seuratObj <- 
    seuratObj %>% 
    SCTransform(assay="Spatial") %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:10)
  
  
  seuratObj <- seuratObj %>% run.SNN.stability(assay="SCT", 
                                               reduction="pca",
                                               dims = 1:30, 
                                               resolution = seq(from=0.1, to=1.5, by=0.05), 
                                               cluster_id="Louvain", algorithm=1)
  
  seuratObj <- seuratObj %>% run.SNN.stability(assay="SCT", reduction="pca",
                                               resolution = seq(from=0.1, to=1.5, by=0.05), 
                                               dims = 1:30, cluster_id="SLM",algorithm=3)
  
  if(dir.exists(plot.folder)==F){dir.create(plot.folder)}
  
  
  
  ggsave(Seurat::DimPlot(seuratObj, group.by = "Louvain"),
         filename=paste0(plot.folder, "/DimPlot_cluster_State1_Louvain.png"))
  
  ggsave(Seurat::DimPlot(seuratObj, group.by = "SLM"),
         filename=paste0(plot.folder, "/DimPlot_cluster_State1_SLM.png"))
  
  
  
  ##Add data to matrix
  #Create DF for all clusters
  
  # PAM Cluster
  
  k=length(unique(seuratObj@meta.data$Louvain))
  
  pca <- seuratObj %>% Seurat::Reductions("pca")
  pca <- pca@cell.embeddings %>% as.data.frame()
  pca <- pca[,1:30]
  pam_PCA <- cluster::pam(pca, k=k)
  
  pca <- seuratObj %>% Seurat::Reductions("pca")
  pca <- pca@cell.embeddings %>% as.data.frame()
  pca <- pca[,1:30]
  HC_PCA <- hclust(dist(pca))
  HC_PCA <- cutree(HC_PCA, k=k)
  
  
  
  cluster_summary=data.frame(barcodes=seuratObj@meta.data %>% rownames(), 
                             PCA_1=c(as.numeric(seuratObj@meta.data$Louvain)+1),
                             PCA_2=c(as.numeric(seuratObj@meta.data$SLM)+1),
                             PCA_4=pam_PCA$clustering %>% as.numeric(),
                             PCA_5=HC_PCA %>% as.numeric())
  
  write.csv(cluster_summary, file=paste0(plot.folder, "/cluster_summary.csv"))
  
  #ggsave(run.fscore(cluster_summary[-c(1)]),
  #       filename=paste0(plot.folder, "/Fscore.png"))
  
  write.csv(run.fscore.mat(cluster_summary[-c(1)]), file=paste0(plot.folder, "/Fscore.csv"))
  
  
  prefix = "PCA_"
  out <- clustree::clustree(cluster_summary, 
                            prefix = prefix)
  
  out.plot <- out$data %>% 
    select(!!sym(prefix), sc3_stability) %>% 
    group_by(!!sym(prefix)) %>% 
    summarise(mean(sc3_stability), sd(sc3_stability)) %>% 
    as.data.frame()
  
  
  ## Validate stabolity
  best.cluster <- which.max(out.plot$`mean(sc3_stability)`)
  
  cluster_summary_UMAP=data.frame(barcodes=seuratObj@meta.data %>% rownames(), 
                                  PCA_1=c(as.numeric(seuratObj@meta.data$Louvain)+1),
                                  PCA_2=c(as.numeric(seuratObj@meta.data$SLM)+1),
                                  PCA_3=pam_PCA$clustering %>% as.numeric(),
                                  UMAP1=seuratObj@reductions$umap@cell.embeddings[,1],
                                  UMAP2=seuratObj@reductions$umap@cell.embeddings[,2])
  
  ggsave(clustree::clustree_overlay(cluster_summary_UMAP, 
                                    prefix = prefix, 
                                    x_value = "UMAP1", 
                                    y_value = "UMAP2"),
         filename=paste0(plot.folder, "/ClusterTree_UMAP.png"))
  
  
  
  
  
  # Merge Clusters together
  stable <- 
    cluster_summary_UMAP %>% 
    mutate(stable = PCA_1+PCA_2+PCA_3 ) %>% 
    count(stable) %>% 
    filter(n>200) %>% 
    pull(stable)
  
  # Select stable clusters
  cluster_summary_UMAP <- 
    cluster_summary_UMAP %>% 
    mutate(stable = PCA_1+PCA_2+PCA_3 ) %>% 
    filter(stable %in% {{stable}})
  
  plot <- clustree::clustree_overlay(cluster_summary_UMAP, 
                                     prefix = prefix, 
                                     x_value = "UMAP1", 
                                     y_value = "UMAP2")
  
  ggsave(plot,
         filename=paste0(plot.folder, "/ClusterTree_UMAP_stable.png"))
  
  
  
  for(i in 1:length(unique(cluster_summary_UMAP$stable))){cluster_summary_UMAP$stable[cluster_summary_UMAP$stable==unique(cluster_summary_UMAP$stable)[i]] <- i }
  
  #ggplot(data=cluster_summary_UMAP, aes(x=UMAP1, y=UMAP2, color=as.factor(stable)))+geom_point()+theme_classic()
  
  ## Add consensus clusers to seurat
  
  seuratObj <- subset(seuratObj, cells=cluster_summary_UMAP$barcodes)
  seuratObj <- 
    seuratObj %>% 
    SCTransform(assay="Spatial") %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:10)
  
  seuratObj@meta.data <- 
    seuratObj@meta.data[cluster_summary_UMAP$barcodes, ] %>% 
    mutate(consensus=cluster_summary_UMAP$stable)
  
  ggsave(Seurat::DimPlot(seuratObj, reduction = "pca", group.by = "consensus" ),
         filename=paste0(plot.folder, "/consensus_PCA.png"))
  ggsave(Seurat::DimPlot(seuratObj, reduction = "umap", group.by = "consensus" ),
         filename=paste0(plot.folder, "/consensus_UMAP.png"))
  
  
  seuratObj <- SetIdent(seuratObj, value=seuratObj@meta.data$consensus)
  
  markers <- FindAllMarkers(seuratObj)
  
  final.clusters <- markers %>% count(cluster) %>% filter(n>50) %>% pull(cluster)
  markers.final <- markers %>% filter(cluster %in% final.clusters)
  
  # Remove duplicated genes 
  markers.final$rows <- 1:nrow(markers.final)
  dup <- markers.final[duplicated(markers.final$gene), ]$gene
  remove <- purrr::map(.x=dup, .f=function(i){
    keep <- markers.final %>% 
      filter(gene=={{i}}) %>% 
      arrange(desc(avg_logFC)) %>% 
      head(1) %>% 
      pull(rows)
    
    remove <- markers.final %>% 
      filter(gene=={{i}}) %>% 
      filter(rows!={{keep}}) %>% 
      pull(rows)
    
    return(remove)
    
  }) %>% unlist()
  markers.final <- markers.final[-c(remove), ]
  final.clusters <- markers.final %>% count(cluster) %>% filter(n>20) %>% pull(cluster)
  markers.final <- markers.final %>% filter(cluster %in% final.clusters)
  rownames(markers.final) <- markers.final$gene
  
  
  ## Add new cluster to seurat
  cluster <- as.numeric(unique(markers.final$cluster))
  bc <- 
    seuratObj@meta.data[seuratObj@meta.data$consensus %in% cluster, ] %>% rownames()
  
  seuratObj <- subset(seuratObj, cells=bc)
  seuratObj <- 
    seuratObj %>% 
    SCTransform(assay="SCT") %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:10)
  
  ggsave(Seurat::DimPlot(seuratObj, reduction = "pca", group.by = "consensus" ),
         filename=paste0(plot.folder, "/consensus_PCA_final.png"))
  ggsave(Seurat::DimPlot(seuratObj, reduction = "umap", group.by = "consensus" ),
         filename=paste0(plot.folder, "/consensus_UMAP_final.png"))
  
  
  
  return(list(seuratObj,markers.final ))
}


#' @title  run.full.cluster.validation (for Spatial Paper Beta Version!!!!!)
#' @author Dieter Henrik Heiland
#' @description run.full.cluster.validation
#' @inherit 
#' @return 
#' @examples 
#' @export
runPCAprojection <- function(PCA,
                             matrix.input, 
                             matrix.transfer, 
                             Dim=2, 
                             var=400, 
                             epochs=50){
  
  note=function(text){return(print(message(paste0(text))))}
  
  
  scaled_org <- matrix.input
  scales_target <- matrix.transfer
  
  #######
  note(c("Start Analysis using ", Dim, " Dimensions"))
  #######
  
  #Train Model for PCA coordinates using a bidirectional RNN
  
  require(keras)
  require(tensorflow)
  require(tfdatasets)
  require(tidyverse)
  
  
  #dim(scaled_org); dim(scaled_target)
  
  # 1. Set sup the data set 
  inter <- generics::intersect(rownames(scaled_org), base::rownames(scaled_target))
  scaled_org <- scaled_org[inter, ]
  scaled_target <- scaled_target[inter, ]
  
  
  
  
  # 2. Select most var Features
  variance <- 
    base::data.frame(var=apply(scaled_org,1, var)) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::arrange(desc(var)) %>% 
    utils::head(var) %>% 
    dplyr::pull(gene)
  
  # 3. Prepare data 
  x_train <- base::scale(t(scaled_org[variance, ]))
  y_train <- PCA[,1:Dim]
  dim(y_train); dim(x_train)
  
  #plot(x_train[,2])
  
  x_test <- t(scaled_target[variance, ])
  dim(x_test)
  
  # 4. Run prediction 
  
  colnames=base::colnames(x_train)
  
  test_df <- x_test %>% 
    tibble::as_tibble(.name_repair = "minimal") %>% 
    stats::setNames(paste0("Gene_", 1:var))
  
  
  Predictions <- purrr::map(.x=1:Dim, .f=function(x){
    
    message(paste0(Sys.time(),"   Start Predict Dimension: ", "PC_", x ))
    
    train_df <- x_train %>% 
      tibble::as_tibble(.name_repair = "minimal") %>% 
      stats::setNames(paste0("Gene_", 1:var)) %>%
      dplyr::mutate(label = y_train[,x])
    
    spec <- feature_spec(train_df, label ~ . ) %>% 
      tfdatasets::step_numeric_column(all_numeric(), normalizer_fn = tfdatasets::scaler_standard()) %>% 
      tfdatasets::fit()
    
    input <- tfdatasets::layer_input_from_dataset(train_df %>% dplyr::select(-label))
    
    output <- input %>% 
      keras::layer_dense_features(tfdatasets::dense_features(spec)) %>% 
      keras::layer_dense(units = 64, activation = "relu") %>%
      keras::layer_dense(units = 32, activation = "relu") %>%
      keras::layer_dense(units = 1) 
    model <- 
      keras::keras_model(input, output) %>% 
      keras::compile(
        loss = "mse",
        optimizer = optimizer_rmsprop(),
        metrics = list("mean_absolute_error")
        )
    
    
    history <- 
      model %>% 
      keras::fit(
        x = train_df %>% dplyr::select(-label),
        y = train_df$label,
        epochs = epochs,
        validation_split = 0.2,
        verbose = 0
    )
    
    test_predictions <- model %>% predict(test_df)
    
    out <- data.frame(PC=test_predictions)
    names(out) <- paste0("PC_",x)
    
    return(out)
    
  })
  
  df_plot <- as.data.frame(do.call(cbind, Predictions))
  
  return(df_plot)
  
  
  
}



#' @title  runSpatialRegression
#' @author Dieter Henrik Heiland
#' @description runSpatialRegression
#' @param object SPATA2 object
#' @param features Character of all factors to be correlated (features/genes/GeneSets)
#' @param model The model to be used: classical:Pearson Correlation; lmSLX, lagsarlm, errorsarlm, CCA or dist
#' @param smooth Logical, if TRUE spatial smooth values
#' @param normalize Logical, if TRUE normalize values
#' @inherit 
#' @return 
#' @examples 
#' @export
runSpatialRegression <- function (object, features, model, smooth = F, normalize = F) {
  color_var <- SPATA2::hlpr_join_with_aes(object, df = SPATA2::getCoordsDf(object), 
                                          color_by = features, normalize = normalize, smooth = smooth)
  coords <- color_var[, c("x", "y")] %>% as.data.frame()
  data <- color_var[, features] %>% as.data.frame()
  rownames(coords) = rownames(data) <- color_var$barcodes
  sp_obj <- sp::SpatialPointsDataFrame(coords = coords, data = data)
  nc <- spdep::knearneigh(sp_obj, k = 6, longlat = T)
  nc.nb <- spdep::knn2nb(nc)
  nb.list <- spdep::nb2listw(nc.nb)
  nb.mat <- spdep::nb2mat(nc.nb)
  runCCADist <- function(object, color_var, a, b) {
    sample <- SPATA2::getSampleNames(object)
    data <- color_var[, c(a, b)] %>% as.data.frame()
    names(data) <- c("a", "b")
    rownames(data) <- color_var$barcodes
    coord_spots <- object %>% SPATA2::getCoordsDf()
    coord_spots <- coord_spots[, c("row", "col")] %>% as.data.frame()
    rownames(coord_spots) <- object %>% SPATA2::getCoordsDf() %>% 
      pull(barcodes)
    coord_spots <- cbind(coord_spots, data[rownames(coord_spots), 
    ])
    coord_spots2 <- coord_spots
    coord_spots2$a <- coord_spots2$b * c(-1)
    getKernel <- function(coord_spots, var) {
      mat <- reshape2::acast(row ~ col, data = coord_spots, 
                             value.var = var)
      mat[is.na(mat)] <- 0
      mat %>% scales::rescale(., c(0, 1)) %>% oce::matrixSmooth()
    }
    real <- proxy::dist(x = getKernel(coord_spots, "a"), 
                        y = getKernel(coord_spots, "b"))
    X <- getKernel(coord_spots, "a")
    Y <- getKernel(coord_spots, "b")
    cca <- CCA::matcor(X, Y)
    return(c(cca$XYcor[!is.na(cca$XYcor)] %>% mean(), mean(na.omit(as.vector(real)))))
  }
  if (model == "classical") {
    mat.cor <- cor(data)
    return(mat.cor)
  }
  if (model == "lmSLX") {
    mat.cor <- matrix(NA, length(features), length(features))
    rownames(mat.cor) = colnames(mat.cor) <- features
    for (i in 1:length(features)) {
      for (j in 1:length(features)) {
        first_feat <- features[i]
        second_feat <- features[j]
        if (first_feat == second_feat) {
          mat.cor[first_feat, second_feat] = 0
        }
        else {
          formula <- as.formula(paste(first_feat, "~", 
                                      second_feat))
          reg1 <- lm(formula, data = sp_obj)
          reg2 = spatialreg::lmSLX(reg1, data = sp_obj, 
                                   nb.list)
          sum_mod_2 <- summary(reg2)
          cor_lag2 <- sum_mod_2$coefficients[3, 1]
          mat.cor[first_feat, second_feat] = as.numeric(cor_lag2)
        }
      }
    }
    return(mat.cor)
  }
  if (model == "lagsarlm") {
    if (length(features) >= 2) 
      stop("Run time to long, reduce number of features")
    first_feat <- features[1]
    second_feat <- paste(features[2:length(features)], collapse = "+")
    formula <- as.formula(paste(first_feat, "~", second_feat))
    reg1 <- lm(formula, data = sp_obj)
    sum_mod_1 <- summary(reg1)
    cor_lag1 <- sum_mod_1$coefficients[2, ]
    reg3 = spatialreg::lagsarlm(reg1, data = sp_obj, nb.list)
    sum_mod_3 <- summary(reg3)
    mat.cor <- data.frame(estimate = as.numeric(sum_mod_3$coefficients[2]), 
                          p.value = as.numeric(sum_mod_3$Wald1$p.value))
  }
  if (model == "errorsarlm") {
    if (length(features) >= 2) 
      stop("Run time to long, reduce number of features")
    first_feat <- features[1]
    second_feat <- paste(features[2:length(features)], collapse = "+")
    formula <- as.formula(paste(first_feat, "~", second_feat))
    reg1 <- lm(formula, data = sp_obj)
    sum_mod_1 <- summary(reg1)
    cor_lag1 <- sum_mod_1$coefficients[2, ]
    reg4 = spatialreg::errorsarlm(reg1, data = sp_obj, nb.list)
    sum_mod_4 <- summary(reg4)
    mat.cor <- data.frame(estimate = as.numeric(sum_mod_4$coefficients[2]), 
                          p.value = as.numeric(sum_mod_4$Wald1$p.value))
  }
  if (model == "CCA") {
    mat.cor <- matrix(NA, length(features), length(features))
    rownames(mat.cor) = colnames(mat.cor) <- features
    for (i in features) {
      for (j in features) {
        if (i == j) {
          mat.cor[i, j] <- 1
        }
        else {
          mat.cor[i, j] <- runCCADist(object, a = i, 
                                      b = j, color_var = color_var)[1]
        }
      }
    }
    return(mat.cor)
  }
  if (model == "dist") {
    mat.cor <- matrix(NA, length(features), length(features))
    rownames(mat.cor) = colnames(mat.cor) <- features
    for (i in features) {
      for (j in features) {
        if (i == j) {
          mat.cor[i, j] <- 0
        }
        else {
          mat.cor[i, j] <- runCCADist(object, a = i, 
                                      b = j, color_var = color_var)[2]
        }
      }
    }
    return(mat.cor)
  }
  return(mat.cor)
}










