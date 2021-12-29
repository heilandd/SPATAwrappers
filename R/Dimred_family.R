
#' @title  runNMF
#' @author Dieter Henrik Heiland adopted from STutility
#' @description runNMF adopted from STutility
#' @param object Seurat object
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export


runNMF <- function (
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  nfactors = 20,
  rescale = TRUE,
  reduction.name = "NMF",
  reduction.key = "factor_",
  n.cores = NULL,
  order.by.spcor = FALSE,
  sort.spcor.by.var = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  var.genes <- features %||% VariableFeatures(object)
  norm.counts <- GetAssayData(object, slot = slot, assay = assay)
  if (rescale) {
    norm.counts <- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }
  if (min(norm.counts) < 0) stop("Negative values are not allowed")
  nmf.results <- rnmf(A = norm.counts[var.genes, ], k = nfactors)
  #nmf.results$W <- swne::ProjectFeatures(norm.counts, nmf.results$H, n.cores = n.cores)
  feature.loadings <- nmf.results$W
  cell.embeddings <- t(nmf.results$H)
  
  # Set cores
  n.cores <- n.cores %||% {detectCores() - 1}
  
  # Order factors based on spatial correlation
  if (order.by.spcor) {
    CN <- do.call(rbind, GetSpatNet(object = object, nNeighbours = NULL, maxdist = NULL))
    resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
    resCN[resCN > 0] <- 1
    empty.CN <- matrix(0, nrow = nrow(cell.embeddings), ncol = nrow(cell.embeddings), dimnames = list(rownames(cell.embeddings), rownames(cell.embeddings)))
    colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
    colnames(resCN) <- gsub(pattern = "^X", replacement = "", x = colnames(resCN))
    empty.CN[rownames(resCN), colnames(resCN)] <- resCN
    listw <- mat2listw(empty.CN)
    fun <- function (x) lag.listw(listw, x, TRUE)
    
    # Calculate the lag matrix from the network
    tablag <- apply(cell.embeddings, 2, fun)
    
    # Split sp.cor by sample
    if (sort.spcor.by.var) {
      sp.cor.split <- do.call(rbind, lapply(unique(GetStaffli(object)@meta.data$sample), function(s) {
        tablag.split <- tablag[GetStaffli(object)@meta.data$sample == s, ]
        cell.embeddings.split <- cell.embeddings[GetStaffli(object)@meta.data$sample == s, ]
        unlist(lapply(1:ncol(cell.embeddings.split), function(i) {
          cor(tablag.split[, i], cell.embeddings.split[, i])
        }))
      }))
      order.vec <- order(apply(sp.cor.split, 2, var))
    } else {
      sp.cor <- unlist(lapply(1:ncol(cell.embeddings), function(i) {
        cor(cell.embeddings[, i], tablag[, i])
      }))
      order.vec <- order(sp.cor, decreasing = TRUE)
    }
    
    cell.embeddings <- cell.embeddings[, order.vec]
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
  }
  
  rownames(x = feature.loadings) <- var.genes
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nfactors)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject (
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}



#' @title  rnmf
#' @author Dieter Henrik Heiland adopted from STutility
#' @description rnmf
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#' 
rnmf <- function (
  A,
  k,
  alpha = 0,
  init = "ica",
  n.cores = 1,
  loss = "mse",
  max.iter = 500,
  ica.fast = F
) {
  if (any(A < 0))
    stop("The input matrix contains negative elements !")
  if (k < 3)
    stop("k must be greater than or equal to 3 to create a viable SWNE plot")
  if (!init %in% c("ica", "nnsvd", "random")) {
    stop("Invalid initialization method")
  }
  A <- as.matrix(A)
  if (any(A < 0)) {
    stop("Input matrix has negative values")
  }
  if (init == "ica") {
    nmf.init <- ica_init(A, k, ica.fast = ica.fast)
  }
  else if (init == "nnsvd") {
    nmf.init <- nnsvd_init(A, k, LINPACK = T)
  }
  else {
    nmf.init <- NULL
  }
  if (is.null(nmf.init)) {
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, n.threads = n.cores,
                          loss = loss, max.iter = max.iter)
  }
  else {
    A.mean <- mean(A)
    zero.eps <- 1e-06
    nmf.init$W[nmf.init$W < zero.eps] <- 0
    nmf.init$H[nmf.init$H < zero.eps] <- 0
    zero.idx.w <- which(nmf.init$W == 0)
    zero.idx.h <- which(nmf.init$H == 0)
    nmf.init$W[zero.idx.w] <- runif(length(zero.idx.w), 0,
                                    A.mean/100)
    nmf.init$H[zero.idx.h] <- runif(length(zero.idx.h), 0,
                                    A.mean/100)
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, init = nmf.init,
                          n.threads = n.cores, loss = loss, max.iter = max.iter)
  }
  colnames(nmf.res$W) <- rownames(nmf.res$H) <- sapply(1:ncol(nmf.res$W),
                                                       function(i) paste("factor", i, sep = "_"))
  return(nmf.res)
}


#' @title  intern
#' @author Dieter Henrik Heiland adopted from STutility
#' @description intern
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export
#'
#'
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

