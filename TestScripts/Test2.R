#IMC test



object <- readRDS("integrated.object.RDS")

SPATA2::getGenes(object)



TAM <- list(c("IBA1","CD68","CD163", "C3", "HLADR", "SPP1"),
            c("EGFR","NeuN", "OLIG1"))
Tcells <- list(c("CD3", "CD4", "FOXP3", "CD8A"),
            c("CD68", "IBA1", "CD11B", "EGFR", "OLIG1"))
Tumor <- list(c("EGFR", "GFAP", "CD44", "HOPX", "CD24"),
              c("CD4", "CD8A", "CD3", "NeuN"))
Neuron <- list(c("NeuN", "THY1", "CALM2"),
               c("EGFR", "OLIG1"))
OD <- list(c("OLIG1"),
           c("EGFR","THY1", "CALM2"))

marker.list <- list(TAM, Tcells, Tumor, Neuron, OD)
names(marker.list) <- c("TAM", "Tcells", "Tumor", "Neuron", "OD")


object.ref <- SPATAwrappers::createReferenceDIM(object, marker.list, booster = 0.3)

object.ref <- Seurat::FindNeighbors(object.ref, reduction = "ref.umap", dims = 1:2)
object.ref <- Seurat::FindClusters(object.ref, reduction = "ref.umap", dims = 1:2, resolution=0.001)

object.ref %>% Seurat::DimPlot(reduction = "ref.umap")
object.ref %>% Seurat::DimPlot(reduction = "ref.umap", split.by = "Region")

Seurat::FeaturePlot(object.ref, reduction = "ref.umap", features = "EGFR", order=T)+scale_color_viridis_c(option="C", limit=c(4,5.5), oob = scales::squish)


NFCN2::Compare_Barplot(object.ref, "seurat_clusters", "Region")


saveRDS(object.ref, "object.ref.imc.RDS")









