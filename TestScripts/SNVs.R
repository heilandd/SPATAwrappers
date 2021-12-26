



# Mt Genotyping Script ----------------------------------------------------

#Load genotype file

path <- "~/Desktop/SpatialTranscriptomics/Visium/Visium/313_T/outs/cellsnp"
setwd(path)
library(R.utils)
library(tidyverse)
library(vcfR)
library("ggtree")
library("dendextend")
library("dplyr")
library("ggplot2")
library("RColorBrewer")



# Import Data And preprocessing -------------------------------------------------------------


setwd(path)
gunzip("cellSNP.cells.vcf.gz", remove=FALSE)
vcf <- read.vcfR("cellSNP.cells.vcf", verbose = FALSE)
bc <- read.csv("cellSNP.samples.tsv", header=F) %>% pull(V1)


#Quality Check and transformation into SNP Matrix

mat <- vcf@gt %>% as.matrix()
dim(mat)
#Get variants
variants <- vcf@fix %>% as.data.frame()
variants.string <- paste0(variants$POS,":", variants$REF, "<", variants$ALT )

# Run Quality fitting
pb <- progress_estimated(ncol(mat)-1)
call.var <- map(.x=2:ncol(mat), .f=function(i){
  pb$tick()$print()
  df <- mat[,i] %>% as.data.frame()
  rownames(df) <- variants.string
  names(df) <- c("V")
  df <- 
    df %>% 
    filter(.!=".:.:.:.:.:.") %>% 
    separate(V, into=c("GT","AD","DP","OTH","PL","ALL"), sep=":") %>% 
    filter(GT != "0/0") %>% 
    separate(PL, into=c("0/0","0/1","1/1"), sep="[,]")
  for(r in c("0/0","0/1","1/1")){df[,r] <- as.numeric(df[,r])}
  
  var <- 
    df %>% 
    mutate(QC=case_when(
      GT=="1/0" & c(df$`0/0`- df$`0/1`) > 10 ~ "ok",
      GT=="1/1" & c(rowMeans(df[,c("0/0","0/1")])-df$`1/1`) > 10 ~"ok",
      TRUE ~ "fail"
    )) %>% 
    filter(DP>1) %>% 
    filter(QC=="ok") %>% 
    rownames()
  return(var)
})
names(call.var) <- bc



# Check frequency of variants

var.freq <- call.var %>% unlist %>% table() %>% as.data.frame()
names(var.freq) <- c("Var", "freq")
var.freq <- var.freq %>% arrange(desc(freq))

var.freq <- 
  var.freq %>% 
  tidyr::separate(Var, into=c("POS", "Mut"), sep=":", remove=F)

unique.var <- unique(var.freq$Mut)

library(patchwork)
wrap_plots(
map(.x=unique.var, function(i){
  print(i)
  ggplot(var.freq %>% filter(Mut==i), aes(x=1:nrow(var.freq %>% filter(Mut==i)), y=freq))+geom_point()+theme_classic()+xlab(i)
})
)



# Clustering --------------------------------------------------------------



#Create a GenoType matrix
mat.geno <- matrix(0, nrow(var.freq), length(call.var))
rownames(mat.geno) <- var.freq$Var
for(i in 1:ncol(mat.geno)){mat.geno[call.var[[i]],i]=1}
colnames(mat.geno) <- bc




#For SNP based spot deconvolution 
#First predict the lineage related sibtypes of SNPs

dim(mat.geno)
filtert.snp <- mat.geno[var.freq$Var[1:1000] %>% as.character(), ]

#measure the distance (russel)
dist.mat <- parallelDist::parallelDist(filtert.snp)
hc <- hclust(dist.mat, method = "ward.D")


# Phylogenic tree Analysis ------------------------------------------------




phylo_tree <- ape::as.phylo(hc)
group <- size <- dendextend::cutree(phylo_tree, 10)
group[] <- LETTERS[1:10][group]
size[] <- sample(size)
group.df <- data.frame(label=names(group), group=group, size=size)
phylo_tree <- dplyr::full_join(phylo_tree, group.df, by='label')

#--- Generate a ggtree with 'daylight' layout
ggt <- ggtree(phylo_tree, layout = 'daylight')
#--- Plot the ggtree
ggt + geom_tippoint(aes(color=group), size=3) + scale_y_reverse()


circular <- ggtree(phylo_tree, layout = 'radial', branch.length = "none")
#--- Plot the ggtree
circular + geom_tippoint(aes(color=group), size=3) + scale_y_reverse()



# Group Spots by heteroeneity ---------------------------------------------

group.snp <- purrr::map(.x=group.df$group %>% unique(), .f=function(g){group.df %>% filter(group==g) %>% pull(label)})
names(group.snp) <- group.df$group %>% unique()


#Connect Spots to group
group.snp <- purrr::map(.x=group.df$group %>% unique(), 
                        .f=function(g){
                          lable <- group.df %>% filter(group==g) %>% pull(label)
                          #get reads
                          filtert.snp.e <- filtert.snp[lable, ]
                          #dim(filtert.snp.e)
                          df.snp <- colSums(filtert.snp.e) %>% as.data.frame()
                          names(df.snp) <- paste0("Group_",g)
                          #df.snp <- df.snp %>% tibble::rownames_to_column("barcodes")
                          #rownames(df.snp)=NULL
                          #print(head(df.snp))
                          #if(g!="A"){df.snp$barcodes==NULL}
                          return(df.snp)
                          })

group.snp <- as.data.frame(do.call(cbind, group.snp))
group.snp <- group.snp %>% tibble::rownames_to_column("barcodes")

# Normalize counts by total counts
t.g <- names(group.snp)[-1]
rs <- rowSums(group.snp[,t.g])

for(i in 1:nrow(group.snp)){
  group.snp[i,t.g] <- rs[i]/(group.snp[i,t.g]+1)
}
group.snp[,t.g] <- scale(group.snp[,t.g])


setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions")
TX <- readRDS("313_T_SPATA_CNV_Pred.RDS")
TX@used_genesets <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/GS_new.RDS")

TX <- TX %>% SPATA2::addFeatures(., group.snp, overwrite = T)
SPATA2::plotSurfaceInteractive(TX)


T.cells <- c("CD8A", "CD4", "CD3D", "CD2")

cor.test <- joinWith(TX, genes = T.cells, features = names(group.snp)[-1], average_genes=T, smooth = T)

cor.test <- cor.test %>% mutate(Group = map(1:nrow(.), function(i){which.max(cor.test[i, t.g]) %>% names()}) %>% unlist())



confuns::plot_statistics_interactive(cor.test)

library(patchwork)
wrap_plots(
  map(.x=names(group.snp)[-1], function(i){
    print(i)
    ggplot(cor.test, aes(x=!!sym(i), y=mean_genes))+geom_point()+theme_classic()+xlab(i)
  })
)


#Function
spatial.plot <- function(object,feature, plot=T){
  
  
  #Run analysis 
  message(paste0(Sys.time(), "  Run MC Simulation for Spatial  Correlations analysis"))
  
  cor.mat <- matrix(NA, length(feature),length(feature));colnames(cor.mat) = rownames(cor.mat) <- feature
  cor.mat <- reshape2::melt(cor.mat)
  
  # fill data 
  cor.mat$value <- pbmcapply::pbmclapply(1:nrow(cor.mat), function(i){
    cor.out <- spatial.mc(P1=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,1]),
                          P2=SPATA2::getFeatureDf(object) %>% pull(cor.mat[i,2]),
                          n=200)
    return(cor.out$Cor)
    
  }, mc.cores = 8) %>% unlist()
  
  message(paste0(Sys.time(), "  Run Auto Corelation analysis"))
  
  cor.Auto <- spatial.ac(object, feature)
  cor.Auto <- as.data.frame(do.call(rbind,cor.Auto)) %>% 
    tibble::rownames_to_column("Var2") %>% 
    mutate(Var1="Autocorelation",
           value=V1) %>%
    dplyr::select(Var1,Var2,value)
  
  cor.mat <- rbind(cor.mat, cor.Auto)
  
  
  message(paste0(Sys.time(), " Plotting "))
  
  if(plot==T){
    corrplot::corrplot.mixed(t(reshape2::acast(cor.mat, Var1~Var2, value.var="value")))
  }
  
  
  return(cor.mat)
  
}
spatial.mc <- function(P1,P2, n = 599){
  # Monte-Carlo Simulation of spatial correlation
  if(length(P2)!=length(P1)) stop("Unequal Inputs")
  message(paste0(Sys.time(), "Start Model"))
  M <- glm(P1 ~ P2)
  #coef(M)[2]
  message(paste0(Sys.time(), "Start MC Simulation"))
  #MC
  I.r <- purrr::map_dbl(.x=1:n, .f=function(i){
    x <- sample(P1, replace=FALSE)
    y <- P2
    # Compute new set of lagged values
    #x.lag <- lag.listw(lw, x)
    # Compute the regression slope and store its value
    M.r    <- glm(y ~ x)
    I.r <- coef(M.r)[2]
    return(I.r)
  })
  
  #hist(I.r, main=NULL, xlab="Spatial-Cor-MC", las=1, xlim=c(-0.5,0.5))
  #abline(v=coef(M)[2], col="red")
  
  #p-value
  message(paste0(Sys.time(), "P-Value Extraction"))
  N.greater <- sum(coef(M)[2] > I.r)
  p <- min(N.greater + 1, n + 1 - N.greater) / (n + 1)
  p
  out <- list(coef(M)[2], p)
  names(out) <- c("Cor", "p")
  return(out)
}
spatial.ac <- function(object,feature){
  
  message(paste0(Sys.time(), "  Start Autocorrelation"))
  
  tissue.pos <- SPATA2::getCoordsDf(object)
  plot.new()
  message(paste0(Sys.time(), "  Transforme in spatial dataset"))
  segments <- pbmcapply::pbmclapply(1:nrow(tissue.pos), function(xx){
    
    segment_c <-plotrix::draw.circle(x=tissue.pos$x[xx], 
                                     y=tissue.pos$y[xx],
                                     radius=10*0.4) 
    
    segment <- as.data.frame(do.call(cbind, segment_c))
    
    segment <- sp::Polygon(cbind(segment$x, segment$y))
    segment <- sp::Polygons(list(segment), tissue.pos$barcodes[xx])
    return(segment)
  }, mc.cores = 8)
  
  SpP = sp::SpatialPolygons(segments, 1:length(segments))
  
  attr <- SPATA2::getFeatureDf(object) %>% as.data.frame()
  rownames(attr) <- attr$barcodes
  attr <- attr[,feature]
  
  SrDf = sp::SpatialPolygonsDataFrame(SpP, attr)
  
  message(paste0(Sys.time(), "  Define neighboring Spots"))
  
  nb <- spdep::poly2nb(SrDf, queen=F, snap=25)
  lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
  
  message(paste0(Sys.time(), "  Computing the Moran’s I statistic"))
  
  stat <- purrr::map(.x=feature, .f=function(x){
    model <- spdep::moran.test(SrDf@data[,x],lw, zero.policy=TRUE)
    model$estimate[[1]]
  })
  names(stat) <- feature
  
  return(stat)
  
  
}




#Spatial Correlation Plot
T.cells <- c("CD8A", "CD4", "CD3D", "CD2")
TX <- 
  joinWith(TX, genes = T.cells, average_genes=T) %>% 
  select(-sample, -x, -y) %>% 
  rename("Tcells":=mean_genes) %>% 
  addFeatures(TX, feature_df=., overwrite = T)

myeloid <- c("CD68", "AIF1", "TMEM119")
TX <- 
    joinWith(TX, genes = myeloid, average_genes=T) %>% 
    select(-sample, -x, -y) %>% 
    rename("myeloid":=mean_genes) %>% 
  addFeatures(TX, feature_df=., overwrite = T)

tumor <- c("EGFR", "CHI3L1")
TX <- 
    joinWith(TX, genes = tumor, average_genes=T) %>% 
    select(-sample, -x, -y) %>% 
    rename("Tumor":=mean_genes) %>% 
  addFeatures(TX, feature_df=.,overwrite = T)


TX@fdata$`313_T` <- TX@fdata$`313_T` %>% mutate(Group = map(1:nrow(.), function(i){which.max(cor.test[i, t.g]) %>% names()}) %>% unlist()) 
TX@fdata$`313_T` <- TX@fdata$`313_T` %>% mutate(Group_var = apply(cor.test[, t.g],1,var) )

data.cor <- joinWith(TX, features = c("pred.tumor", "Group_var"), smooth_span = 0.1)

plot(data.cor$pred.tumor, data.cor$Group_var)


SPATA2::plotStatisticsInteractive(SPATA2::getFeatureDf(TX))


spCor <- spatial.plot(TX, feature = c(names(group.snp)[-1], "Tcells", "myeloid", "Tumor"), plot=F)


spCor %>% 
  filter(Var2 %in% c("Tcells", "myeloid", "Tumor")) %>% 
  #mutate(value=scales::rescale(value, c(-1,1))) %>% 
  reshape2::acast(., Var1~Var2, value.var =  "value") %>% 
  corrplot::corrplot(pch = 5, col=viridis::plasma(50), is.corr=F)






#Add SNP Matrix
TX <- SPATA2::addExpressionMatrix(TX, expr_mtr=filtert.snp, mtr_name="SNP")
TX <- SPATA2::setActiveExpressionMatrix(TX, mtr_name="SNP")

TX <- SPATA2::addFeatures(TX, cluster.var, overwrite = T)

plotSurfaceInteractiveSingleCell(TX)


# Read in the st Data 

spata.obj <- readRDS()







# VCF Data WG -------------------------------------------------------------

setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/mtDNA/Test")
vcf <- read.vcfR("269_T.vcf", verbose = FALSE)


#Quality Check and transformation into SNP Matrix
mat <- vcf@gt %>% as.matrix()

#Get variants
variants <- vcf@fix %>% as.data.frame()
names(variants)
variants.string <- data.frame(SNP =paste0("Chr_",variants$CHROM,"_", variants$POS,":", variants$REF, ">", variants$ALT),
                              Pos= paste0(variants$CHROM,":",variants$POS),
                              Mut=paste0(variants$REF, "<", variants$ALT)
                              )

# Run Quality fitting
pb <- progress_estimated(ncol(mat)-1)
call.var <- map(.x=2:ncol(mat), .f=function(i){
  pb$tick()$print()
  df <- mat[,i] %>% as.data.frame()
  rownames(df) <- variants.string
  names(df) <- c("V")
  df <- 
    df %>% 
    filter(.!=".:.:.:.:.:.") %>% 
    separate(V, into=c("GT","AD","DP","OTH","PL","ALL"), sep=":") %>% 
    filter(GT != "0/0") %>% 
    separate(PL, into=c("0/0","0/1","1/1"), sep="[,]")
  
  for(r in c("0/0","0/1","1/1")){df[,r] <- as.numeric(df[,r])}
  
  var <- 
    df %>% 
    mutate(QC=case_when(
      GT=="1/0" & c(df$`0/0`- df$`0/1`) > 10 ~ "ok",
      GT=="1/1" & c(rowMeans(df[,c("0/0","0/1")])-df$`1/1`) > 10 ~"ok",
      TRUE ~ "fail"
    )) %>% 
    filter(DP>1) %>% 
    filter(QC=="ok") %>% 
    rownames()
  return(var)
})
names(call.var) <- bc

# Read in Count table 
counts.alt <- read.csv("alt_filtered.csv")
#counts.ref <- read.csv("ref_filtered.csv")

filter.alt <- rowSums(counts.alt[,-1]) != 0
counts.alt <- counts.alt[filter.alt, ]
dim(counts.alt)

#filter SNPs 


re.counts.alt <- 
  reshape2::melt(counts.alt) %>% 
  rename("Count":=value) %>% 
  rename("Pos":=SNV) %>% 
  rename("barcodes":=variable) %>% 
  filter(Count!=0) %>% 
  left_join(.,variants.string, by="Pos" )


#Get Mutations per spot

mutations <- {re.counts.alt}$Mut %>% unique()

snp.var <- 
  re.counts.alt %>% 
  mutate(barcodes=str_replace_all(as.character(barcodes), pattern="[.]", "-"))

mut.spot <- 
  snp.var %>% 
  pivot_wider(names_from = Mut, values_from = Count, values_fill = 0) %>% 
  group_by(barcodes) %>% 
  summarise_at(.vars={{mutations}},.funs=sum) %>% 
  ungroup() %>% 
  mutate(snp.var= map_dbl(.x=barcodes, .f=function(bc){{{snp.var}} %>% filter(barcodes %in% bc) %>% pull(SNP) %>% unique() %>% length() })) %>% 
  left_join(.,getFeatureDf(TX) %>% select(barcodes, nCount_Spatial), by="barcodes") %>% 
  as.data.frame() %>% 
  mutate(Total.snp=rowSums(.[,mutations]))

mut.spot <- mut.spot %>% mutate(var.score=c((Total.snp)/(snp.var+1)))


ggplot(data=mut.spot, mapping = aes(x=var.score, y=nCount_Spatial, color=var.score))+
  geom_point(size=mut.spot$var.score+1)+theme_classic()+
  xlim(1,1.5)+
  SPATA2::scale_color_add_on(clrsp = "Reds", limit=c(1,1.5))

ggplot(re.counts.alt %>% arrange(Count), aes(fill=Mut, y=Count, x=barcodes)) + 
  geom_bar(position="stack", stat="identity")+theme_void()+
  scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(13))



setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions")
TX <- readRDS("269_T_SPATA_CNV_Pred.RDS")
TX@used_genesets <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/GS_new.RDS")
TX <- TX %>% updateSpataObject()


TX <- TX %>% SPATA2::addFeatures(., mut.spot, feature_names = c(mutations, "Total.snp", "snp.var", "var.score"), overwrite = T)
TX@fdata$`269_T`[is.na(TX@fdata$`269_T`)]=0
TX %>% getFeatureDf() %>% pull(var.score) %>% is.infinite() %>% table()
TX %>% plotSurfaceInteractive()



# SNP Allel mutation type by cluster

data.plot <- joinWith(TX, features = c(mutations, "seurat_clusters" )) %>% select(mutations, "seurat_clusters") %>% pivot_longer(., cols={mutations})

ggplot(data.plot, aes(fill=name, y=value, x=seurat_clusters)) + 
  geom_bar(position="stack", stat="identity")+theme_classic()+
  scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(13))




# Check for Oncogenes -----------------------------------------------------

oncoGenes <- read.table("~/Desktop/SpatialTranscriptomics/Visium/Visium/mtDNA/OncoGenes.txt", header=T)


#Anotate gene


library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("hgnc_symbol")
filters <- c("chromosome_name","start","end")
values <- list(chromosome= c(10),
               start=c("102476888"),
               end=c("102476888"))


clean.pos <- 
  snp.var %>% 
  filter(Count>5) %>% 
  tidyr::separate(.,col=SNP, into=c("X","chromosome", "start"), sep="[_]") %>% 
  tidyr::separate(.,start, into=c("start", "XX"), sep="[:]") %>% 
  dplyr::select(-X,-XX) %>% 
  mutate(end=start)


library(future)
library(furrr)
plan("multiprocess", workers = 7)
options(future.globals.maxSize = 3000 * 1024^2)
clean.pos.genes <- furrr::future_map_dfr(.x=1:nrow(clean.pos), .f=function(i){
    pb$tick()$print()
    res <- getBM(attributes =  attributes, filters =  filters,
               values = as.list(clean.pos[i,c("chromosome", "start", "end")]) , mart = mart)
    if(nrow(res)>1){res <- cbind(res, clean.pos[i,])}else{
      res <- cbind(data.frame(hgnc_symbol=NA), clean.pos[i,])
    }
    
    
    return(res)
  }, .progress = T)

clean.pos.genes.filter <- clean.pos.genes %>% filter(hgnc_symbol %in% oncoGenes$Gene) %>% arrange(hgnc_symbol)
clean.pos.genes %>% dplyr::count(hgnc_symbol) %>% arrange(desc(n))

mut.df <- 
  SPATA2::getFeatureDf(TX) %>% 
  mutate(CD151_mut=ifelse(barcodes %in% c(
                                        clean.pos.genes %>% 
                                          filter(hgnc_symbol=="DNAJC25-GNG10") %>% 
                                          pull(barcodes) ),
                                        "1", "0")) %>% 
  dplyr::select(barcodes, CD151_mut)

TX <- SPATA2::addFeatures(TX, mut.df, overwrite = T)

TX %>% plotSurfaceInteractive()



# Mutation Data 275_T -----------------------------------------------------



setwd("~/Desktop/SpatialTranscriptomics/Visium/Visium/All_SPATA_Revisions")
spata.obj <- readRDS("275_T_SPATA_CNV_Pred.RDS")
spata.obj@used_genesets <- readRDS("~/Desktop/SpatialTranscriptomics/Visium/Visium/GS_new.RDS")
spata.obj <- spata.obj %>% updateSpataObject()

path.snv="~/Desktop/SpatialTranscriptomics/Visium/Visium/275_T/275_SNP/alt_filtered.csv"
path.vcf="~/Desktop/SpatialTranscriptomics/Visium/Visium/275_T/275_SNP/275_T.vcf"

spata.obj <- SPATAwrappers::importSNVmatrix(spata.obj,path.snv=path.snv, path.vcf=path.vcf)
spata.obj <- SPATAwrappers::runSpotwiseSNV(spata.obj, as.features = T)

spata.obj %>% plotSurfaceInteractive()
spata.obj <- SPATA2::createSegmentation(spata.obj)
library(patchwork)
plotSurface(spata.obj, color_by = "seurat_clusters")+
plotSNVtype(spata.obj, feature = "seurat_clusters")

#Integrate information in the Subtype heatmap

gs <- c("Module_consensus_1","Module_consensus_2","Module_consensus_4","Module_consensus_3","Module_consensus_5")
dev.off()
SPATAwrappers::plotHeatmap(spata.obj, gs=gs, thr=c(0.67,0.7,0.7,0.7,0.65))[[1]]
out <- SPATAwrappers::plotHeatmap(spata.obj, gs=gs,thr=c(0.67,0.7,0.7,0.7,0.65), plot.type = "point", feature = "var.score", sample.type = "seurat_clusters")
plotSNVtype(spata.obj, order=out[[3]]$barcodes, type="fill")
out[[1]]

all.com.plots <- 
map(.x=spata.obj@data$`275_T`$SNV$Conversion.type, .f=function(con){
  SPATAwrappers::plotHeatmap(spata.obj, gs=gs,thr=c(0.67,0.7,0.7,0.7,0.65), 
                             plot.type = "point", feature = con, sample.type = "seurat_clusters")[[1]]+Seurat::NoLegend()
})
all.com.plots %>% patchwork::wrap_plots()

SPATA2::plotStatisticsInteractive(out[[3]] %>% mutate(module.2=str_remove(module,pattern="Module_consensus_")))












