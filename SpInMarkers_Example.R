library(gridExtra)
library(patchwork)
library(dplyr)
library(hdf5r)
library(EnhancedVolcano)
library(ggpubr)
library(stringr)
library(viridis)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(grid)
library(gridtext)
load10XExpr<- function(visiumDir = NULL,h5filename= 'filtered_feature_bc_matrix.h5'){
  h5FilePath <- dir(path = visiumDir,pattern = h5filename,full.names = T)
  
  hf <- hdf5r::h5file(filename = h5FilePath, mode='r')
  mat <- names(hf)
  
  
  counts <- hf[[paste0(mat, '/data')]]
  indices <- hf[[paste0(mat, '/indices')]]
  indptr <- hf[[paste0(mat, '/indptr')]]
  shp <- hf[[paste0(mat, '/shape')]]
  features <- hf[[paste0(mat, '/features/name')]][]
  barcodes <- hf[[paste0(mat, '/barcodes')]][]
  spMat <- Matrix::sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )
  spMat <- log2(1+spMat)
  features <- make.unique(names = features)
  rownames(spMat) <- features
  colnames(spMat) <- barcodes
  hf$close_all()
  return(spMat)
}


patient_id <- 'PATIENT_ID'
visiumDir <- "PATH_TO_VISIUM_DIR"
cogapsFilePath <- 'PATH_TO_COGAPS_DIR'
OUTPUT_DIR <- 'PATH_TO_OUTPUT_DIR'
## Loading data
pVisiumPath <- paste0(visiumDir,patient_id)
fullMat <- load10XExpr(pVisiumPath)
good_gene_threshold <- 3
goodGenes <- rownames(fullMat)[apply(fullMat,1,function(x) sum(x>0)>=good_gene_threshold)]
fullMat <- fullMat[goodGenes,]
st_hcc <- Load10X_Spatial(pVisiumPath)

SpInMarkers <- readRDS(gsub(x =  cogapsFilePath, pattern = '.rds', replacement = '_SpInMarkers.rds'))

hotspot <- SpInMarkers$hotspotRegions

expMar <- as.matrix(fullMat)
hotspot <- hotspot[colnames(expMar),]

patternFirst <- 1
patternSecond <- 8
region <- hotspot[,patternFirst]
region <- ifelse(!is.na(region) & !is.na(hotspot[,patternSecond]),"Interacting",ifelse(!is.na(region),region,hotspot[,patternSecond]))
##setting the rownames for region to combine with seurat metadata
region_df <- data.frame(region)
region_df <- region_df %>% replace(is.na(.), "Neither")
rownames(region_df) <- rownames(hotspot)
region <- factor(region, levels = c("Pattern_1", "Interacting", "Pattern_8"))

## Directing outputs to output folder
setwd(OUTPUT_DIR)

st_hcc@meta.data <- cbind(st_hcc@meta.data, region_df)
SpatialDimPlot(st_hcc, group.by = "region") & scale_fill_viridis(discrete = TRUE, alpha=0.6) &
  theme(legend.position='right', legend.title = element_text(size=16), legend.text = element_text(size=12)) &
  guides(fill = guide_legend(override.aes = list(size=5)))

mplot2 <- data.frame(t(expMar[,!is.na(region)]))
region <- region[!is.na(region)]
mplot2 <- cbind(mplot2, region)


#geneList <- c( 'CDH5', 'CD34')
#geneList <- c('ALB', 'ALDH1A1', 'IGLC2', 'CCL19')
geneList <- c('TGFB1', 'PECAM1')
plist <- list()

my_comparisons <- list( c("Interacting", "Pattern_1"), c("Pattern_1", "Pattern_8"),c("Interacting", "Pattern_8") )


for (ii in 1:length(geneList)){
  plist[[ii]]<- mplot2 %>% ggplot( aes_string(x="region", y=geneList[ii], fill="region")) +
    geom_violin() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme(
      legend.position="none",
      plot.title = element_text(size=14),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 14)
    ) +
    ggtitle(paste0(geneList[ii]," Expression (Log)")) +
    xlab("") + 
    scale_x_discrete(labels=c('Cancer', 'Interacting', 'Immune')) + 
    stat_compare_means(comparisons = my_comparisons, size = 4)
}
n <- length(plist)
nCol <- floor(sqrt(n))


title <- textGrob("HCC-1 (R)", gp = gpar(fontsize = 20))

grid.arrange(grobs = plist, ncol = 2, top = title)
