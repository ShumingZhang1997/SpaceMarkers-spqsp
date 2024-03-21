## author: Atul Deshpande
## email: adeshpande@jhu.edu
rm(list = ls())

getSpatialParameters_spQSP <- function(spatialPatterns){
  good_gene_threshold <- 3;  
  sigmaRes <- max(floor(min(diff(range(spatialPatterns$x)),diff(range(spatialPatterns$y)))/250),1)
  sigmaVec <- seq(2,40*sigmaRes,sigmaRes)
  threshVec <- seq(1,2,0.1)
  
  #goodGenes <- rownames(fullMat)[apply(fullMat,1,function(x) sum(x>0)>=good_gene_threshold)]
  patternList <- colnames(spatialPatterns)[startsWith(colnames(spatialPatterns),"Pattern_")]
  optParams<-sapply(patternList, function(i) unlist(getOptimalSigmaThresh(pattern = spatialPatterns[,i], 
                                                                          locs = data.frame(x = spatialPatterns$x, y = spatialPatterns$y),
                                                                          sigmaVec, threshVec)))
  return(optParams)
}


library(ggplot2)
library(dplyr)
library(tidyr)
library(matrixStats)
library(caret)
library(gridExtra)
library(patchwork)
library(hdf5r)
library(ggpubr)
library(stringr)
library(viridis)

source("./R/preprocessing.R")
source("./R/getSpatialParameters.R")
source("./R/getInteractingGenes.R")
source('./R/find_genes_of_interest_nonparametric_fast.R')

working_dir <- 'PATH_TO_WORKING_DIR'


#spQSP Simulation results
input_folder <- './sample_simulated_result'
#Desired output folder
output_folder <- './out_sample_simulated_result'
setwd(working_dir)
# SpInMarkersMode: defaut mode is "residual". You can also set "DE" mode for Differential Expression mode.
SpInMarkersMode = "DE"  
# SpinMarkersRefPattern is the pattern whose "interaction" with every other pattern we want to study. 
# Pattern_1: Cancer Cell
# Pattern_2 CD8+ T cell
SpinMarkersRefPattern = "Pattern_1" 

patient_sample <- commandArgs(trailingOnly = TRUE)[1]


y_slice <- 5
cellPath <- sprintf('%s/snapShots/cell_280.csv', input_folder)
cellData_raw <- read.csv(cellPath)

#filter T_EFF and TEXH
cellData <- cellData_raw[!(cellData_raw$State %in% c(3,5)),]

#preprocessing and organize cytokine data
cytokinePath <- sprintf('%s/snapShots/grid_core_280.csv', input_folder)
cytokineData <- read.csv(cytokinePath)

#preprocessing and organize cytokine data
vasPath <- sprintf('%s/snapShots/vas_280.csv', input_folder)
vasData <- read.csv(vasPath)

##Organize Type Data
cellData$Pattern_ <- as.character(cellData$Type)
dmy <- dummyVars("~ Pattern_", cellData)
cellData <- cbind(cellData,  data.frame(predict(dmy, newdata=cellData)))

cellData$State_ <- as.character(cellData$State)
dmyState <- dummyVars("~ State_", cellData)
cellData <- cbind(cellData,  data.frame(predict(dmyState, newdata=cellData)))

##Addup all cells in each coordinate
cellData_2D <- cellData[cellData$y == y_slice,]
#typeData_2D <- subset (cellData_2D, select = c(x,z,Pattern_1, Pattern_2, Pattern_3, Pattern_4, Pattern_5))
typeData_2D <- subset (cellData_2D, select = c(x,z,Pattern_1, Pattern_2))
names(typeData_2D)[names(typeData_2D) == 'z'] <- 'y'

typeData_2D_agg <- aggregate(cbind(typeData_2D$Pattern_1, typeData_2D$Pattern_2),
                             list(typeData_2D$x, typeData_2D$y), sum)

#typeData_2D_agg <- aggregate(cbind(typeData_2D$Pattern_1, typeData_2D$Pattern_2, typeData_2D$Pattern_3, cellData_2D$Pattern_4, typeData_2D$Pattern_5 ),
#list(typeData_2D$x, typeData_2D$y), sum)

names(typeData_2D_agg) <- c('x', 'y', 'Pattern_1', 'Pattern_2')
#names(typeData_2D_agg) <- c('x', 'y', 'Pattern_1', 'Pattern_2', 'Pattern_3', 'Pattern_4', 'Pattern_5')
spPatterns_spQSP <- cbind(barcodes = paste(typeData_2D_agg$x, "_", typeData_2D_agg$y), typeData_2D_agg)
rownames(spPatterns_spQSP) <- spPatterns_spQSP$barcodes


expressionData = merge(cytokineData, vasData, by.x=c('x', 'y', 'z'), by.y=c('x', 'y', 'z'))

expressionData_2D <- expressionData[expressionData$y == 5,]
expressionData_2D <- subset (expressionData_2D, select = -c(y))
names(expressionData_2D)[names(expressionData_2D) == 'z'] <- 'y'
rownames(expressionData_2D) <- paste(expressionData_2D$x, "_", expressionData_2D$y)
filtered_expressionData_2D <- expressionData_2D[rownames(spPatterns_spQSP),]
colnames(filtered_expressionData_2D)[13] <- 'Vasculature'
cgMat_spQSP <- t(subset (filtered_expressionData_2D, select = -c(x,y)))
fullMat_spQSP <- as(cgMat_spQSP, "sparseMatrix") 

#spPatterns_spQSP <- subset(spPatterns_spQSP, select = -c(Pattern_3, Pattern_4, Pattern_5) )

optParams_spQSP <- getSpatialParameters_spQSP(spPatterns_spQSP)
optParams <- optParams_spQSP

SpInMarkers_spQSP <- getInteractingGenes(data = fullMat_spQSP, reconstruction = NULL, 
                                         spatialPatterns = spPatterns_spQSP, 
                                         refPattern = SpinMarkersRefPattern, mode = SpInMarkersMode)

#replace Pattern_1 value to either NA(0) or "Pattern_1" (1)
hotspot_matrix <- subset(spPatterns_spQSP, select = -c(x, y, barcodes))
rownames(hotspot_matrix) <- rownames(SpInMarkers_spQSP[["hotspotRegions"]])
# 
hotspot_matrix$Pattern_1[hotspot_matrix$Pattern_1 == 0] <- NA
hotspot_matrix$Pattern_1[hotspot_matrix$Pattern_1 == 1] <- "Pattern_1"
#hotspot_matrix$Pattern_1 <- SpInMarkers_spQSP[["hotspotRegions"]][,'Pattern_1']
#hotspot_matrix$Pattern_4[hotspot_matrix$Pattern_4 == 0] <- NA
#hotspot_matrix$Pattern_4[hotspot_matrix$Pattern_4 == 1] <- "Pattern_4"
#hotspot_matrix$Pattern_5[hotspot_matrix$Pattern_5 == 0] <- NA
#hotspot_matrix$Pattern_5[hotspot_matrix$Pattern_5 == 1] <- "Pattern_5"
# ###
hotspot_matrix$Pattern_2 <- SpInMarkers_spQSP[["hotspotRegions"]][,'Pattern_2']
#hotspot_matrix$Pattern_3 <- SpInMarkers_spQSP[["hotspotRegions"]][,'Pattern_3']
hotspot_matrix <- as.matrix(hotspot_matrix)

SpInMarkers_spQSP[["hotspotRegions"]] <- hotspot_matrix
SpInMarkers_spQSP <- getInteractingGenes(data = fullMat_spQSP, reconstruction = NULL, 
                                         spatialPatterns = spPatterns_spQSP, refPattern = SpinMarkersRefPattern, 
                                         mode = SpInMarkersMode, hotspotRegions = hotspot_matrix)


dir.create(file.path(output_folder), recursive = TRUE)
setwd(output_folder)
for (i in seq(1,length(SpInMarkers_spQSP$interacting_genes)))
{
  filename <- paste0(gsub(x = 'spQSP_SPIN',pattern = '.rds',replacement = '_SpInMarkers_'),names(SpInMarkers_spQSP$interacting_genes[[i]])[1],".txt")
  write.table(rownames(SpInMarkers_spQSP$interacting_genes[[i]]),file = filename,row.names = F,col.names = F, quote = F)
}  


write.table((!is.na(SpInMarkers_spQSP$hotspotRegions))*1,file = 'hotspot.txt',row.names = T,col.names = T, quote = F)
#write.table((!is.na(hotspot_matrix))*1,file = 'hotspot.txt',row.names = T,col.names = T, quote = F)
write.table(subset(spPatterns_spQSP, select = c(x, y)),file = 'spatial_coord.txt',row.names = T,col.names = T, quote = F)

saveRDS(SpInMarkers_spQSP, file = 'spQSP_SPIN.rds')
################################################################################################################################################################

SpInMarkers_spQSP <- readRDS('spQSP_SPIN.rds')

geneList <- c("IL_2","Vasculature","VEGFA","TGFB" )

region <- SpInMarkers_spQSP$hotspotRegions[,1]
region <- ifelse(!is.na(region) & !is.na(SpInMarkers_spQSP$hotspotRegions[,2]),"Interacting",ifelse(!is.na(region),region,SpInMarkers_spQSP$hotspotRegions[,2]))
region <- factor(region, levels = c("Pattern_1", "Interacting", "Pattern_2"))
my_comparisons <- list( c("Interacting", "Pattern_1"),c("Interacting", "Pattern_2"), c("Pattern_1", "Pattern_2") )

png(filename = "Expression_boxplots_spQSP.png", units="in",width = 9,height = 8, res=600)
plist <- list()
mplot2 <- filtered_expressionData_2D[!is.na(region),]
region <- region[!is.na(region)]
mplot2 <- cbind(mplot2, region)
for (ii in 1:length(geneList)){
  if(geneList[ii] == 'Vasculature'){
    title <- ''
    ylabel = 'Vascular Density'
  }
  else if(geneList[ii] == 'IL_2'){
    title <- ''
    ylabel = 'IL2 Concentration (ng/mL)'
  } 
  else{
    title = ''
    unit <- ' Concentration (ng/mL)'
    ylabel = paste0(geneList[ii], unit)
  }
  plist[[ii]]<- mplot2 %>% ggplot( aes_string(x="region", y=geneList[ii], fill="region")) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    #geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_bw() +
    theme(
      legend.position="none",
      plot.title = element_text(size=14),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    ggtitle(title) +
    xlab("") + 
    ylab(ylabel)+
    scale_x_discrete(labels=c('Cancer', 'Interacting', 'CD8')) +
    stat_compare_means(comparisons = my_comparisons, size = 4)
}
n <- length(plist)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plist, ncol=2))

dev.off()













