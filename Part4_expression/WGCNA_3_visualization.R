#Visualization of gene correlation network characteristics and modules. 
# Based on Peter Langfelder tutorials
suppressPackageStartupMessages({
  library(dplyr, quietly = T)
  library(ggplot2, quietly = T)
  library(WGCNA, quietly = T)
  library(RColorBrewer, quietly =T)
})
options(stringsAsFactors = F)
enableWGCNAThreads()

#Paths to files
filesD <- 'out_45nup/'
inputData <- paste0(filesD,'data_step1_2025_11_26.RData')
inputNet <- paste0(filesD, 'manual_net_2025_11_26.RData')
load(inputData)
load(inputNet)

suffix <- gsub("[^0-9]", "_", as.character(Sys.Date()))
plotTOMF <- paste0(filesD, 'TOM_heatmap', suffix, '.tiff')
dendroEigenF <- paste0(filesD, 'Eigengene_network_dendrogram', suffix, '.tiff')
heatEigenF <- paste0(filesD, 'Eigengene_net_heatmap', suffix, '.tiff')
sft_power <- 12
enhance_power <- 10

#Visualize gene network
plotTOM <- (1-TOMsimilarityFromExpr(datExpr, power = sft_power))^enhance_power
# Transform dissTOM with a power to make moderately strong connections more visible
plotTOM <- plotTOM^enhance_power
# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA

tiff(file = plotTOMF)
title <- "Network heatmap plot, all genes"
TOMplot(plotTOM, geneTree, moduleColors, main = title)
invisible(dev.off())
rm(plotTOM)

#Vizualize eigengene network
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
df <- MEs
MET <- orderMEs(df)
#Eigengene network dendrogram
tiff(filename = dendroEigenF)
title <- "Eigengene dendrogram"
plotEigengeneNetworks(MET, title, marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
invisible(dev.off())

#Eigengene network heatmap
cols <- brewer.pal(n = 7, name = "YlGnBu")
tiff(filename = heatEigenF)
title <- "Eigengene correlation heatmap"
plotEigengeneNetworks(MET, title, marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90, signed= T,
                      heatmapColors = cols, plotAdjacency = F)
invisible(dev.off())
