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
rm(list = ls())
#Paths to files
filesD <- 'WGCNA_blood/'
inputData <- paste0(filesD,'data_step1_2026_05_26.RData')
inputNet <- "blood_manual_net_minMod30_meDissThr0.1_2026_05_27"
inputNetF <- paste0(filesD, inputNet, '.RData')
load(inputData)
load(inputNetF)

suffix <- gsub("[^0-9]", "_", as.character(Sys.Date()))
plotTOMF <- paste0(filesD, 'TOM_heatmap', inputNet, "_", suffix, '.tiff')
dendroEigenF <- paste0(filesD, 'Eigengene_network_dendrogram', inputNet, "_", suffix, '.tiff')
heatEigenF <- paste0(filesD, 'Eigengene_net_heatmap', inputNet, "_", suffix, '.tiff')
sft_power <- 12
enhance_power <- 10

#Visualize gene network
plotTOM <- (1-TOMsimilarityFromExpr(datExpr, power = sft_power))
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
