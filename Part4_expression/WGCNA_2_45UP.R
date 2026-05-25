#Weighted Gene Coexpression Network Analysis part 2, for network construction
#Part 2: Network construction and module detection
library(dplyr)
library(ggplot2)
library(WGCNA)

options(stringsAsFactors = F)
enableWGCNAThreads()

input <- paste0('./out_45nup/data_step1_2025_11_26.RData')
load(file = input)

powers <- c(c(1:10), seq(from = 12, to=20, by=2)) #possible thresholds
sign <- T
#Parameters: inspection of scale free topology by network type
if (sign){
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed', 
                           verbose = 5)
  filename <- "out_45nup/Soft_threshold_signed_net.tiff"} else {
  sft <- pickSoftThreshold(datExpr, powerVector = powers, 
                           networkType = 'unsigned', verbose = 5)
  filename <- "out_HR/Soft_threshold_unsigned_net.tiff"}

#Plot how does data fit the scale free network model and mean connectivity
#based on the threshold used.
x <- sft$fitIndices[,1]
xlab <- "Soft Threshold (power)"
y1 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
ylab1 <- "Scale Free Topology Model Fit, signed R^2"
y2 <- sft$fitIndices[,5]
ylab2 <- "Mean Connectivity"

cex1 <- 0.9
tiff(file = filename, units="in", width=15, height=5, res=300)
par(mfrow = c(1,2))

plot(x, y1, xlab= xlab, ylab= ylab1, type="n", main = paste("Scale independence"))
text(x, y1, labels=powers, cex=cex1, col="red")
abline(h=0.85, col="orange")
abline(h=0.8, col="red")

plot(x, y2, xlab=xlab, ylab=ylab2, type="n", main = ylab2)
text(x, y2, labels=powers, cex=cex1, col="red")
invisible(dev.off())

#Set the soft threshold parameter at 12 as it is the lowest at which
#we obtain a R^2 of 0.85 though ideally it should reach 0.9
sft_power <- 12
minMod <- 25

#Automatic network construction for different module sizes
net <- blockwiseModules(datExpr, networkType = 'signed', minModuleSize = minMod,
                        power = sft_power, randomSeed = 1,
                        numericLabels = F, verbose = 5)
moduleColors <- labels2colors(net$colors)
geneTree <- net$dendrograms[[1]]

#Plot
jpeg(file = "out_45nup/autonet_constr_minMod25.jpeg", units="in", width=15, height=5, res=300)
blocks <- moduleColors[net$blockGenes[[1]]]
plotDendroAndColors(geneTree, blocks, "Module colors",dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
invisible(dev.off())

#save results
moduleLabels <- net$colors
MEs <- net$MEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "out_45nup/automatic_network_minMod25.RData")
sizes <- sort(table(moduleColors), decreasing=T)

jpeg(filename="out_45nup/distribution_size_minMod25_auto.jpeg")
hist(as.numeric(sizes),
     breaks = 20,
     main = "Distribution module size",
     xlab = "Number genes by module")
dev.off()

minMod <- 25
MEDissThr <- 0.15 #threshold to which merge modules in manual construction
#height cut of 0.25 which corresponds to a correlation of 75% among modules

##STEP BY STEP NETWORK CONSTRUCTION

#co expression similarity and adjacency
adjac <- adjacency(datExpr, power = sft_power, type = 'signed')
#topological overlap matrix: compute how much directly and indirectly two genes are
tom <- TOMsimilarity(adjac, TOMType = 'signed', verbose = 5)
diss_tom <- 1-tom
#Gene tree plot
geneTree <- hclust(as.dist(diss_tom), method = 'average')
title <- 'Gene clustering on TOM-based dissimilarity'
plot(geneTree, xlab = '', sub = '', main = title, labels = F, hang = 0.04)

#Save this objects, as running previous commands is very computationally expensive
#save(tom, diss_tom, geneTree, file = "tempRdata.RData")
#load(file="tempRdata.RData")

#Cut the dendrogram
dynamicMods <- cutreeDynamic(geneTree, distM = diss_tom, deepSplit =2, 
                             pamRespectsDendro =F, minClusterSize = minMod)
table(dynamicMods)
dynColors <- labels2colors(dynamicMods)
title <- 'Gene dendrogram and module colors'
pdf("out_45nup/geneTree_manualnet.pdf")
plotDendroAndColors(geneTree, dynColors, 'Dynamic Tree Cut', dendroLabels = F, 
                    hang=0.03, addGuide = T, guideHang = 0.05, main = title)
invisible(dev.off())

#Merge modules with similar expression profile
MEList <- moduleEigengenes(datExpr, colors = dynColors)
MEs <-MEList$eigengenes
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method = 'average')

#Plot
title <- 'Clustering of Module Eigengenes'
pdf(file = "out_45nup/module_eigengene_cluster_final.pdf", width = 20, height = 8)
plot(METree, main = title)
abline(h=MEDissThr, col = 'red')
invisible(dev.off())

#Merge modules with high correlation among them
merge <- mergeCloseModules(datExpr, dynColors, cutHeight = MEDissThr, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

#Plot the modules with and without the merging
pdf(file="out_45nup/gene_dendro_tree_manual_final.pdf", wi = 20, he = 6)
plotDendroAndColors(geneTree, cbind(dynColors, mergedColors), 
                    c('Dynamic Tree Cut','Merged dynamic'),
                    dendroLabels = F, hang = 0.03, guideHang = 0.05, addGuide = T)
invisible(dev.off())

#Save the results
moduleColors <- mergedColors
colorOrder <- c('grey', standardColors(n=NULL))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

#Save module colors and labels for subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file="out_45nup/manual_net_2025_11_26.RData")

sizes <- sort(table(moduleColors), decreasing=T)
jpeg(filename="out_45nup/distribution_size_manual.jpeg")
hist(as.numeric(sizes),
     breaks = 20,
     main = "Distribution module size",
     xlab = "Number genes by module")
dev.off()
