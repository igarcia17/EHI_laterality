library(dplyr)
library(WGCNA)
library(ggplot2)
library(ggrepel)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

#Functions
run_step_1 <- F
run_step_2 <- F
#Parameters
tissue <- "Amygdala"

if(run_step_1){
  
  cutoffs <- c("0", "60")
  oldies <- F
  suffix <- gsub("[^0-9]", "_", as.character(Sys.Date()))
  outF <- paste0("WGCNA_brain/", tissue, "/")
  samtreeF <- paste0(outF, "sampleTree_", suffix, ".jpeg")
  contreeF <- paste0(outF, "sample_connectivity_", suffix, ".jpeg")
  res_step1F <- paste0(outF, "data_step1_", tissue, "_", suffix, ".RData")
  filesF <- paste0("Imputed_brain/Brain_", tissue, "/")
  
config_df <- as.data.frame(readxl::read_xlsx("All_data_available_genetic_exp.xlsx"
                                             , sheet = 1))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>% mutate(Diagnostic = factor(Diagnostic),
                Cutoff0 = factor(Cutoff0),
                Cutoff40= factor(Cutoff40),
                Cutoff60 = factor(Cutoff60),
                Cutoff80 = factor(Cutoff80),Cutoff90 = factor(Cutoff90),
                Age = as.numeric(Age),
                Score10items =as.numeric(Score10items),
                PRS_handedness= as.numeric(PRS_handedness),
                PRS_amb = as.numeric(PRS_amb),
                PRS_BD = as.numeric(PRS_BD),
                PRS_SCZ = as.numeric(PRS_SCZ),
                Age_group =factor(Age_group),
                WAVE = factor(WAVE)
  )%>%
  mutate(across(
    c(Cutoff0, Cutoff40, Cutoff60, Cutoff80, Cutoff90),
    ~ factor(.x, levels = c("RIGHT", "MIXED", "LEFT"))))%>%
  mutate(across(
    paste0("Cutoff", cutoffs),
    ~ factor(ifelse(.x == "MIXED",
                    "NON-LATERAL",
                    "LATERAL"),
             levels = c("LATERAL", "NON-LATERAL")),
    .names = "Lateral_Cutoff{str_remove(.col, 'Cutoff')}"
  ))%>%
  mutate(across(
    paste0("Cutoff", cutoffs),
    ~ factor(ifelse(.x == "RIGHT",
                    "RH",
                    "NRH"),
             levels = c("RH", "NRH")),
    .names = "NRH_{str_remove(.col, 'Cutoff')}"
  ))


if(oldies){
  config_df <- config_df %>%
    filter(Oldies_dataset == "YES") %>% arrange(ID) %>% 
    filter(has_RNA=="YES")%>%droplevels()
}else{
  config_df <- config_df%>%filter(has_RNA=="YES")%>%arrange(ID) %>% droplevels()
}

sapply(config_df, class)

#Load genetic expression data
files <- sub("\\.tsv$", "", config_df$File)
# IDs desde input$File, quitando .tsv, porque tienen esa R
sample_ids <- sub("\\.tsv$", "", config_df$File)

expr_list <- lapply(sample_ids, function(sample_id) {
  matched_file <- list.files(
    path = filesF,
    pattern = paste0("^", sample_id, "_.*\\.tsv$"),
    full.names = TRUE)
  
  if (length(matched_file) == 0) {
    warning(paste("No se encontró archivo para:", sample_id))
    return(NULL)
  }
  
  if (length(matched_file) > 1) {
    warning(paste("Más de un archivo para:", sample_id, "- usando el primero"))
    matched_file <- matched_file[1]
  }
  
  df <- read.table(
    matched_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  counts <- df[, 2]
  names(counts) <- df[, 1]
  
  return(counts)
})
expr_matrix <- do.call(cbind, expr_list)
colnames(expr_matrix) <- config_df$ID
rm(expr_list)
gc()

print(paste0("Loaded counts of ", tissue))

#Ensrue metadata is correctly ordered
all(colnames(expr_matrix) %in% config_df$ID)
config_df <- config_df[match(colnames(expr_matrix), config_df$ID), ]
config_df <- config_df%>%
  select(-c("has_RNA",
            'Oldies_dataset', 'Cutoff0', 'Cutoff40', 'Cutoff60', 'Cutoff80', 
            'Cutoff90', "PC1_oldies", "PC2_oldies", "PC3_oldies","PC4_oldies",
            "PC1_all","PC2_all","PC3_all","PC4_all"))

datExpr <- as.data.frame(t(expr_matrix))

#QC
gsg <- goodSamplesGenes(datExpr, verbose=3)
#Cluster to detect outliers in samples
sampleTree <- hclust(dist(datExpr), method = 'average')
sample_names <- sampleTree$labels
clust <- cutreeStatic(sampleTree, cutHeight = 50, minSize = 1)
tab <- table(clust)
split(sampleTree$labels, clust)

jpeg(file = samtreeF, width = 1700, height = 1700)
par(mar = c(0,4,2,0))
title <- "Sample clustering to detect outliers"
plot(sampleTree, main = title, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
invisible(dev.off())

#Eliminate outliers
outliers <- c("UAM_0031", "UAM_0321")
datExpr <- datExpr[rownames(datExpr) != outliers,]
config_df <- config_df %>% filter(ID != outliers)

save(datExpr, config_df, file = res_step1F)

rm(list = setdiff(ls(), c("res_step1F", "tissue")))
}

###########STEP 2###############################################################

suffix <- gsub("[^0-9]", "_", as.character(Sys.Date()))
outF <- paste0("WGCNA_brain/", tissue, "/")
res_step1F <- paste0(outF, "data_step1_", tissue, "_", suffix, ".RData")
load(res_step1F)

sign <- T
run_scale_free_topology <- T

#Parameters: inspection of scale free topology by network type
if(run_scale_free_topology){
if (sign){
  powers <- c(c(1:10), seq(from = 12, to=20, by=2)) #possible thresholds
  sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed', 
                           verbose = 5)
  filename <- paste0(outF, "Soft_threshold_signed_net.tiff")} else {
    sft <- pickSoftThreshold(datExpr, powerVector = powers, 
                             networkType = 'unsigned', verbose = 5)
    filename <- paste0(outF, "Soft_threshold_unsigned_net.tiff")}

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
}

run_step_2<-F
#Set the soft threshold parameter at 12 as it is the lowest at which
#we obtain a R^2 of 0.85 though ideally it should reach 0.9
sft_power <- 20

if(runt_step_2){
auto_net_build <- T
if(auto_net_build){
  minMod <- 15
#Automatic network construction for different module sizes
net <- blockwiseModules(datExpr, networkType = 'signed', minModuleSize = minMod,
                        power = sft_power, randomSeed = 1,
                        numericLabels = F, verbose = 5)
moduleColors <- labels2colors(net$colors)
geneTree <- net$dendrograms[[1]]

#Plot
jpeg(file = paste0(outF, "autonet_constr_minMod",minMod,".jpeg"), units="in", width=15, height=5, res=300)
blocks <- moduleColors[net$blockGenes[[1]]]
plotDendroAndColors(geneTree, blocks, "Module colors",dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)
invisible(dev.off())

#save results
moduleLabels <- net$colors
MEs <- net$MEs
save(MEs, moduleLabels, moduleColors, geneTree, file = paste0(outF, "automatic_network_minMod",minMod,".RData"))
sizes <- sort(table(moduleColors), decreasing=T)

jpeg(filename=paste0(outF, "distribution_size_minMod",minMod,"_auto.jpeg"))
hist(as.numeric(sizes),
     breaks = 20,
     main = "Distribution module size",
     xlab = "Number genes by module")
dev.off()

}


#STEP BY STEP NETWORK CONSTRUCTION
run_manual_net <- F

if(run_manual_net){
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
save(tom, diss_tom, geneTree, file = paste0(outF,"blood_manual_net_input_matrixes_", suffix, ".RData"))
}
load(file=paste0(outF,"blood_manual_net_input_matrixes_2026_05_26.RData"))

minMod <- 30
MEDissThr <- 0.10 #threshold to which merge modules in manual construction
#height cut of 0.25 which corresponds to a correlation of 75% among modules

#Cut the dendrogram
dynamicMods <- cutreeDynamic(geneTree, distM = diss_tom, deepSplit =2, 
                             pamRespectsDendro =F, minClusterSize = minMod)
table(dynamicMods)
dynColors <- labels2colors(dynamicMods)
title <- 'Gene dendrogram and module colors'
filename <- paste0(outF, "blood_geneTree_manualnet_minMod",minMod,"_meDissThr",MEDissThr,"_",suffix,".pdf")
pdf(filename)
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
filename <- paste0(outF, "module_eigengene_cluster_final_minMod",minMod, "_meDissThr", MEDissThr, "_", suffix, ".pdf" )
pdf(file = filename, width = 20, height = 8)
plot(METree, main = title)
abline(h=MEDissThr, col = 'red')
invisible(dev.off())

#Merge modules with high correlation among them
merge <- mergeCloseModules(datExpr, dynColors, cutHeight = MEDissThr, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

#Plot the modules with and without the merging
filename <- paste0(outF, "gene_dendro_tree_manual_final_minMod",minMod, "_meDissThr", MEDissThr, "_", suffix, ".pdf" )
pdf(file=filename, wi = 20, he = 6)
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
res_step_2F <- paste0(outF, "blood_manual_net_minMod",minMod, "_meDissThr", MEDissThr, "_", suffix, ".RData" )
save(MEs, moduleLabels, moduleColors, geneTree, file=res_step_2F)

sizes <- sort(table(moduleColors), decreasing=T)
filename <- paste0(outF, "distribution_size_manual_minMod",minMod, "_meDissThr", MEDissThr, "_", suffix, ".jpeg" )
jpeg(filename=filename)
hist(as.numeric(sizes),
     breaks = 20,
     main = "Distribution module size",
     xlab = "Number genes by module")
dev.off()

}




