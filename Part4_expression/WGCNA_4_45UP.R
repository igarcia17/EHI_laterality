#Weighted Gene Coexpression Network Analysis part 3, to relate modules to external
#clinical traits and gene significance
#Based on Peter Langfelder tutorials
#Annotation based on:
#https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html#4_Over-Representation_Analysis_with_clusterProfiler_-_RNA-seq

library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(WGCNA, quietly = T)
library(org.Hs.eg.db, quietly = T)
library(clusterProfiler)
library(msigdbr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

filesD <- 'out_45nup/'
inputData <- paste0(filesD,'data_step1_2025_11_26.RData')
inputNet <- paste0(filesD, 'manual_net_2025_11_26.RData')
load(inputData)
load(inputNet)

#For output
subD <- "day_26_11_2025/"
modtraitF <- paste0(filesD, subD, "module-trait_relation.jpeg")
intramodFprefix <- paste0(filesD, subD, 'MM_and_GS/intramodular_')
geneInfoF <- paste0(filesD, subD, 'MM_and_GS/genes_info.tsv')
finalresR <- paste0(filesD, subD, 'MM_and_GS/summary_sigMods.RData')
gNumberF <- paste0(filesD, subD, 'MM_and_GS/genes_in_significant_modules.txt')
annotationFprefix <- paste0(filesD, subD, 'annotation/enrich_results_')
genelistsFprefix <- paste0(filesD, subD, "MM_and_GS/genes-ENSGID")

#Parameters
nSamples <- nrow(datExpr)
alpha <- 0.05

category <- 'C5'
subcategory <- NULL
nomenclature_is_ENSGID <- TRUE
#A) Quantify module trait associations

#Transform qualitative/factor traits to binary code
sapply(config_df, class)

levels(config_df$Tube_type) <- c(0,1)
colnames(config_df)[3] <- 'is_PXG'
levels(config_df$Gender) <- c(1,0)
colnames(config_df)[4] <- 'is_female'
levels(config_df$PHQ9_Category) <- c(3,4,1,2,0)
colnames(config_df)[12] <- 'Depression'
levels(config_df$FamHxANY_BD_SCZ_related) <- c(0,1)
colnames(config_df)[10] <- 'FamHx'
levels(config_df$Q26_4BDMeds) <- c(2,0,1,2) #possibleMeds considered as "Current"
colnames(config_df)[15] <- 'anyBDmedEver'
levels(config_df$LiTherapeuticRange) <- c(4,2,3,1)
levels(config_df$CRP_as_afactor) <- c(1,0) # when 1 is more abnormal
levels(config_df$eGFR_as_factor) <- c(1,0) # when 1 is more abnormal
levels(config_df$ISS_MoodState) <- c(-1, 0, 1, 1)
config_df$RNAprocessingBatch <- NULL
config_df$Collection2Arrival_minutes <-NULL
config_df$TimeAtRT_3hoursbin <- NULL
config_df$RNAprocessing2Storage <- NULL
config_df$Arrival2RNAproc_hours_new <- NULL
config_df$is_PXG <- NULL
config_df$General4BDresponse <- NULL
config_df$SUI_Thoughts <- NULL
config_df$LiTherapeuticRange <- NULL
##Exploratory analysis
moduleTrait <- corAndPvalue(MEs, config_df[,3:23])
moduleTraitCor <- moduleTrait$cor
moduleTraitPvalue <- moduleTrait$p

#Get all together:
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #give shape
#Graphical representation:
jpeg(file = modtraitF, width = 1050, height = 1300, quality = 100)
par(mar = c(6, 8.5, 3, 3))
title <- "Module-trait relationships"
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(config_df)[3:23], xLabelsAngle = 45,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, main = title)
invisible(dev.off())

moduleNames <- substring(names(MEs), 3) #3 to take away prefix
#correlation of each gene in each sample to each module
modMemb <- corAndPvalue(datExpr, MEs)
geneModuleMembership <- as.data.frame(modMemb$cor)
names(geneModuleMembership) <- paste("MM", moduleNames, sep="_")

MMPvalue <- as.data.frame(modMemb$p)
names(MMPvalue) <- paste("p_MM", moduleNames, sep="_")

genes <- names(datExpr)
symbol <- mapIds(org.Hs.eg.db, keys = genes, column = 'SYMBOL', 
                 keytype = 'ENSEMBL', multiVals = 'first')

#-----------------Mood state----------------------------------------------------
#correlation of each gene to each condition, thus, gene significance
Mood <- as.data.frame(config_df$ISS_MoodState)
names(Mood) <- 'condition'
geneSignif <- corAndPvalue(datExpr, Mood)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(Mood), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(Mood), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$ISS_MoodState < alpha)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())
}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  if(nomenclature_is_ENSGID){
    write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F
                , col.names = F
                , row.names=F)
  } else {write.table(as.data.frame(modGenes), file = filename, quote = F,
                row.names = FALSE, col.names = FALSE)}
  
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, Mood, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------Depression-------------------------------------------------------
#correlation of each gene to each condition, thus, gene significance
depression <- as.data.frame(config_df$Depression)
names(depression) <- 'condition'
geneSignif <- corAndPvalue(datExpr, depression)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(depression), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(depression), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$Depression < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())
}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
                , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)
  }

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, depression, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------Li_concentration--------------------------------------------
#correlation of each gene to each condition, thus, gene significance
li_range <- as.data.frame(config_df$LiRESULT)
names(li_range) <- 'condition'
geneSignif <- corAndPvalue(datExpr, li_range)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(li_range), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(li_range), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$LiRESULT < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
geneNumber <- length(modGenes)
line <- paste0('Module ', mod, ': ', geneNumber)
write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, li_range, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------LiALDA-------------------------------------------------------
#correlation of each gene to each condition, thus, gene significance
li_range <- as.data.frame(config_df$LiALDA)
names(li_range) <- 'condition'
geneSignif <- corAndPvalue(datExpr, li_range)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(li_range), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(li_range), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$LiALDA < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, li_range, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------CRP_as_afactor-------------------------------------------------------
#correlation of each gene to each condition, thus, gene significance
CRP_as_afactor <- as.data.frame(config_df$CRP_as_afactor)
names(CRP_as_afactor) <- 'condition'
geneSignif <- corAndPvalue(datExpr, CRP_as_afactor)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(CRP_as_afactor), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(CRP_as_afactor), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$CRP_as_afactor < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, CRP_as_afactor, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------eGFR_as_factor------------------------------------------------
#correlation of each gene to each condition, thus, gene significance
eGFR_as_factor <- as.data.frame(config_df$eGFR_as_factor)
names(eGFR_as_factor) <- 'condition'
geneSignif <- corAndPvalue(datExpr, eGFR_as_factor)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(eGFR_as_factor), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(eGFR_as_factor), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$eGFR_as_factor < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, eGFR_as_factor, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------anyBDmedEver--------------------------------------------------


#correlation of each gene to each condition, thus, gene significance
anyBDmedEver <- as.data.frame(config_df$anyBDmedEver)
names(anyBDmedEver) <- 'condition'
geneSignif <- corAndPvalue(datExpr, anyBDmedEver)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(anyBDmedEver), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(anyBDmedEver), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$anyBDmedEver < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, CRP_as_afactor, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------SUI_Attempts--------------------------------------------------

#correlation of each gene to each condition, thus, gene significance
SUI_Attempts <- as.data.frame(config_df$SUI_Attempts)
names(SUI_Attempts) <- 'condition'
geneSignif <- corAndPvalue(datExpr, SUI_Attempts)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(SUI_Attempts), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(SUI_Attempts), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$SUI_Attempts < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, SUI_Attempts, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)

#-----------------MARS--------------------------------------------------

#correlation of each gene to each condition, thus, gene significance
MARS <- as.data.frame(config_df$MARS)
names(MARS) <- 'condition'
geneSignif <- corAndPvalue(datExpr, MARS)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(MARS), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(MARS), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$MARS < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, MARS, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)


#-----------------FamHx--------------------------------------------------

#correlation of each gene to each condition, thus, gene significance
FamHx <- as.data.frame(config_df$FamHx)
names(FamHx) <- 'condition'
geneSignif <- corAndPvalue(datExpr, FamHx)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(FamHx), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(FamHx), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$FamHx < 0.05)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
for (mod in sigMods){
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors==mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  xlab <- paste("Module Membership in", mod, "module")
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(filename = paste0(intramodFprefix, mod, '.tiff'), units="in", width=5, 
       height=5, res=300)
  verboseScatterplot(x, y, xlab = xlab, ylab = ylab, main = title,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = mod)
  invisible(dev.off())}

#Save the genes present in each significant module
for (mod in sigMods){
  modGenes <- (moduleColors==mod)
  modGenes <- symbol[modGenes]
  filename <- paste0(genelistsFprefix, mod, ".txt")
  write.table(rownames(as.data.frame(modGenes)), file = filename, quote =F, col.names = F
              , row.names=F)
  geneNumber <- length(modGenes)
  line <- paste0('Module ', mod, ': ', geneNumber)
  write(line, file = gNumberF, append = T)}

#Gene information data frame
geneInfo <- data.frame(geneSymbol = symbol, moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)
# Order modules by their significance for the condition
eigengeneTraitCorrelation <- WGCNA::cor(MEs, FamHx, method = "p"
                                        , use = "pairwise.complete.obs")
modOrder <- order(-abs(eigengeneTraitCorrelation))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo)
  geneInfo <- data.frame(geneInfo, 
                         geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  modMem <- paste0("MM.", moduleNames[modOrder[mod]])
  pModMem <- paste0("p.MM.", moduleNames[modOrder[mod]])
  names(geneInfo) <- c(oldNames, modMem,pModMem)
}

write.table(geneInfo, file = geneInfoF, sep = "\t", row.names = T, col.names = NA,
            quote = F)
#Save general summary results for farther experiments
save(geneInfo, sigMods, file = finalresR)















### E) ANNOTATION: Over Representation Analysis
#Make a data frame with each gene related to their assigned module
gene_module_df <- data.frame(rownames(geneInfo), geneInfo$moduleColor)
colnames(gene_module_df) <- c('gene', 'module')
#Retrieve data base
hs_msigdb_df <- msigdbr(species = 'Homo sapiens', collection = category, 
                        subcategory = subcategory)%>% 
  dplyr::select(gs_name, ensembl_gene)
#Make a background gene set with all possible genes, to compare later on
background_set <- unique(as.character(gene_module_df$gene))

#For each significant module make a over representation analysis
set.seed(13)
for (mod in moduleNames){
  #Make a list of genes in the module
  active_genes <- gene_module_df %>%
    dplyr::filter(module == mod) %>%
    dplyr::pull("gene")
  #Make analysis
  enriching <- enricher(gene = active_genes,
                        pvalueCutoff = 0.1,
                        pAdjustMethod = "BH",
                        universe = background_set,
                        TERM2GENE = hs_msigdb_df)
  enrich_res <- data.frame(enriching@result)
  write.table(enrich_res, file = paste0(annotationFprefix, mod, '.txt'),
              quote = F, row.names=F, col.names=T)
  enrich_res_sig <- enrich_res %>%
    dplyr::filter(p.adjust < alpha)
  #Just if there is any significant enriched term
  if(nrow(enrich_res_sig) != 0){
    write.table(enrich_res_sig, file = paste0(annotationFprefix, 'sig_', mod, '.txt'),
                quote = F, row.names=F, col.names=T)
    #Only make plots if there are more than 1 enriched terms
    if (nrow(enrich_res_sig) > 1){
      enrich_plot <- enrichplot::dotplot(enriching)
      #upset_plot <- enrichplot::upsetplot(enriching)
      
      filename <- paste0(annotationFprefix, 'dotplot_', mod, '.jpeg')
      ggsave(enrich_plot, file=filename, device = "jpeg", units= "in", 
             height = 10, width = 15)
      
      filename <- paste0(annotationFprefix, 'upsetplot_', mod, '.jpeg')
      #ggsave(upset_plot, file=filename, device = "jpeg", units= "in", 
      # height = 15, width = 20)
    }
  } else {
    message <- paste0('The module ', mod, 
                      ' doesnt have significantly enriched terms')
    print(message)
  }
  #To check that process is going OK
  message <- paste0('Finished with module ', mod)
  print(message)
}

