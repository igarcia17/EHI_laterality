library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(WGCNA, quietly = T)
library(org.Hs.eg.db, quietly = T)
library(clusterProfiler)
library(msigdbr)
library(openxlsx)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

tissue <- "Amygdala"
filesD <- paste0("WGCNA_brain/", tissue, "/")

inputData <- paste0(filesD,'data_step1_', tissue, "_2026_05_28.RData")
inputNet <- "automatic_network_minMod15"
inputNetF <- paste0(filesD, inputNet, '.RData')
load(inputData)
load(inputNetF)

#For output
subD <- paste0(inputNet, "/")
modtraitF <- paste0(filesD, subD, "module-trait_relation.jpeg")
intramodFprefix <- paste0(filesD, subD, 'MM_and_GS/intramodular_')
geneInfoF <- paste0(filesD, subD, 'genes_info.tsv')
finalresR <- paste0(filesD, subD, 'MM_and_GS/summary_sigMods.RData')
gNumberF <- paste0(filesD, subD, 'MM_and_GS/genes_in_significant_modules.txt')
annotationFprefix <- paste0(filesD, subD, 'enrich_results_')
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
levels(config_df$SEX) <- c(0, 1)
levels(config_df$STATUS) <- c(0,1)
colnames(config_df)[4] <- 'is_BD'
levels(config_df$Age_group) <- c(6, 1,2,3,4,5)
levels(config_df$Ethnicity) <- c(0,1)
colnames(config_df)[13] <- 'is_nonEU'
levels(config_df$Diagnostic) <- c(1, 2, 0,3,1) #cOMPLETELY MADE UP
levels(config_df$Lateral_Cutoff0)<- c(0,1)
levels(config_df$Lateral_Cutoff60)<- c(0,1)
levels(config_df$NRH_0)<- c(0,1)
levels(config_df$NRH_60)<- c(0,1)
config_df$WAVE <- NULL
config_df$RIN <- NULL
config_df$PRS_amb <- NULL
config_df$PRS_SCZ <- NULL
config_df$is_nonEU <- NULL
config_df$Age_group <- NULL
config_df$Lateral_Cutoff0 <- NULL
config_df$Lateral_Cutoff60 <- NULL

colnames(config_df)

##Exploratory analysis##ExplNULLoratory analysis
moduleTrait <- corAndPvalue(MEs, config_df[,3:dim(config_df)[2]])
moduleTraitCor <- moduleTrait$cor
moduleTraitPvalue <- moduleTrait$p

#Get all together:
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #give shape
#Graphical representation:
jpeg(file = modtraitF, width = 1150, height = 1400, quality = 100)
par(mar = c(6, 8.5, 3, 3))
title <- "Module-trait relationships"
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(config_df)[3:dim(config_df)[2]], xLabelsAngle = 45,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, main = title)
invisible(dev.off())

moduleNames <- substring(names(MEs), 3) #3 to take away prefix

wb <- createWorkbook()

for (name in names(moduleTrait)) {
  obj <- moduleTrait[[name]]
  if (is.vector(obj)) {
    obj <- data.frame(value = obj)
    rownames(obj) <- names(moduleTrait[[name]])
  }
  addWorksheet(wb, substr(name, 1, 31))
  writeData(wb,
            sheet = substr(name, 1, 31),
            x = obj,
            rowNames = TRUE)
}

saveWorkbook(wb,
             file = paste0(filesD, "moduleTrait.xlsx"),
             overwrite = TRUE)
#correlation of each gene in each sample to each module
modMemb <- corAndPvalue(datExpr, MEs)
geneModuleMembership <- as.data.frame(modMemb$cor)
names(geneModuleMembership) <- paste("MM", moduleNames, sep="_")

MMPvalue <- as.data.frame(modMemb$p)
names(MMPvalue) <- paste("p_MM", moduleNames, sep="_")

genes <- names(datExpr)
symbol <- mapIds(org.Hs.eg.db, keys = genes, column = 'SYMBOL', 
                 keytype = 'ENSEMBL', multiVals = 'first')



#correlation of each gene to each condition, thus, gene significance
STATUS <- as.data.frame(config_df$is_BD)
names(STATUS) <- 'condition'
geneSignif <- corAndPvalue(datExpr, STATUS)
geneTraitSignificance <- geneSignif$cor
geneTraitSignificance <- as.data.frame(geneTraitSignificance)
names(geneTraitSignificance) <- paste("GS.", names(STATUS), sep="")

GSPvalue <- geneSignif$p
GSPvalue <- as.data.frame(GSPvalue)
names(GSPvalue) <- paste("p.GS.", names(STATUS), sep="")

moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$is_BD < alpha)
sigMods <- substring(rownames(sigMods), 3)

#Plot significance of each gene in every significant module
# against their correlation to the module.
# Plot significance of each gene in every significant module
# against their correlation to the module.

for (mod in sigMods){
  
  column <- match(mod, moduleNames)
  moduleGenes <- moduleColors == mod
  
  x <- abs(geneModuleMembership[moduleGenes, column])
  y <- abs(geneTraitSignificance[moduleGenes, 1])
  
  # Skip problematic modules
  if (length(x) == 0 ||
      length(y) == 0 ||
      all(is.na(x)) ||
      all(is.na(y))) {
    
    message("Skipping module: ", mod,
            " (empty or invalid data)")
    next
  }
  
  xlab <- paste("Module Membership in", mod, "module")
  ylab <- "Gene significance for condition"
  title <- "Module membership vs. gene significance\n"
  
  tiff(
    filename = paste0(intramodFprefix, mod, '.tiff'),
    units = "in",
    width = 5,
    height = 5,
    res = 300
  )
  
  tryCatch({
    
    verboseScatterplot(
      x,
      y,
      xlab = xlab,
      ylab = ylab,
      main = title,
      cex.main = 1.2,
      cex.lab = 1.2,
      cex.axis = 1.2,
      col = mod
    )
    
  }, error = function(e) {
    
    message("Error plotting module ", mod, ": ", e$message)
    
  })
  
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
eigengeneTraitCorrelation <- WGCNA::cor(MEs, STATUS, method = "p"
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

wb <- createWorkbook()
wb_sig <- createWorkbook()

for (mod in moduleNames){
  active_genes <- gene_module_df %>%
    dplyr::filter(module == mod) %>%
    dplyr::pull("gene")
  set.seed(13)
  enriching <- enricher(gene = active_genes,
                        pvalueCutoff = 0.1,
                        pAdjustMethod = "BH",
                        universe = background_set,
                        TERM2GENE = hs_msigdb_df)
  sheet_name <- substr(as.character(mod), 1, 31)
  
  if (is.null(enriching)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name,
              x = data.frame(message = paste0("No gene can be mapped for module ", mod)))
    
    print(paste0("No enrichment result for module ", mod, ". Continuing..."))
    next
  }
  
  enrich_res <- data.frame(enriching@result)
  addWorksheet(wb, sheet_name)
  addWorksheet(wb_sig, sheet_name)
  enrich_res_sig <- enrich_res %>%
    dplyr::filter(p.adjust < alpha)
  writeData(wb, sheet = sheet_name, x = enrich_res)
  writeData(wb_sig, sheet = sheet_name, x = enrich_res_sig)
  if(nrow(enrich_res_sig) != 0){
    
    if (nrow(enrich_res_sig) > 1){
      enrich_plot <- enrichplot::dotplot(enriching)
      
      filename <- paste0(annotationFprefix, 'dotplot_', mod, '.jpeg')
      ggsave(enrich_plot, file = filename, device = "jpeg",
             units = "in", height = 10, width = 15)
    }
    
  } else {
    print(paste0('The module ', mod,
                 ' doesnt have significantly enriched terms'))
  }
  
  print(paste0('Finished with module ', mod))
}
saveWorkbook(wb,
             file = paste0(filesD, subD, "ORA_enrichment_all_modules.xlsx"),
             overwrite = TRUE)
saveWorkbook(wb_sig,
             file = paste0(filesD, subD, "ORA_sig_enrichment_all_modules.xlsx"),
             overwrite = TRUE)
