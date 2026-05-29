library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(WGCNA, quietly = T)
library(org.Hs.eg.db, quietly = T)
library(clusterProfiler)
library(msigdbr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

filesD <- 'WGCNA_blood/'
inputData <- paste0(filesD,'data_step1_2026_05_27.RData')
inputNet <- "automatic_network_minMod25"
inputNetF <- paste0(filesD, inputNet, '.RData')
load(inputData)
load(inputNetF)

#For output
subD <- paste0(inputNet, "_27_05_2026/")
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

##Exploratory analysis
moduleTrait <- corAndPvalue(MEs, config_df[,3:22])
moduleTraitCor <- moduleTrait$cor
moduleTraitPvalue <- moduleTrait$p

#Get all together:
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor) #give shape
#Graphical representation:
jpeg(file = modtraitF, width = 1050, height = 1300, quality = 100)
# par(mar = c(6, 8.5, 3, 3))
title <- "Module-trait relationships"
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(config_df)[3:22], xLabelsAngle = 45,
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
sigMods <- filter(moduleTraitPvalue, moduleTraitPvalue$NRH_0 < alpha)
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

