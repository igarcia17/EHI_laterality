#Brain-Regional Gene Expression Imputed from the Blood Transcriptome by BrainGENIE Recapitulates Dysregulation Observed in the Postmortem Brain in Alzheimer’s Disease
#https://www.medrxiv.org/content/10.1101/2025.10.13.25337831v1
library(dplyr)
library(limma)
library(splines)
library(edgeR)
library('org.Hs.eg.db')
library(EnhancedVolcano, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(msigdbr, quietly = T)
library(tidyverse)
library(ggplot2)
library(tidyr)
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

cutoffs <- c("0", "60", "90")
oldies <- F
tissues <- c("Spinal_cord_cervical_c-1")
db_sets <- msigdbr(species = 'Homo sapiens', collection = "C5", 
                   subcategory = NULL)%>% 
  dplyr::select(gs_name, ensembl_gene)

input <- as.data.frame(readxl::read_xlsx("All_data_available_genetic_exp.xlsx"
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
  input <- input %>%
    filter(Oldies_dataset == "YES") %>% 
    filter(has_RNA=="YES")%>%droplevels()
}else{
  input <- input%>%filter(has_RNA=="YES")%>%droplevels()
}

#Take cellular proportions into account
cell_prop <- as.data.frame(data.table::fread("CIBERSORTx_725_indiv_RNAseq_MadManic.txt"))
cells <- colnames(cell_prop)[!colnames(cell_prop) %in% c("Mixture", "RMSE", "P-value", "Correlation")]

input <- input %>%
  mutate(row_order_original = row_number())

input <- input %>%
  left_join(cell_prop, by = c("File" = "Mixture"))
stopifnot(all(input$row_order_original == seq_len(nrow(input))))

cell_mat <- input[, cells, drop = FALSE]
cell_mat <- as.data.frame(lapply(cell_mat, as.numeric))
cell_mat <- as.data.frame(lapply(cell_mat, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

set.seed(1)
pca <- prcomp(cell_mat,
              center = TRUE,
              scale. = TRUE)
cell_pcs <- as.data.frame(pca$x[, 1:2])
colnames(cell_pcs) <- c("CellPC1", "CellPC2")
input <- bind_cols(input, cell_pcs)


for(tissue in tissues){
  print(paste0("Starting with region ", tissue))
  if (!dir.exists(paste0("Brain_results_DGE_GSEA/Brain_BD_BDvscontrols_",tissue))) {
    dir.create(paste0("Brain_results_DGE_GSEA/Brain_BD_BDvscontrols_",tissue))
  }

  #To load counts
inputF <- paste0("Imputed_brain/Brain_", tissue, "/")
files <- sub("\\.tsv$", "", input$File)
# IDs desde input$File, quitando .tsv, porque tienen esa R
sample_ids <- sub("\\.tsv$", "", input$File)

expr_list <- lapply(sample_ids, function(sample_id) {
  matched_file <- list.files(
    path = inputF,
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
colnames(expr_matrix) <- input$ID
rm(expr_list)
gc()

print(paste0("Loaded counts of ", tissue))
pathF <- paste0("Brain_results_DGE_GSEA/Brain_BD_BDvscontrols_",tissue)
outF <- paste0(pathF, "/all_results_", tissue, ".txt")
out_005_F <- paste0(pathF,"/degs_adjP_005_", tissue, ".txt")
gsea_out_F <- paste0(pathF,"/gsea_", tissue, ".txt")


#Model design
design <- model.matrix(~ NRH_0 + SEX + ns(Age, df=3) # permite modelizacion polinomica de la edad
                       + RIN + Ethnicity + WAVE + CellPC1 + CellPC2 + STATUS
                       , data = input)

#Fit model
set.seed(1)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "STATUSYES", number = Inf)
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(results), column="SYMBOL", 
                 keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Hs.eg.db'), keys=row.names(results),
                      column="GENENAME", keytype="ENSEMBL", 
                      multiVals="first") #to obtain description
results <- cbind(symbol, results)
results$description <- description

deg_005 <- results %>%
  dplyr::filter(adj.P.Val < 0.05)
write.table(results, file=outF, quote=FALSE, sep="\t", col.names=NA)
write.table(deg_005, file=out_005_F, quote=FALSE, sep="\t", col.names=NA)

jpeg(filename = paste0(pathF,"/Volcano_",tissue,".jpeg"), units="in", width=8, height=12, res=300)
EnhancedVolcano(results, lab = results$symbol, x = 'logFC', y = 'P.Value',
                pCutoff = 0.0000001, FCcutoff= 0.1, 
                #pCutOff is p value for last significant acc to adjp value
                ylim = c(0, 11), xlim = c(-0.6, 0.6), labSize = 3,
                legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 80, title = '', arrowheads = FALSE,
                subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
invisible(dev.off())

#Perform GSEA
results <- results[complete.cases(results), ]
geneList <- results$t
names(geneList) <- rownames(results)
geneList <- sort(geneList, decreasing = TRUE)

set.seed(1)
egs <- GSEA(geneList = geneList, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

write.table(egs_df, file = gsea_out_F, sep= "\t", quote = F, row.names = T)
}


###Only in BD

input <- as.data.frame(readxl::read_xlsx("All_data_available_genetic_exp.xlsx"
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
  input <- input %>%
    filter(Oldies_dataset == "YES") %>%filter(STATUS =="YES")%>% 
    filter(has_RNA=="YES")%>%droplevels()
}else{
  input <- input%>%filter(has_RNA=="YES")%>%filter(STATUS =="YES")%>%droplevels()
}

#Take cellular proportions into account
cell_prop <- as.data.frame(data.table::fread("CIBERSORTx_725_indiv_RNAseq_MadManic.txt"))
cells <- colnames(cell_prop)[!colnames(cell_prop) %in% c("Mixture", "RMSE", "P-value", "Correlation")]

input <- input %>%
  mutate(row_order_original = row_number())

input <- input %>%
  left_join(cell_prop, by = c("File" = "Mixture"))
stopifnot(all(input$row_order_original == seq_len(nrow(input))))

cell_mat <- input[, cells, drop = FALSE]
cell_mat <- as.data.frame(lapply(cell_mat, as.numeric))
cell_mat <- as.data.frame(lapply(cell_mat, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

set.seed(1)
pca <- prcomp(cell_mat,
              center = TRUE,
              scale. = TRUE)
cell_pcs <- as.data.frame(pca$x[, 1:2])
colnames(cell_pcs) <- c("CellPC1", "CellPC2")
input <- bind_cols(input, cell_pcs)


for(tissue in tissues){
  
  if (!dir.exists(paste0("Brain_results_DGE_GSEA/Brain_NRH60_onlyBD_",tissue))) {
    dir.create(paste0("Brain_results_DGE_GSEA/Brain_NRH60_onlyBD_",tissue))
  }
  
  #To load counts
  inputF <- paste0("Imputed_brain/Brain_", tissue, "/")
  files <- sub("\\.tsv$", "", input$File)
  # IDs desde input$File, quitando .tsv, porque tienen esa R
  sample_ids <- sub("\\.tsv$", "", input$File)
  
  expr_list <- lapply(sample_ids, function(sample_id) {
    matched_file <- list.files(
      path = inputF,
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
  colnames(expr_matrix) <- input$ID
  rm(expr_list)
  gc()
  
  
  pathF <- paste0("Brain_results_DGE_GSEA/Brain_NRH60_onlyBD_",tissue)
  outF <- paste0(pathF, "/all_onlyBD_results_", tissue, ".txt")
  out_005_F <- paste0(pathF,"/degs_onlyBD_adjP_005_", tissue, ".txt")
  gsea_out_F <- paste0(pathF,"/gsea_onlyBD_", tissue, ".txt")
  
  
  #Model design
  design <- model.matrix(~ NRH_60 + SEX + ns(Age, df=3) # permite modelizacion polinomica de la edad
                         + RIN + Ethnicity + WAVE +CellPC1+CellPC2
                         , data = input)
  
  #Fit model
  set.seed(1)
  fit <- lmFit(expr_matrix, design)
  fit <- eBayes(fit)
  
  results <- topTable(fit, coef = "NRH_60NRH", number = Inf)
  symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(results), column="SYMBOL", 
                   keytype="ENSEMBL", multiVals="first") #to obtain gene symbols
  
  description <- mapIds(get('org.Hs.eg.db'), keys=row.names(results),
                        column="GENENAME", keytype="ENSEMBL", 
                        multiVals="first") #to obtain description
  results <- cbind(symbol, results)
  results$description <- description
  
  deg_005 <- results %>%
    dplyr::filter(adj.P.Val < 0.05)
  write.table(results, file=outF, quote=FALSE, sep="\t", col.names=NA)
  write.table(deg_005, file=out_005_F, quote=FALSE, sep="\t", col.names=NA)
  
  jpeg(filename = paste0(pathF,"/Volcano_onlyBD_",tissue,".jpeg"), units="in", width=8, height=12, res=300)
  EnhancedVolcano(results, lab = results$symbol, x = 'logFC', y = 'P.Value',
                  pCutoff = 0.000002, FCcutoff= 0.3, 
                  #pCutOff is p value for last significant acc to adjp value
                  ylim = c(0, 11), xlim = c(-0.6, 0.6), labSize = 3,
                  legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                  widthConnectors = 0.5, max.overlaps = 80, title = '', arrowheads = FALSE,
                  subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
  invisible(dev.off())
  
  #Perform GSEA
  results <- results[complete.cases(results), ]
  geneList <- results$t
  names(geneList) <- rownames(results)
  geneList <- sort(geneList, decreasing = TRUE)
  
  set.seed(1)
  egs <- GSEA(geneList = geneList, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
              seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
  egs_df <- data.frame(egs@result)
  egs_df <- egs_df[, -c(1,2)]
  
  write.table(egs_df, file = gsea_out_F, sep= "\t", quote = F, row.names = T)
}






