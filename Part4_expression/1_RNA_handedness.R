library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(edgeR)
library(limma)
library(splines)
library('org.Hs.eg.db')
library(EnhancedVolcano, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(msigdbr, quietly = T)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())
########################Parameters
cutoffs <- c("0", "60", "90")
oldies <- T
filesF <- "RNAseq-counts_blood/"
##########Functions

########## Input

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

sapply(input, class)

if(oldies){
  input <- input %>%
    filter(Oldies_dataset == "YES") %>% 
    filter(has_RNA=="YES")%>%droplevels()
}else{
  input <- input%>%filter(has_RNA=="YES")%>%droplevels()
}

#Load genetic expression data
files <- input$File
files <- paste0(filesF, input$File)

expr_list <- lapply(files, function(f) {
  df <- read.table(f,
                   header = F,
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
  
  counts <- df[, 2]
  names(counts) <- df[, 1]
  
  return(counts)
})

expr_matrix <- do.call(cbind, expr_list)
colnames(expr_matrix) <- input$ID
rm(expr_list)
gc()

expr_matrix <- expr_matrix[-c(63141:63145),]
#Ensrue metadata is correctly ordered
all(colnames(expr_matrix) %in% input$ID)
meta <- input[match(colnames(expr_matrix), input$ID), ]
meta <- meta%>%
  select(-c("has_RNA", 'PRS_handedness', 'PRS_amb', 'PRS_BD', 'PRS_SCZ',
            'Oldies_dataset', 'Cutoff0', 'Cutoff40', 'Cutoff60', 'Cutoff80', 
            'Cutoff90', 'PC1_all', 'PC1_oldies', 'PC2_all', 'PC2_oldies', 
            'PC3_all', 'PC3_oldies','PC4_all', 'PC4_oldies'))

#Make a dge object of EdgeR + filter
dge <- DGEList(counts = expr_matrix)

keep <- filterByExpr(dge, group = input$NRH_0, min.count=20, min.total.count=25,
                     min.prop=0.5) #Minimum 20 counts in 50% of NRH
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge, method = "TMM")
#Model design
design <- model.matrix(~ NRH_0 + SEX + ns(Age, df=3) # permite modelizacion polinomica de la edad
                       + RIN + Ethnicity + WAVE, data = meta)
v <- voom(dge, design, plot = TRUE)

#Fit model
set.seed(13)
fit <- lmFit(v, design)
fit <- eBayes(fit)

#To check what comparison to adress
colnames(fit$coefficients)

results <- topTable(fit, coef = "NRH_0NRH", number = Inf)
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(results), column="SYMBOL", 
                 keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Hs.eg.db'), keys=row.names(results),
                      column="GENENAME", keytype="ENSEMBL", 
                      multiVals="first") #to obtain description
results <- cbind(symbol, results)
results$description <- description
results$Gene <- rownames(results)
results <- results[, c("Gene", setdiff(colnames(results), "Gene"))]

writexl::write_xlsx(results, path="DEGs_blood_RNA_NRH_0.xlsx",col_names = T)


jpeg(filename = paste0("Volcano_NRH.jpeg"), units="in", width=8, height=12, res=300)
EnhancedVolcano(results, lab = results$symbol, x = 'logFC', y = 'P.Value',
                pCutoff = 0.05, FCcutoff= 0.3, 
                #pCutOff is p value for last significant acc to adjp value
                ylim = c(0, 11), xlim = c(-0.6, 0.6), labSize = 3,
                legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 80, title = '', arrowheads = FALSE,
                subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
invisible(dev.off())

#Perform GSEA
db_sets <- msigdbr(species = 'Homo sapiens', collection = "C5", 
                   subcategory = NULL)%>% 
  dplyr::select(gs_name, ensembl_gene)
results <- results[complete.cases(results), ]
geneList <- results$t
names(geneList) <- rownames(results)
geneList <- sort(geneList, decreasing = TRUE)

set.seed(13)
egs <- GSEA(geneList = geneList, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

egs_df$col <- rownames(egs_df)
egs_df <- egs_df[, c("col", setdiff(colnames(egs_df), "col"))]

writexl::write_xlsx(egs_df, path="DEGs_blood_RNA_NRH_0_GSEA.xlsx",col_names = T)

#####Now only for BD
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

sapply(input, class)

#Filter only BD
if(oldies){
  input <- input %>%
    filter(Oldies_dataset == "YES") %>% 
    filter(has_RNA=="YES")%>%filter(STATUS=="YES")%>%droplevels()
}else{
  input <- input%>%filter(has_RNA=="YES")%>%filter(STATUS=="YES")%>%droplevels()
}

#Load genetic expression data
files <- input$File
files <- paste0(filesF, input$File)

expr_list <- lapply(files, function(f) {
  df <- read.table(f,
                   header = F,
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
  
  counts <- df[, 2]
  names(counts) <- df[, 1]
  
  return(counts)
})

expr_matrix <- do.call(cbind, expr_list)
colnames(expr_matrix) <- input$ID
rm(expr_list)
gc()

expr_matrix <- expr_matrix[-c(63141:63145),]

#Ensrue metadata is correctly ordered
all(colnames(expr_matrix) %in% input$ID)
meta <- input[match(colnames(expr_matrix), input$ID), ]
meta <- meta%>%
  select(-c("has_RNA", 'PRS_handedness', 'PRS_amb', 'PRS_BD', 'PRS_SCZ',
            'Oldies_dataset', 'Cutoff0', 'Cutoff40', 'Cutoff60', 'Cutoff80', 
            'Cutoff90', 'PC1_all', 'PC1_oldies', 'PC2_all', 'PC2_oldies', 
            'PC3_all', 'PC3_oldies','PC4_all', 'PC4_oldies'))

#Make a dge object of EdgeR + filter
dge <- DGEList(counts = expr_matrix)

keep <- filterByExpr(dge, group = input$NRH_60, min.count=20, min.total.count=25,
                     min.prop=0.5) #Minimum 20 counts in 50% of NRH ###17694
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge, method = "TMM")
#Model design
design <- model.matrix(~ NRH_60 + SEX + ns(Age, df=3) # permite modelizacion polinomica de la edad
                       + RIN + Ethnicity, data = meta)
v <- voom(dge, design, plot = TRUE)

#Fit model
set.seed(13)
fit <- lmFit(v, design)
fit <- eBayes(fit)

#To check what comparison to adress
colnames(fit$coefficients)

results <- topTable(fit, coef = "NRH_60NRH", number = Inf)
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(results), column="SYMBOL", 
                 keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Hs.eg.db'), keys=row.names(results),
                      column="GENENAME", keytype="ENSEMBL", 
                      multiVals="first") #to obtain description
results <- cbind(symbol, results)
results$description <- description
results$Gene <- rownames(results)
results <- results[, c("Gene", setdiff(colnames(results), "Gene"))]

writexl::write_xlsx(results, path="DEGs_blood_RNA_NRH_60_onlyBD.xlsx",col_names = T)


jpeg(filename = paste0("Volcano_NRH60_onlyBD.jpeg"), units="in", width=8, height=12, res=300)
EnhancedVolcano(results, lab = results$symbol, x = 'logFC', y = 'P.Value',
                pCutoff = 0.05, FCcutoff= 0.3, 
                #pCutOff is p value for last significant acc to adjp value
                ylim = c(0, 11), xlim = c(-0.6, 0.6), labSize = 3,
                legendLabSize = 9, legendIconSize = 5, drawConnectors = TRUE,
                widthConnectors = 0.5, max.overlaps = 80, title = '', arrowheads = FALSE,
                subtitle= '', gridlines.major = FALSE, gridlines.minor = FALSE)
invisible(dev.off())

#Perform GSEA
db_sets <- msigdbr(species = 'Homo sapiens', collection = "C5", 
                   subcategory = NULL)%>% 
  dplyr::select(gs_name, ensembl_gene)
results <- results[complete.cases(results), ]
geneList <- results$t
names(geneList) <- rownames(results)
geneList <- sort(geneList, decreasing = TRUE)

set.seed(13)
egs <- GSEA(geneList = geneList, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

egs_df$col <- rownames(egs_df)
egs_df <- egs_df[, c("col", setdiff(colnames(egs_df), "col"))]

writexl::write_xlsx(egs_df, path="DEGs_blood_RNA_NRH_60_GSEA_onlyBD.xlsx",col_names = T)


#
#
#
##
##
##
##
##
##
##
##
#
