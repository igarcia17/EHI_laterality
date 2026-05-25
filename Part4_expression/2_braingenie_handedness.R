workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())
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

filter_genes_blood <- T
tissues <- c("Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia",
             "Cerebellar_Hemisphere", "Cerebellum", "Cortex", "Frontal_Cortex_BA9",
             "Hippocampus", "Hypothalamus", 
             "Nucleus_accumbens_basal_ganglia",
             "Putamen_basal_ganglia","Spinal_cord_cervical_c-1",
             "Substantia_nigra")

genes_to_include <- as.data.frame(data.table::fread("all_genes.tsv"))[,1]
if (!file.exists("genes_to_include.RData")) {
  save(genes_to_include, file = "genes_to_include.RData")}

for(tissue in tissues){
  
  if (!dir.exists(paste0("SA_out_",tissue))) {
    dir.create(paste0("SA_out_",tissue))
  }
  
inputF <- paste0("Brain_", tissue, "_merged_counts.txt")
outF <- paste0("SA_out_",tissue,"/all_results_", tissue, ".txt")
out_005_F <- paste0("SA_out_",tissue,"/degs_adjP_005_", tissue, ".txt")
out_001_F <- paste0("SA_out_",tissue,"/degs_adjP_001_", tissue, ".txt")
gsea_out_F <- paste0("SA_out_",tissue,"/gsea_", tissue, ".txt")

counts <- as.data.frame(data.table::fread(inputF))
meta <- as.data.frame(readxl::read_xlsx("config_file_DGE_MadManic.xlsx"))%>%
  mutate(across(c(Sex, Case_Ctrl_SA, Ethnicity,
                  Case_Ctrl_BD,  Wave, lt_exposed, est_exposed, 
                  bzd_exposed, ap_exposed,ad_exposed), as.factor))%>%
  mutate(Sex = relevel(Sex, ref = "Male"),
         Ethnicity =relevel(Ethnicity, ref="EUR"),
         Case_Ctrl_BD=relevel(Case_Ctrl_BD, ref = "Control"),
         Case_Ctrl_SA=relevel(Case_Ctrl_SA, ref= "Control")
                  )
cell_prop <- as.data.frame(data.table::fread("CIBERSORTx_725_indiv_RNAseq_MadManic.txt"))
cell_prop$Mixture <- cell_prop$Mixture |>
  gsub(".tsv", "", x= _)

if(filter_genes_blood){
counts <- counts %>%
  dplyr::filter(Gene %in% genes_to_include)}

meta$Sample_name <- meta$Sample_name |>
  gsub("CSI", "CSIC", x = _) |>
  gsub("_", "R_", x = _)

samples_selected <- meta$Sample_name
meta <- meta %>%
  left_join(cell_prop, by = c("Sample_name" = "Mixture"))

if (!file.exists("meta.RData")) {
  save(meta, file = "meta.RData")}

counts <- counts%>%
  dplyr::select(Gene, all_of(samples_selected))

# Set counts as matrix
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- as.matrix(counts)

#Order metadata
all(colnames(counts) %in% meta$Sample_name)
meta <- meta[match(colnames(counts), meta$Sample_name), ]

#Model design
design <- model.matrix(~ Case_Ctrl_SA + Sex + ns(Age, df=3) # permite modelizacion polinomica de la edad
                       + RIN + Ethnicity 
                       #+ PC1_CIBERSORT + PC2_CIBERSORT + PC3_CIBERSORT + PC4_CIBERSORT + PC5_CIBERSORT
                       + lt_exposed + bzd_exposed + ap_exposed
                       + ad_exposed + est_exposed, data = meta)
#Tambien añadir cibersortX and SVA

#Fit model
set.seed(1)
fit <- lmFit(counts, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "Case_Ctrl_SACase", number = Inf)
symbol <- mapIds(get('org.Hs.eg.db'), keys=row.names(results), column="SYMBOL", 
                 keytype="ENSEMBL", multiVals="first") #to obtain gene symbols

description <- mapIds(get('org.Hs.eg.db'), keys=row.names(results),
                      column="GENENAME", keytype="ENSEMBL", 
                      multiVals="first") #to obtain description
results <- cbind(symbol, results)
results$description <- description

deg_005 <- results %>%
  dplyr::filter(adj.P.Val < 0.05)
deg_001 <- results %>%
  dplyr::filter(adj.P.Val < 0.01)
write.table(results, file=outF, quote=FALSE, sep="\t", col.names=NA)
write.table(deg_005, file=out_005_F, quote=FALSE, sep="\t", col.names=NA)
write.table(deg_001, file=out_001_F, quote=FALSE, sep="\t", col.names=NA)

jpeg(filename = paste0("SA_out_",tissue,"/MA_plot_",tissue,".jpeg"), units="in", width=8, height=10, res=300)
sig <- results$adj.P.Val < 0.01
plot(results$AveExpr, results$logFC,
     pch = 16, cex = 0.5, col = "grey30",
     xlab = "Average Expression",
     ylab = "Log2 Fold Change",
     title = paste0("MA plot ", tissue))

points(results$AveExpr[sig], results$logFC[sig],
       col = "red", pch = 16, cex = 0.5)
abline(h = 0, col = "red", lwd = 2)
invisible(dev.off())


jpeg(filename = paste0("SA_out_",tissue,"/Volcano_",tissue,".jpeg"), units="in", width=8, height=12, res=300)
EnhancedVolcano(results, lab = results$symbol, x = 'logFC', y = 'P.Value',
                pCutoff = 0.000002, FCcutoff= 0.3, 
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

set.seed(1)
egs <- GSEA(geneList = geneList, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

write.table(egs_df, file = gsea_out_F, sep= "\t", quote = F, row.names = T)
}








################################################################################
### Codigo cuando relacion de una covariable con el cambio no es lineal
design <- model.matrix(~ sexo + ns(edad_num, df = 3), data = meta)
###



jpeg(filename = PCAF, width=900, height=900, quality=300)
pca <- DESeq2::plotPCA(counts)

title <- "Principal Components Plot"

pca + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst)), size=5, point.padding = 0.6)
invisible(dev.off())
