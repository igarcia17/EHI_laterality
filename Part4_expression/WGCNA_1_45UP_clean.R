library(dplyr)
library(readODS)
library(DESeq2)
library(WGCNA)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(variancePartition)
options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
suffix <- gsub("[^0-9]", "_", as.character(Sys.Date()))
outF <- "out_45nup/"
samtreeF <- paste0(outF, "sampleTree_", suffix, ".jpeg")
res_step1F <- paste0(outF, "data_step1_", suffix, ".RData")
configF <- "./cov_45nup/45nup_RNAseq_covar_postI_261125.ods"

config_df <- read_ods(configF, sheet = "Clean_Covar_RNA", col_names = TRUE, na="ND"
                      , strings_as_factors = T, skip=1)%>%
  filter(!grepl("CONT", Sample))%>%
  filter(!grepl("-EDTA", Sample))

sapply(config_df, class)
config_df$Filename <- as.character(config_df$Filename)
config_df$Sample <- as.character(config_df$Sample)
config_df$Notes <- NULL
levels(config_df$PHQ9_Category) <- c(3, 4, 1, 2, 0)
config_df$Notes <- NULL
config_df$Diagnosis_confidence <- NULL
config_df$Age_years <- NULL
config_df$eGFR_score <- NULL
config_df$`CRP_HNL(ref<=4.9)` <- NULL
config_df$AnyMeds_AntiInflam <- NULL
config_df$AnyMeds_Gastrointestinal <- NULL
config_df$AnyMeds_Hormone <- NULL
config_df$AnyMeds_Systemic <- NULL
config_df$AnyMeds_Psych <- NULL
config_df$Li_Response <- NULL
config_df$Val_Response <- NULL
config_df$Olanz_Response <- NULL
config_df$Quet_Response <- NULL
config_df$custom_erase <- NULL

#get gene names from first file
first_counts <- read.table(as.character(config_df$Filename[1]), header = TRUE, row.names = 1)
combined_counts <- data.frame(row.names = rownames(first_counts))

for (i in seq_len(nrow(config_df))) {
  sample_name <- as.character(config_df$Sample[i])
  file_path <- as.character(config_df$Filename[i])
  counts <- read.table(file_path, header = TRUE, row.names = 1)
  combined_counts[[sample_name]] <- counts[, 1]
}

#Filter low count genes, with less than 10 counts in 90% of features
#Originally 78903 genes

keep <- rowSums(combined_counts > 10) > (ncol(combined_counts) * 0.9)
df <- combined_counts[keep,] #drops to 16133
df <- head(df, -4)
rm(keep,combined_counts,counts, first_counts) #Remove because big size

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = config_df,
                              design = ~ 1) #to avoid design
dds <- estimateSizeFactors(dds)
vst <- vst(dds, blind = FALSE)
mat <- assay(vst)
datExpr <- as.data.frame(t(mat))

#Quality control
#First check possible sources of confounding
pca <- prcomp(datExpr, scale. = TRUE)
pc_scores <- as.data.frame(pca$x)
#For explained variance:
ve <- pca$sdev^2 / sum(pca$sdev^2)
vars <- colnames(config_df)[3:32]

cor_matrix <- sapply(vars, function(vname) {
  v <- config_df[[vname]]
  if (is.numeric(v)) apply(pc_scores, 2, function(pc) cor(pc, v, use="pairwise.complete.obs"))
  else apply(pc_scores, 2, function(pc) summary(lm(pc ~ v, data=config_df))$r.squared)
})

cor_df <- melt(t(cor_matrix))
colnames(cor_df) <- c("Variable", "PC", "Value")
cor_df <- cor_df %>%
  filter(grepl("^PC", PC)) %>%
  mutate(PCnum = as.numeric(gsub("PC", "", PC))) %>%
  filter(PCnum <= 10)

# Heatmap
ggplot(cor_df, aes(x = PC, y = Variable, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(abs(Value) > 0.2, sprintf("%.2f", Value), "")),
            color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "Variables vs PCs - no batch correction",
       x = "PC",
       y = "Feature",
       fill = "Cor/R²")

# Correlation between variables
association_matrix <- function(df) {
  vars <- names(df)
  n <- length(vars)
  mat <- matrix(NA, n, n, dimnames = list(vars, vars))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      x <- df[[i]]
      y <- df[[j]]
      print(x)
      # Ambos numéricos → correlación de Pearson
      if (is.numeric(x) && is.numeric(y)) {
        mat[i,j] <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
        
        # Uno numérico y otro factor → R² de un modelo lineal
      } else if (is.numeric(x) && is.factor(y)) {
        model <- lm(x ~ y)
        mat[i,j] <- summary(model)$r.squared
      } else if (is.factor(x) && is.numeric(y)) {
        model <- lm(y ~ x)
        mat[i,j] <- summary(model)$r.squared
        
        # Ambos factores → Cramér’s V
      } else if (is.factor(x) && is.factor(y)) {
        x <- droplevels(x)
        y <- droplevels(y)
        tbl <- table(x, y)
        mat[i,j] <- vcd::assocstats(tbl)$cramer
      }
    }
  }
  return(mat)
}
M <- association_matrix(config_df[,-c(1,2)])
cols <- brewer.pal(n = 7, name = "YlGnBu")
pheatmap(abs(M),
         display_numbers = T,
         cluster_rows = F,
         cluster_cols = F,
         fontsize_number = 10,
         color = cols,
         angle_col = 315,
         title = "Absolute correlation between variables")
#Cluster to detect outliers in samples
sampleTree <- hclust(dist(datExpr), method = 'average')
jpeg(file = samtreeF, width = 1700, height = 1700)
par(mar = c(0,4,2,0))
title <- "Sample clustering to detect outliers - no batch correction"
plot(sampleTree, main = title, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
invisible(dev.off())

#QC plots from dds object
var <- "Tube_type"
fname <- paste0("out_45nup/PCA_by_",var,"_",suffix,".jpeg")
jpeg(filename = fname, width=900, height=900)
pca <- DESeq2::plotPCA(vst, intgroup = var, pcsToUse=1:2)
title <- "Principal Components Plot"
pca + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst)), size=5, point.padding = 0.6)
invisible(dev.off())


#Reload data to remove batch effects
first_counts <- read.table(as.character(config_df$Filename[1]), header = TRUE, row.names = 1)
combined_counts <- data.frame(row.names = rownames(first_counts))
for (i in seq_len(nrow(config_df))) {
  sample_name <- as.character(config_df$Sample[i])
  file_path <- as.character(config_df$Filename[i])
  counts <- read.table(file_path, header = TRUE, row.names = 1)
  combined_counts[[sample_name]] <- counts[, 1]
}
keep <- rowSums(combined_counts > 10) > (ncol(combined_counts) * 0.9)
df <- combined_counts[keep,]
rm(keep,combined_counts,counts, first_counts) #Remove because big size
dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = config_df,
                              design = ~ 1) #to avoid design
dds <- estimateSizeFactors(dds)
vst <- vst(dds, blind = FALSE)
mat <- assay(vst)

#To remove batch effect of those variables with missing values, the NAs have to be imputed.
des <- matrix(1, nrow = ncol(mat), ncol = 1)
config_df <- config_df %>%
  mutate(across(c(Arrival2RNAproc_hours_new, Collection2Arrival_minutes, TimeAtRT_3hoursbin),
                ~ if_else(is.na(.), median(., na.rm = TRUE), .)))
cov_mat <- as.matrix(config_df[,c(6,8,5)])
mat <- limma::removeBatchEffect(mat, batch = vst$RNAprocessingBatch
                                , covariates = cov_mat
                                , design = des)
assay(vst) <- mat
datExpr <- as.data.frame(t(mat))

#Quality control
pca <- prcomp(datExpr, scale. = TRUE)
pc_scores <- as.data.frame(pca$x)
ve <- pca$sdev^2 / sum(pca$sdev^2)
vars <- colnames(config_df)[3:31]
cor_matrix <- sapply(vars, function(vname) {
  v <- config_df[[vname]]
  if (is.numeric(v)) apply(pc_scores, 2, function(pc) cor(pc, v, use="pairwise.complete.obs"))
  else apply(pc_scores, 2, function(pc) summary(lm(pc ~ v, data=config_df))$r.squared)
})

cor_df <- melt(t(cor_matrix))
colnames(cor_df) <- c("Variable", "PC", "Value")
cor_df <- cor_df %>%
  filter(grepl("^PC", PC)) %>%
  mutate(PCnum = as.numeric(gsub("PC", "", PC))) %>%
  filter(PCnum <= 10)

# Heatmap
ggplot(cor_df, aes(x = PC, y = Variable, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(abs(Value) > 0.2, sprintf("%.2f", Value), "")),
            color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "Variables vs PCs",
       x = "PC",
       y = "Feature",
       fill = "Cor/R²")

form <- ~  Gender + TimeAtRT_3hoursbin + Tube_type  +
  Arrival2RNAproc_hours_new + Collection2Arrival_minutes
vobj <- fitExtractVarPartModel(t(datExpr), form, config_df)
vp <- sort(colMeans(vobj))
df_qc <- data.frame(
  variable = names(vp),
  media = round(vp, 3)
)

ggplot(df_qc, aes(x = variable, y = 1, fill = media)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = media), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) 

#Cluster to detect outliers in samples
sampleTree <- hclust(dist(datExpr), method = 'average')
jpeg(file = samtreeF, width = 1700, height = 1700)
par(mar = c(0,4,2,0))
title <- "Sample clustering to detect outliers"
plot(sampleTree, main = title, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
invisible(dev.off())

#QC plots from dds object
var <- "Tube_type"
fname <- paste0("out_45nup/PCA_by_",var,"_",suffix,".jpeg")
jpeg(filename = fname, width=900, height=900)
pca <- DESeq2::plotPCA(vst, intgroup = var, pcsToUse = c(1,2))
title <- "Principal Components Plot"
pca + ggtitle(title) + 
  geom_point(size = 6) +
  theme(plot.title = element_text(size=40, hjust = 0.5, face = "bold"), axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=15)) +
  geom_text_repel(aes(label=colnames(vst)), size=5, point.padding = 0.6)
invisible(dev.off())
#-------------Remove outliers and save definitive data
save(datExpr, config_df, file = res_step1F)
