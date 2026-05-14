library(dplyr)
library(ggplot2)
library(tableone)
library(caret)
library(pheatmap)
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

association_matrix <- function(df) {
  vars <- names(df)
  n <- length(vars)
  mat <- matrix(NA, n, n, dimnames = list(vars, vars))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      x <- df[[i]]
      y <- df[[j]]
      
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
cols <- RColorBrewer::brewer.pal(n =9, name = "YlGnBu")

#Load input 15,646 indiv
input <- as.data.frame(readxl::read_xlsx("dawba_hand.xlsx", sheet = 1, 
                                         na="#N/A"))

mental_health <- colnames(input)[4:28]
handedness_items <- colnames(input)[29:35]
mental_factorial <- mental_health[-c(20, 1)]
handedness_factorial <- handedness_items[-7]
df <- input%>%
  mutate(across(-c(1, 2), as.numeric))
df <- df[!apply(df[, 4:35], 1, function(x) all(is.na(x))), ] ##11,142 samples after removing those with no info at all
df <- df[!apply(df[, 29:35], 1, function(x) all(is.na(x))), ] ##9,963 samples after removing those with no info for handedness
df <- df[!apply(df[, 4:28], 1, function(x) all(is.na(x))), ] ##8,820 samples after removing those with no info for mental health

df_factorial <- df %>%
  mutate(
    across(
      all_of(mental_factorial),
      ~ factor(
        .,
        levels = c(2, 1),
        labels = c("NO", "YES")
      )))%>%
  mutate(across( all_of(handedness_factorial), ~factor(
    ., levels = c(2,3,1),
    labels = c("Right", "Either", "Left")
  ) ))%>%
  select(where(~ {
    if (is.factor(.)) {
      all(table(.) >= 5)
    } else {
      TRUE
    }
  })) %>%mutate(
    Sex =factor(Sex, levels = c("Male", "Female")),
    Age = as.numeric(Age)
  )%>%droplevels()

#PCA of handedness
X <- df[,29:35]
X_clean <- X ##impute missing values by mean
for (i in seq_along(X_clean)) {
  X_clean[[i]][is.na(X_clean[[i]])] <- mean(X_clean[[i]], na.rm = TRUE)
}

pca <- prcomp(X_clean, center = TRUE, scale. = TRUE)
plot(pca, type = "l") #plot percentage of variances
pca$rotation 
scores <- as.data.frame(pca$x)
scores$group <- df$Sex

ggplot(scores, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_minimal()

cor_mat <- cor(X_clean, use = "pairwise.complete.obs", method="pearson")
pheatmap(
  cor_mat,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  display_numbers = TRUE,
  border_color = NA
)
handedness_cor <- association_matrix(df_factorial[,28:34])
pheatmap(handedness_cor*100,
         display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 12,
         color = cols[1:6], angle_col = 315,
         title = "Correlation (%)")
#dev.off()

mental_cor <- association_matrix(df_factorial[,4:27])
pheatmap(mental_cor*100,
         display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 12,
         color = cols[1:6], angle_col = 315,
         title = "Correlation (%)")
#dev.off()

all_cor <- association_matrix(df_factorial)
pheatmap(all_cor[2:27, 28:34]*100,
         display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 12,
         color = cols[1:6], angle_col = 315,
         title = "Correlation (%)")

temp <- mental_health[-c(1, 20)]

df_summary <- df_factorial %>%
  
  # 1. at least 85% of data
  filter(
    rowMeans(!is.na(.)) >= 0.85
  ) %>%
  
  # 2. Fno more than 3 NAs in mental health
  filter(
    rowSums(is.na(across(any_of(mental_health)))) <= 3
  ) %>%
  
  # 3. create any mental
  mutate(
    ANY_MENTAL = case_when(
      if_any(any_of(temp),
             ~ !is.na(.) & toupper(trimws(as.character(.))) == "YES") ~ "YES",
      TRUE ~ "NO"
    ),
    ANY_MENTAL = factor(ANY_MENTAL, levels = c("NO", "YES"))
  ) %>%
  select(-any_of(temp))
all_cor_superclean <- association_matrix(df_summary)
pheatmap(all_cor_superclean[c(2,3,4,5,13),c(6,7,8,9,10,11,12)]*100,
         display_numbers = T,
         cluster_rows = F, cluster_cols = F, fontsize_number = 12,
         color = cols[1:6], angle_col = 315,
         title = "Correlation (%)")


vars <- colnames(df_summary)[-c(1,13)]
categorical <- vars[-c(2,3,4,11)]
num <- vars[c(2:4,11)]

tabla1 <- CreateTableOne(
  vars = vars, 
  strata = "ANY_MENTAL", 
  data = df_summary, 
  factorVars = categorical,
  addOverall = TRUE  # Para que incluya la columna "Total"
)

# Mostrar la tabla con un formato limpio
t <- print(tabla1, 
           showAllLevels = TRUE, 
           nonnormal = num,
           formatOptions = list(big.mark = ","), pDigits=16
) 
write.csv(t, "Summary-ALSPAC.csv")
















