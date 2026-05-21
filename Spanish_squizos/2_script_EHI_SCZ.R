library(dplyr)
library(tidyr)
library(effsize)
library(psych)
library(psychTools)
library(ggplot2)
library(forcats)
library(mvtnorm)
library(class)
library(VIM)
library(pheatmap)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

##### Functions
calculate_LQ_EHI <- function(df){
  df_temp <- df %>%
    mutate(
      Right10 = rowSums(
        across(Item1:Item10, ~ case_when(
          . == 1 ~ 2,
          . %in% c(2,3) ~ 1,
          . %in% c(4, 5) ~ 0,
          TRUE ~ NA_real_
        )),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      Left10 = rowSums(
        across(Item1:Item10, ~ case_when(
          . == 5 ~ 2,
          . %in% c(4,3) ~ 1,
          . %in% c(1,2) ~ 0,
          TRUE ~ NA_real_
        )),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      Score10items = round(((Right10 - Left10) / (Right10 + Left10) *100),0)
    ) %>%select(ID, Score10items)
  return(df_temp)
}
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
###### Load input



casesF <- "../Part1_input_prep/cases_EHI_190526.xlsx"
controlsF <- "../Part1_input_prep/controls_clean_EHI_190526.xlsx"
sczF <- "BNADN_EHI_SCZ_210526.xlsx"

cases <- as.data.frame(readxl::read_xlsx(casesF, sheet = 1, 
                                              col_types="text", na="#N/A"))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  ))

ehi_cases <-  as.data.frame(readxl::read_xlsx(casesF, sheet = 2, 
                                              col_types="text", na="#N/A"))[,1:13]

controls <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 1, col_types="text", na="#N/A"))%>%
  filter(!(ID %in% cases$ID))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  ))%>%
  distinct(ID, .keep_all = TRUE)
ehi_controls <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 2, na="#N/A"))[,1:13]%>%
  distinct(ID, .keep_all = TRUE)

scz_df <- as.data.frame(readxl::read_xlsx(sczF, sheet = 2, col_types="text", na="#N/A"))%>%
  filter(!(ID %in% cases$ID))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  ))%>%
  distinct(ID, .keep_all = TRUE)
ehi_scz <- as.data.frame(readxl::read_xlsx(sczF, sheet = 4))%>%
  distinct(ID, .keep_all = TRUE)

cols <- c("ID", "Item1", "Item2", "Item3", "Item4", "Item5", "Item6", "Item7",
          "Item8", "Item9", "Item10", "Item11", "Item12")
colnames(ehi_cases) <- cols
colnames(ehi_controls) <- cols
colnames(ehi_scz) <- cols[1:11]


all_ehi <- rbind(ehi_cases,ehi_controls)%>%
  distinct(ID, .keep_all = TRUE)%>%
  mutate(across(starts_with("Item"), as.numeric))
all_ehi<-all_ehi[,1:11]
all_ehi <- rbind(all_ehi,ehi_scz)%>%
  distinct(ID, .keep_all = TRUE)%>%
  mutate(across(starts_with("Item"), as.numeric))

#Missing values imputation with KNN
complete_data <- all_ehi[complete.cases(all_ehi),]
set.seed(13)
prop_missing <- 0.2  
mask <- matrix(runif(nrow(complete_data)*ncol(complete_data))<prop_missing,nrow=nrow(complete_data))
data_missing <- complete_data
data_missing[mask] <- NA

ks <- seq(2,25)
rmse <- numeric(length(ks))

for (i in seq_along(ks)){
  k_val <- ks[i]
  imp <- kNN(data_missing,k=k_val,imp_var=FALSE)
  
  mat_complete <- as.matrix(complete_data)
  mat_imp <- as.matrix(imp)
  true_vals <- as.numeric(mat_complete[mask])
  pred_vals <- as.numeric(mat_imp[mask])
  rmse[i] <- sqrt(mean((true_vals-pred_vals)^2,na.rm=TRUE))
}

plot(ks,rmse,type="b",pch=1,
     xlab="k",ylab="Imputation RMSE",
     main = "Selecting K via RMSE")
best_index <- which.min(rmse)
best_k <- ks[best_index] #4 es la mejor k con 0.2 missing values  

#imputation
all_ehi_imputed <- kNN(all_ehi,k=best_k,imp_var=FALSE)

#### Data integration and LQ calculus
all_samples <- rbind(cases, controls)%>%
  distinct(ID, .keep_all = TRUE)%>%
  mutate(BD_patient = BD_patient %>% 
           fct_recode("YES" = "SI") %>% 
           fct_relevel("NO"))

all_samples <- rbind(all_samples, scz_df)%>%
  distinct(ID, .keep_all = TRUE)%>%
  mutate(BD_patient = BD_patient %>% 
           fct_relevel("NO"))

all_LQ <- calculate_LQ_EHI(all_ehi_imputed)%>%
  left_join(all_samples %>% select(ID, Sex, Age, BD_patient, Diagnostic), by = "ID")%>%
  mutate(Age = as.numeric(Age),
         Sex = factor(
           if_else(Sex == "HOMBRE", "HOMBRE",
                   if_else(Sex == "MUJER", "MUJER", NA_character_)),
           levels = c("HOMBRE", "MUJER")),
         Diagnostic = factor(Diagnostic)
         )%>%
  mutate(
    Age = if_else(Age < 10 | Age > 1000, NA_real_, Age)
  )
writexl::write_xlsx(all_LQ, path= "all_samples_all_data_LQ_included_SCZ.xlsx")

import_df <- all_LQ%>%
  left_join(all_ehi_imputed, by ="ID")
writexl::write_xlsx(import_df, path= "all_samples_all_data_LQ_PLUS_items_included_SCZ.xlsx")


#All LQ scores
ggplot(all_LQ, aes(x = Score10items)) +
  geom_density(alpha = 0.5) + labs(title = "EHI LQ distribution") +
  theme_minimal()
#Differing BD status
ggplot(all_LQ, aes(x = Score10items, fill = BD_patient)) +
  geom_density(alpha = 0.5) +labs(title = "EHI LQ distribution by BD/no BD") +
  theme_minimal()
#Histogram
ggplot(all_LQ, aes(x = Score10items)) +
  geom_histogram(
    binwidth = 10,
    boundary = -100,closed = "right",
    fill = "cadetblue",
    color = "black"
  ) +
  geom_text(
    stat = "bin",
    binwidth = 10,
    boundary = -100,closed = "right",
    aes(label = after_stat(count)),
    vjust = -0.5
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 10)) +
  labs(
    title = "EHI LQ histogram distribution",
    x = "Score",
    y = "Frecuencia"
  ) +
  theme_minimal()

#Separated histogram
ggplot(all_LQ, aes(x = Score10items, fill = BD_patient)) +
  geom_histogram(
    binwidth = 10,
    boundary = -100,
    position = "dodge",
    color = "black"
  ) +
  geom_text(
    stat = "bin",
    binwidth = 10,
    boundary = -100,
    position = position_dodge(width = 10),
    aes(label = after_stat(count), group = BD_patient),
    vjust = -0.5,
    size = 3
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 10)) +
  scale_fill_manual(values = c("NO" = "lightsalmon", "YES" = "slateblue2", "SCZ" = "cadetblue")) +
  labs(
    title = "LQ distribution by BD-nonBD",
    x = "Score",
    y = "N",
    fill = "BD_patient"
  ) +
  theme_minimal()

