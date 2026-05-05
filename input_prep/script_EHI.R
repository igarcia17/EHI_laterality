library(dplyr)
library(tidyr)
library(effsize)
library(psych)
library(psychTools)
library(ggplot2)
library(mvtnorm)
library(class)
library(VIM)
library(pheatmap)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

casesF <- "cases_EHI_300426_CHECKED_DIAGNOSTICS.xlsx"
controlsF <- "controls_EHI_300426_DIAGNOSTICSUNIFORM.xlsx"

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

cols <- c("ID", "Item1", "Item2", "Item3", "Item4", "Item5", "Item6", "Item7",
          "Item8", "Item9", "Item10", "Item11", "Item12")
colnames(ehi_cases) <- cols
colnames(ehi_controls) <- cols

all_ehi <- rbind(ehi_cases,ehi_controls)%>%
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

ehi_items_corr <- cor(all_ehi_imputed[,2:11], use = "complete.obs")
pheatmap(ehi_items_corr, 
         display_numbers = TRUE,       # Mostrar los valores de correlación
         number_format = "%.2f",       # Formato con 2 decimales
         number_color = "black",       # Color de los números
         fontsize_number = 10,         # Tamaño de los números
         color = colorRampPalette(c("white", "lightyellow", "lightsalmon"))(100), 
         main = "EHI items correlation",
         border_color = "white",       
         treeheight_row = 30,          
         treeheight_col = 30          
)

calculate_LQ_EHI <- function(df){
  df_temp <- df %>%
    mutate(
      Right10 = rowSums(
        across(Item1:Item10, ~ case_when(
          . == 1 ~ 2,
          . %in% c(2) ~ 1,
          . %in% c(3,4, 5) ~ 0,
          TRUE ~ NA_real_
        )),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      Left10 = rowSums(
        across(Item1:Item10, ~ case_when(
          . == 5 ~ 2,
          . %in% c(4) ~ 1,
          . %in% c(1,2,3) ~ 0,
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

all_samples <- rbind(cases, controls)%>%
  distinct(ID, .keep_all = TRUE)

all_LQ <- calculate_LQ_EHI(all_ehi_imputed)%>%
  left_join(all_samples %>% select(ID, BD_patient), by = "ID")


cases_LQ <- calculate_LQ_EHI(ehi_cases)
controls_LQ <- calculate_LQ_EHI(ehi_controls)

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
  scale_x_continuous(breaks = seq(-100, 100, 10)) +scale_fill_manual(values = c("NO" = "lightsalmon", "YES" = "slateblue2")) +
  labs(
    title = "LQ distribution by BD-nonBD",
    x = "Score",
    y = "N",
    fill = "BD_patient"
  ) +
  theme_minimal()




#Variables categoricas segun umbrales
data_cases <- data_cases %>%
  mutate(
    Score10_binary = cut(
      Score10items,
      breaks = c(-101, 0, 101),
      labels = c("Left", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_binary = cut(
      Score12items,
      breaks = c(-101, 0, 101),
      labels = c("Left", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score10_thr50 = cut(
      Score10items,
      breaks = c(-101, -50, 50, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_thr50 = cut(
      Score12items,
      breaks = c(-101, -50, 50, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score10_thr60 = cut(
      Score10items,
      breaks = c(-101, -60, 60, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_thr60 = cut(
      Score12items,
      breaks = c(-101, -60, 60, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score10_thr40 = cut(
      Score10items,
      breaks = c(-101, -40, 40, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_thr40 = cut(
      Score12items,
      breaks = c(-101, -40, 40, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

#Setting up controls
data_controls <- data_controls %>%
  mutate(
    Right10 = rowSums(
      across(Item1:Item10, ~ case_when(
        . == 1 ~ 2,
        . %in% c(2, 3) ~ 1,
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
  )%>%
  mutate(
    Right12 = rowSums(
      across(Item1:Item12, ~ case_when(
        . == 1 ~ 2,
        . %in% c(2, 3) ~ 1,
        . %in% c(4, 5) ~ 0,
        TRUE ~ NA_real_
      )),
      na.rm = TRUE
    )
  ) %>%
  mutate(
    Left12 = rowSums(
      across(Item1:Item12, ~ case_when(
        . == 5 ~ 2,
        . %in% c(4,3) ~ 1,
        . %in% c(1,2) ~ 0,
        TRUE ~ NA_real_
      )),
      na.rm = TRUE
    )
  ) %>%
  mutate(
    Score12items = round(((Right12 - Left12) / (Right12 + Left12) *100),0)
  )
#Density plot of scores
ggplot(data_controls) +
  geom_density(aes(x = Score10items), color = "lightblue", fill = "lightblue", alpha = 0.9) +
  geom_density(aes(x = Score12items), color = "coral", fill = "coral", alpha = 0.3) +
  labs(
    title = "EHI in 10 items or 12 items",
    x = "Score",
    y = "Density"
  ) +
  theme_minimal()

#Variables categoricas segun umbrales
data_controls <- data_controls %>%
  mutate(
    Score10_binary = cut(
      Score10items,
      breaks = c(-101, 0, 101),
      labels = c("Left", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_binary = cut(
      Score12items,
      breaks = c(-101, 0, 101),
      labels = c("Left", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score10_thr50 = cut(
      Score10items,
      breaks = c(-101, -50, 50, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_thr50 = cut(
      Score12items,
      breaks = c(-101, -50, 50, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score10_thr60 = cut(
      Score10items,
      breaks = c(-101, -60, 60, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_thr60 = cut(
      Score12items,
      breaks = c(-101, -60, 60, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score10_thr40 = cut(
      Score10items,
      breaks = c(-101, -40, 40, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    ),
    Score12_thr40 = cut(
      Score12items,
      breaks = c(-101, -40, 40, 101),
      labels = c("Left", "Ambidextrous", "Right"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

