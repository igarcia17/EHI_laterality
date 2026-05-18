library(dplyr)
library(tidyverse)
library(tableone)
library(ggplot2)
library(tidyr)
library(openxlsx)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())
########################Parameters
oldies <- TRUE
variables <- c("Age","PRS_handedness","PRS_amb","PRS_BD","PRS_SCZ")
##########Functions

create_comparison_table <- function(df, Strata_var){
  # Crear la tabla base
  table1 <- CreateTableOne(
    vars = c("SEX",variables), 
    strata = Strata_var, 
    data = df, 
    factorVars = c("SEX"),
    addOverall = TRUE
  )
  
  t1 <- print(table1, 
              showAllLevels = TRUE, 
              nonnormal = variables,
              formatOptions = list(big.mark = ","), pDigits=16)
  t1 <- as.data.frame(t1)
  return(t1)
}
create_combined_comparison_table <- function(df, var1, var2){
  # 1. Creamos una nueva columna que combine ambas variables
  # El uso de interaction() crea grupos como "Yes.Right", "No.Right", etc.
  df_temp <- df %>%
    mutate(Combined_Strata = interaction(!!sym(var1), !!sym(var2), sep = " - "))
  
  # 2. Crear la tabla usando la nueva variable combinada
  table1 <- CreateTableOne(
    vars = c("SEX", variables), 
    strata = "Combined_Strata", 
    data = df_temp, 
    factorVars = c("SEX"),
    addOverall = TRUE
  )
  
  # 3. Formatear y exportar
  t1 <- print(table1, 
              showAllLevels = TRUE, 
              nonnormal = variables,
              formatOptions = list(big.mark = ","), 
              pDigits = 16,
              printToggle = FALSE) # Evita que se imprima doble en la consola
  
  t1 <- as.data.frame(t1)
  return(t1)
}

########## Input

input <- as.data.frame(readxl::read_xlsx("All_data_available_PRS.xlsx", sheet = 2, 
                                         col_types="text", na="#N/A"))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>% mutate(Diagnostic = factor(Diagnostic),
                Cutoff0 = factor(Cutoff0),
                Cutoff40= factor(Cutoff40),
                Cutoff60 = factor(Cutoff60),
                Cutoff80 = factor(Cutoff80),Cutoff90 = factor(Cutoff90),
                Age = as.numeric(Age),
                PRS_handedness= as.numeric(PRS_handedness),
                PRS_amb = as.numeric(PRS_amb),
                PRS_BD = as.numeric(PRS_BD),
                PRS_SCZ = as.numeric(PRS_SCZ)
                )

sapply(input, class)

if(oldies){
  input <- input %>%
    filter(Oldies_dataset == "YES") %>% droplevels()
}

df <- input%>%
  drop_na(FID)

df_long <- df %>%
  pivot_longer(cols = starts_with("PRS"), names_to = "PRS_Type", values_to = "Value") %>%
  group_by(PRS_Type) %>%
  mutate(
    # Centrar en 0 (Restar la media)
    Value_Centered = Value - mean(Value, na.rm = TRUE),
    # Escalar de 0 a 1: (x - min) / (max - min)
    Value_Scaled = (Value_Centered - min(Value_Centered)) / (max(Value_Centered) - min(Value_Centered))
  ) %>%
  ungroup()

# 2. Visualización con ggplot2
ggplot(df_long, aes(x = Value_Scaled, fill = STATUS)) +
  geom_density(alpha = 0.5) + 
  facet_wrap(~ PRS_Type, scales = "free_y", ncol = 2) + 
  scale_fill_manual(values = c("YES" = "#00BFC4", "NO" = "#F8766D")) +
  theme_minimal() +
  labs(
    title = "PRS distribution",
    x = "Scaled",
    y = "Densidad",
    fill = "BD patient"
  )
###Getting into the comparisons


pgs_BD <- create_comparison_table(df, "STATUS")
pgs_hand_cut0 <- create_comparison_table(df, "Cutoff0")
pgs_hand_cut40<- create_comparison_table(df, "Cutoff40")
pgs_hand_cut60<- create_comparison_table(df, "Cutoff60")
pgs_hand_cut80<- create_comparison_table(df, "Cutoff80")
pgs_hand_cut90<- create_comparison_table(df, "Cutoff90")

wb <- createWorkbook()
addWorksheet(wb, "all_PGS_all_indiv_byBD")
addWorksheet(wb, "all_indiv_by_cutoff0")
addWorksheet(wb, "all_indiv_by_cutoff40")
addWorksheet(wb, "all_indiv_by_cutoff60")
addWorksheet(wb, "all_indiv_by_cutoff80")
addWorksheet(wb, "all_indiv_by_cutoff90")

writeData(wb, "all_PGS_all_indiv_byBD", pgs_BD, rowNames = TRUE)
writeData(wb, "all_indiv_by_cutoff0", pgs_hand_cut0, rowNames = TRUE)
writeData(wb, "all_indiv_by_cutoff40", pgs_hand_cut40, rowNames = TRUE)
writeData(wb, "all_indiv_by_cutoff60", pgs_hand_cut60, rowNames = TRUE)
writeData(wb, "all_indiv_by_cutoff80", pgs_hand_cut80, rowNames = TRUE)
writeData(wb, "all_indiv_by_cutoff90", pgs_hand_cut90, rowNames = TRUE)

#PGS only in BD, strata by cutoffs
df_BD <- df%>% filter(STATUS=="YES") %>%droplevels()

pgs_BD_hand_cut0 <- create_comparison_table(df_BD, "Cutoff0")
pgs_BD_hand_cut40<- create_comparison_table(df_BD, "Cutoff40")
pgs_BD_hand_cut60<- create_comparison_table(df_BD, "Cutoff60")
pgs_BD_hand_cut80<- create_comparison_table(df_BD, "Cutoff80")
pgs_BD_hand_cut90<- create_comparison_table(df_BD, "Cutoff90")

addWorksheet(wb, "BD_indiv_by_cutoff0")
addWorksheet(wb, "BD_indiv_by_cutoff40")
addWorksheet(wb, "BD_indiv_by_cutoff60")
addWorksheet(wb, "BD_indiv_by_cutoff80")
addWorksheet(wb, "BD_indiv_by_cutoff90")

writeData(wb, "BD_indiv_by_cutoff0", pgs_BD_hand_cut0, rowNames = TRUE)
writeData(wb, "BD_indiv_by_cutoff40", pgs_BD_hand_cut40, rowNames = TRUE)
writeData(wb, "BD_indiv_by_cutoff60", pgs_BD_hand_cut60, rowNames = TRUE)
writeData(wb, "BD_indiv_by_cutoff80", pgs_BD_hand_cut80, rowNames = TRUE)
writeData(wb, "BD_indiv_by_cutoff90", pgs_BD_hand_cut90, rowNames = TRUE)

##PGS by seggregating BD and hand type, by cutoff
comb_0 <- create_combined_comparison_table(df, "STATUS", "Cutoff0")
comb_40 <- create_combined_comparison_table(df, "STATUS", "Cutoff40")
comb_60 <- create_combined_comparison_table(df, "STATUS", "Cutoff60")
comb_80 <- create_combined_comparison_table(df, "STATUS", "Cutoff80")
comb_90 <- create_combined_comparison_table(df, "STATUS", "Cutoff90")

addWorksheet(wb, "Comb_cutoff0")
addWorksheet(wb, "Comb_cutoff40")
addWorksheet(wb, "Comb_cutoff60")
addWorksheet(wb, "Comb_cutoff80")
addWorksheet(wb, "Comb_cutoff90")

writeData(wb, "Comb_cutoff0", comb_0, rowNames = TRUE)
writeData(wb, "Comb_cutoff40", comb_40, rowNames = TRUE)
writeData(wb, "Comb_cutoff60", comb_60, rowNames = TRUE)
writeData(wb, "Comb_cutoff80", comb_80, rowNames = TRUE)
writeData(wb, "Comb_cutoff90", comb_90, rowNames = TRUE)



saveWorkbook(wb, "PGS_in_hand_groups_and_BD.xlsx", overwrite = TRUE)
#
#
#
#

