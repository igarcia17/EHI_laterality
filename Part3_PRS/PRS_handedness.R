library(dplyr)
library(tidyverse)
library(tableone)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(forcats)
library(pscl)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())
########################Parameters
cutoffs <- c("0", "60", "90")
oldies <- T
write_count_comparisons <- FALSE
variables <- c("PRS_handedness","PRS_amb","PRS_BD","PRS_SCZ")
##########Functions

create_comparison_table <- function(df, Strata_var){
  table1 <- CreateTableOne(vars = c("SEX",variables), strata = Strata_var, 
    data = df, factorVars = c("SEX"),
    addOverall = TRUE)
  
  t1 <- print(table1, showAllLevels = TRUE, nonnormal = c("Age"), 
              formatOptions = list(big.mark = ","), pDigits=16)
  t1 <- as.data.frame(t1)
  return(t1)
}
create_combined_comparison_table <- function(df, var1, var2){
  # 1. Creamos una nueva columna que combine ambas variables
  # El uso de interaction() crea grupos como "Yes.Right", "No.Right", etc.
  df_temp <- df %>%
    mutate(Combined_Strata = interaction(!!sym(var1), !!sym(var2), sep = " - "))
  table1 <- CreateTableOne(
    vars = c("SEX", variables), strata = "Combined_Strata", 
    data = df_temp, factorVars = c("SEX"),addOverall = TRUE)
  
  t1 <- print(table1, 
              showAllLevels = TRUE, nonnormal = c("Age"),
              formatOptions = list(big.mark = ","), pDigits = 16,
              printToggle = FALSE) 
    t1 <- as.data.frame(t1)
  return(t1)
}

########## Input

input <- as.data.frame(readxl::read_xlsx("All_data_available_PRS.xlsx", sheet = 2, 
                                         na="#N/A"))%>%
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
                PRS_SCZ = as.numeric(PRS_SCZ)
                )

sapply(input, class)

if(oldies){
  input <- input %>%
    filter(Oldies_dataset == "YES") %>% droplevels()
}

#Set up variables of interest
df <- input%>%
  drop_na(FID)%>%
  mutate(Score10items = as.numeric(Score10items))%>%
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
# Normalize PRS values
#before
colMeans(df[variables], na.rm = TRUE)
apply(df[variables], 2, sd, na.rm = TRUE)

df[variables] <- scale(df[variables])

#after
colMeans(df[variables], na.rm = TRUE)
apply(df[variables], 2, sd, na.rm = TRUE)

#Visualize
df_long <- df %>%
  pivot_longer(cols = starts_with("PRS"), names_to = "PRS_Type", values_to = "Value") 

# Plot of the 4 PRS in BD patients and controls
ggplot(df_long, aes(x = Value, fill = NRH_0)) +
  geom_density(alpha = 0.5) + 
  facet_wrap(~ PRS_Type, scales = "free_y", ncol = 2) + 
  #scale_fill_manual(values = c("RH" = "#F8766D","NRH" = "#00BFC4" )) +
  theme_minimal() +
  labs(
    #title = "PRS distribution between left and right handers",
    x = "Scaled",
    y = "Densidad"
    #, fill = "Left (NRH) or right (RH)"
  )
dev.off()

df_only_BD <- df%>%filter((STATUS=="YES"))

#Run logistic regression
m1<- glm(NRH_0 ~ SEX + Age + PC1_oldies+ PC2_oldies + PRS_handedness, data=df,
         family=binomial
         )
m_b<- glm(NRH_0 ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df, family=binomial)
summary(m1)
r2_full <-pR2(m1)["r2CU"] ###Full r² (%) Nagelkarke
r2_base <- pR2(m_b)["r2CU"]
prs_r2 <- r2_full - r2_base ###PRS r² (%) Nagelkarke
prs_r2
r2_full

m1<- lm(Score10items ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_handedness, data=df)
m_b<- lm(Score10items ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

m1<- glm(STATUS ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_handedness, data=df, family=binomial)
m_b<- glm(STATUS ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df, family=binomial)
summary(m1)
r2_full <-pR2(m1)["r2CU"] ###Full r² (%) Nagelkarke
r2_base <- pR2(m_b)["r2CU"]
prs_r2 <- r2_full - r2_base ###PRS r² (%) Nagelkarke
prs_r2
r2_full

##Severity
m1<- lm(WHODAS ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_handedness, data=df_only_BD)
m_b<- lm(WHODAS ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared
summary(m1)$r.squared - summary(m_b)$r.squared

m1<- lm(CGI ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_handedness, data=df_only_BD)
m_b<- lm(CGI ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

m1<- lm(INES ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_handedness, data=df_only_BD)
m_b<- lm(INES ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

m1<- lm(GAF ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_handedness, data=df_only_BD)
m_b<- lm(GAF ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

###PRS of BD
m1<- glm(STATUS ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df, family=binomial)
m_b<- glm(STATUS ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df, family=binomial)
summary(m1)
r2_full <-pR2(m1)["r2CU"] ###Full r² (%) Nagelkarke
r2_base <- pR2(m_b)["r2CU"]
prs_r2 <- r2_full - r2_base ###PRS r² (%) Nagelkarke
prs_r2
r2_full

##Severity using PRS of BD
m1<- lm(WHODAS ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df_only_BD)
m_b<- lm(WHODAS ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared
summary(m1)$r.squared - summary(m_b)$r.squared

m1<- lm(CGI ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df_only_BD)
m_b<- lm(CGI ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

m1<- lm(INES ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df_only_BD)
m_b<- lm(INES ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

m1<- lm(GAF ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df_only_BD)
m_b<- lm(GAF ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

#Left vs right predicted by PRS of BD in all samples and only in cases
m1<- glm(NRH_0 ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df, family=binomial)
m_b<- glm(NRH_0 ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df, family=binomial)
summary(m1)
r2_full <-pR2(m1)["r2CU"] ###Full r² (%) Nagelkarke
r2_base <- pR2(m_b)["r2CU"]
prs_r2 <- r2_full - r2_base ###PRS r² (%) Nagelkarke
prs_r2
r2_full

m1<- glm(NRH_0 ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df_only_BD, family=binomial)
m_b<- glm(NRH_0 ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD, family=binomial)
summary(m1)
r2_full <-pR2(m1)["r2CU"] ###Full r² (%) Nagelkarke
r2_base <- pR2(m_b)["r2CU"]
prs_r2 <- r2_full - r2_base ###PRS r² (%) Nagelkarke
prs_r2
r2_full

m1<- lm(Score10items ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df)
m_b<- lm(Score10items ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke

m1<- lm(Score10items ~ SEX + Age+ PC1_oldies+ PC2_oldies + PRS_BD, data=df_only_BD)
m_b<- lm(Score10items ~ SEX + Age + PC1_oldies+ PC2_oldies, data=df_only_BD)
summary(m1)
summary(m1)$r.squared###Full r² (%) Nagelkarke
summary(m1)$r.squared - summary(m_b)$r.squared###PRS r² (%) Nagelkarke









#
#
#

###Getting into the comparisons
if(write_count_comparisons){
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
}
#
#
#