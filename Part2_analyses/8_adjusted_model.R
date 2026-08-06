# Calcular diferencias ajustadas por covariables en uso de zurdera
library(dplyr)
library(ggplot2)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

thresholds <- c(0, 40, 60, 80, 90)
input <- as.data.frame(readxl::read_xlsx("../Part1_input_prep/all_samples_all_data_LQ_PLUS_items.xlsx", sheet = 1, 
                                         col_types="text", na="#N/A"))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>%
  mutate(across(
    -c(ID, Diagnostic),                    
    ~ if (!is.factor(.)) as.numeric(.) else .
  )) %>% mutate(Diagnostic = factor(Diagnostic))

input_clean <- subset(input, is.finite(Age))
set.seed(13)
match_obj <- MatchIt::matchit(BD_patient ~ Age, 
                              data = input_clean, 
                              method = "nearest", 
                              ratio = 1)
balanced_age <- MatchIt::match.data(match_obj)

df <- balanced_age %>%mutate(
  
  # --- NRH columns ---
  NRH_0  = factor(if_else(Score10items <= 0,  "NRH", "RH"),
                  levels = c("RH", "NRH")),
  NRH_40 = factor(if_else(Score10items <= 40, "NRH", "RH"),
                  levels = c("RH", "NRH")),
  NRH_60 = factor(if_else(Score10items <= 60, "NRH", "RH"),
                  levels = c("RH", "NRH")),
  NRH_80 = factor(if_else(Score10items <= 80, "NRH", "RH"),
                  levels = c("RH", "NRH")),
  NRH_90 = factor(if_else(Score10items <= 90, "NRH", "RH"),
                  levels = c("RH", "NRH")),
  
  # --- Lateralization columns ---
  LAT_0  = factor(if_else(
    Score10items >= 0  & Score10items <= 0,"Non-lateral","Lateral")),
  
  LAT_40 = factor(if_else(
    Score10items >= -40 & Score10items <= 40,"Non-lateral","Lateral")),
  
  LAT_60 = factor(if_else(
    Score10items >= -60 & Score10items <= 60,"Non-lateral","Lateral")),
  
  LAT_80 = factor(if_else(
    Score10items >= -80 & Score10items <= 80,"Non-lateral", "Lateral" )),
  
  LAT_90 = factor(if_else(
    Score10items >= -90 & Score10items <= 90,"Non-lateral","Lateral")))

controls <- df %>% filter(BD_patient =="NO") %>% 
  mutate(group_Age = if_else(Age >= 40, "OLD", "YOUNG"),group_Age = factor(group_Age, levels = c("YOUNG", "OLD")))%>%
  mutate(NRH_0  = factor(if_else(Score10items <= 0,  "NRH", "RH"), levels = c("RH", "NRH")),
         NRH_40 = factor(if_else(Score10items <= 40, "NRH", "RH"),levels = c("RH", "NRH")),
         NRH_60 = factor(if_else(Score10items <= 60, "NRH", "RH"),levels = c("RH", "NRH")),
         NRH_80 = factor(if_else(Score10items <= 80, "NRH", "RH"),levels = c("RH", "NRH")),
         NRH_90 = factor(if_else(Score10items <= 90, "NRH", "RH"),levels = c("RH", "NRH")),
         LAT_0  = factor(if_else(Score10items >= 0  & Score10items <= 0,"Non-lateral","Lateral")),
         LAT_40 = factor(if_else(Score10items >= -40 & Score10items <= 40,"Non-lateral","Lateral")),
         LAT_60 = factor(if_else(Score10items >= -60 & Score10items <= 60,"Non-lateral","Lateral")),
         LAT_80 = factor(if_else(Score10items >= -80 & Score10items <= 80,"Non-lateral", "Lateral" )),
         LAT_90 = factor(if_else(Score10items >= -90 & Score10items <= 90,"Non-lateral","Lateral")))

scales <- as.data.frame(readxl::read_xlsx("Escalas_clinicas_MadManic.xlsx", sheet = 1, 
                                          na=c("#N/A", "#N/D")))
scales[[1]] <- gsub("IE", "IEV", scales[[1]])
colnames(scales) <- c("ID", "INES", "GAF", "CGI", "WHODAS")

cases <- df %>% filter(BD_patient =="YES") %>% 
  mutate(group_Age = if_else(Age >= 40, "OLD", "YOUNG"),group_Age = factor(group_Age, levels = c("YOUNG", "OLD")))%>%
  mutate(NRH_0  = factor(if_else(Score10items <= 0,  "NRH", "RH"), levels = c("RH", "NRH")),
         NRH_40 = factor(if_else(Score10items <= 40, "NRH", "RH"),levels = c("RH", "NRH")),
         NRH_60 = factor(if_else(Score10items <= 60, "NRH", "RH"),levels = c("RH", "NRH")),
         NRH_80 = factor(if_else(Score10items <= 80, "NRH", "RH"),levels = c("RH", "NRH")),
         NRH_90 = factor(if_else(Score10items <= 90, "NRH", "RH"),levels = c("RH", "NRH")),
         LAT_0  = factor(if_else(Score10items >= 0  & Score10items <= 0,"Non-lateral","Lateral")),
         LAT_40 = factor(if_else(Score10items >= -40 & Score10items <= 40,"Non-lateral","Lateral")),
         LAT_60 = factor(if_else(Score10items >= -60 & Score10items <= 60,"Non-lateral","Lateral")),
         LAT_80 = factor(if_else(Score10items >= -80 & Score10items <= 80,"Non-lateral", "Lateral" )),
         LAT_90 = factor(if_else(Score10items >= -90 & Score10items <= 90,"Non-lateral","Lateral"))) %>%
  left_join(scales, by = "ID")


###Models

#Differences in control group
summary(glm(NRH_0 ~  Age + Sex, data=controls, family = "binomial"))
summary(glm(NRH_40 ~  Age + Sex, data=controls, family = "binomial"))
summary(glm(NRH_60 ~  Age + Sex, data=controls, family = "binomial")) ##Sex
summary(glm(NRH_80 ~  Age + Sex, data=controls, family = "binomial"))##Sex
summary(glm(NRH_90 ~  Age + Sex, data=controls, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~  Age + Sex, data=controls, family = "binomial"))
summary(glm(LAT_40 ~  Age + Sex, data=controls, family = "binomial"))
summary(glm(LAT_60 ~  Age + Sex, data=controls, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ Age + Sex, data=controls, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ Age + Sex, data=controls, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#regarding NRH or lat
summary(glm(NRH_0 ~ BD_patient + Age + Sex, data=df, family = "binomial"))
summary(glm(NRH_40 ~ BD_patient + Age + Sex, data=df, family = "binomial"))
summary(glm(NRH_60 ~ BD_patient + Age + Sex, data=df, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ BD_patient + Age + Sex, data=df, family = "binomial"))##Sex
summary(glm(NRH_90 ~ BD_patient + Age + Sex, data=df, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ BD_patient + Age + Sex, data=df, family = "binomial"))
summary(glm(LAT_40 ~ BD_patient + Age + Sex, data=df, family = "binomial"))
summary(glm(LAT_60 ~ BD_patient + Age + Sex, data=df, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ BD_patient + Age + Sex, data=df, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ BD_patient + Age + Sex, data=df, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#Within BD
summary(glm(NRH_0 ~  Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_40 ~  Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_60 ~  Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(NRH_80 ~  Age + Sex, data=cases, family = "binomial"))##Sex
summary(glm(NRH_90 ~  Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~  Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_40 ~  Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_60 ~  Age + Sex, data=cases, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#Only females
fems <- df %>% filter(Sex =="MUJER") %>% droplevels()
summary(glm(NRH_0 ~ BD_patient + Age, data=fems, family = "binomial"))
summary(glm(NRH_40 ~ BD_patient + Age, data=fems, family = "binomial"))
summary(glm(NRH_60 ~ BD_patient + Age, data=fems, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ BD_patient + Age, data=fems, family = "binomial"))##Sex
summary(glm(NRH_90 ~ BD_patient + Age, data=fems, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ BD_patient + Age, data=fems, family = "binomial"))
summary(glm(LAT_40 ~ BD_patient + Age, data=fems, family = "binomial"))
summary(glm(LAT_60 ~ BD_patient + Age, data=fems, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ BD_patient + Age, data=fems, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ BD_patient + Age, data=fems, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#Only males
males <- df %>% filter(Sex=="HOMBRE") %>% droplevels()
summary(glm(NRH_0 ~ BD_patient + Age, data=males, family = "binomial"))
summary(glm(NRH_40 ~ BD_patient + Age, data=males, family = "binomial"))
summary(glm(NRH_60 ~ BD_patient + Age, data=males, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ BD_patient + Age, data=males, family = "binomial"))##Sex
summary(glm(NRH_90 ~ BD_patient + Age, data=males, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ BD_patient + Age, data=males, family = "binomial"))
summary(glm(LAT_40 ~ BD_patient + Age, data=males, family = "binomial"))
summary(glm(LAT_60 ~ BD_patient + Age, data=males, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ BD_patient + Age, data=males, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ BD_patient + Age, data=males, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#Severity scales

#WHODAS
summary(glm(NRH_0 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_40 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_60 ~ WHODAS + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))##Sex
summary(glm(NRH_90 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_40 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_60 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ WHODAS + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ WHODAS + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#CGI
summary(glm(NRH_0 ~ CGI + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_40 ~ CGI + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_60 ~ CGI + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ CGI + Age + Sex, data=cases, family = "binomial"))##Sex
summary(glm(NRH_90 ~ CGI + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ CGI + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_40 ~ CGI + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_60 ~ CGI + Age + Sex, data=cases, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ CGI + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ CGI + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#INES
summary(glm(NRH_0 ~ INES + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_40 ~ INES + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_60 ~ INES + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ INES + Age + Sex, data=cases, family = "binomial"))##Sex
summary(glm(NRH_90 ~ INES + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ INES + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_40 ~ INES + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_60 ~ INES + Age + Sex, data=cases, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ INES + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ INES + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

#GAF
summary(glm(NRH_0 ~ GAF + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_40 ~ GAF + Age + Sex, data=cases, family = "binomial"))
summary(glm(NRH_60 ~ GAF + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ GAF + Age + Sex, data=cases, family = "binomial"))##Sex
summary(glm(NRH_90 ~ GAF + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ GAF + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_40 ~ GAF + Age + Sex, data=cases, family = "binomial"))
summary(glm(LAT_60 ~ GAF + Age + Sex, data=cases, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ GAF + Age + Sex, data=cases, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ GAF + Age + Sex, data=cases, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

# heatmap correlatoins between scales and sex or age
library(dplyr)
library(pheatmap)

var <- c("Sex", "Age")
sca <- c("WHODAS", "INES", "CGI", "GAF")

# Preparar las variables
df_cor <- cases %>%
  select(all_of(c(var, sca))) %>%
  mutate(
    Sex = case_when(
      Sex == "HOMBRE"   ~ 0,
      Sex == "MUJER" ~ 1,
      TRUE ~ NA_real_
    )
  )
# Matriz de correlaciones
cor_matrix <- cor(
  df_cor,
  use = "pairwise.complete.obs",
  method = "spearman"
)
cor_heatmap <- cor_matrix[var, sca, drop = FALSE]

# Heatmap
pheatmap(
  cor_heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Correlations with clinical scales",
  border_color = NA
)




summary(glm(BD_patient ~ Age + Sex + NRH_0, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + NRH_40, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + NRH_60, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + NRH_80, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + NRH_90, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + LAT_0, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + LAT_40, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + LAT_60, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + LAT_80, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex + LAT_90, data=df, family = "binomial"))

###iNTERACTION

summary(glm(BD_patient ~ Age + Sex*NRH_0, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*NRH_40, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*NRH_60, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*NRH_80, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*NRH_90, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*LAT_0, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*LAT_40, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*LAT_60, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*LAT_80, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex*LAT_90, data=df, family = "binomial"))

##wITH LQ
summary(glm(BD_patient ~ Age + Sex + Score10items, data=df, family = "binomial"))
summary(glm(BD_patient ~ Age + Sex * Score10items, data=df, family = "binomial"))




females_clean <- df%>%
  filter(Sex=="MUJER")%>%droplevels()
males_clean <- df%>%filter(Sex=="HOMBRE")%>%droplevels()

#On females
###Models

summary(glm(BD_patient ~ Age  + NRH_0, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_40, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_60, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_80, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_90, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_0, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_40, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_60, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_80, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_90, data=females_clean, family = "binomial"))

###iNTERACTION

summary(glm(BD_patient ~ Age *NRH_0, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_40, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_60, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_80, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_90, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_0, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_40, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_60, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_80, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_90, data=females_clean, family = "binomial"))

##wITH LQ
summary(glm(BD_patient ~ Age  + Score10items, data=females_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  * Score10items, data=females_clean, family = "binomial"))


#regarding NRH or lat
summary(glm(NRH_0 ~ BD_patient + Age , data=females_clean, family = "binomial"))
summary(glm(NRH_40 ~ BD_patient + Age , data=females_clean, family = "binomial"))
summary(glm(NRH_60 ~ BD_patient + Age , data=females_clean, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ BD_patient + Age , data=females_clean, family = "binomial"))##Sex
summary(glm(NRH_90 ~ BD_patient + Age , data=females_clean, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ BD_patient + Age , data=females_clean, family = "binomial"))
summary(glm(LAT_40 ~ BD_patient + Age , data=females_clean, family = "binomial"))
summary(glm(LAT_60 ~ BD_patient + Age , data=females_clean, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~ BD_patient + Age , data=females_clean, family = "binomial")) ##Sex
summary(glm(LAT_90 ~ BD_patient + Age , data=females_clean, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status



#On males
###Models

summary(glm(BD_patient ~ Age  + NRH_0, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_40, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_60, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_80, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + NRH_90, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_0, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_40, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_60, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_80, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  + LAT_90, data=males_clean, family = "binomial"))

###iNTERACTION

summary(glm(BD_patient ~ Age *NRH_0, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_40, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_60, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_80, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *NRH_90, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_0, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_40, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_60, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_80, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age *LAT_90, data=males_clean, family = "binomial"))

##wITH LQ
summary(glm(BD_patient ~ Age  + Score10items, data=males_clean, family = "binomial"))
summary(glm(BD_patient ~ Age  * Score10items, data=males_clean, family = "binomial"))


#regarding NRH or lat
summary(glm(NRH_0 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(NRH_40 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(NRH_60 ~ BD_patient + Age , data=males_clean, family = "binomial")) 
summary(glm(NRH_80 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(NRH_90 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(LAT_0 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(LAT_40 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(LAT_60 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(LAT_80 ~ BD_patient + Age , data=males_clean, family = "binomial"))
summary(glm(LAT_90 ~ BD_patient + Age , data=males_clean, family = "binomial"))


BD_clean <- df%>%
  filter(BD_patient=="YES")%>%droplevels()
controls_clean <- df%>%filter(BD_patient=="NO")%>%droplevels()

summary(glm(NRH_0 ~ Age + Sex, data=BD_clean, family = "binomial"))
summary(glm(NRH_40 ~  Age + Sex, data=BD_clean, family = "binomial"))
summary(glm(NRH_60 ~  Age + Sex, data=BD_clean, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ Age + Sex, data=BD_clean, family = "binomial"))##Sex
summary(glm(NRH_90 ~  Age + Sex, data=BD_clean, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ Age + Sex, data=BD_clean, family = "binomial"))
summary(glm(LAT_40 ~ Age + Sex, data=BD_clean, family = "binomial"))
summary(glm(LAT_60 ~ Age + Sex, data=BD_clean, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~  Age + Sex, data=BD_clean, family = "binomial")) ##Sex
summary(glm(LAT_90 ~  Age + Sex, data=BD_clean, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status

summary(glm(NRH_0 ~ Age + Sex, data=controls_clean, family = "binomial"))
summary(glm(NRH_40 ~  Age + Sex, data=controls_clean, family = "binomial"))
summary(glm(NRH_60 ~  Age + Sex, data=controls_clean, family = "binomial")) ##Sex
summary(glm(NRH_80 ~ Age + Sex, data=controls_clean, family = "binomial"))##Sex
summary(glm(NRH_90 ~  Age + Sex, data=controls_clean, family = "binomial"))##Sex + BD + Age -> there are strong changes in NRH with age sex and BD status
summary(glm(LAT_0 ~ Age + Sex, data=controls_clean, family = "binomial"))
summary(glm(LAT_40 ~ Age + Sex, data=controls_clean, family = "binomial"))
summary(glm(LAT_60 ~ Age + Sex, data=controls_clean, family = "binomial"))#Some in BD, some in Sex
summary(glm(LAT_80 ~  Age + Sex, data=controls_clean, family = "binomial")) ##Sex
summary(glm(LAT_90 ~  Age + Sex, data=controls_clean, family = "binomial"))##Sex + BD + Age -> there are strong changes in lateralization with age sex and BD status










