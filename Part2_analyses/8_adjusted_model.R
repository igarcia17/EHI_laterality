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

###Models

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










