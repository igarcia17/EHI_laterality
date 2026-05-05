library(dplyr)
library(tidyr)
library(effsize)
library(psych)
library(psychTools)
library(ggplot2)
library(stringr)
library(openxlsx)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

controlsF <- "datos_cuestionarios_UAM-CSIC_16_10_2025.xlsx"

c_1 <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 1, 
                                         col_types="text", na="#N/A")) %>%
  select(ID, SEXO, Edad, `Diagnóstico psiquiátrico propio`, `Detalles diagnóstico`)
c_2  <- as.data.frame(readxl::read_xlsx(controlsF, sheet =2, 
                                        col_types="text", na="#N/A")) %>%
  select(ID, SEXO, Edad, `Diagnóstico psiquiátrico propio`, `Detalles diagnóstico`)
c_3  <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 3, 
                                        col_types="text", na="#N/A")) %>%
  select(ID, SEXO, Edad, `Diagnóstico psiquiátrico propio`, `Detalles diagnóstico`)
c_4 <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 4, 
                                   col_types="text", na="#N/A")) %>%
  select(ID, SEXO, Edad, `Diagnóstico psiquiátrico propio`, `Detalles diagnóstico`)
c_5 <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 5, 
                                       col_types="text", na="#N/A")) %>%
  select(ID, SEXO, Edad, `Diagnóstico psiquiátrico propio`, `Detalles diagnóstico`)
c_6 <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 6, 
                                       col_types="text", na="#N/A")) %>%
  select(ID, SEXO, Edad, `Diagnóstico psiquiátrico propio`, `Detalles diagnóstico`)

ControlsCBMSO <- rbind(c_1, c_2, c_3, c_3, c_4, c_5, c_6)%>%
  mutate(across(everything(), ~ ifelse(is.na(.), "NO", .)))%>%
  mutate(across(where(is.character),
                ~ toupper(iconv(., from = "UTF-8", to = "ASCII//TRANSLIT"))))%>%
  mutate(ID = paste0(
    str_sub(ID, 1, 3), "_",
    str_pad(str_extract(ID, "\\d+"), width = 4, pad = "0")
  ))

cases <- ControlsCBMSO%>%
  filter(str_detect(`Detalles diagnóstico`, "BIPOLAR|TB|TAB|BD|BIPOILAR|ESQUIZOAFECTIVO"))
controls <- ControlsCBMSO%>%
  filter(!str_detect(`Detalles diagnóstico`, "BIPOLAR|TB|TAB|BD|BIPOILAR|ESQUIZOAFECTIVO"))%>%
  filter(!str_detect(`Detalles diagnóstico`, "AUTISMO|OBSESIVO|AUTISTA|ESQUIZOF|PERSONALIDAD|TEA|TLP|GRAVE|PSICOT|MAYOR|BOURCAUT"))%>%
  distinct(ID, .keep_all = TRUE)

edi_sheet <- as.data.frame(readxl::read_xlsx(controlsF, sheet = 7, 
                                             col_types="text", na="#N/A"))%>%
  mutate(ID = paste0(
    str_sub(ID, 1, 3), "_",
    str_pad(str_extract(ID, "\\d+"), width = 4, pad = "0")
  ))

ehi_cases <- edi_sheet %>%
  filter(ID %in% cases$ID)
ehi_controls <- edi_sheet %>%
  filter(ID %in% controls$ID)%>%
  distinct(ID, .keep_all = TRUE)

controls <- controls%>%
  filter(ID %in% ehi_controls$ID)

###Load MadManic data

casesF <- "Raw_MadManic_InesG.xlsx"
cases_madmanic <- as.data.frame(readxl::read_xlsx(casesF, sheet = 1, 
                                         col_types="text", na="#N/A"))%>%
  mutate(across(where(is.character),
                ~ toupper(iconv(., from = "UTF-8", to = "ASCII//TRANSLIT"))))%>%
  mutate(Genetic_code = paste0(
    str_sub(Genetic_code, 1, 3), "_",
    str_pad(str_extract(Genetic_code, "\\d+"), width = 4, pad = "0")
  ))%>%
  mutate(Year = as.numeric(Year),
         Edad = 2024 - Year)

meta_cases <- cases_madmanic %>%
  mutate(Cases = "YES") %>%
  select(Genetic_code, Sex, Edad, Cases, Diagnostic)%>%
  distinct(Genetic_code, .keep_all = TRUE)
  
data_cases <- cases_madmanic %>%
  select(-user_id, -Sex, -Edad,-Year, -Diagnostic)%>%
  distinct(Genetic_code, .keep_all = TRUE)

meta_cols <- c("ID", "Sex", "Age", "BD_patient", "Diagnostic")
colnames(meta_cases) <- meta_cols
colnames(cases) <- meta_cols
colnames(controls) <- meta_cols
controls$BD_patient <- "NO"

Cases_all <- rbind(meta_cases, cases)
colnames(ehi_cases)<- c(colnames(data_cases), "Total")
Cases_all_ehi<-rbind(ehi_cases[,1:13], data_cases)


###SAVE an EXCEL file for controls
wb <- createWorkbook()
addWorksheet(wb, "meta_controls")
writeData(wb, "meta_controls", controls)
addWorksheet(wb, "EHI_controls")
writeData(wb, "EHI_controls", ehi_controls)
saveWorkbook(wb, "controls_EHI_300426.xlsx", overwrite = TRUE)

#SAVE an EXCEL file for cases
wb <- createWorkbook()
addWorksheet(wb, "meta_cases")
writeData(wb, "meta_cases", Cases_all)
addWorksheet(wb, "EHI_cases")
writeData(wb, "EHI_cases", Cases_all_ehi)
saveWorkbook(wb, "cases_EHI_300426.xlsx", overwrite = TRUE)




