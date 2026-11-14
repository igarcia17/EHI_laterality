#
#Severity index phormula: (rapid_cycle + n_hosp + n_sui + n_manic + n_depressive + n_mixed + psychosis + n_hipomaniacos)/(Age - onset)
#

library(dplyr)
library(stringr)
library(tidyr)
library(effsize)
library(psych)
library(psychTools)
library(ggplot2)
library(openxlsx)
library(mvtnorm)
library(class)
library(VIM)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))


inputF <- "Enero_2026/madmanic-1_20251219084233234383.xlsx" #de este excel sacamos todos los demás datos

df0 <- as.data.frame(readxl::read_xlsx(inputF, sheet = 1)) %>%
  select(LitioSocio94, BirthDay, LitioSocio01, LUNDBECK12, LitioSocio21, LitioSocio41,
         LUNDBECK17, LUNDBECK18, LUNDBECK20, LitioSocio57, LitioSocio59,
         LUNDBECK23, LUNDBECK25, LitioSocio45, LitioSocio48, LitioSocio51, LitioSocio54)%>%
  filter(!is.na(LitioSocio94), LitioSocio94 != "NO ANALÍTICA")%>%
  mutate(across(c(BirthDay, LitioSocio01, LUNDBECK12),
                ~ str_extract(.x, "\\d{4}$")))%>%
  mutate(BirthDay = as.numeric(BirthDay),
         LitioSocio01= as.numeric(LitioSocio01),
         LUNDBECK12= as.numeric(LUNDBECK12))

colnames(df0) <- c("ID", "birth_year", "pheno_year", "initial_disease_year", "psych_traits","n_hosp",
                   "n_manic", "n_depressive", "n_mixed", "ever_psych", "actual_psych",
                   "n_sui", "rapid_cycle", "age_first_manic", "age_first_depressive", "age_first_hipo", "age_first_mixed")

#Manejar edad de comienzo
#LUNDBECK12 -> fecha inicio enfermedad
#LitioSocio45 - EDAD PRIMER EPISODIO MANIACO
#LitioSocio48 - EDAD PRIMER EPISODIO DEPRESIVO
#LitioSocio51 - EDAD PRIMER EPISODIO HIPOMANIACO
#LitioSocio54 - EDAD PRIMER EPISODIO MIXTO

df0$initial_disease_year #hay que convertirlo a edad
df0$age_first_manic #algunos son numeros y otros años y otros texto: los que son 0 pasarlos a Na y hay que ponerlo como numero
df0$age_first_depressive #algunos son numeros y otros años y otros texto: los que son 0 pasarlos a Na y hay que ponerlo como numero


clean_age <- function(x, birth_year) {
  # Convertir a character
  x <- as.character(x)
  # Convertir textos tipo "no sabe", "nunca", "lo desconoce", etc. a NA
  x[x == "no sabe, alrededor 35"] <- "35" #caso específico para hipomanic, porque lleva no sabe delante
  x[x == "2 por mes desde los 35"] <- "35" #caso específico para hipomanic
  x[grepl("\\bno\\b|nunca|desconoce|sabe|recuerda|2 al año", x, ignore.case = TRUE)] <- NA #lo de NO lo ponemos asi para que sea esa palabra solo, que no coja diagNOstico
  x[x == "marzo 2019 y con 26 años"] <- "26" #solo para depressive episodes, es el único donde sale esto y quiero que coja 26 años y no 2019
  # Extraer el primer número que aparezca
  num <- suppressWarnings(as.numeric(str_extract(x, "\\d+"))) #sacar el primer número
  # Si no hay número → NA
  num[is.na(num)] <- NA
  # Convertir 0 a NA
  num[num == 0] <- NA #si es 0, es NA
  # Detectar si es un año (>= 1900 y <= 2026)
  is_year <- num >= 1900 & num <= 2026
  is_year[is.na(is_year)] <- FALSE #si es NA, es FALSE
  # Convertir años a edad
  num[is_year] <- num[is_year] - birth_year[is_year]
  return(num)
}

df0 <- df0 %>% mutate(age_first_manic_clean = clean_age(age_first_manic, birth_year))
df0 <- df0 %>% mutate(age_first_depressive_clean = clean_age(age_first_depressive, birth_year))
df0 <- df0 %>% mutate(age_first_hipo_clean = clean_age(age_first_hipo, birth_year))
df0 <- df0 %>% mutate(age_first_mixed_clean = clean_age(age_first_mixed, birth_year))
df0 <- df0 %>% mutate(initial_disease_age = ifelse( #transformamos el año de inicio de enfermedad en edad
    !is.na(initial_disease_year) & initial_disease_year > 1900 & birth_year < 2025, #he quitado los dos que nacían en 2025 así
    initial_disease_year - birth_year,
    NA
  ))



#Obtener la edad de onset como el mínimo de primer episodio maniaco, mixto, hipomaniaco, depresivo y initial_disease_year
df0 <- df0 %>%
  mutate(age_onset = pmin( #nos quedamos con el mínimo de todas
    initial_disease_age,
    age_first_manic_clean,
    age_first_depressive_clean,
    age_first_hipo_clean,
    age_first_mixed_clean,
    na.rm = TRUE
  )) #la edad de onset va de 5 a 73 exito

df0 <- df0 %>% mutate(pheno_age = ifelse( #transformamos el año de inicio de enfermedad en edad
  !is.na(pheno_year) & pheno_year > 1900 & birth_year < 2025 & birth_year != pheno_year, #he quitado los dos que nacían en 2025 así
  pheno_year - birth_year, #la edad de febo
  NA
))

df0 <- df0 %>% mutate(disease_years = pheno_age - age_onset) #son NA los que son NA en alguna de las dos categorías

df0 <- df0 %>% select(-age_first_manic, -age_first_depressive, -age_first_hipo, -age_first_mixed, -initial_disease_age, -pheno_year)


#Seguir leyendo inputs

df0 <- df0 %>% #poner los identificadores con 4 digitospara mapear con los del otro excel
  mutate(
    ID_complete = str_replace(ID, "(.*_)(\\d+)$", function(x) {
      prefix <- str_match(x, "(.*_)(\\d+)$")[,2]
      num    <- str_match(x, "(.*_)(\\d+)$")[,3]
      paste0(prefix, str_pad(num, width = 4, pad = "0"))
    })
  )


inputF2 <- "Z:/BioBANCO/Genotyping_Santiago/Datos Crudos Genotyping CEL/FAM_BD-Controls_Wave1_MadManic.xlsx" #de este Excel sacamos: hospitalizaciones, psicosis, episodios hipomaniacos, episodios maniacos
df_F2 <- readxl::read_xlsx(inputF2, sheet = 1, col_types = "text") #se lee todo como texto para que lo coja

df_psico <- df_F2 %>% dplyr::select( #nos quedamos solo con las columnas que necesitamos
  ID, Edad,
  `Psychosis (0: Unkown; 1: without psychosis; 2: with psychosis)`,
  `hospitalizaciones (mania y depresion)`,
  N_episodios_maniacos,
  N_hipomania) %>% #cambiamos nombres de columnas
  dplyr::rename(
    psychosis = `Psychosis (0: Unkown; 1: without psychosis; 2: with psychosis)`,
    n_hosp = `hospitalizaciones (mania y depresion)`,
    n_manic = N_episodios_maniacos,
    n_hipomanic = N_hipomania)  %>% #cambiamos los IE para poder mapearlos
  mutate(    ID = gsub("^IEV_", "IE_", ID)) %>% #filtramos las muestras que son las que están en df0
  dplyr::filter(ID %in% df0$ID_complete)

setdiff(df0$ID_complete, df_psico$ID) #hay 20 individuos que están en el excel de madmanic pero no en el de genotipado, como no se genotiparon no los tenemos en cuenta

#filtramos las muestras que están genotipadas
df0 <- df0 %>% 
  dplyr::filter(ID_complete %in% df_psico$ID) #nos quedan 368 individuos

#Juntar los dos dataframes en 1 solo
df0_clean <- df0 %>% #se eliminan las columnas comunes
  select(
    -n_hosp,
    -n_manic,
    -ever_psych,
    -actual_psych,
    -psych_traits
  )

df_complete <- df0_clean %>% #se unen los dos df
  left_join(df_psico, by = c("ID_complete" = "ID"))

df_complete <- df_complete %>%
  mutate(across(
    c(n_manic, n_hipomanic),
    ~ na_if(na_if(trimws(.), ""), "NA")
  ))


#Variables dicotomicas: psicosis (sale de df_psico) y ciclacion rapida (sale de MadManic)
df_dich <- df_complete %>%
  mutate(rapid_cycle = case_when(
    rapid_cycle == 1 ~ "NO",
    rapid_cycle ==2 ~ "YES",
    rapid_cycle ==3 ~ "YES",
    TRUE ~ NA_character_
  ))%>% 
  mutate(rapid_cycl_count = if_else((rapid_cycle == "YES"), 1, 0))%>%
  select(-rapid_cycle) %>%
  mutate(
    psychosis = case_when(
      psychosis == 0 ~ NA,#0 es unknown
      psychosis == 1 ~ 0, #1 es no psicosis, asi que ponemos 0
      psychosis == 2 ~ 1, #2 es psicosis, asi que ponemos 1
      TRUE ~ NA_real_
    ))


#Variables dicotómicas

df_multi <- df_dich %>%
  mutate(
    n_hosp_clean = case_when( #detectar la frase especial ANTES de extraer números
      str_detect(n_hosp, "1ro con 17 años") ~ 6,   #caso especial
      TRUE ~ NA_real_
    )
  ) %>%
  mutate( #extraer número SOLO si no es el caso especial
    n_hosp = case_when(
      !is.na(n_hosp_clean) ~ n_hosp_clean,  # mantener el 6
      TRUE ~ str_extract(n_hosp, "\\d+") |> as.numeric()
    )
  ) %>%
  mutate( #establecer categorias
    n_hosp_cat = case_when(
      n_hosp == 0 ~ 0,
      n_hosp %in% 1:2 ~ 1,
      n_hosp %in% 3:5 ~ 2,
      n_hosp > 5 ~ 3,
      TRUE ~ NA_real_
    )
  )%>%
  mutate( #pasamos a intentos de suicidio
  n_sui = case_when(
    str_detect(n_sui, regex("varios", ignore_case = TRUE)) ~ 6, #si pone varios, es más de 5
    TRUE ~ str_extract(n_sui, "\\d+") |> as.numeric()
  )) %>%
  mutate(
  n_sui_cat = case_when(
    n_sui == 0 ~ 0,
    n_sui %in% 1:2 ~ 1,
    n_sui %in% 3:5 ~ 2,
    n_sui > 5 ~ 3,
  TRUE ~ NA_real_
  )) %>% 
  mutate(
    n_manic = str_extract(n_manic, "\\d+") |> as.numeric()
  ) %>%
  mutate (
    n_manic_cat = case_when(
      n_manic == 0 ~ 0,
      n_manic %in% 1:2 ~ 1,
      n_manic %in% 3:5 ~ 2,
      n_manic > 5 ~ 3,
      TRUE ~ NA_real_
    )) %>% 
  mutate(
    n_hipomanic = str_extract(n_hipomanic, "\\d+") |> as.numeric()
  ) %>%
  mutate (
    n_hipomanic_cat = case_when(
      n_hipomanic == 0 ~ 0,
      n_hipomanic %in% 1:2 ~ 1,
      n_hipomanic %in% 3:5 ~ 2,
      n_hipomanic > 5 ~ 3,
      TRUE ~ NA_real_
    )) %>%
  select(-n_manic, -n_hipomanic)

colnames(df_multi)

df_multi_final <- df_multi %>% 
  mutate(
    n_dep_clean = case_when( #Detectar casos especiales que implican >5 episodios
      str_detect(n_depressive, regex("muchisimos|múltiples|varias fases|imposible|crónica|mas de|describe|al año|continuado|mas de 3|5-6|pequeños", ignore_case = TRUE)) ~ 6,
      str_detect(n_depressive, regex("primer brote", ignore_case = TRUE)) ~ 2,
      str_detect(n_depressive, regex("no recuerda haber tenido ninguno", ignore_case = TRUE)) ~ 0,
      str_detect(n_depressive, regex("no recuerda|mirar|dice no recordar", ignore_case = TRUE)) ~ NA_real_,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    n_dep_num = case_when( #Extraer número solo si no es caso especial
      !is.na(n_dep_clean) ~ n_dep_clean,
      TRUE ~ str_extract(n_depressive, "\\d+") |> as.numeric() #comprobar 37, 97, 113, 124, 130, 151, 155, 164, 171, 204, 278, 313 es NA, 360
    )
  ) %>%
  mutate(
    # Recodificar en categorías
    n_depressive_cat = case_when(
      is.na(n_dep_num) ~ NA_real_,
      n_dep_num == 0 ~ 0,
      n_dep_num %in% 1:2 ~ 1,
      n_dep_num %in% 3:5 ~ 2,
      n_dep_num > 5 ~ 3
    )
  ) %>%
  mutate(
    n_mixed_clean = case_when( #Detectar casos especiales que implican >5 episodios
      str_detect(n_mixed, regex("varios", ignore_case = TRUE)) ~ 6,
      str_detect(n_mixed, regex("no recuerda", ignore_case = TRUE)) ~ NA_real_,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    n_mixed_num = case_when( #Extraer número solo si no es caso especial
      !is.na(n_mixed_clean) ~ n_mixed_clean,
      TRUE ~ str_extract(n_mixed, "\\d+") |> as.numeric()
    )
  ) %>%
  mutate(
    n_mixed_cat = case_when( #Recodificar en categorías
      is.na(n_mixed_num) ~ NA_real_,
      n_mixed_num == 0 ~ 0,
      n_mixed_num %in% 1:2 ~ 1,
      n_mixed_num %in% 3:5 ~ 2,
      n_mixed_num > 5 ~ 3
    )
  ) %>%
  select(-n_mixed_clean, -n_mixed_num, -n_dep_num, -n_dep_clean, -n_depressive, -n_mixed, -n_hosp_clean, -n_hosp, -n_sui)

colnames(df_multi_final)

#Extract the columns we need
df_severity_index <- df_multi_final %>% select(ID, ID_complete, disease_years, psychosis, rapid_cycl_count, n_hosp_cat, n_sui_cat, n_manic_cat, n_hipomanic_cat, n_depressive_cat, n_mixed_cat)

#Remove individuals with more than 3 NAs
df_severity_index <- df_severity_index %>% filter(rowSums(is.na(.)) <= 3)

#remove individuals with NA in disease_years
df_severity_index <- df_severity_index %>% filter(!is.na(disease_years) & disease_years > 1) #se quitan los que tienen 0 o 1 año de enfermedad



############
#Imputation#
############


complete_data <- df_severity_index[complete.cases(df_severity_index),] #51 con todos los datos
set.seed(1234)
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
#points(6,rmse[5],col="red",pch=20,cex=2)
#abline(h=rmse[5],lty=2,col="red")

results <- data.frame(ks,rmse)

best_index <- which.min(rmse)
best_k <- ks[best_index] #16 es la mejor k con 0.2 missing values  

#imputation
df_imputed_severity_index <- kNN(df_severity_index,k=16,imp_var=FALSE)

########################################
#Index of Number of Events and Severity#
########################################

df_imputed_severity_index <- df_imputed_severity_index %>%
  mutate(
    severity_index =
      (rapid_cycl_count +
         n_hosp_cat +
         n_sui_cat +
         n_manic_cat +
         n_depressive_cat +
         n_mixed_cat +
         psychosis +
         n_hipomanic_cat) / log2(disease_years + 1) #ELIMINAR LOS QUE TIENEN DISEASE YEARS = 0 PARA QUE NO DE INFINITO
  )

#Represent the distribution
df_imputed_severity_index %>%
  ggplot(aes(x = severity_index)) +
  geom_histogram(bins = 30, fill = "#4C72B0", color = "white", alpha = 0.8) +
  labs(
    title = "Distribución del Index of Number of Events and Severity (INES)",
    x = "Severity Index",
    y = "Frecuencia"
  ) +
  theme_minimal(base_size = 14)

df_imputed_severity_index %>%
  ggplot(aes(x = severity_index)) +
  geom_density(fill = "#4C72B0", color = "white", alpha = 0.8) +
  labs(
    title = "Distribución del Index of Number of Events and Severity (INES)",
    x = "Severity Index",
    y = "Frecuencia"
  ) +
  theme_minimal(base_size = 14)


#######################################
#Extract EPIC individuals and imputate#
#######################################

input_EPIC <- "Z:/BioBANCO/EPIC_V2/Selection_239_10_EPIC.xlsx" #de este Excel sacamos: hospitalizaciones, psicosis, episodios hipomaniacos, episodios maniacos
EPIC_individuals <- readxl::read_xlsx(input_EPIC, sheet = 1)[,3] #se lee todo como texto para que lo coja

df_EPIC_final <- df_multi_final %>%
  filter(ID_complete %in% EPIC_individuals)

colSums(is.na(df_EPIC_final)) #numero de NA en cada cateogoria

NA_por_individuo <- df_EPIC_final %>% #numero de NA en cada individuo
  mutate(NA_count = rowSums(is.na(.))) %>%
  select(ID_complete, NA_count)

sum(NA_por_individuo$NA_count != 0) #79 de los 130 tienen al menos 1 NA

df_epic_for_imputation <- df_EPIC_final %>%
  mutate(NA_count = rowSums(is.na(.))) %>%
  filter(NA_count <= 3)

# import data -----

data <- df_epic_for_imputation %>%
  mutate(disease_years = as.numeric(pheno_year) - as.numeric(initial_disease_year) + 1) %>% 
  dplyr::select(-pheno_year, -initial_disease_year, -ID, -birth_year, -NA_count)

data <- data %>% #se transofrman todos los 
  mutate(across(-ID_complete, as.numeric))
  

# cross-validation -----
complete_data <- data[complete.cases(data),] #51 con todos los datos
set.seed(1234)
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
points(6,rmse[5],col="red",pch=20,cex=2)
abline(h=rmse[5],lty=2,col="red")

results <- data.frame(ks,rmse)

best_index <- which.min(rmse)
best_k <- ks[best_index]  
# imputation for EPIC

imputed_epic_data_k6 <- kNN(data,k=6,imp_var=FALSE)

#Obtain the severiry index
epic_with_severity_index <- imputed_epic_data_k6 %>%
  mutate(
    severity_index =
      (rapid_cycl_count +
         n_hosp_cat +
         n_sui_cat +
         n_manic_cat +
         n_depressive_cat +
         n_mixed_cat +
         psychosis +
         n_hipomanic_cat) / log2(disease_years + 1)
  )

#Represent the distribution
epic_with_severity_index %>%
  ggplot(aes(x = severity_index)) +
  geom_histogram(bins = 30, fill = "#4C72B0", color = "white", alpha = 0.8) +
  labs(
    title = "Distribución del Severity Index",
    x = "Severity Index",
    y = "Frecuencia"
  ) +
  theme_minimal(base_size = 14)

#Save EPIC severity indexes
write.xlsx(epic_with_severity_index, "D:/severity_index_EPIC_samples.xlsx")



#Anotaciones de texto (no ejecutar)

#df_multi$n_depressive 
#2 (pequeños unos 20) es más de 5
#no recuerda es más de 5
#"No sabe definir, episodios en concreto, es algo más continuado en el tiempo" es más de 5 
#"no recuerda haber tenido ninguno, no aparece nada en su historial" es 0
#"Diagnosticada de depresión crónica, no se puede especificar" es más de 5
#"dice no recordar uno como tal" es NA
#"Imposible definir" es más de 5
#"no especificado en historía clínica, donde consta: \"varias fases depresivas\"" es más de 5
#df_multi$n_mixed 
#"en la historia aparecen varios pero ella no los reconoce" es más de 5
#"no recuerda" es NA


#Problema: los NA reemplazarlos con knn 
#Si algun indiviudo tiene más de 3 missing (datos en menos de 5 escalas) se elimina

