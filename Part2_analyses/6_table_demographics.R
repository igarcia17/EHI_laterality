#Hacer tablas demograficas descriptivas de la muestra
library(dplyr)
library(ggplot2)
library(tableone)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())


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

#Aqui se ve que efectivamente si usamos todos los individuos los grupos de BD y control no están iguales en edad
shapiro.test(input$Age)
shapiro.test(input$Score10items)
shapiro.test(input$Item11)
shapiro.test(input$Item12)

# 1. Definimos las variables que queremos en la tabla
variables <- c("Sex", "Age", "Score10items", "Item11", "Item12")

# 2. Especificamos cuáles de ellas son categóricas (factores)
categoricas <- c("Sex")

# Crear la tabla base
tabla1 <- CreateTableOne(
  vars = variables, 
  strata = "BD_patient", 
  data = input, 
  factorVars = categoricas,
  addOverall = TRUE  # Para que incluya la columna "Total"
)

# Mostrar la tabla con un formato limpio
 t <- print(tabla1, 
      showAllLevels = TRUE, 
      nonnormal = c("Age", "Score10items", "Item11", "Item12"),
      formatOptions = list(big.mark = ","), pDigits=16
      ) 
write.csv(t, "Tabla1_Demograficos.csv")

#Esto es equivalenete a usar la tabla que ya se guardó porque usa la misma seed
input_clean <- subset(input, is.finite(Age))
set.seed(13)
match_obj <- MatchIt::matchit(BD_patient ~ Age, 
                              data = input_clean, 
                              method = "nearest", 
                              ratio = 1)
balanced_age <- MatchIt::match.data(match_obj)

#Se rehace la tabla demografica con age-matched indivs
tabla2 <- CreateTableOne(
  vars = variables, 
  strata = "BD_patient", 
  data = balanced_age, 
  factorVars = categoricas,
  addOverall = TRUE  # Para que incluya la columna "Total"
)

t2 <- print(tabla2, 
           showAllLevels = TRUE, 
           nonnormal = c("Age", "Score10items", "Item11", "Item12"),
           formatOptions = list(big.mark = ","), pDigits=16
) 
write.csv(t2, "Tabla1_Demograficos_matchedAge.csv")

#Se explora de donde vienen las diferencias
set_BDI <- input_clean[input_clean$Diagnostic =="BD-I",]%>% droplevels()
set_BDII <- input_clean[input_clean$Diagnostic =="BD-II",]%>% droplevels()
set_cont <- input_clean[input_clean$Diagnostic =="NO",] %>% droplevels()
set_cases <- input_clean[input_clean$BD_patient =="YES",]%>% droplevels()
set_cont_balanced <- balanced_age[balanced_age$Diagnostic == "NO",]%>% droplevels()
summary(set_BDI$Age)
summary(set_BDII$Age) #Sobre todo por BD-II
summary(set_cases$Age)
summary(set_cont_balanced$Age)
summary(set_cont$Age)


#Se rehace la tabla demografica solo con BD
tabla3 <- CreateTableOne(
  vars = variables, 
  strata = "Diagnostic", 
  data = set_cases, 
  factorVars = categoricas,
  addOverall = TRUE  # Para que incluya la columna "Total"
)

t3 <- print(tabla3, 
            showAllLevels = TRUE, 
            nonnormal = c("Age", "Score10items", "Item11", "Item12"),
            formatOptions = list(big.mark = ","), pDigits=16
) 
write.csv(t3, "Tabla1_Demograficos_onlyCases.csv")

#Plots
jpeg(filename='Age_distribution_by_diagnsotic_including_controls.jpeg')
ggplot(input_clean, aes(x = Age, fill = Diagnostic)) +
  geom_density(alpha = 0.5) +  # alpha da transparencia para ver dónde se solapan
  scale_fill_manual(values = c("BD-I" = "slateblue2", "BD-II" = "pink",
                               "NO" = "darkorange", "SZAFF" = "red", "UNKNOWN" = "grey")) + # Colores personalizados
  labs(
    title = "Distribución de Edad por Diagnóstico",
    x = "Edad",
    y = "Densidad",
    fill = "Grupo"
  ) +
  theme_minimal()
dev.off()

jpeg(filename='Age_distribution_by_diagnsotic_including_balanced_controls.jpeg')
ggplot(balanced_age, aes(x = Age, fill = Diagnostic)) +
  geom_density(alpha = 0.5) +  # alpha da transparencia para ver dónde se solapan
  scale_fill_manual(values = c("BD-I" = "slateblue2", "BD-II" = "pink",
                               "NO" = "darkorange", "SZAFF" = "red", "UNKNOWN" = "grey")) + # Colores personalizados
  labs(
    title = "Distribución de Edad por Diagnóstico",
    x = "Edad",
    y = "Densidad",
    fill = "Grupo"
  ) +
  theme_minimal()
dev.off()

jpeg(filename='Age_distribution_by_BD_vs_controls_including_controls.jpeg')
ggplot(input_clean, aes(x = Age, fill = BD_patient)) +
  geom_density(alpha = 0.5) +  # alpha da transparencia para ver dónde se solapan
  scale_fill_manual(values = c("YES" = "slateblue2", 
                               "NO" = "darkorange")) + # Colores personalizados
  labs(
    title = "Distribución de Edad por status",
    x = "Edad",
    y = "Densidad",
    fill = "BD patient"
  ) +
  theme_minimal()
dev.off()

jpeg(filename='Age_distribution_by_BD_vs_controls_including_balanced_controls.jpeg')
ggplot(balanced_age, aes(x = Age, fill = BD_patient)) +
  geom_density(alpha = 0.5) +  # alpha da transparencia para ver dónde se solapan
  scale_fill_manual(values = c("YES" = "slateblue2",
                               "NO" = "darkorange")) + # Colores personalizados
  labs(
    title = "Distribución de Edad por status",
    x = "Edad",
    y = "Densidad",
    fill = "BD patient"
  ) +
  theme_minimal()
dev.off()
