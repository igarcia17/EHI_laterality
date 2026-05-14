library(dplyr)
library(ggplot2)
library(tableone)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())


input <- as.data.frame(readxl::read_xlsx("../input_prep/all_samples_all_data_LQ_PLUS_items.xlsx", sheet = 1, 
                                         col_types="text", na="#N/A"))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>%
  mutate(across(
    -c(ID, Diagnostic),                    
    ~ if (!is.factor(.)) as.numeric(.) else .
  )) %>% mutate(Diagnostic = factor(Diagnostic))


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
write.csv(t, "Tabla1_Demograficose.csv")

input_clean <- subset(input, is.finite(Age))
set.seed(13)
match_obj <- MatchIt::matchit(BD_patient ~ Age, 
                              data = input_clean, 
                              method = "nearest", 
                              ratio = 1)
balanced_age <- MatchIt::match.data(match_obj)
tabla2 <- CreateTableOne(
  vars = variables, 
  strata = "BD_patient", 
  data = balanced_age, 
  factorVars = categoricas,
  addOverall = TRUE  # Para que incluya la columna "Total"
)

# Mostrar la tabla con un formato limpio
t2 <- print(tabla2, 
           showAllLevels = TRUE, 
           nonnormal = c("Age", "Score10items", "Item11", "Item12"),
           formatOptions = list(big.mark = ","), pDigits=16
) 
write.csv(t2, "Tabla1_Demograficose_matchedAge.csv")

set_BDI <- input[input$Diagnostic =="BD-I",]
set_BDII <- input[input$Diagnostic =="BD-II",]
set_cont <- input[input$Diagnostic =="NO",]
set_cases <- input[input$BD_patient =="YES",]
summary(set_BDI$Age)
summary(set_BDII$Age)
summary(set_cases$Age)
summary(set_cont$Age)


ggplot(input, aes(x = Age, fill = Diagnostic)) +
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




