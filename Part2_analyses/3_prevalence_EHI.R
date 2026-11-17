library(dplyr)
library(tidyr)
library(ggplot2)
#Para funciones basicas de prevalencia y hacer un subset de individuos age-matched
#Para graficos para la presentacion del congreso
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())
source("../functions.R")

#Load input data
input <- as.data.frame(readxl::read_xlsx("../Part1_input_prep/all_samples_all_data_LQ_PLUS_items.xlsx", sheet = 1, 
                                         col_types="text", na="#N/A"))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>%
  mutate(across(
    -c(ID, Diagnostic),                    
    ~ if (!is.factor(.)) as.numeric(.) else .)) %>% mutate(Diagnostic = factor(Diagnostic))

input <- subset(input, is.finite(Age))

#Contar uso de manos en general, comparando entre casos y controles, sin match por edad
count_by_group(input, cutoff = 40, var = "Sex")

###Contar uso de manos en hombres y en mujeres, comparando entre casos y controles, sin match por edad
females <- input%>%  filter(Sex=="MUJER")%>%droplevels()
males <- input%>%filter(Sex=="HOMBRE")%>%droplevels()

count_by_group(females, cutoff=60, var="BD_patient")
count_by_group(males, cutoff=60, var="BD_patient")  

#Proportion tests
prop_test_by_group(input, cutoff = 0, var = "Sex")

#Taking Age into play
jpeg(filename='../Age_in_BD_patients.jpeg')
ggplot(
  input, aes(x = Age, fill = BD_patient)
) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Age distribution in BD",
    x = "Age",
    y = "Density"
  ) +
  theme_minimal()
dev.off()

####Match controls and cases
input_clean <- subset(input, is.finite(Age))
set.seed(13)
match_obj <- MatchIt::matchit(BD_patient ~ Age, 
                     data = input_clean, 
                     method = "nearest", 
                     ratio = 1)
balanced_age <- MatchIt::match.data(match_obj)

jpeg(filename='../After_balance_age_in_BD_patients.jpeg')
ggplot(
  balanced_age, aes(x = Age, fill = BD_patient)
) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Age distribution in BD",
    x = "Age",
    y = "Density"
  ) +
  theme_minimal()
dev.off()

writexl::write_xlsx(balanced_age, path= "balanced_samples_all_data_LQ.xlsx")
count_by_group(balanced_age, cutoff=80, var="BD_patient") 

#Balanced age LQ distrib
jpeg("All_samples_balanced_age_LQ_dist.jpeg")
ggplot(balanced_age, aes(x = Score10items)) +
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
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  labs(
    title = "EHI LQ histogram distribution",
    x = "Score",
    y = "Frecuencia"
  ) +
  theme_minimal()
dev.off()

#Separated histogram
jpeg("Before_balanced_By_BD_all_samples_balanced_age_LQ_dist.jpeg")
ggplot(input, aes(x = Score10items, fill = BD_patient)) +
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
  scale_x_continuous(breaks = seq(-100, 100, 20)) +
  scale_fill_manual(values = c("NO" = "lightsalmon", "YES" = "slateblue2")) +
  labs(
    title = "LQ distribution by BD-nonBD",
    x = "Score",
    y = "N",
    fill = "BD_patient"
  ) +
  theme_minimal()
dev.off()

#Boxplot of hand and foot use
jpeg("Foot_and_eye_preference.jpeg")
density_df <- balanced_age %>%
    filter(
          (is.na(Item11) | Item11 != 3.5) &
                (is.na(Item12) | Item12 != 3.5)
      ) %>%
    select(BD_patient, Item11, Item12) %>%
  tidyr::pivot_longer(
          cols = c(Item11, Item12),
          names_to = "Item",
          values_to = "Score"
      ) %>%
    filter(!is.na(Score), !is.na(BD_patient))
ggplot(density_df,
            aes(x = Score,
                             fill = factor(BD_patient),
                             color = factor(BD_patient))) +
    geom_density(alpha = 0.3, linewidth = 1, adjust =1.5) +
    facet_wrap(~ Item) +
    theme_minimal() +
    labs(
          x = "Score",
          y = "Density",
          fill = "Group",
          color = "Group"
      )
dev.off()

input <- create_laterality_vars(input)
balanced_age <- create_laterality_vars(balanced_age)
m_nrh <- glm(NRH_0 ~ BD_patient+Age+Sex, data = balanced_age, family=binomial)
m_lat <- glm(LAT_90 ~BD_patient+Age+ Sex, data= balanced_age, family=binomial)
summary(m_lat)

#Tabla 3
females_clean <- balanced_age%>%
  filter(Sex=="MUJER")%>%droplevels()
males_clean <- balanced_age%>%filter(Sex=="HOMBRE")%>%droplevels()

count_by_group(females_clean, cutoff=0)
count_by_group(males_clean, cutoff=0)

count_by_group(females_clean, cutoff=0, var="BD_patient")
count_by_group(males_clean, cutoff=0, var="BD_patient")  

prop_test_by_group(females_clean, cutoff = 0, var="BD_patient")

summary(glm(NRH_0 ~ BD_patient+Age, data = females_clean, family=binomial))
summary(glm(LAT_90 ~BD_patient+Age, data= males_clean, family=binomial))

#Age on controls
all_controls <- input_clean %>% filter(BD_patient == "NO") %>%droplevels()%>%
  mutate(
    lateralized_90 = ifelse(abs(Score10items) > 90, 1, 0),
    NRH_60 = ifelse(Score10items <= 60, 1, 0)
  )

age_prop <- all_controls %>%
  filter(
    !is.na(Age),
    !is.na(Score10items)
  ) %>%
  mutate(
    Age_group = cut(
      Age,
      breaks = c(-Inf, 20, 30, 40, 50, 60, 70, 80, Inf),
      labels = c("≤20", "21–30", "31–40", "41–50",
                 "51–60", "61–70", "71–80", "≥81")
    )
  ) %>%
  group_by(Age_group) %>%
  summarise(
    n_total = n(),
    n_90 = sum(abs(Score10items) > 90),
    prop_90 = n_90 / n_total,
    .groups = "drop"
  )

ggplot(age_prop, aes(x = Age_group, y = prop_90)) +
  geom_col(width = 0.8, fill = "#E76F51") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    x = "Age group",
    y = "Lateralized individuals (%)"
  )+theme(
  text = element_text(size = 20)
)

#Aqui se ve que la edad influencia mucho en la lateralizacion cuando 90,
# y realmente en NRH no muchisimo 0.086

m <- glm(lateralized_90 ~ Age+Sex, data = all_controls)
m2 <- glm(NRH_60 ~Age +Sex, data = all_controls)
summary(m)
summary(m2)


