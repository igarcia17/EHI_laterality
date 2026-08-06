#Igual que el script 3_prevalence.R, funcion para contar por grupos

library(dplyr)
library(ggplot2)

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
input_clean <- subset(input, is.finite(Age))
set.seed(13)
match_obj <- MatchIt::matchit(BD_patient ~ Age, 
                              data = input_clean, 
                              method = "nearest", 
                              ratio = 1)
balanced_age <- MatchIt::match.data(match_obj)
count_by_group <- function(df, cutoff, var = NULL) {
  
  allowed_vars <- c("Sex", "BD_patient", "Diagnostic", "group_Age")
  # función interna para crear grupos
  make_groups <- function(data) {
    data %>%
      mutate(
        group = case_when(
          Score10items < -cutoff ~ "left",
          Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
          Score10items > cutoff ~ "right"))}
  
  if (is.null(var) || var == "None") {
    df2 <- make_groups(df)
    return(
      df2 %>%
        group_by(group) %>%
        summarise(n = n(), .groups = "drop") %>%
        # añadir lateralized
        bind_rows(
          df2 %>%
            filter(group %in% c("left", "right")) %>%
            summarise(n = n()) %>%
            mutate(group = "lateralized")) %>%
        mutate(
          total = sum(n[group != "lateralized"]),
          perc = 100 * n / total,
          label = sprintf("%d (%.2f%%)", n, perc)
        ) %>%
        select(group, label) %>%
        arrange(match(group, c("left", "mixed", "right", "lateralized"))))}
  
  if (!(var %in% allowed_vars)) {
    stop("var debe ser una de: Sex, BD_patient, Diagnostic o NULL/None")
  }
  df2 <- make_groups(df)
  res <- df2 %>%
    group_by(group, .data[[var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[var]]) %>%
    mutate(
      total_var = sum(n),
      perc = 100 * n / total_var,
      label = sprintf("%d (%.2f%%)", n, perc)) %>% ungroup()
  # añadir lateralized por cada nivel de var
  lateralized <- df2 %>%
    filter(group %in% c("left", "right")) %>%
    group_by(.data[[var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(.data[[var]]) %>%
    mutate(
      total_var = sum(res$n[res[[var]] == .data[[var]]], na.rm = TRUE),
      perc = 100 * n / total_var,
      group = "lateralized",
      label = sprintf("%d (%.2f%%)", n, perc)
    ) %>% select(group, all_of(var), label)
  bind_rows(res, lateralized) %>%
    arrange(.data[[var]], match(group, c("left", "mixed", "right", "lateralized")))
}

#With all individuals
count_by_group(input, 90, "BD_patient")

#With matched-age controls
count_by_group(balanced_age, 90)
count_by_group(balanced_age, 90, "BD_patient")

#Only for controls
controls <- input %>% filter(BD_patient =="NO") %>% mutate(group_Age = if_else(Age >= 40, "OLD", "YOUNG"),
                                                           group_Age = factor(group_Age, levels = c("YOUNG", "OLD")))
summary(controls)
count_by_group(controls, 0, "Sex")
count_by_group(controls, 90, "group_Age")

#Withing BD
onlyfem <- balanced_age %>% filter(Sex=="MUJER")
onlymale <- balanced_age%>% filter(Sex=="HOMBRE")
count_by_group(onlyfem, 90)
count_by_group(onlyfem, 90, "BD_patient")

count_by_group(onlymale, 90)
count_by_group(onlymale, 90, "BD_patient")







