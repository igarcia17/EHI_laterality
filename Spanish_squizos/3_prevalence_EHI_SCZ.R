library(dplyr)
library(ggplot2)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())
###Functions
count_by_group <- function(df, cutoff, var = NULL) {
  allowed_vars <- c("Sex", "BD_patient", "Diagnostic")
  # función interna para crear grupos
  make_groups <- function(data) {
    data %>%
      mutate(
        group = case_when(
          Score10items < - cutoff ~ "left",
          Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
          Score10items > cutoff ~ "right"
        )
      )
  }
  
  # =========================
  # CASO 1: SIN VARIABLE
  # =========================
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
            mutate(group = "lateralized")
        ) %>%
        
        mutate(
          total = sum(n[group != "lateralized"]),
          perc = 100 * n / total,
          label = sprintf("%d (%.2f%%)", n, perc)
        ) %>%
        select(group, label) %>%
        arrange(match(group, c("left", "mixed", "right", "lateralized")))
    )
  }
  
  # =========================
  # CASO 2: CON VARIABLE
  # =========================
  
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
      label = sprintf("%d (%.2f%%)", n, perc)
    ) %>%
    ungroup()
  
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
    ) %>%
    select(group, all_of(var), label)
  
  bind_rows(res, lateralized) %>%
    arrange(.data[[var]], match(group, c("left", "mixed", "right", "lateralized")))
}
#####Load input
input <- as.data.frame(readxl::read_xlsx("all_samples_all_data_LQ_PLUS_items_included_SCZ.xlsx", sheet = 1))%>%
  mutate(Score10items =as.numeric(Score10items), 
         Sex = factor(Sex),
         BD_patient = factor(BD_patient),
         Diagnostic = factor(Diagnostic),
         Age = as.numeric(Age)) %>%
  mutate(across(
    -c(ID),                    
    ~ if (!is.factor(.)) as.numeric(.) else .
  )) %>% mutate(Diagnostic = factor(Diagnostic))
  
levels(input$BD_patient) <- c("CONTROL", "SCZ", "BD")


count_by_group(input, cutoff = 0, var = "BD_patient")
count_by_group(input, cutoff = 60, var = "BD_patient")
count_by_group(input, cutoff = 90, var = "BD_patient")



### Analysis by sex

females <- input%>%  filter(Sex=="MUJER")%>%droplevels()
males <- input%>%filter(Sex=="HOMBRE")%>%droplevels()

count_by_group(females, cutoff=90, var="BD_patient")
count_by_group(males, cutoff=90, var="BD_patient")  
  
#Taking Age into play
jpeg(filename='Age_in_BD_patients.jpeg')
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

out_BD <- input%>%
  filter(!(BD_patient == "BD"))%>%droplevels()
####Match controls and cases
input_clean <- subset(out_BD, is.finite(Age))
set.seed(13)
match_obj <- MatchIt::matchit(BD_patient ~ Age, 
                     data = input_clean, 
                     method = "nearest", 
                     ratio = 1)
balanced_age <- MatchIt::match.data(match_obj)

jpeg(filename='balance_age_in_SCZ_patients.jpeg')
ggplot(
  balanced_age, aes(x = Age, fill = BD_patient, color=BD_patient)
) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Age distribution in SCZ",
    x = "Age",
    y = "Density"
  ) +
  theme_minimal()
dev.off()

writexl::write_xlsx(balanced_age, path= "balanced_samples_all_data_LQ_SCZ.xlsx")
count_by_group(balanced_age, cutoff=0, var="BD_patient") 

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

females_clean <- balanced_age%>%
  filter(Sex=="MUJER")%>%droplevels()
males_clean <- balanced_age%>%filter(Sex=="HOMBRE")%>%droplevels()

count_by_group(females_clean, cutoff=0)
count_by_group(males_clean, cutoff=0)  

count_by_group(females_clean, cutoff=0, var="BD_patient")
count_by_group(males_clean, cutoff=0, var="BD_patient")  





