#Demographics of SCZ and comparison
library(dplyr)
library(ggplot2)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

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

#Before everything, demographics of SCZ
#Remove items with NAs in sex and age
scz_alone <- input%>%filter(BD_patient=="SCZ",   !is.na(Sex),
                            !is.na(Age))%>%droplevels()
summary(scz_alone)

scz_demographics <- scz_alone %>%
  summarise(
    N = n(),
    `Female (%)` = sprintf(
      "%d (%.1f%%)",
      sum(Sex == "MUJER"),
      100 * mean(Sex == "MUJER")
    ),
    `Age (median [IQR])` = sprintf(
      "%.1f [%.1f–%.1f]",
      median(Age),
      quantile(Age, 0.25),
      quantile(Age, 0.75)
    ),
    `LQ (median [IQR])` = sprintf(
      "%.1f [%.1f–%.1f]",
      median(Score10items),
      quantile(Score10items, 0.25),
      quantile(Score10items, 0.75)
    )
  )

print(scz_demographics)

# Remove individuals with missing Sex or Age
df_comparison <- input %>%
  filter(
    BD_patient %in% c("CONTROL", "SCZ", "BD"),
    !is.na(Sex),
    !is.na(Age)
  ) %>%
  droplevels()


# Function to compare two groups
compare_groups <- function(data, group1, group2) {
  
  df <- data %>%
    filter(BD_patient %in% c(group1, group2)) %>%
    droplevels()
  
  # Sex: Chi-squared test
  p_sex <- chisq.test(
    table(df$Sex, df$BD_patient)
  )$p.value
  
  # Age: Wilcoxon rank-sum test
  p_age <- wilcox.test(
    Age ~ BD_patient,
    data = df,
    exact = FALSE
  )$p.value
  
  # LQ: Wilcoxon rank-sum test
  p_lq <- wilcox.test(
    Score10items ~ BD_patient,
    data = df,
    exact = FALSE
  )$p.value
  
  return(c(
    Sex = p_sex,
    Age = p_age,
    LQ = p_lq
  ))
}

# SCZ vs Controls
p_scz_controls <- compare_groups(
  df_comparison, "SCZ", "CONTROL"
)

# SCZ vs BD
p_scz_bd <- compare_groups(
  df_comparison, "SCZ", "BD"
)

# Results table
pvalue_table <- data.frame(
  Variable = c("Sex", "Age", "LQ"),
  `SCZ vs Controls` = format.pval(
    p_scz_controls,
    digits = 3,
    eps = 1e-16
  ),
  `SCZ vs BD` = format.pval(
    p_scz_bd,
    digits = 3,
    eps = 1e-16
  ),
  check.names = FALSE
)

print(pvalue_table)

#With matched controls to BD
#For all data, not matched: ../Part1_input_prep/all_samples_all_data_LQ_PLUS_items.xlsx
#For matched data: ../Part2_analyses/balanced_samples_all_data_LQ.xlsx
balanced_controls <- as.data.frame(readxl::read_xlsx("../Part1_input_prep/all_samples_all_data_LQ_PLUS_items.xlsx", sheet = 1))%>%
  mutate(Score10items =as.numeric(Score10items), 
         Sex = factor(Sex),
         BD_patient = factor(BD_patient),
         Diagnostic = factor(Diagnostic),
         Age = as.numeric(Age)) %>%
  mutate(across(
    -c(ID),                    
    ~ if (!is.factor(.)) as.numeric(.) else .
  )) %>% mutate(Diagnostic = factor(Diagnostic))

levels(balanced_controls$BD_patient) <- c("CONTROL", "BD")

# ==========================================================
# SCZ vs BD and Controls (matched dataset)
# ==========================================================

# Check group labels before merging
table(balanced_controls$BD_patient, useNA = "ifany")
table(scz_alone$BD_patient, useNA = "ifany")

# Make sure both datasets have the same group labels
balanced_controls <- balanced_controls %>%
  mutate(BD_patient = as.character(BD_patient))

scz_alone <- scz_alone %>%
  mutate(BD_patient = as.character(BD_patient))

# Merge SCZ with matched BD and controls
df_matched_scz <- bind_rows(
  balanced_controls,
  scz_alone
) %>%
  filter(
    BD_patient %in% c("CONTROL", "BD", "SCZ"),
    !is.na(Sex),
    !is.na(Age)
  ) %>%
  mutate(
    BD_patient = factor(
      BD_patient,
      levels = c("CONTROL", "BD", "SCZ")
    )
  ) %>%
  droplevels()

# Check sample sizes
table(df_matched_scz$BD_patient)

# ==========================================================
# Statistical comparisons
# ==========================================================

# SCZ vs matched controls
p_scz_matched_controls <- compare_groups(
  df_matched_scz, "SCZ", "CONTROL"
)

# SCZ vs BD
p_scz_matched_bd <- compare_groups(
  df_matched_scz, "SCZ", "BD"
)

# Results table
pvalue_matched_table <- data.frame(
  Variable = c("Sex", "Age", "LQ"),
  `SCZ vs matched Controls` = format.pval(
    p_scz_matched_controls,
    digits = 3,
    eps = 1e-16
  ),
  `SCZ vs BD` = format.pval(
    p_scz_matched_bd,
    digits = 3,
    eps = 1e-16
  ),
  check.names = FALSE
)

print(pvalue_matched_table)

####
# ==========================================================
# Descriptive statistics: matched BD, SCZ and controls
# ==========================================================

library(dplyr)
library(tidyr)

# Function: median [Q1, Q3]
median_iqr <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(NA_character_)
  
  sprintf(
    "%.2f [%.2f, %.2f]",
    median(x),
    quantile(x, 0.25),
    quantile(x, 0.75)
  )
}

# Function: descriptive statistics
get_demographics <- function(df, group_name) {
  
  tibble(
    Group = group_name,
    N_total = nrow(df),
    
    `Sex = Female (%)` = sprintf(
      "%d (%.1f%%)",
      sum(df$Sex == "MUJER", na.rm = TRUE),
      100 * mean(df$Sex == "MUJER", na.rm = TRUE)
    ),
    
    `Age (median [IQR])` =
      median_iqr(df$Age),
    
    `LQ (median [IQR])` =
      median_iqr(df$Score10items),
    
    `Foot (median [IQR])` =
      median_iqr(df$Item11),
    
    `Eye (median [IQR])` =
      median_iqr(df$Item12),
    
    N_Sex = sum(!is.na(df$Sex)),
    N_Age = sum(!is.na(df$Age)),
    N_LQ = sum(!is.na(df$Score10items)),
    N_Foot = sum(!is.na(df$Item11)),
    N_Eye = sum(!is.na(df$Item12))
  )
}

# ==========================================================
# Define datasets
# ==========================================================

# a) All: BD + SCZ + controls
df_all <- df_matched_scz

# b) BD + controls
df_bd_controls <- df_matched_scz %>%
  filter(BD_patient %in% c("BD", "CONTROL"))

# c) BD only
df_bd <- df_matched_scz %>%
  filter(BD_patient == "BD")

# d) Controls only
df_controls <- df_matched_scz %>%
  filter(BD_patient == "CONTROL")

# ==========================================================
# Calculate descriptive statistics
# ==========================================================

demographics <- bind_rows(
  get_demographics(df_all, "All (BD + SCZ + Controls)"),
  get_demographics(df_bd_controls, "BD + Controls"),
  get_demographics(df_bd, "BD"),
  get_demographics(df_controls, "Controls")
)

# Table for manuscript
demographics_table <- demographics %>%
  select(
    Group,
    `Sex = Female (%)`,
    `Age (median [IQR])`,
    `LQ (median [IQR])`,
    `Foot (median [IQR])`,
    `Eye (median [IQR])`
  ) %>%
  pivot_longer(
    cols = -Group,
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = Value
  )

print(demographics_table, width = Inf)

# Sample sizes for each variable
sample_sizes <- demographics %>%
  select(
    Group, N_total,
    N_Sex, N_Age, N_LQ, N_Foot, N_Eye
  )

print(sample_sizes, width = Inf)

# ==========================================================
# BD vs age-matched controls: updated p-values
# ==========================================================

df_bd_controls <- df_matched_scz %>%
  filter(BD_patient %in% c("BD", "CONTROL")) %>%
  droplevels()

# Check sample sizes
table(df_bd_controls$BD_patient)

# Chi-squared: Sex
p_sex <- chisq.test(
  table(df_bd_controls$Sex, df_bd_controls$BD_patient)
)$p.value

# Wilcoxon: Age
p_age <- wilcox.test(
  Age ~ BD_patient,
  data = df_bd_controls,
  exact = FALSE
)$p.value

# Wilcoxon: LQ
p_lq <- wilcox.test(
  Score10items ~ BD_patient,
  data = df_bd_controls,
  exact = FALSE
)$p.value

# Wilcoxon: Foot
p_foot <- wilcox.test(
  Item11 ~ BD_patient,
  data = df_bd_controls,
  exact = FALSE
)$p.value

# Wilcoxon: Eye
p_eye <- wilcox.test(
  Item12 ~ BD_patient,
  data = df_bd_controls,
  exact = FALSE
)$p.value

# Results table
pvalue_bd_controls <- data.frame(
  Variable = c("Sex", "Age", "LQ", "Foot", "Eye"),
  `BD vs age-matched Controls` = format.pval(
    c(p_sex, p_age, p_lq, p_foot, p_eye),
    digits = 3,
    eps = 1e-16
  ),
  check.names = FALSE
)

print(pvalue_bd_controls)
