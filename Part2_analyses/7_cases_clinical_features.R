library(tibble)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(forcats)
library(mvtnorm)
library(class)
library(VIM)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

### Functions of interest  #####################################
calculate_LQ_EHI <- function(df){
  df_temp <- df %>%
    mutate(
      Right10 = rowSums(
        across(Item1:Item10, ~ case_when(
          . == 1 ~ 2,
          . %in% c(2,3) ~ 1,
          . %in% c(4, 5) ~ 0,
          TRUE ~ NA_real_
        )),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      Left10 = rowSums(
        across(Item1:Item10, ~ case_when(
          . == 5 ~ 2,
          . %in% c(4,3) ~ 1,
          . %in% c(1,2) ~ 0,
          TRUE ~ NA_real_
        )),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      Score10items = round(((Right10 - Left10) / (Right10 + Left10) *100),0)
    ) %>%select(ID, Score10items)
  return(df_temp)
}
association_matrix <- function(df) {
  vars <- names(df)
  n <- length(vars)
  mat <- matrix(NA, n, n, dimnames = list(vars, vars))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      x <- df[[i]]
      y <- df[[j]]
      
      # Ambos numéricos → correlación de Pearson
      if (is.numeric(x) && is.numeric(y)) {
        mat[i,j] <- cor(x, y, use = "pairwise.complete.obs", method = "pearson")
        
        # Uno numérico y otro factor → R² de un modelo lineal
      } else if (is.numeric(x) && is.factor(y)) {
        model <- lm(x ~ y)
        mat[i,j] <- summary(model)$r.squared
      } else if (is.factor(x) && is.numeric(y)) {
        model <- lm(y ~ x)
        mat[i,j] <- summary(model)$r.squared
        
        # Ambos factores → Cramér’s V
      } else if (is.factor(x) && is.factor(y)) {
        x <- droplevels(x)
        y <- droplevels(y)
        tbl <- table(x, y)
        mat[i,j] <- vcd::assocstats(tbl)$cramer
      }
    }
  }
  return(mat)
}
count_by_group <- function(df, cutoff, var = NULL) {
  
  allowed_vars <- c("Sex", "BD_patient", "Diagnostic")
  # función interna para crear grupos
  make_groups <- function(data) {
    data %>%
      mutate(
        group = case_when(
          Score10items < -cutoff ~ "left",
          Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
          Score10items > cutoff ~ "right"
        )
      )
  }
  
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
make_groups_summary <- function(data, cutoff, lateral = FALSE) {
  
  library(dplyr)
  library(tidyr)
  
  # crear grupos base
  data_grouped <- data %>%
    mutate(
      group = case_when(
        Score10items < -cutoff ~ "left",
        Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
        Score10items > cutoff ~ "right"
      )
    )
  
  # modo lateral
  if (lateral) {
    data_grouped <- data_grouped %>%
      mutate(
        group = ifelse(group == "mixed", "mixed", "lateralized")
      )
    
    group_levels <- c("lateralized", "mixed")
    
  } else {
    group_levels <- c("left", "mixed", "right")
  }
  
  vars <- c("WHODAS_total32", "IcgGe", "INES", "ScoreEeag")
  
  # n por grupo
  n_groups <- data_grouped %>%
    count(group)
  
  n_groups <- setNames(n_groups$n, n_groups$group)
  
  # stats
  get_stats <- function(df, var) {
    m <- mean(df[[var]], na.rm = TRUE)
    s <- sd(df[[var]], na.rm = TRUE)
    sprintf("%.2f (%.2f)", m, s)
  }
  
  # p-values
  get_p <- function(df, var) {
    
    df <- df %>% filter(!is.na(.data[[var]]), !is.na(group))
    
    if (lateral) {
      
      wilcox.test(
        df[[var]][df$group == "lateralized"],
        df[[var]][df$group == "mixed"]
      )$p.value
      
    } else {
      
      kruskal.test(df[[var]] ~ df$group)$p.value
    }
  }
  
  summary_table <- lapply(vars, function(v) {
    
    stats_df <- data_grouped %>%
      group_by(group) %>%
      summarise(
        value = get_stats(cur_data(), v),
        .groups = "drop"
      ) %>%
      pivot_wider(names_from = group, values_from = value)
    
    # evitar columnas fantasma
    stats_df <- stats_df %>%
      select(any_of(c("variable", "left", "mixed", "right", "lateralized")))
    
    p_val <- get_p(data_grouped, v)
    
    stats_df %>%
      mutate(
        variable = v,
        p_value = sprintf("p=%.3g", p_val)
      ) %>%
      select(variable, any_of(group_levels), p_value)
  }) %>%
    bind_rows()
  
  # n seguro
  n_groups <- n_groups[group_levels]
  
  # nombres finales
  new_names <- c(
    "variable",
    paste0(group_levels, " (n=", n_groups, ")"),
    if (lateral) "Wilcoxon p" else "Kruskal-Wallis p"
  )
  
  colnames(summary_table) <- new_names
  
  return(summary_table)
}
make_cutoff_boxplot_facet <- function(data, cutoffs, variable) {
  
  library(dplyr)
  library(ggplot2)
  
  df_all <- lapply(cutoffs, function(cutoff) {
    
    data %>%
      mutate(
        group = case_when(
          Score10items < -cutoff ~ "left",
          Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
          Score10items > cutoff ~ "right"
        ),
        cutoff = factor(cutoff)
      )
  }) %>%
    bind_rows()
  
  df_all <- df_all %>%
    filter(!is.na(.data[[variable]]), !is.na(group))
  
  # calcular n y mediana por grupo/cutoff
  summary_df <- df_all %>%
    group_by(cutoff, group) %>%
    summarise(
      n = n(),
      med = median(.data[[variable]], na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(df_all, aes(x = group, y = .data[[variable]], fill = group)) +
    
    # puntos individuales
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    
    # boxplot
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    
    # línea de mediana
    stat_summary(fun = median, geom = "crossbar", width = 0.6, color = "black") +
    
    # n dentro del plot
    geom_text(
      data = summary_df,
      aes(
        x = group,
        y = max(df_all[[variable]], na.rm = TRUE),
        label = paste0("n=", n)
      ),
      vjust = -0.5,
      size = 3,
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~cutoff) +
    
    labs(
      title = variable,
      x = "Group",
      y = variable
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "none"
    )
  
  return(p)
}
################################################################################


#Load and clean clinical variables
scales <- as.data.frame(readxl::read_xlsx("Escalas_clinicas_MadManic.xlsx", sheet = 1, 
                                          na=c("#N/A", "#N/D")))
scales[[1]] <- gsub("IE", "IEV", scales[[1]])

meta_cases <- as.data.frame(readxl::read_xlsx("../input_prep/cases_EHI_070526.xlsx", sheet = 1, 
                                         col_types="text", na=c("#N/A", "#N/D")))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>% mutate(Diagnostic = factor(Diagnostic))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  ))%>%
  left_join(scales, by = c("ID" = "ID_complete"))%>%
  mutate(across(starts_with("NUMBER"), 
                ~ { x <- as.numeric(.)
                case_when(
                  x == 0 ~ 0,
                  x %in% c(1,2) ~ 1,
                  x >=3 & x <=5 ~ 2,
                  x >=6 & x <=10 ~ 3,
                  x >=11 ~4,
                  TRUE ~ NA_real_
                ) }))

meta_clean <- meta_cases%>%
  filter(rowSums(is.na(.)) < 7) %>% mutate(Age = as.numeric(Age))%>% select(-BD_patient)%>%
  mutate(across(where(is.factor), ~ {  # set all columns that have a level "No" as their reference
    if ("NO" %in% levels(.x)) relevel(.x, ref = "NO") else .x
  }))

###Load EHI and compute LQ

ehi_cases <-  as.data.frame(readxl::read_xlsx("../input_prep/cases_EHI_070526.xlsx", sheet = 2, 
                                              col_types="text", na="#N/A"))[,1:13]%>%
  mutate(across(-1, as.numeric))

cols <- c("ID", "Item1", "Item2", "Item3", "Item4", "Item5", "Item6", "Item7",
          "Item8", "Item9", "Item10", "Item11", "Item12")
colnames(ehi_cases) <- cols
LQ_cases <- calculate_LQ_EHI(ehi_cases)
all_df <- meta_clean %>% left_join(LQ_cases, by="ID") 

cor <- association_matrix(all_df[,-1])
pheatmap(cor, 
         display_numbers = TRUE,       
         number_format = "%.2f",       
         number_color = "black",       
         fontsize_number = 10,         
         color = colorRampPalette(c("slateblue4", "lightyellow", "darkorange"))(50), 
         main = "Clinical features and LQ",
         border_color = "white",       
         treeheight_row = 30,          
         treeheight_col = 30)

all_df_selection <- all_df%>%
  select(ID, Score10items, Sex, Age, Diagnostic, INES,
         IcgGe, WHODAS_total32, ScoreEeag)
cor_sel <- association_matrix(all_df_selection[,-1])
pheatmap(cor_sel, 
         display_numbers = TRUE,       
         number_format = "%.2f",       
         number_color = "black",       
         fontsize_number = 10,         
         color = colorRampPalette(c("slateblue4", "lightyellow", "darkorange"))(50), 
         main = "Clinical features and LQ",
         border_color = "white",       
         treeheight_row = 30,          
         treeheight_col = 30)

make_groups_summary(all_df_selection, 0)
whodas_plots <- make_cutoff_boxplot_facet(all_df_selection, c(0, 40, 60, 80, 90), "WHODAS_total32")
whodas_plots
cgi_plots <- make_cutoff_boxplot_facet(all_df_selection, c(0, 40, 60, 80, 90), "IcgGe")
cgi_plots
gaf_plots <- make_cutoff_boxplot_facet(all_df_selection, c(0, 40, 60, 80, 90), "ScoreEeag")
gaf_plots
ines_plots<- make_cutoff_boxplot_facet(all_df_selection, c(0, 40, 60, 80, 90), "INES")
ines_plots




