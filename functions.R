##Funciones para analisis del proyecto de lateralidad

count_by_group <- function(df, cutoff, var = NULL) {
  
  allowed_vars <- c("Sex", "BD_patient", "Diagnostic", "Group_age")
  group_order <- c("left", "mixed", "right", "lateralized", "NRH")
  
  # Crear grupos de lateralidad
  make_groups <- function(data) {
    data %>%
      mutate(
        group = case_when(
          Score10items < -cutoff ~ "left",
          Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
          Score10items > cutoff ~ "right",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(group))
  }
  
  if (!is.null(var) && var != "None" &&
      !(var %in% allowed_vars)) {
    stop("var debe ser Sex, BD_patient, Diagnostic o NULL/None")
  }
  
  df2 <- make_groups(df)
  
  # Función auxiliar para calcular los recuentos
  count_groups <- function(data) {
    
    total <- nrow(data)
    
    counts <- data %>%
      count(group, name = "n")
    
    # Añadir lateralized (left + right) y NRH (left + mixed)
    extra <- tibble(
      group = c("lateralized", "NRH"),
      n = c(
        sum(data$group %in% c("left", "right")),
        sum(data$group %in% c("left", "mixed"))
      )
    )
    
    bind_rows(counts, extra) %>%
      mutate(
        perc = if (total > 0) 100 * n / total else NA_real_,
        label = sprintf("%d (%.2f%%)", n, perc),
        group = factor(group, levels = group_order)
      ) %>%
      arrange(group) %>%
      mutate(group = as.character(group)) %>%
      select(group, label)
  }
  
  # Sin estratificación
  if (is.null(var) || var == "None") {
    return(count_groups(df2))
  }
  
  # Estratificado por variable
  df2 %>%
    group_by(.data[[var]]) %>%
    group_modify(~ count_groups(.x)) %>%
    ungroup() %>%
    select(group, all_of(var), label) %>%
    arrange(.data[[var]], match(group, group_order))
}

prop_test_by_group <- function(df, cutoff, var = "Sex") {
  allowed_vars <- c("Sex", "Group_age", "Diagnostic", "BD_patient")
  if (!(var %in% allowed_vars)) {
    stop("var debe ser 'Sex' o 'Group_age'")
  }
  
  # Clasificar la lateralidad según el umbral
  df2 <- df %>%
    mutate(
      group = case_when(
        Score10items < -cutoff ~ "left",
        Score10items >= -cutoff & Score10items <= cutoff ~ "mixed",
        Score10items > cutoff ~ "right",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group), !is.na(.data[[var]]))
  
  # Comprobar que existen exactamente dos grupos
  levels_var <- unique(as.character(df2[[var]]))
  
  if (length(levels_var) != 2) {
    stop("La variable de comparación debe tener exactamente dos grupos")
  }
  
  # Número total de individuos en cada grupo
  totals <- df2 %>%
    count(.data[[var]], name = "total")
  
  # Calcular categorías binarias
  df2 <- df2 %>%
    mutate(
      left = group == "left",
      mixed = group == "mixed",
      NRH = group %in% c("left", "mixed"),
      lateralized = group %in% c("left", "right")
    )
  
  # Recuento de individuos en cada categoría
  counts <- df2 %>%
    group_by(.data[[var]]) %>%
    summarise(
      across(
        c(left, mixed, NRH, lateralized),
        ~ sum(.x)
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(left, mixed, NRH, lateralized),
      names_to = "group",
      values_to = "n"
    ) %>%
    left_join(totals, by = var) %>%
    mutate(
      perc = 100 * n / total,
      label = sprintf("%d (%.2f%%)", n, perc)
    )
  
  # Test de proporciones para cada categoría
  results <- counts %>%
    group_by(group) %>%
    summarise(
      p_value = prop.test(
        x = n,
        n = total,
        alternative = "two.sided",
        correct = TRUE
      )$p.value,
      .groups = "drop"
    )
  
  # Tabla de resultados
  counts %>%
    select(all_of(var), group, label) %>%
    pivot_wider(
      names_from = all_of(var),
      values_from = label
    ) %>%
    left_join(results, by = "group") %>%
    mutate(
      group = factor(
        group,
        levels = c("left", "mixed", "NRH", "lateralized")
      )
    ) %>%
    arrange(group) %>%
    mutate(group = as.character(group))
}

create_laterality_vars <- function(df, cutoffs = c(0, 40, 60, 80, 90)) {
  
  library(dplyr)
  
  # Comprobar que existe Score10items
  if (!"Score10items" %in% names(df)) {
    stop("El dataset debe contener la columna 'Score10items'")
  }
  
  for (cutoff in cutoffs) {
    
    # Non-right-handedness (Left + Mixed)
    df[[paste0("NRH_", cutoff)]] <- factor(
      if_else( 
        df$Score10items <= cutoff,
        "NRH", "RH"
      ),
      levels = c("RH", "NRH")
    )
    
    # Lateralization (Left + Right)
    df[[paste0("LAT_", cutoff)]] <- factor(
      if_else(
        df$Score10items >= -cutoff &
          df$Score10items <= cutoff,
        "Non-lateral", "Lateral"
      ),
      levels = c("Lateral", "Non-lateral")
    )
  }
  
  return(df)
}
