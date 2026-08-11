workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tibble)

prs_results <- tribble(
  ~Target, ~Base, ~Sample, ~beta, ~SE, ~PRS_r2, ~Full_r2, ~P_value,
  "NRH_0",      "PRS_handedness", "All",           0.228,  0.115, 0.009, 0.018, 0.0465,
  "LQ",         "PRS_handedness", "All",          -2.854,  1.340, 0.004, 0.014, 0.0333,
  "BD_patient", "PRS_handedness", "All",          -0.015,  0.081, 0.000, 0.430, 0.852,
  "WHODAS",     "PRS_handedness", "Only BD", 0.977,  1.074, 0.003, 0.008, 0.363,
  "INES",       "PRS_handedness", "Only BD",       0.006,  0.005, 0.003, 0.019, 0.291,
  "GAF",        "PRS_handedness", "Only BD",      -1.152,  0.892, 0.006, 0.018, 0.1974,
  "CGI",        "PRS_handedness", "Only BD",       0.166,  0.072, 0.017, 0.025, 0.022,
  "BD_patient", "PRS_BD",         "All",            0.677,  0.091, 0.052, 0.482, 1.17e-13,
  "WHODAS",     "PRS_BD",         "Only BD",        1.225,  1.191, 0.004, 0.009, 0.305,
  "INES",       "PRS_BD",         "Only BD",        0.006,  0.006, 0.003, 0.018, 0.328,
  "GAF",        "PRS_BD",         "Only BD",       -0.625,  0.981, 0.001, 0.013, 0.5247,
  "CGI",        "PRS_BD",         "Only BD",        0.031,  0.080, 0.001, 0.008, 0.695,
  "NRH_0",      "PRS_BD",         "All",           -0.010,  0.117, 0.000, 0.009, 0.920,
  "NRH_0",      "PRS_BD",         "Only BD",       -0.529,  0.231, 0.037, 0.053, 0.022,
  "LQ",         "PRS_BD",         "All",            0.177,  1.372, 0.000, 0.010, 0.898,
  "LQ",         "PRS_BD",         "Only BD",        5.083,  2.612, 0.011, 0.017, 0.0525
)

prs_results <- prs_results %>%
  mutate(
    lower = beta - 1.96 * SE,
    upper = beta + 1.96 * SE,
    significant = P_value < 0.05
  )

ggplot(
  prs_results,
  aes(
    x = beta,
    y = reorder(Target, beta),
    xmin = lower,
    xmax = upper
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_errorbarh(
    height = 0.2
  ) +
  geom_point(
    aes(shape = Sample),
    size = 3
  ) +
  facet_wrap(
    ~ Base,
    scales = "free_x"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = expression(beta~"(95% CI)"),
    y = NULL,
    shape = "Sample"
  )

prs_hand <- prs_results %>%
  filter(Base == "PRS_handedness") %>%
  mutate(
    lower = beta - 1.96 * SE,
    upper = beta + 1.96 * SE,
    significant = P_value < 0.05,
    Target = factor(
      Target,
      levels = rev(c(
        "NRH_0", "LQ", "BD_patient",
        "WHODAS", "INES", "GAF", "CGI"
      ))
    )
  )
ggplot(prs_hand,
       aes(
         x = beta,
         y = Target,
         xmin = lower,
         xmax = upper
       )) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  geom_errorbarh(
    height = 0.15,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(color = significant, shape = Sample),
    size = 3.5
  ) +
  
  scale_color_manual(
    values = c(
      "FALSE" = "grey40",
      "TRUE" = "#E76F51"
    ),
    labels = c(
      "FALSE" = "P ≥ 0.05",
      "TRUE" = "P < 0.05"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "All" = 16,       # círculo sólido
      "Only BD" = 17    # triángulo sólido
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  labs(
    x = expression(beta~"(95% CI)"),
    y = NULL,
    color = NULL,
    shape = "Sample"
  ) +
  
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14)
  )



prs_bd <- prs_results %>%
  filter(Base == "PRS_BD") %>%
  mutate(
    Sample = recode(
      Sample,
      "Only BD (265)" = "Only BD"
    ),
    lower = beta - 1.96 * SE,
    upper = beta + 1.96 * SE,
    significant = P_value < 0.05,
    Target = factor(
      Target,
      levels = rev(c(
        "NRH_0",
        "LQ",
        "BD_patient",
        "WHODAS",
        "INES",
        "GAF",
        "CGI"
      ))
    )
  )
pd <- position_dodge(width = 0.5)

ggplot(prs_bd,
       aes(
         x = beta,
         y = Target,
         xmin = lower,
         xmax = upper,
         group = Sample
       )) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  geom_errorbarh(
    aes(color = significant),
    height = 0.15,
    linewidth = 0.8,
    position = pd
  ) +
  
  geom_point(
    aes(
      color = significant,
      shape = Sample
    ),
    size = 3.5,
    position = pd
  ) +
  
  scale_color_manual(
    values = c(
      "FALSE" = "grey40",
      "TRUE" = "#E76F51"
    ),
    labels = c(
      "FALSE" = "P ≥ 0.05",
      "TRUE" = "P < 0.05"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "All" = 16,
      "Only BD" = 17
    )
  ) +
  
  theme_minimal(base_size = 15) +
  
  labs(
    x = expression(beta~"(95% CI)"),
    y = NULL,
    color = NULL,
    shape = "Sample"
  ) +
  
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
