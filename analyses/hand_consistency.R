library(dplyr)
library(ggplot2)
library(ggrepel)

workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))
rm(list = ls())


input <- as.data.frame(readxl::read_xlsx("../input_prep/all_samples_all_data_LQ_PLUS_items.xlsx", sheet = 1, 
                                         col_types="text", na="#N/A"))%>%
  mutate(across(
    where(~ n_distinct(., na.rm = TRUE) == 2),
    as.factor
  )) %>%
  mutate(across(-c(ID, Diagnostic),~ if (!is.factor(.)) as.numeric(.) else .)) %>% 
  mutate(Diagnostic = factor(Diagnostic))

lat_df <- input %>%
  rowwise() %>%
  mutate(consistency_sd = 
           (sd(c_across(7:16))),
         hand_cons_inverse =
           1/(consistency_sd + 1),
         counts_of_three= sum(c_across(7:16) == 3, na.rm = TRUE),
         ) %>%
  ungroup()
n_above <- sum(lat_df$hand_cons_inverse > 0.9, na.rm = TRUE)
n_below <- sum(lat_df$hand_cons_inverse < 0.9, na.rm = TRUE)

ggplot(lat_df, aes(x = Score10items,
               y = hand_cons_inverse,
               color = BD_patient)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ID),
                  size = 4,
                  max.overlaps = 30,
                  box.padding = 0.40) +
  theme_minimal()

writexl::write_xlsx(lat_df, path= "account_for_consistent_hand_use.xlsx")
poi <- lat_df%>%
  filter(counts_of_three >= 5)

writexl::write_xlsx(poi, path= "../People_Of_Interest_a_lot_of_threes.xlsx")


