# Title: Table 1 for CC-Oxygen
# Author: Ed Palmer

# Set up
library(tidyverse)
library(tableone)
library(knitr)
# Run 01-preliminaries.R first

setwd("####/oxygen/")
load("####/cohorts_final.RData")

table_one_data <- bind_rows(
  cohort_1 %>% mutate(subset = "day 1 subset"),
  cohort_2 %>% mutate(subset = "day 3 subset"),
  cohort_3 %>% mutate(subset = "day 5 subset"),
  cohort_4 %>% mutate(subset = "day 7 subset")) %>%
  mutate(exposed = if_else(hyperoxia_13 == 1, "yes", "no"))

table_one_data <- table_one_data %>%
  mutate(ethnicity = case_when(
    ethnicity == "white British" ~ "white British",
    ethnicity == "white other" ~ "white other",
    ethnicity == "black/black British African" ~ "black/black British African",
    ethnicity == "Asian/Asian British other" ~ "Asian/Asian British other",
    ethnicity == "black/black British Caribbean" ~ "black/black British Caribbean",
    ethnicity == "Asian/Asian British Indian" ~ "Asian/Asian British Indian",
    is.na(ethnicity) ~ "other or not stated",
    !is.null(ethnicity) ~ "other or not stated"
  ))

table_one_data <- table_one_data %>%
  select(tw_hyperoxia_13,
         cumulative_hyperoxia_13,
         exposed,
         age,
         weight,
         sex,
         apache_score,
         system,
         location_in,
         prior_dependency,
         is_medical,
         surgical_classification,
         ethnicity,
         spell_los,
         unit_mortality,
         subset) %>%
  rename(
    `time weighted mean hyperoxaemia` = tw_hyperoxia_13,
    `cumulative hyperoxaemia` = cumulative_hyperoxia_13,
    `exposed to hyperoxaemia` = exposed,
    `age (years)` = age,
    `APACHE II score` = apache_score,
    `unit mortality` = unit_mortality,
    `primary organ system` = system,
    `prior location` = location_in,
    `prior dependency` = prior_dependency,
    `medical` = is_medical,
    `surgical classification` = surgical_classification,
    `length of stay (days)` = spell_los,
    `weight (kg)` = weight
    )


table_1_full <- CreateTableOne(
  vars = names(table_one_data)[-length(names(table_one_data))],
  data = table_one_data,
  includeNA = TRUE,
  strata = "subset",
  test = FALSE)

table_1_full <- print(table_1_full,
                      nonnormal = c("apache II score",
                                    "length of stay (days)",
                                    "age (years)",
                                    "time weighted mean hyperoxaemia",
                                    "cumulative hyperoxaemia"),
                      printToggle = FALSE, showAllLevels = FALSE)

table_1_out <- as.tibble(table_1_full) %>%
  mutate(characteristic = rownames(table_1_full)) %>%
  select(characteristic, `day 1 subset`, `day 3 subset`, `day 5 subset`, `day 7 subset`)

any_exposure <- bind_rows(
  cohort_1 %>% mutate(subset = "day 1 subset"),
  cohort_2 %>% mutate(subset = "day 3 subset"),
  cohort_3 %>% mutate(subset = "day 5 subset"),
  cohort_4 %>% mutate(subset = "day 7 subset")) %>%
  group_by(subset, hyperoxia_13 == 1) %>%
  tally() %>%
  rename(exposed = `hyperoxia_13 == 1`) %>%
  group_by(subset) %>%
  mutate(perc = n/sum(n)*100)

save(table_1_out, any_exposure, file = "./manuscript/table_one.RData")

## Move to 03-logisitc-regression
