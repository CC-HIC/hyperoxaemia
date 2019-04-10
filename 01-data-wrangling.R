## Data Wrangling Code
## Primary Author: Dr. Ed Palmer (MRC PhD Student (Clinical Data Science) / Anaesthetic Intensivist SpR)
## Please contact me on github @Doc_Ed if you find any mistakes in this code

setwd("####/oxygen")

library(tidyverse)
library(inspectEHR) # note extract is masked by tidyverse
library(lubridate)
library(cowplot)
source("####/00-cpp-functions.R")

# Establish a DB connection and retrive tables
ctn <- connect(username = "####",
               password = "####",
               database = "####",
               system = "####")
tbls <- retrieve_tables(ctn)

# Grab non-longitudinal data ----------------
demo_codes <- tribble(
  ~hic_codes, ~short_name,
  "NIHR_HIC_ICU_0399", "primary_admission_reason",
  "NIHR_HIC_ICU_0088", "secondary_admission_reason",
  "NIHR_HIC_ICU_0074", "pmhx",
  "NIHR_HIC_ICU_0068", "location_in",
  "NIHR_HIC_ICU_0398", "admission_type",
  "NIHR_HIC_ICU_0004", "treatment_code",
  "NIHR_HIC_ICU_0011", "pre_surgical_prep",
  "NIHR_HIC_ICU_0027", "surgical_classification",
  "NIHR_HIC_ICU_0409", "apache_score",
  "NIHR_HIC_ICU_0055", "prior_dependency",
  "NIHR_HIC_ICU_0058", "ethnicity",
  "NIHR_HIC_ICU_0033", "dob",
  "NIHR_HIC_ICU_0019", "weight",
  "NIHR_HIC_ICU_0021", "cpr",
  "NIHR_HIC_ICU_0093", "sex",
  "NIHR_HIC_ICU_0094", "organ_donor",
  "NIHR_HIC_ICU_0103", "withdrawal",
  "NIHR_HIC_ICU_0097", "unit_mortality",
  "NIHR_HIC_ICU_0095", "hospital_mortality",
  "NIHR_HIC_ICU_0100", "transferring_unit_number",
  "NIHR_HIC_ICU_0101", "transferring_unit_id",
  "NIHR_HIC_ICU_0013", "ccu_transfer_group",
  "NIHR_HIC_ICU_0930", "ultimate_unit_mortality",
  "NIHR_HIC_ICU_0098", "utlimate_hospital_mortality"
)

dtb <- extract_demographics(metadata = tbls[["variables"]],
                            events = tbls[["events"]],
                            code_names = demo_codes$hic_codes,
                            rename = demo_codes$short_name)

# Grab longitudinal data ---------------------
# Set cut offs for LOCF procedure here
long_codes <- tribble(
  ~hic_codes, ~short_name, ~locf_limit, ~lower_bound_off,
  "NIHR_HIC_ICU_0411", "start_dttm", NA, NA,
  "NIHR_HIC_ICU_0149", "peak_ap", 6, NA,
  "NIHR_HIC_ICU_0152", "mean_ap", 6, NA,
  "NIHR_HIC_ICU_0151", "peep", 6, NA,
  "NIHR_HIC_ICU_0550", "tidal_vol", 6, NA,
  "NIHR_HIC_ICU_0129", "spo2", 6, NA,
  "NIHR_HIC_ICU_0549", "rr_spont", 6, NA,
  "NIHR_HIC_ICU_0144", "ventilation", 12, NA,
  "NIHR_HIC_ICU_0130", "sxo2", 6, NA,
  "NIHR_HIC_ICU_0132", "pxo2", 6, NA,
  "NIHR_HIC_ICU_0150", "fio2", 6, NA,
  "NIHR_HIC_ICU_0126", "airway", 12, NA,
  "NIHR_HIC_ICU_0913", "pf_ratio", 12, NA)

ltb <- extract_timevarying(
  events = tbls[["events"]],
  metadata = collect(tbls[["variables"]]),
  code_names = long_codes$hic_codes,
  rename = long_codes$short_name,
  chunk_size = 10000,
  cadance = 0.5)

save(dtb, ltb, demo_codes, long_codes, file = "####/end_data_extraction.RData")
load("####/cam_additional.RData")

## ltb_imp = imported data from cambridge
dtb %>%
  filter(episode_id %in% unique(ltb_imp$episode_id))

# General set up
core <- make_core(ctn)
reference <- make_reference(ctn)
episode_length <- epi_length(core, reference, tbls[["events"]])
episodes <- collect(tbls[["episodes"]])
spells <- identify_spells(episode_length, episodes, minutes = 360)

# remove extreme negative times
# this likely represents events that occur outside the ICU epsiode but
# have become attached inadvertently the to current episode

# remove back by 24 hours (48 steps)
ltb <- ltb %>%
  filter(time > -48)

## Lets weave in the cambridge data here
## Are there ANY Cambridhe patients with PaO2 data? We don't want to overwrite

ltb %>%
  filter(episode_id %in% unique(ltb_imp$episode_id)) %>%
  filter(!is.na(pxo2)) %>%
  group_by(episode_id) %>%
  tally()

ltb <- ltb %>%
  filter(episode_id != "*****") #(removing single overlapping patient)

# Set a row ID
ltb <- ltb %>%
  arrange(episode_id, time) %>%
  mutate(ltb_id = seq(1, n()))

cam <- ltb %>%
  filter(episode_id %in% unique(ltb_imp$episode_id))

others <- ltb %>%
  filter(!(episode_id %in% unique(ltb_imp$episode_id)))

cam <- cam %>%
  select(-pxo2, -NIHR_HIC_ICU_0132.meta.1) %>%
  full_join(ltb_imp, by = c("episode_id", "time"))

ltb_backup <- ltb

ltb <- bind_rows(cam, others) %>%
  arrange(episode_id, time) %>%
  group_by(episode_id)

# Expand this dense tb (add NAs where appropraite)
ltb <- ltb %>%
  select(episode_id, time) %>%
  split(., .$episode_id) %>% 
  imap(function(base_table, epi_id) {
    tibble(episode_id = as.integer(epi_id),
           time = seq(min(base_table$time, 0L),
                      max(base_table$time, 0L),
                      by = 0.5))
  }) %>%
  bind_rows() %>%
  left_join(ltb, by = c("episode_id" = "episode_id", "time" = "time"))

## Resequence Primary Key
ltb <- ltb %>%
  arrange(episode_id, time) %>%
  mutate(ltb_id = seq(1, n()))

# Check that time coverage meets expectations
# hist(ltb$time)

## Make a site level key
site_key <- tbls[["episodes"]] %>%
  left_join(tbls[["provenance"]], by = c("provenance" = "file_id")) %>%
  select(episode_id, site) %>%
  collect()

## Sanity Check
ltb %>%
  left_join(site_key, by = "episode_id") %>%
  group_by(site) %>%
  summarise(sum(!is.na(pxo2)))

## episode validation was previously completed and loaded into the DB
validated_episodes <- tbls[["episode_validation"]] %>%
  filter(validity == 0) %>%
  collect() %>%
  select(episode_id) %>%
  pull()

# Remove any of the episodes we aren't going to be using
ltb <- ltb %>%
  filter(episode_id %in% validated_episodes)

## DATA MODIFICATIONS ====
## PaO2 ====

# Determine the value at which there is only a small risk of making
# A mistake on re-classifying venous o2s
venous <- na.omit(ltb$pxo2[ltb$NIHR_HIC_ICU_0132.meta.1 == 2])
venous_cutoff <- quantile(venous, probs = c(0.98))
rm(venous)

## There is "some" semantic mislabelling here
## Using the semantic labelling, we can define a cut off point
## where we are only calling PvO2s PaO2s 2% of the time in the unlabelled data
## Since we are only looking at values above 13.3 anyway,
## I doubt this would
## be enough to introduce bias.

ltb %>%
  filter(!is.na(pxo2)) %>%
  filter(pxo2 < 100) %>%
  mutate(NIHR_HIC_ICU_0132.meta.1 = factor(NIHR_HIC_ICU_0132.meta.1,
                                           levels = c(1, 2, 3),
                                           labels = c("arterial", "venous", "unknown"))) %>%
  ggplot(aes(x = pxo2, fill = NIHR_HIC_ICU_0132.meta.1)) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = venous_cutoff))

# Reclassify ABG based on our strict criteria
ltb <- ltb %>%
  mutate(NIHR_HIC_ICU_0132.meta.1 = factor(NIHR_HIC_ICU_0132.meta.1,
                                           levels = c(1, 2, 3),
                                           labels = c("arterial", "venous", "unknown"))) %>%
  mutate(pao2 = if_else(
    ((NIHR_HIC_ICU_0132.meta.1 == "arterial") |
       ((is.na(NIHR_HIC_ICU_0132.meta.1)) & (pxo2 > venous_cutoff)) |
       ((NIHR_HIC_ICU_0132.meta.1 == "unknown") & (pxo2 > venous_cutoff))),
    pxo2, as.numeric(NA))) %>%
  select(-pxo2, -NIHR_HIC_ICU_0132.meta.1) %>%
  mutate(pao2 = if_else(
    pao2 > 100, as.numeric(NA), pao2
  ))

rm(venous_cutoff)

pao2_global_final <- ltb %>%
  filter(!is.na(pao2)) %>%
  ggplot(aes(x = pao2)) +
  geom_density() +
  geom_vline(xintercept = 13.3, linetype = 2) +
  xlab("PaO2 (kPa)") +
  ylab("Sample density")

pao2_global_final %>%
  ggsave(plot = ., filename = "####/pao2_global_final.png",
         width = 10, height = 8, dpi = 300, units = "in")

## FIO2 ====
ggplot(ltb, aes(fio2)) + geom_density()

# Fix incorrect units
# Simple transformation

ltb <- ltb %>%
  mutate(
    fio2 = if_else(
      fio2 >= 21 | fio2 <= 100, fio2/100, fio2)) %>%
  mutate(
    fio2 = if_else(
      fio2 < 0.21 | fio2 > 1, as.numeric(NA), fio2))

fio2_global_final <- ltb %>%
  filter(!is.na(fio2)) %>%
  ggplot(aes(x = fio2)) +
  geom_histogram(bins = 20) +
  xlim(0.2, 1) +
  xlab("FiO2") +
  ylab("Number of Samples")

fio2_global_final %>%
  ggsave(plot = ., filename = "####/fio2_global_final.png",
         width = 10, height = 8, dpi = 300, units = "in")

## PF Ratio

ltb %>%
  left_join(site_key, by = "episode_id") %>%
  filter(!is.na(pf_ratio)) %>%
  ggplot(aes(x = pf_ratio, fill = site)) + geom_density(alpha = 0.5)

## GSTT are given in mmHg, so lets transform them

ltb <- ltb %>%
  left_join(site_key, by = "episode_id") %>%
  mutate(pf_ratio = if_else(
    site == "GSTT" & !is.na(pf_ratio), pf_ratio/7.6, pf_ratio
  ))

ltb %>%
  filter(!is.na(pf_ratio)) %>%
  filter(pf_ratio <= 100) %>%
  ggplot(aes(x = pf_ratio, fill = site)) + geom_density(alpha = 0.5)

## Let's null out impossible values

ltb <- ltb %>%
  mutate(pf_ratio = if_else(pf_ratio > 100 | pf_ratio == 0, as.numeric(NA), pf_ratio))

ltb %>%
  filter(!is.na(pf_ratio)) %>%
  filter(pf_ratio <= 100) %>%
  ggplot(aes(x = pf_ratio, fill = site)) + geom_density(alpha = 0.5)

## Those all look correct now.
# Use only validated spells. Note, 0000000000 validates, though is
# evidently erroneous, hence the exception here
## Cambridge APACHE II are incorrect. This is a confounder that is central to the
## analysis so we will drop them from the data.

valid_spells <- spells %>%
  mutate(is_val_nhs = validate_nhs(nhs_number)) %>%
  mutate(is_val_nhs = if_else(nhs_number == "0000000000", FALSE, is_val_nhs)) %>%
  filter(is_val_nhs == TRUE) %>%
  filter(site != "RGT") %>%
  select(-is_val_nhs)

## Start a tracker for consort diagram
## Consort tracker will run from valid_spells as we pair cases away

consort_info <- tribble(
  ~stage, ~n,
  "overall episodes", nrow(episode_length %>% distinct(episode_id)),
  "valid episodes", nrow(episode_length %>% filter(validity == 0, site != "RGT") %>% distinct(episode_id)),
  "refactored as spells", nrow(valid_spells %>% distinct(spell_id))
)

## Add a sequence counter and retain only the FIRST admission
valid_spells_f <- valid_spells %>%
  arrange(nhs_number, epi_start_dttm) %>%
  group_by(nhs_number) %>%
  mutate(first_spell = min(spell_id)) %>% # we can do this as identify_spells increments based on start datetime
  mutate(first_spell = if_else(first_spell == spell_id, 1L, 0L)) %>%
  filter(first_spell == 1L) %>%
  group_by(spell_id) %>%
  arrange(spell_id, epi_start_dttm) %>%
  mutate(episode_position = as.integer(1:n())) %>%
  select(-first_spell)

consort_info <- bind_rows(consort_info, tribble(
  ~stage, ~n,
  "index ICU admission only", nrow(valid_spells_f %>% distinct(spell_id))
))

valid_spells_f <- valid_spells_f %>%
  group_by(spell_id) %>%
  mutate(spell_los = sum(los))

length(unique(valid_spells_f$episode_id)) == nrow(valid_spells_f)

# remove data from ltb that comes from invalid episodes
ltb <- inner_join(ltb, valid_spells_f, by = "episode_id")

# negative information for episode position > 1 must be dropped
# This drop a very small number of rows, as expected.

ltb <- ltb %>% 
  filter(episode_position == 1 | time >= 0)

# Remove information from an episode than extends beyond its length of stay
ltb <- ltb %>%
  mutate(los_hours = los*24) %>%
  filter(time <= los_hours)

# Save prior to longitudinal imputation ====
save.image(file = "####/pre_imputation.RData")

ltb <- ltb %>%
  select(-site.x, -validity) %>%
  rename(site = site.y)

ltb %>%
  group_by(site) %>%
  summarise(sum(!is.na(pao2)))

library(zoo)

## prevents dplyr throwing an error if there aren't enough samples to
## impute. 24 as we want to impute by 12 hours (half a shift) and the
## cadance is half hourly

custom_impute <- function(x) {
  if (sum(!is.na(x)) < 2) {
    return(x)
  } else {
    return(na.approx(x, maxgap = 24, na.rm = FALSE))
  }
}

## Linear imputation proceedure

ltb <- ltb %>%
  group_by(spell_id) %>%
  arrange(spell_id, episode_position, time) %>%
  mutate_at(.vars = vars(spo2, pao2, fio2, tidal_vol, peep, mean_ap, rr_spont, peak_ap, pf_ratio),
            .funs = custom_impute)

save.image(file = "####/post_imputation.RData")

# Checking availibility of data that is central to analysis
# PaO2 (and potentially PF, FiO2)

names(ltb)

ltb %>%
  group_by(spell_id) %>%
  summarise(pao2_null = sum(is.na(pao2))/n(),
            fio2_null = sum(is.na(fio2))/n(),
            pf_null = sum(is.na(pf_ratio))/n()) %>%
  filter(pao2_null == 1 & fio2_null == 1 & pf_null == 1) %>%
  select(spell_id) %>%
  left_join(spells, by = "spell_id") %>%
  group_by(site) %>%
  tally

no_data <- ltb %>%
  group_by(spell_id) %>%
  summarise(pao2_null = sum(is.na(pao2))/n()) %>%
  filter(pao2_null == 1) %>%
  select(spell_id) %>%
  pull()

consort_info <- bind_rows(consort_info, tribble(
  ~stage, ~n,
  "at least 1 pao2", nrow(valid_spells_f %>%
                            filter(!(spell_id %in% no_data)) %>% distinct(spell_id))
))

## Add in PF ratios directly when they are missing

ltb <- ltb %>%
  mutate(pf = if_else(!is.na(pao2) & !is.na(fio2), pao2/fio2, as.numeric(NA))) %>%
  mutate(pf = if_else(is.na(pf) & !is.na(pf_ratio), pf_ratio, pf))

ltb %>%
  filter(!is.na(pf)) %>%
  ggplot(aes(pf)) + geom_density()

ltb <- ltb %>%
  mutate(pf = if_else(pf > 100 | pf == 0, as.numeric(NA), pf))

ltb <- select(ltb, -los_hours) %>%
  rename(los_episode = los)

# Recalibrate times over spells
ltb <- ltb %>%
  filter(time >= 0) %>%
  group_by(spell_id) %>%
  arrange(spell_id, episode_position, time) %>%
  mutate(time = seq(0, n()/2 - 0.5, by = 0.5))

ltb <- ltb %>%
  ungroup() %>%
  mutate(day_number = as.integer(time/24))

save.image(file = "####/check_point1.RData")

## Conventional Hyperoxia ----
## AUC by trapezoid rule

## We want to take a vector of x co-ordinates, and a vector of
## y co-ordinates and return the area
## The x co-ordinates are always going to be 30 minutes (0.5 hours)
## So the area under the curve needs to be defined by the difference

## Calcualtes the hourly AUC by trapezoid estimation
## NAs are respected
## a single NA is inserted at the end to account for the fact that areas will
## always produce one less than the number of points

auc <- function(x, y) {
  auc <- apply(cbind(y[-length(y)], y[-1]), 1, mean) * (x[-1] - x[-length(x)])
  return(c(auc, NA))
}

# We want to measure standard hyperoxic dose exposure
# We will use 13.3 because it can only be iatrogenic

# This gives by hour the AUC by trapezoid rule
ltb <- ltb %>%
  group_by(spell_id) %>%
  mutate(pao2_auc = auc(time, pao2))

# This gives us the AUC above the 13.3 line by half hour
ltb <- ltb %>%
  mutate(hyperoxic_13_dose = if_else(pao2_auc - (13.3*0.5) < 0, 0, pao2_auc - (13.3*0.5), missing = 0))

ltb %>%
  select(episode_id, time, pao2, pao2_auc, hyperoxic_13_dose)

testthat::expect_false(any(is.na(ltb$hyperoxic_13_dose)))
## Everything has a value, so we are safe to use cumulative sum

# Now lets just create cumulative oxygen doses
ltb <- ltb %>% 
  group_by(spell_id) %>%
  arrange(spell_id, time) %>%
  mutate(cumulative_hyperoxia_13 = cumsum(hyperoxic_13_dose))

ltb %>%
  filter(!(spell_id %in% no_data)) %>%
  ggplot(aes(cumulative_hyperoxia_13)) + geom_density()

## There is a significant weighting of zeros here
## Note will need to use Royston method of correction.

## Spell LOS ====
valid_spells_f <- valid_spells_f %>%
  rename(episode_los = los)

# Mortality Checks
glimpse(dtb)
table(dtb$unit_mortality, useNA = "always")

# We want just unit mortality, as we cannot account for oxygen
# exposure outside
# of the unit, and therefore it is not a useful outcome measure
# (from a causal inference point of view)

## Fix Apache Score (missing is coded as -1)
dtb <- dtb %>%
  mutate(apache_score = if_else(
    apache_score < 0, as.integer(NA), apache_score, missing = as.integer(NA)))

dtb <- left_join(valid_spells_f, dtb, by = "episode_id")

## Compress Spells ===

## Spells are comprised of more than 1 episode. So we want information from the final
## episode (death) and from the first (apache) where relevent

## Take the LATEST info for death
backwards <- dtb %>%
  group_by(spell_id) %>%
  arrange(spell_id, desc(episode_position)) %>%
  select(spell_id, unit_mortality, hospital_mortality) %>%
  summarise_all(.funs = function(x) unique(na.omit(x))[1])

## Take the EARLIEST info for everything else
forwards <- dtb %>%
  group_by(spell_id) %>%
  arrange(spell_id, episode_position) %>%
  select(-episode_id, -episode_los, -episode_position, -hospital_mortality, -unit_mortality) %>%
  summarise_all(.funs = function(x) unique(na.omit(x))[1])

## Join back together
dtb <- inner_join(backwards, forwards, by = "spell_id")
rm(forwards, backwards)

table(dtb$unit_mortality, useNA = "always")
table(dtb$site, useNA = "always")
# check % missingness

dtb %>%
  ungroup() %>%
  summarise_all(.funs = function(x) (sum(is.na(x))/length(x))*100) %>%
  glimpse()

# Treatment limitation orders

table(dtb$withdrawal, useNA = "always")
table(dtb$organ_donor, useNA = "always")

table(dtb$withdrawal, dtb$organ_donor, useNA = "always")

# can also remove records with no PaO2 data here
# We have done this above already for the consort

dtb %>%
  filter(!(spell_id %in% no_data)) %>%
  group_by(site) %>%
  tally

dtb <- dtb %>%
  filter(!(spell_id %in% no_data)) %>%
  filter(withdrawal == "N" | is.na(withdrawal))

consort_info <- bind_rows(consort_info, tribble(
  ~stage, ~n,
  "without treatment limitation orders", nrow(dtb)
))

## DOB to Age ====

dtb <- dtb %>%
  mutate(age = as.integer(lubridate::year(epi_start_dttm) - lubridate::year(dob)))

dtb %>%
  group_by(site, is.na(apache_score)) %>%
  tally

ox_apache <- read_csv("####/oxford_apacheii_####.csv")

ox_apache <- ox_apache %>%
  rename(apache_score = APACHE_II,
         nhs_number = NHS_number) %>%
  mutate(adm = paste(Unit_admission_date, Unit_admission_time)) %>%
  mutate(adm = lubridate::dmy_hms(adm)) %>%
  select(apache_score, nhs_number, adm) %>%
  mutate(nhs_number = as.character(nhs_number))

## Insert Oxford APACHE scores here
ox <- filter(dtb, site == "OUH")
others <- filter(dtb, site != "OUH")

distinct(ox, spell_id)

ox <- left_join(ox %>% select(-apache_score), ox_apache, by = "nhs_number") %>%
  filter(epi_start_dttm >= adm - hours(6) & epi_start_dttm <= adm + hours(6)) %>%
  # we will relax the admission time by 6 hours either way, as this is the softening factor we use for
  # spell identification
  select(spell_id, apache_score) %>%
  distinct(spell_id, .keep_all = TRUE) %>%
  left_join(ox %>% select(-apache_score), ., by = "spell_id")

dtb <- bind_rows(ox, others)
ggplot(dtb, aes(apache_score, fill = site)) + geom_density(alpha = 0.5)

## Markers of specific disease processes
## COPD

# Fix formatting issues
dtb <- dtb %>%
  separate(col = primary_admission_reason, into = c(paste0("icnarc_1_", LETTERS[1:5])), sep = "\\.") %>%
  separate(col = secondary_admission_reason, into = c(paste0("icnarc_2_", LETTERS[1:5])), sep = "\\.") %>%
  separate(col = pmhx, into = c(paste0("icnarc_3_", LETTERS[1:5])), sep = "\\.") %>%
  mutate_at(.vars = vars(contains("icnarc")), .funs = str_pad, width = 2, side = "left", pad = "0") %>%
  mutate(icn_1 = icnarc_1_A,
         icn_2 = icnarc_1_B) %>%
  unite(col = primary_admission_reason, icnarc_1_A:icnarc_1_E, sep = ".", remove = TRUE) %>%
  unite(col = secondary_admission_reason, icnarc_2_A:icnarc_2_E, sep = ".", remove = TRUE) %>%
  unite(col = pmhx, icnarc_3_A:icnarc_3_E, sep = ".", remove = TRUE)

copd <- c("01.01.02.27.02", "02.01.02.27.02", "01.01.02.30.06", "02.01.02.30.06")

dtb <- dtb %>%
  mutate(is_medical = if_else(icn_1 == "02", 1L, 0L, missing = as.integer(NA)))

dtb <- dtb %>%
  mutate(has_copd = if_else(primary_admission_reason %in% copd |
                              secondary_admission_reason %in% copd |
                              pmhx %in% copd, 1L, 0L, missing = as.integer(0L)))

table(dtb$icn_1, dtb$is_medical, useNA = "always")

dtb <- dtb %>%
  mutate(is_medical = if_else(
    is.na(icn_1) & !is.na(surgical_classification),
    0L, is_medical
  ))

# fix weight

dtb <- dtb %>%
  mutate(weight = if_else(weight < 20, as.numeric(NA), weight))

# Icnarc codes

icnarc_codes <- tribble(
  ~level2, ~system,
  "01", "Respiratory",
  "02", "Cardiovascular",
  "03", "Gastrointestinal",
  "04", "Neurological",
  "06", "Poisoning",
  "07", "Genito, urinary",
  "08", "Endocrine, Metabolic, Thermoregulation and Poisoning",
  "09", "Haematological/Immunological",
  "10", "Musculoskeletal",
  "11", "Dermatological",
  "12", "Psychiatric",
  "13", "Trauma"
)

dtb <- dtb %>%
  rename(system_code = icn_2)

dtb <- left_join(dtb, icnarc_codes, by = c("system_code" = "level2"))

dtb <- dtb %>%
  mutate(unit_mortality = factor(
    unit_mortality,
    levels = c("A", "D"),
    labels = c("alive", "deceased"))) %>%
  mutate(hospital_mortality = factor(
    hospital_mortality,
    levels = c("A", "D"),
    labels = c("alive", "deceased"))) %>%
  mutate(is_medical = factor(
    is_medical,
    levels = c(0L, 1L),
    labels = c("surgical", "medical")
  )) %>%
  mutate(location_in = case_when(
    location_in == "B" ~ "other",
    location_in == "C" ~ "other",
    location_in == "E" ~ "ed",
    location_in == "G" ~ "other",
    location_in == "H" ~ "icu/hdu",
    location_in == "I" ~ "icu/hdu",
    location_in == "M" ~ "other",
    location_in == "R" ~ "theatres",
    location_in == "S" ~ "other",
    location_in == "T" ~ "theatres",
    location_in == "U" ~ "icu/hdu",
    location_in == "W" ~ "ward"
  )) %>%
  mutate(location_in = factor(
    location_in,
    levels=c("ed", "icu/hdu", "theatres", "ward", "other"))) %>%
  mutate(planned_admission = if_else(
    admission_type %in% c("L", "U"), "unplanned", "planned", missing = as.character(NA)
  )) %>%
  mutate(planned_admission = factor(
    planned_admission,
    levels = c("unplanned", "planned"),
    labels = c("unplanned", "planned")
  )) %>%
  mutate(surgical_classification = factor(
    surgical_classification,
    levels = c("L", "S", "U", "M"),
    labels = c("elective", "scheduled", "urgent", "emergency"),
    ordered = TRUE
  )) %>%
  mutate(prior_dependency = if_else(
    prior_dependency %in% c("J", "N", "T"), "any", "none", missing = as.character(NA))) %>%
  mutate(prior_dependency = factor(
    prior_dependency,
    levels = c("any", "none"),
    labels = c("any", "none")
  )) %>%
  mutate(ethnicity = factor(
    ethnicity,
    levels = c("A", "B", "C", "D", "E", "F", "G", "H",
               "J", "K", "L", "M", "N", "P", "R", "S",
               "Z"),
    labels = c("white British",
               "white Irish",
               "white other",
               "mixed white and black Carribean",
               "mixed white and black African",
               "mixed white and Asian",
               "mixed any other",
               "Asian/Asian British Indian",
               "Asian/Asian British Pakistani",
               "Asian/Asian British Bangladeshi",
               "Asian/Asian British other",
               "black/black British Carribean",
               "black/black British African",
               "black/black British other",
               "other ethnic group-Chinese",
               "any other ethnic group",
               "not stated"))
  ) %>%
  mutate(cpr = if_else(cpr == "N", "no", "yes")) %>%
  mutate(cpr = factor(
    cpr,
    levels = c("yes", "no"),
    labels = c("yes", "no")
  )) %>%
  mutate(sex = factor(
    sex,
    levels = c("F", "M"),
    labels = c("female", "male"))) %>%
  mutate(system = factor(
    system,
    levels = icnarc_codes$system,
    labels = icnarc_codes$system)) %>%
  mutate(has_copd = factor(
    has_copd,
    levels = c(0, 1),
    labels = c("no", "yes"))) %>%
  mutate(site = factor(
    site,
    levels = c("GSTT", "RYJ", "UCL", "OUH"),
    labels = c("Guy's and St. Thomas'", "Imperial", "University College London", "Oxford")))

# Lets grab our cohorts
# days 1, 3, 5, 7

## Remember the cumulative hyperoxia is a measure on a
## 30 minute cadance
## So for a tw mean, we need to divide by the number of hours
## We will refer to this as "hyperoxia dose"

cohort_1 <- ltb %>%
  filter(time == 1*24) %>%
  select(spell_id, time, cumulative_hyperoxia_13) %>%
  mutate(tw_hyperoxia_13 = cumulative_hyperoxia_13/time) %>%
  select(spell_id, cumulative_hyperoxia_13, tw_hyperoxia_13)

cohort_2 <- ltb %>%
  filter(time == 3*24) %>%
  select(spell_id, time, cumulative_hyperoxia_13) %>%
  mutate(tw_hyperoxia_13 = cumulative_hyperoxia_13/time) %>%
  select(spell_id, cumulative_hyperoxia_13, tw_hyperoxia_13)

cohort_3 <- ltb %>%
  filter(time == 5*24) %>%
  select(spell_id, time, cumulative_hyperoxia_13) %>%
  mutate(tw_hyperoxia_13 = cumulative_hyperoxia_13/time) %>%
  select(spell_id, cumulative_hyperoxia_13, tw_hyperoxia_13)

cohort_4 <- ltb %>%
  filter(time == 7*24) %>%
  select(spell_id, time, cumulative_hyperoxia_13) %>%
  mutate(tw_hyperoxia_13 = cumulative_hyperoxia_13/time) %>%
  select(spell_id, cumulative_hyperoxia_13, tw_hyperoxia_13)

cohort_1 <- inner_join(dtb, cohort_1, by = "spell_id")
cohort_2 <- inner_join(dtb, cohort_2, by = "spell_id")
cohort_3 <- inner_join(dtb, cohort_3, by = "spell_id")
cohort_4 <- inner_join(dtb, cohort_4, by = "spell_id")

consort_info <- bind_rows(consort_info, tribble(
  ~stage, ~n,
  "los >= 1 day", nrow(cohort_1),
  "los >= 3 days", nrow(cohort_2),
  "los >= 5 days", nrow(cohort_3),
  "los >= 7 days", nrow(cohort_4)
))


# Adding indicator variable for zero hyperoxia exposure (royston)

cohort_1 <- cohort_1 %>%
  mutate(hyperoxia_13 = if_else(cumulative_hyperoxia_13 != 0, 1L, 0L))
cohort_2 <- cohort_2 %>%
  mutate(hyperoxia_13 = if_else(cumulative_hyperoxia_13 != 0, 1L, 0L))
cohort_3 <- cohort_3 %>%
  mutate(hyperoxia_13 = if_else(cumulative_hyperoxia_13 != 0, 1L, 0L))
cohort_4 <- cohort_4 %>%
  mutate(hyperoxia_13 = if_else(cumulative_hyperoxia_13 != 0, 1L, 0L))

cohort_1 %>%
  ggplot(aes(cumulative_hyperoxia_13)) + geom_histogram()

cohort_1 %>%
  arrange(desc(cumulative_hyperoxia_13)) %>%
  select(spell_id, cumulative_hyperoxia_13, tw_hyperoxia_13)

exemplar <- ltb %>%
  filter(spell_id == "#####") %>%
  ungroup() %>%
  select(time, fio2, spo2, pao2) %>%
  mutate(fio2 = fio2 * 100) %>%
  gather(key = monitor, value = value, fio2:pao2)

exemplar %>%
  ggplot(aes(x = time, y = value, group = monitor, colour = monitor)) +
  geom_path() +
  geom_hline(yintercept = 13.3, linetype = 2)

save.image("####/end_data_wrangling.RData")
save(cohort_1, cohort_2, cohort_3, cohort_4, consort_info, file = "####/cohorts_final.RData")
## Moving to table 1 script ====
