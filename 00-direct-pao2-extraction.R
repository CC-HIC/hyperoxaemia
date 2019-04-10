library(tidyverse)
library(inspectEHR)
library(lubridate)

## Location of the supplementary files
camf <- list.files("####")[grepl(pattern = "csv", x = list.files("####"))]

## initialise, read in and bind these tables
for (file in seq_along(camf)) {
  
  if (file == 1) {
    df <- read_csv(paste0("####", camf[file]))
  } else {
    temp <- read_csv(paste0("####", camf[file]))
    df <- bind_rows(df, temp)  
  }
  
}

## data is in long form, separate out what we need
poc_label <- df %>%
  filter(type == "POC BLOOD SPECIMEN TYPE") %>%
  select(nhs, date, value) %>%
  rename(specimen = value)

poc_pco2 <- df %>%
  filter(type == "POC PCO2 TEMP") %>%
  select(nhs, date, value) %>%
  rename(co2 = value)

poc_ph <- df %>%
  filter(type == "POC PH") %>%
  select(nhs, date, value) %>%
  rename(ph = value)

poc_po2 <- df %>%
  filter(type == "POC PO2 TEMP") %>%
  select(nhs, date, value) %>%
  rename(o2 = value)

## Join back together
cdf <- reduce(
  list(poc_label, poc_pco2, poc_ph, poc_po2),
  full_join,
  by = c("nhs" = "nhs", "date" = "date"))

## modify formatting to match that of the database structure
cdf <- cdf %>%
  select(nhs, date, specimen, o2) %>%
  mutate(o2 = as.numeric(o2)) %>%
  rename(nhs_number = nhs, datetime = date, real = o2) %>%
  mutate(integer = case_when(
    specimen == "Arterial blood" ~ as.integer(1),
    specimen == "Venous blood" ~ as.integer(2),
    TRUE ~ as.integer(3)
  )) %>%
  select(-specimen) %>%
  filter(real < 100)

## Check this follows the distribution we are expecting
cdf %>%
  ggplot(aes(real, fill = as.factor(integer))) +
  geom_density(alpha = 0.5)

## The basic work is now done, now we need to get this back into
## the main analysis

## Establish a DB connection
ctn <- connect(username = "####",
               password = "####",
               database = "####",
               system = "####")
tbls <- retrieve_tables(ctn)

## Setup
core <- make_core(ctn)
reference <- make_reference(ctn)
episode_length <- epi_length(core, reference, tbls[["events"]])
episodes <- collect(tbls[["episodes"]])
spells <- identify_spells(episode_length, episodes, minutes = 360)

##  Modify data type
cdf <- mutate(cdf, nhs_number = as.character(nhs_number))

## Sanity Check
cdf %>%
  filter(integer == 1 | is.na(integer) | integer == 3) %>%
  group_by(nhs_number) %>%
  tally()

spells %>% filter(nhs_number %in% unique(cdf$nhs_number))

full_join(spells, cdf, by = "nhs_number") %>%
  filter(datetime >= epi_start_dttm,
         datetime <= epi_end_dttm) %>%
  filter(site == "RGT") %>%
  distinct(episode_id)

## Manipulation work to get things into the long-form data structure
df <- left_join(cdf, left_join(episode_length %>% select(-los, -validity, -site),
                               episodes %>% select(-start_date, -provenance),
                by = "episode_id"), by = "nhs_number") %>%
                  filter(datetime >= epi_start_dttm,
                         datetime <= epi_end_dttm) %>%
  select(real, integer, datetime, episode_id, epi_start_dttm)

cam_start <- select(df, episode_id, epi_start_dttm)
df <- select(df, -epi_start_dttm)
cam_start <- distinct(cam_start, episode_id, .keep_all = TRUE)

df <- df %>%
  mutate(code_name = "NIHR_HIC_ICU_0132")

cam_start <- cam_start %>%
  rename(datetime = epi_start_dttm) %>%
  mutate(real = as.numeric(NA),
         integer = as.integer(NA),
         code_name = "NIHR_HIC_ICU_0411")

## This is now in the form expected by extract_tyimevarying:
## this can now do the heavy lifting for us
add_imp <- bind_rows(cam_start, df) %>%
  arrange(episode_id, datetime) %>%
  mutate(event_id = seq(1, n()))

ltb_imp <- extract_timevarying(
  events = add_imp,
  metadata = collect(tbls[["variables"]]),
  code_names = c("NIHR_HIC_ICU_0411", "NIHR_HIC_ICU_0132"),
  rename = c("start_dttm", "pxo2"),
  chunk_size = 10000,
  cadance = 0.5)

### And APACHE Scores
camf <- list.files("####")

for (file in seq_along(camf)) {
  
  if (file == 1) {
    df <- read_csv(paste0("####", camf[file]))
  } else {
    temp <- read_csv(paste0("####", camf[file]))
    df <- bind_rows(df, temp)  
  }
  
}

## Miscalculation in the PF score

df <- df %>%
  tidyr::separate(col = `fio2/pao2`, into = c("fio2", "pao2"), sep = "/") %>%
  mutate_at(vars(fio2, pao2), as.numeric) %>%
  mutate(oxy_score = case_when(
    pao2 <= 55/7.6 ~ 4,
    pao2 <= 60/7.6 ~ 3,
    pao2 <= 70/7.6 ~ 1,
    pao2 > 70/7.6 ~ 0
  )) %>%
  rowwise() %>%
  mutate(apache_score = `APACHE II` - `oxygenation score` + oxy_score)

apch <- select(df, apache_score, nhs, adm, dis) %>%
  rename(nhs_number = nhs) %>%
  mutate(nhs_number = as.character(nhs_number)) %>%
  ungroup()

## No further work needed here as there is no time info to wrestle with
save(ltb_imp, apch, file = "####/cam_additional.RData")

