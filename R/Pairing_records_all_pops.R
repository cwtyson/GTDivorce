## Process standard format data to get brooding data for all pops

############################
## Housekeeping
############################

library(tidyverse)
library(here)

############################
## Process brood data
############################

## Read in data in standard format
standard_format_data <- readRDS(paste0(here::here(), "/data/standard_format/20210802/standard_format.RDS"))


## Keep only brood data
brood_dat <- standard_format_data[[1]]

## Summarize number of cases where female has more than one social partner in one year
clutch_sum <- brood_dat %>%

  filter(Species == "PARMAJ") %>%
  filter(!is.na(FemaleID) | !is.na(MaleID)) %>%
  select(PopID, BreedingSeason, LocationID, ClutchType_calculated, FemaleID, MaleID, Brood_exp = ExperimentID) %>%
  filter(!is.na(FemaleID)) %>%
  group_by(PopID, FemaleID, BreedingSeason) %>%
  summarise(num_partners = n_distinct(MaleID),
            num_partners_known = n_distinct(MaleID, na.rm = T))

table(clutch_sum$PopID, clutch_sum$num_partners)

## Get brood data for GTs when either male or female is known
brood_dat_ref <- brood_dat %>%

  filter(Species == "PARMAJ") %>%
  filter(!is.na(FemaleID) | !is.na(MaleID)) %>%
  # distinct(PopID, BreedingSeason, paste(FemaleID, MaleID), .keep_all = T) %>%
  mutate(ExperimentID = case_when(ExperimentID == "NA;NA;NA" ~ FALSE,
                                  ExperimentID == "NA" ~ FALSE,
                                  ExperimentID == "FALSE" ~ FALSE,
                                  ExperimentID == "" ~ FALSE,
                                  ExperimentID == "UNKNOWN" ~ FALSE,
                                  is.na(ExperimentID) ~ FALSE,
                                  grepl("control", ExperimentID) ~ FALSE,
                                  !is.na(ExperimentID) ~ TRUE)) %>%
  select(PopID, BreedingSeason, LocationID, FemaleID, MaleID, Brood_exp = ExperimentID, contains("Date_observed"), ClutchType_calculated) %>%

  ## Set 'unknown' IDs to NA: 1000000
  mutate(FemaleID = case_when(FemaleID == "1000000" ~ NA_character_,
                              TRUE ~ FemaleID),
         MaleID = case_when(MaleID %in% c("1000000", "UNKNOWN") ~ NA_character_,
                              TRUE ~ MaleID)) %>%

  mutate_at(vars(FemaleID, MaleID), .funs = toupper)


## For each female, determine if she had a known breeding partner in the year
brood_dat_pairs <- brood_dat_ref %>%
  group_by(PopID, FemaleID, BreedingSeason) %>%
  mutate(Num_male_partners = case_when(!is.na(FemaleID) ~ n_distinct(MaleID, na.rm = T)),
         MaleID_known = purrr::map_chr(.x = list(unique(stats::na.omit(.data$MaleID))),
                                       .f = ~{
                                         if(length(..1) == 0){
                                           return(NA_character_)
                                         } else if(length(..1) == 1){
                                           return(..1)
                                         } else {
                                           return(NA_character_)
                                         }
                                       }),
         MaleID_recode = case_when(Num_male_partners == 2 ~ MaleID,
                                   TRUE ~ MaleID_known)) %>%

  ## For each male, determine if he had a known breeding partner in the year
  group_by(PopID, MaleID, BreedingSeason) %>%
  mutate(Num_female_partners = case_when(!is.na(MaleID) ~ n_distinct(FemaleID, na.rm = T)),
         FemaleID_known = purrr::map_chr(.x = list(unique(stats::na.omit(.data$FemaleID))),
                                       .f = ~{
                                         if(length(..1) == 0){
                                           return(NA_character_)
                                         } else if(length(..1) == 1){
                                           return(..1)
                                         } else {
                                           return(NA_character_)
                                         }
                                       }),
         FemaleID_recode = case_when(Num_female_partners == 2 ~ FemaleID,
                                   TRUE ~ FemaleID_known)) %>%

  ## Within each year, keep only distinct pairs
  distinct(PopID, BreedingSeason, FemaleID, MaleID, .keep_all = T)

  ##

  ## Determine pairing order for males
  brood_dast_pairs_ordered <- brood_dat_pairs %>%
  filter(Num_female_partners == 2) %>%
  # group_by(PopID, MaleID, BreedingSeason) %>%
  arrange(PopID, MaleID, BreedingSeason, ClutchType_calculated)




## Get first year breeding for each individual
first_breed_year <- brood_dat_pairs %>%
  select(PopID, BreedingSeason, IndvID = FemaleID) %>%
  bind_rows(brood_dat_ref %>%
              select(PopID, BreedingSeason, IndvID = MaleID)) %>%
  na.omit() %>%
  group_by(PopID, IndvID) %>%
  summarise(breed_year = min(BreedingSeason)) %>%
  distinct(PopID, IndvID, .keep_all = T)


## Get capture records
cap_dat <- standard_format_data[[2]]

## Create data frame of capture records after the breeding year with information on experiments
cap_dat_exp <- cap_dat %>%
  filter(Species == "PARMAJ") %>%

  ## Keep individuals that are in the breeding records
  filter(paste(CapturePopID, IndvID) %in% paste(first_breed_year$PopID, first_breed_year$IndvID)) %>%

  mutate(ExperimentID = case_when(ExperimentID == "NA;NA;NA" ~ NA_character_,
                                  ExperimentID == "NA" ~ NA_character_,
                                  ExperimentID == "FALSE" ~ NA_character_,
                                  grepl("control", ExperimentID) ~ NA_character_,
                                  TRUE ~ ExperimentID)) %>%
  select(PopID = CapturePopID, BreedingSeason, IndvID, Sex_observed, ExperimentID, CaptureDate) %>%

  ## Join in first breeding year
  left_join(first_breed_year,
            by = c("PopID", "IndvID")) %>%

  ## Remove records from before first breeding year
  filter(BreedingSeason >= breed_year) %>%

  group_by(PopID, IndvID) %>%
  mutate(Cap_exp = case_when(any(!is.na(ExperimentID)) ~ TRUE,
                             TRUE ~ FALSE),
         years = n_distinct(BreedingSeason),
         capture = "yes") %>%

  distinct(PopID, IndvID, BreedingSeason, .keep_all = T) %>%
  select(-ExperimentID, -breed_year)


## Join missing records from the capture data of individuals in the brood data
pair_dat <- brood_dat_pairs %>%
  mutate(brood = "yes") %>%

  ## Join capture records
  full_join(cap_dat_exp %>%
              filter(Sex_observed == "F"),
            by = c("PopID", "BreedingSeason", "FemaleID" = "IndvID")) %>%
  select(-Sex_observed) %>%
  rename(Record_f = capture,
         Female_exp = Cap_exp,
         Female_cap_dat = CaptureDate) %>%

  ## Join records from capture data
  full_join(cap_dat_exp %>%
              filter(Sex_observed == "M"),
            by = c("PopID", "BreedingSeason", "MaleID" = "IndvID")) %>%
  select(-Sex_observed) %>%
  rename(Record_m = capture,
         Male_exp = Cap_exp,
         Male_cap_dat = CaptureDate) %>%
  distinct(PopID, BreedingSeason, FemaleID, MaleID, .keep_all = T) %>%

  select(PopID, BreedingSeason, LocationID, ClutchType_calculated, contains("Date_observed"), Brood_exp, brood,
         FemaleID,Female_cap_dat, Record_f, Num_male_partners, Female_exp,
         MaleID, Male_cap_dat, Record_m, Num_female_partners, Male_exp) %>%
  mutate(across(c(Brood_exp, Female_exp, Male_exp), ~ replace_na(., FALSE))) %>%

  ## Label records that were added from captures
  group_by(PopID, FemaleID) %>%
  mutate(extra_f = any(is.na(brood))) %>%
  group_by(PopID, MaleID) %>%
  mutate(extra_m = any(is.na(brood))) %>%

  ## If an individual's nest was experimental, set individual experiment to TRUE and recalculate whether an individual was ever experimented on
  mutate(Female_exp = case_when(Brood_exp == TRUE ~ TRUE,
                                TRUE ~ Female_exp),
         Male_exp = case_when(Brood_exp == TRUE ~ TRUE,
                              TRUE ~ Male_exp)) %>%
  group_by(PopID, FemaleID) %>%
  mutate(Female_exp = case_when(any(Female_exp) ~ TRUE,
                                  TRUE ~ FALSE),
         Male_exp = case_when(any(Male_exp) ~ TRUE,
                                TRUE ~ FALSE))

## Create order within each pair
pair_dat_ref <- pair_dat %>%
  mutate(order = case_when(ClutchType_calculated == "first" ~ 1,
                           ClutchType_calculated == "replacement" ~ 2,
                           ClutchType_calculated == "second" ~ 3,
                           is.na(ClutchType_calculated) ~ 4))  %>%
  arrange(PopID, BreedingSeason, FemaleID, order) %>%
  group_by(PopID, BreedingSeason, FemaleID) %>%
  mutate(no_order_f = ifelse(n_distinct(order) == 1 & Num_male_partners == 2, "yes", "no")) %>%
  group_by(PopID, BreedingSeason, MaleID) %>%
  mutate(no_order_m = ifelse(n_distinct(order) == 1 & Num_female_partners == 2, "yes", "no"))

## Summarize experimental birds
brood_exp_sum <- pair_dat_ref %>%
  group_by(PopID) %>%
  summarise(Females_exp = n_distinct(FemaleID[Female_exp == TRUE], na.rm = T),
            Females_exp_y = sum(Female_exp == TRUE),
            Males_exp = n_distinct(MaleID[Male_exp == TRUE], na.rm = T),
            Males_exp_y = sum(Male_exp == TRUE),
            Brood_exp = sum(Brood_exp == TRUE))


## Filter out experimental birds
pair_dat_filter <- pair_dat_ref %>%
  filter(Brood_exp == FALSE,
         Female_exp == FALSE,
         Male_exp == FALSE) %>%
  select(pop = PopID,
         nest = LocationID,
         year = BreedingSeason,
         female = FemaleID,
         f_pair_num = Num_male_partners,
         male = MaleID,
         m_pair_num = Num_female_partners,
         order,
         no_order_f,
         no_order_m)

## Save full file
readr::write_csv(pair_dat_filter, file = paste0(here(),"/data/pop_pairs/all_pop_pairs.csv"))
