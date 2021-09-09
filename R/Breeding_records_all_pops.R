## Breeding records with a known individual for all populations

############################
## Housekeeping
############################

library(tidyverse)
library(here)

############################
## Process pair data
############################

## Read in data in standard format
pair_dat <- read.csv(paste0(here(), "/data/standard_format/20210531/standard_format_Brood_data.csv"))

## Reformat
pair_dat_ref <- pair_dat %>%
  select(Species, PopID, BreedingSeason, LocationID, FemaleID, MaleID, LayDate_observed, ClutchType_calculated, ExperimentID) %>%
  filter(Species == "PARMAJ" &

           ## Keep records with either male or female ID
           (!is.na(FemaleID) | !is.na(MaleID)) &

           ## Drop suspicious IDs - these are clearly not unique within the Pop and do not refer to different individuals
           !(FemaleID %in% c(0,1,"1e+06")) &
           !(MaleID %in% c(0,1))) %>%

  ## Recode ExperimentID
  mutate(ExperimentID = case_when(ExperimentID == "NA;NA;NA" ~ NA_character_,
                                  ExperimentID == "NA" ~ NA_character_,
                                  ExperimentID == "FALSE" ~ NA_character_,
                                  grepl("control", ExperimentID) ~ NA_character_,
                                  TRUE ~ ExperimentID),
         Experiment = case_when(is.na(ExperimentID) ~ "n",
                                TRUE ~ "y")) %>%

  distinct(PopID, BreedingSeason, LocationID, FemaleID, MaleID, LayDate_observed, ClutchType_calculated, Experiment) %>%
  arrange(PopID, BreedingSeason, FemaleID, LayDate_observed)

## Save
write.csv(pair_dat_ref, "./data/processed/GT_breeding_records.csv")


## Process pairs
pair_dat_proc <- pair_dat_ref %>%
  group_by(PopID, FemaleID) %>%
  mutate(num_yrs = n_distinct(BreedingSeason[is.na(FemaleID) == F]),
         partner = case_when(is.na(FemaleID) == F & is.na(MaleID) == F &
                               lag(MaleID) == MaleID &
                               lag(BreedingSeason) != BreedingSeason ~ "same",
                             is.na(FemaleID) == F & is.na(MaleID) == F &
                               lag(MaleID) != MaleID &
                               lag(BreedingSeason) != BreedingSeason ~ "switched",
                             TRUE ~ NA_character_),
         ncharID_F = nchar(FemaleID),
         ncharID_M = nchar(MaleID))

## Visualize proportion of pairs that are reuniting each year
(repair_fig <- pair_dat_proc %>%
    filter(ClutchType_calculated == "first") %>%
    group_by(PopID, BreedingSeason) %>%
    mutate(repairs = n_distinct(FemaleID[partner == "same" & !is.na(partner)]),
           tot_pairs = n_distinct(FemaleID[is.na(partner) == F]),
           repair_prop = repairs/tot_pairs) %>%
    group_by(PopID) %>%
    mutate(tot_repairs = sum(repairs)) %>%
    filter(tot_repairs > 0) %>%
    ggplot() +
    geom_line(aes(x = BreedingSeason, y = repair_prop, color = tot_pairs)) +
    theme_minimal() +
    facet_wrap(.~PopID, scales = "free"))





