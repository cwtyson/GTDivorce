## Capture data for all populations

############################
## Housekeeping
############################

library(tidyverse)
library(here)

############################
##
############################

## Read in capture data in standard format
cap_dat <- read.csv(paste0(here(), "/data/standard_format/20210531/standard_format_Capture_data.csv"))

## Filter and reformat
cap_dat_ref <- cap_dat %>%
  filter(Species == "PARMAJ") %>%
  group_by(BreedingSeason, IndvID) %>%

  ## TODO: No experiments in PARMAJ - suspicious
  mutate(Experiment = case_when(any(!is.na(ExperimentID)) ~ "y",
                                TRUE ~ "n")) %>%
  distinct(IndvID, BreedingSeason, .keep_all = T) %>%
  select(CapturePopID, IndvID, BreedingSeason, Age_calculated, Experiment) %>%
  arrange(CapturePopID, IndvID, BreedingSeason) %>%
  na.omit()

## Save
write.csv(cap_dat_ref, "./data/processed/GT_capture_records.csv")

