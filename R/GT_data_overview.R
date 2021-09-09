## Script to summarize data for each population

############################
## Housekeeping
############################

library(tidyverse)
library(here)

############################
## Compare processed data to available data
############################

## Read in summarized data
gt_pop_sum <- read.csv("data/processed/GT_pop_RS_summary.csv")

## Read in population overview sheet
gt_pop_all <- readxl::read_xlsx(path = "data/Population_for_env_data.xlsx")

## Join processed data to see which are still missing
gt_pop_pro <- gt_pop_all %>%
  rename(PopID = pop_code) %>%
  left_join(gt_pop_sum %>% distinct(PopID) %>% mutate(done = "yes") , by = "PopID")

## Get list of populations still missing
gt_pops$pop_code[!(gt_pops$pop_code %in% unique(gt_pop_sum$PopID))]

############################
## Summary info for each pop
############################

## Data from standard format
bro_dat <- read.csv(paste0(here(), "/data/standard_format/standard_format_Brood_data.csv"))

## Summary
pair_year_sum <- bro_dat %>%
  filter(!is.na(MaleID) & !is.na(FemaleID)) %>%
  group_by(PopID) %>%
  summarise(num_pairs = n_distinct(paste(MaleID, FemaleID)),
            num_years = max(BreedingSeason) + 1 - min(BreedingSeason)) %>%
  left_join(gt_pop_pro %>% select(PopID, done), by = "PopID")

## Plot
ggplot(pair_year_sum) +
  geom_col(aes(x= as.factor(BreedingSeason), y = num_pairs)) +
  facet_wrap(~PopID, shrink = T, scales = "free") +
  theme_minimal()

