## Script to process standard format Brood data and create a metric of the average number of chicks fledged per nest

############################
## Housekeeping
############################

library(tidyverse)
library(here)
library(wesanderson)

## Function for when all NA
sumna <- function(x) {
  sumna <- NULL
  return(ifelse(all(is.na(x)), NA, sum(na.omit(x))))
}

############################
## Process fledging counts to get measure of environmental data
############################

## Read in data in standard format
bro_dat <- read.csv(paste0(here::here(), "/data/standard_format/20210531/standard_format_Brood_data.csv"))


## Filter full data
bro_dat_f <- bro_dat %>%

  ## Merge columns across protocol versions
  mutate(NumberFledged = case_when(is.na(NumberFledged_observed) ~  NumberFledged,
                                   TRUE ~ NumberFledged_observed),
         ClutchSize = case_when(is.na(ClutchSize_observed) ~  ClutchSize,
                                TRUE ~ ClutchSize_observed),
         BroodSize = case_when(is.na(BroodSize_observed) ~  BroodSize,
                               TRUE ~ BroodSize_observed),
         LayDate_observed = case_when(is.na(LayDate_observed) ~  LayDate,
                                      TRUE ~ LayDate_observed)) %>%

  ## Select relevant columns
  select(Species, PopID, BreedingSeason, LocationID, FemaleID, MaleID,
         ClutchType_calculated, LayDate_observed, ClutchSize, BroodSize, NumberFledged,
         ExperimentID) %>%

  ## Recode ExperimentID
  mutate(ExperimentID = case_when(ExperimentID == "NA;NA;NA" ~ NA_character_,
                                  ExperimentID == "NA" ~ NA_character_,
                                  ExperimentID == "FALSE" ~ NA_character_,
                                  grepl("control", ExperimentID) ~ NA_character_,
                                  TRUE ~ ExperimentID)) %>%

  ## Filter
  filter(Species == "PARMAJ",
         is.na(ExperimentID)) %>%
  select(-ExperimentID) %>%
  arrange(PopID)


## Calculate fledge metrics
bro_dat_sum <- bro_dat_f %>%

  ## Get total fledged in first or replacement clutches
  group_by(PopID, BreedingSeason) %>%
  mutate(first_clutch_avg = mean(NumberFledged[ClutchType_calculated == "first"], na.rm = T),
         second_clutch_avg = mean(NumberFledged[ ClutchType_calculated == "second" | ClutchType_calculated == "replacement" ], na.rm = T),
         prop_second_clutch = n_distinct(LocationID[ClutchType_calculated == "second"])/
           max(n()),
         clutches = max(n())) %>%
  rowwise() %>%
  mutate(avg_sum_clutch = sum(first_clutch_avg, second_clutch_avg, na.rm = T)) %>%

  ## Total fledged by female
  group_by(PopID, BreedingSeason, FemaleID) %>%
  mutate(tot_fledged_fem = case_when(!is.na(FemaleID) & !(FemaleID %in% c(0,1,"1e+06")) ~ as.integer(sumna(NumberFledged)),
                                     TRUE ~ NA_integer_)) %>%

  ## Average totals each year
  group_by(PopID, BreedingSeason) %>%
  mutate(females = n_distinct(FemaleID, na.rm = T),
         fledge_fem_avg = mean(tot_fledged_fem, na.rm = T)) %>%

  ## Total by nest box
  group_by(PopID, BreedingSeason, LocationID) %>%
  mutate(tot_fledged_nest = sumna(NumberFledged)) %>%

  ## Average totals each year
  group_by(PopID, BreedingSeason) %>%
  mutate(nests = n_distinct(LocationID, na.rm = T),
         fledge_nest_avg = mean(tot_fledged_nest, na.rm = T)) %>%

  select(PopID, BreedingSeason,
         clutches, prop_second_clutch,
         avg_sum_clutch,
         females, fledge_fem_avg,
         nests, fledge_nest_avg) %>%
  distinct(PopID, BreedingSeason, .keep_all = T) %>%

  ## Get difference in avg fledge counts
  mutate(fem_clutch_ratio = females / clutches,
         fledge_ratio = fledge_fem_avg/fledge_nest_avg) %>%
  group_by(PopID) %>%
  mutate(fledge_ratio_avg = round(mean(fledge_ratio, na.rm = T), 2))


## Adjust based on relationship between fledge metric ratio and proportion of second clutches (by pop)
bro_dat_sum_adj <- bro_dat_sum %>%
  group_by(PopID) %>%

  ## Create PopID for populations with much information on
  mutate(prop_NA = sum(is.na(fledge_ratio))/length(fledge_ratio),
         PopID_NA = ifelse(prop_NA > 0.50, "NA", PopID)) %>%

  ## Fit model for each PopID_NA
  group_by(PopID_NA) %>%
  nest() %>%
  mutate(fit = map(data, ~ lm(fledge_ratio ~ prop_second_clutch, data = .)),
         fledge_ratio_pred = map2(.x = fit, .y = data, ~predict(.x, newdata = .y))) %>%
  unnest(cols = c(data, fledge_ratio_pred)) %>%
  select(-fit) %>%

  ## Adjust fledge_nest_avg based on prediction for each population
  group_by(PopID) %>%
  mutate(fledge_nest_avg_adj = fledge_nest_avg*fledge_ratio_pred)


## Correlations
bro_dat_sum_cor <- bro_dat_sum_adj %>%
  filter(!is.na(fledge_fem_avg) & !is.na(fledge_nest_avg_adj)) %>%
  group_by(PopID) %>%
  mutate(cor = round(as.numeric(cor(fledge_fem_avg, fledge_nest_avg_adj, use = "complete.obs", method = "spearman")),2),
         count = n_distinct(BreedingSeason))

weighted.mean(x = as.numeric(names(bro_dat_sum_cor %>% distinct(PopID, .keep_all = T) %>% pull(name = cor))),
              w = bro_dat_sum_cor %>% distinct(PopID, .keep_all = T) %>% pull(name = count),
              na.rm = T)


## Visualize: Scatterplot of 3 metrics
pal <- wes_palette("Zissou1", 10, type = "continuous")
bro_dat_sum_adj %>%
  pivot_longer(cols = c(fledge_fem_avg, fledge_nest_avg, fledge_nest_avg_adj), names_to = "metric") %>%
  # filter(PopID == "MIS") %>%
  ggplot() +
  geom_point(aes(x = BreedingSeason, y = value, color = metric), alpha = 0.7) +
  labs(x = "Year", y = "# Fledged", color = "Metric") +
  theme_minimal() +
  facet_wrap(~PopID, scales = "free", drop = T) +
  scale_color_manual(breaks = c("fledge_fem_avg", "fledge_nest_avg", "fledge_nest_avg_adj"),
                     labels = c("Per female", "Per nest", "Per nest adjusted"),
                     values = c(pal[6],pal[1],pal[10]))+
  theme(axis.text.x = element_text(angle = 90, size = 6),
        legend.position = c(0.75,0.05),
        legend.direction = "horizontal")
ggsave(filename = "plots/Fledge_metrics.png", width = 15, units = "in", height = 10)


## Standardize
bro_dat_std <- bro_dat_sum_adj %>%

  ## Rename
  rename(pop_year_fledge_avg = fledge_nest_avg_adj) %>%

  ## Calculate population avg and sd
  group_by(PopID) %>%
  mutate(pop_tot_fledge_avg = mean(pop_year_fledge_avg, na.rm = T),
         pop_tot_fledge_sd = sd(pop_year_fledge_avg, na.rm = T)) %>%

  ## Standardize for each pop (could be done with the function scale, but this makes it explicit)
  mutate(envStd_fledge_num = (pop_year_fledge_avg - pop_tot_fledge_avg) / pop_tot_fledge_sd) %>%
  distinct(PopID, BreedingSeason, .keep_all = T) %>%
  select(PopID, BreedingSeason, pop_year_fledge_avg, pop_tot_fledge_avg, pop_tot_fledge_sd,  envStd_fledge_num)


## Calculate lag
bro_dat_std_lag <- bro_dat_std %>%

  ## Calculate absolute change in envStd values
  arrange(PopID, BreedingSeason) %>%
  group_by(PopID) %>%
  mutate(envStd_fledge_num_lag = abs(envStd_fledge_num - lag(envStd_fledge_num))) %>%

  ## Standardize envStd values
  group_by(PopID) %>%
  mutate(gm_envStd_fledge_num_lag = mean(envStd_fledge_num_lag, na.rm = T),
         gsd_envStd_fledge_num_lag = sd(envStd_fledge_num_lag, na.rm = T)) %>%
  group_by(PopID, BreedingSeason) %>%
  mutate(Ztrans_envStd_fledge_num_lag = (envStd_fledge_num_lag - gm_envStd_fledge_num_lag) / gsd_envStd_fledge_num_lag) %>%

  select(PopID, BreedingSeason, envStd_fledge_num, Ztrans_envStd_fledge_num_lag)


## Save
write.csv(bro_dat_std_lag, file = "data/processed/GT_fledge_counts_std.csv", row.names = F)

############################
## Visualize
############################

## Visualize standardized environmental measures
bro_dat_std_lag <- read.csv("data/processed/GT_fledge_counts_std.csv")


## Visualize change between years in standardized environmental
bro_dat_std_lag %>%
  ggplot() +

  ## Fledge counts
  geom_line(aes(x = as.numeric(BreedingSeason), y = envStd_fledge_num)) +
  geom_point(aes(x = as.numeric(BreedingSeason), y = envStd_fledge_num, color = Ztrans_envStd_fledge_num_lag)) +

  facet_wrap(~PopID, scales= "free", shrink = T, drop = T) +
  labs(x = "Breeding Season",
       y = "Standardized average fledging count \n([annual mean - pop mean] / pop sd)",
       color = "Standardized difference \nin average fledge count \nfrom previous year") +
  theme_minimal() +
  scale_color_gradient2(low = "blue",
                        mid = "grey",
                        high = "red",
                        na.value = NA) +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        legend.position = c(0.75,0.05),
        legend.direction = "horizontal")
ggsave(filename = "plots/GT_flege_counts_plot.png", width = 15, units = "in", height = 10)

