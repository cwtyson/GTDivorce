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

## Filter
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


# ## Plot of hatch dates by populations
# bro_dat_f %>%
#   mutate(jd = lubridate::yday(lubridate::ymd(LayDate_observed)),
#          min_lay = min(jd, na.rm = T) + 30) %>%
#   # filter(PopID == "MTV") %>%
#   ggplot() +
#   geom_density_ridges(aes(x = jd, y = ClutchType_calculated)) +
#   facet_wrap(~PopID, scales = "free") +
#   theme_ridges()

# ## Proportion of clutch types per pop
# (female_id_tab <- bro_dat_f %>%
#     mutate(id = ifelse(is.na(FemaleID),0,1)) %>%
#     group_by(PopID, ClutchType_calculated) %>%
#     summarise(tot_obs = n(),
#               female_id = sum(id),
#               perc_id = female_id/tot_obs) %>%
#     ggplot() +
#     geom_col(aes(x = ClutchType_calculated, y = perc_id, fill = ClutchType_calculated)) +
#     facet_wrap(~PopID) +
#     theme_minimal() +
#     theme(axis.text.x = element_blank()) +
#     labs(x = NULL, y = "% records with Female ID"))

# ## Summary of Pops
# female_id_tab <- bro_dat_f %>%
#   group_by(PopID, BreedingSeason) %>%
#   summarise(nests = n_distinct(LocationID, na.rm = T),
#             fems = n_distinct(FemaleID, na.rm = T))

## Calculate fledge metrics
bro_dat_sum <- bro_dat_f %>%

  ## Get total fledged in first or replacement clutches
  group_by(PopID, BreedingSeason) %>%
  mutate(first_clutch_avg = mean(NumberFledged[ClutchType_calculated == "first"], na.rm = T),
         # replacement_clutch_avg = mean(NumberFledged[ClutchType_calculated == "replacement"], na.rm = T),
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

## Relationship between fledge metric ratio and proportion of second clutches
bro_dat_sum %>%
  # filter(PopID == "LIE") %>%
  ggplot() +
  geom_text(aes(x = prop_second_clutch, y = fledge_ratio, label = paste(BreedingSeason, sep = "-"), color = PopID), size = 3) +
  stat_smooth(method = "lm", aes(x = prop_second_clutch, y = fledge_ratio)) +
  facet_wrap(~PopID, scales = "free",drop = F) +
  labs(x = "Proportion second clutches",
       y = "Fledge ratio") +
  theme_minimal()

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

## Scatterplot of 3 metrics
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

## Difference between adjusted and raw values
bro_dat_sum_adj %>%
  mutate(raw_diff = fledge_fem_avg - fledge_nest_avg,
         adj_diff = fledge_fem_avg - fledge_nest_avg_adj,
         diff_diff = adj_diff - raw_diff) %>%
  pivot_longer(cols = c(raw_diff, adj_diff), names_to = "metric") %>%
  ggplot() +
  geom_col(aes(x = BreedingSeason, y = value, fill = metric), position = "dodge") +
  theme_minimal() +
  facet_wrap(~PopID, scales = "free")









## Scatterplot of fledge female vs fledge nest adjusted for all Pops
ggplot(bro_dat_sum_cor) +
  geom_point(aes(x = fledge_fem_avg, y = fledge_nest_avg_adj, color = fem_clutch_ratio)) +
  geom_text(aes(x = 2, y = 15, label = cor), data = . %>% distinct(PopID, .keep_all = T)) +
  geom_abline(slope = 1) +
  scale_x_continuous(limits= c(0,15)) +
  scale_y_continuous(limits= c(0,15)) +
  facet_wrap(~PopID, scales = "free") +
  theme_minimal() +
  labs(y = "Avg fledge per nest", x = "Avg fledge per female") +
  theme(legend.position = c(0.75,0.05),
        legend.direction = "horizontal")
ggsave(filename = "plots/Nest_vs_female_fledge_avg.png", width = 15, units = "in", height = 10)

## Timeseries of fledge nest adjusted for all Pops
ggplot(bro_dat_sum_cor) +
  geom_point(aes(x = BreedingSeason, y = fledge_nest_avg_adj, color = fledge_fem_avg)) +
  facet_wrap(~PopID, scales = "free") +
  theme_minimal() +
  labs(y = "Avg fledge per nest (adjusted)", x = "Year")


## Ratio of fledge per female and per nest as a function of the ratio of IDd females to breeding records
ggplot(bro_dat_sum_cor) +
  geom_point(aes(x = fem_clutch_ratio, y = fledge_fem_avg, size = females, color = PopID)) +
  stat_smooth(method = "gam", aes(x = fem_clutch_ratio, y = fledge_fem_avg)) +
  labs(x = "Ratio of females to number of breeding attempts",
       y = "Average fledglings per female",
       size = "# females") +
  theme_minimal()

ggplot(bro_dat_sum_cor) +
  geom_point(aes(x = females, y = fledge_fem_avg, size = females, color = PopID)) +
  stat_smooth(method = "gam", aes(x = females, y = fledge_fem_avg)) +
  facet_wrap(~PopID, scales = "free") +
  labs(x = "Number nests",
       y = "Average fledglings per female - Average fledglings per nest",
       size = "# females") +
  theme_minimal()


## Relationship between ratio of fledge metrics and number of second clutches
ggplot(bro_dat_sum_cor) +
  geom_point(aes(x = prop_second_clutch, y = fledge_ratio, size = females, color = PopID)) +
  stat_smooth(method = "lm", aes(x = prop_second_clutch, y = fledge_ratio)) +
  facet_wrap(~PopID, scales = "free") +
  labs(x = "Proportion second clutches",
       y = "Fledge ratio",
       size = "# females") +
  theme_minimal()



## Visualize populations without information on female ID
bro_dat_sum_cor %>%
  filter(fem_clutch_ratio < 0.5) %>%
  ggplot() +
  geom_point(aes(x = BreedingSeason, y = fledge_nest_avg), alpha = 0.6, color = "steelblue") +
  geom_point(aes(x = BreedingSeason, y = fledge_fem_avg), alpha = 0.6) +
  facet_wrap(~PopID, scales = "free_x") +
  labs(x = "Year", y = "Mean total chicks fledged per female") +
  theme_minimal()


## Correlation between total fledged by location and female depending on clutch type
bro_dat_sum %>%
  filter(!is.na(tot_fledged_fem) & !is.na(tot_fledged_nest)) %>%
  distinct(PopID, BreedingSeason,LocationID, ClutchType_calculated,.keep_all = T) %>%
  mutate(shape_type = as.factor(case_when(is.na(FemaleID) ~ 1,
                                          TRUE ~ 2))) %>%
  ggplot() +
  geom_jitter(aes(x = tot_fledged_fem, y = tot_fledged_loc, shape = shape_type), alpha = 0.6, color = "steelblue", size = 1) +
  facet_grid(ClutchType_calculated~.) +
  labs(x = "Total fledge (per female)", y = "Total fledge (per location)") +
  theme_minimal()
ggsave(filename = "plots/Fledge_counts_by_clutchtype_plot.png", width = 15, units = "in", height = 10)



## Visualize populations without information on female ID
ggplot(bro_dat_sum_cor) +
  geom_jitter(aes(x = BreedingSeason, y = fledge_nest_avg), alpha = 0.6, color = "steelblue", size = 1) +
  geom_jitter(aes(x = BreedingSeason, y = fledge_fem_avg), alpha = 0.6, size = 1) +
  facet_wrap(~PopID, scales = "free_x") +
  labs(x = "Year", y = "Mean total chicks fledged per female") +
  theme_minimal()

## Visualize number of females in the year relative to fledging numbers
bro_dat_sum %>%
  distinct(PopID, BreedingSeason, name, .keep_all = T) %>%
  ggplot() +
  geom_col(aes(x = BreedingSeason, y = value, fill = name), position = "dodge") +
  facet_wrap(~PopID, scales = "free") +
  labs(x = "Year", y = "Banded females") +
  theme_minimal()

## Visualize relationship between average fledging number and number of females IDd
bro_dat_sum %>%
  distinct(PopID, BreedingSeason, .keep_all = T) %>%
  filter(name == "counts_f") %>%
  ggplot() +
  geom_point(aes(x = tot_fledged, y = value)) +
  stat_smooth(method = "lm", aes(x = tot_fledged, y = value)) +
  facet_wrap(~PopID, scales = "free") +
  # labs(x = "Year", y = "Banded females") +
  theme_minimal()

## Plot average values of total fledged
bro_dat_sum %>%
  distinct(PopID, BreedingSeason, .keep_all = T) %>%
  ggplot() +
  geom_line(aes(x = BreedingSeason, y = pop_year_fledge_avg), alpha = 0.6) +
  geom_hline(aes(yintercept = pop_tot_fledge_avg), color = "red", alpha = 0.6) +
  # scale_color_gradient(low = "blue", high = "yellow") +
  facet_wrap(~PopID, scales = "free_x") +
  labs(x = "Year", y = "Mean total chicks fledged per female") +
  theme_minimal()

## Standardize
bro_dat_std <- bro_dat_sum %>%

  ## Calculate summary stats for each year
  group_by(PopID, BreedingSeason) %>%
  mutate(envStd_fledge_num = (pop_year_fledge_avg - pop_tot_fledge_avg) / pop_tot_fledge_sd) %>%
  distinct(PopID, BreedingSeason, .keep_all = T) %>%
  select(PopID, BreedingSeason, pop_year_fledge_avg, pop_tot_fledge_avg, pop_tot_fledge_sd,  envStd_fledge_num) %>%
  na.omit()

## Plot standardized value
bro_dat_std %>%
  ggplot() +
  geom_line(aes(x = BreedingSeason, y = envStd_fledge_num), alpha = 0.6) +
  facet_wrap(~PopID, scales = "free_x") +
  theme_minimal()

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


