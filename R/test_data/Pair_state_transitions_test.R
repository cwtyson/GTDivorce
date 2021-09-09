## Test state transition script with mock data

library(tidyverse)
library(here)

pair_test_dat <- readxl::read_xlsx("R/data/pair_state_test_data.xlsx",
                                   sheet = "sheet1")

########## Females
pair_dat_complete_test <- pair_test_dat %>%

  filter(pop %in% c("HOG")) %>%
  mutate(captured = ifelse(!is.na(female), TRUE, NA)) %>%
  group_by(pop) %>%
  tidyr::complete(female, year) %>%
  mutate(captured = ifelse(is.na(captured), FALSE, captured)) %>%

  ## Get previous partner for the female
  arrange(pop, female, year, order) %>%
  group_by(pop, female) %>%
  mutate(multi_partner = case_when(any(f_pair_num == 2) ~ 1,
                                   TRUE ~ 0),
         prev_partner = lag(male)) %>%

  ## In the current year, is the previous partner breeding?
  group_by(pop, year) %>%
  mutate(prev_partner_breeding = case_when((prev_partner %in% male ) & (!is.na(prev_partner)) ~ TRUE,
                                           !(prev_partner %in% male ) & (!is.na(prev_partner)) ~ FALSE)) %>%

  ## In the current year, is the previous partner breeding with a known individual?
  group_by(pop, year) %>%
  mutate(partner_known = case_when(prev_partner %in% male[captured == TRUE] & !is.na(prev_partner)  ~ TRUE,
                                   !(prev_partner %in% male[captured == TRUE] ) & !is.na(prev_partner) ~ FALSE)) %>%

  arrange(pop, female, year, order) %>%
  group_by(pop, female) %>%
  mutate(caps = sum(captured, na.rm = T),
         event_code = case_when(

           ## 0: not captured
           captured == FALSE &
             (is.na(lag(male)) | is.na(male)) &
             (isTRUE(prev_partner_breeding) | is.na(prev_partner_breeding)) ~ "0",

           ## 1: Focal captured, partner captured, same partner as previous partner
           captured == "TRUE" &
             (male == lag(male) & !is.na(male)) ~ "1",

           ## 2: Focal captured, current partner captured, different from previous partner
           captured == "TRUE" &
             !is.na(male) &
             (male != lag(male)) ~ "2",

           ## 3: Focal captured, partner captured, relationship with previous partner unknown
           captured == "TRUE" &
             is.na(lag(male)) &
             !is.na(male) ~ "3",

           ## 4: Focal captured, current partner not captured, previous partner breeding elsewhere in current year
           captured == "TRUE" &
             is.na(male) &
             !is.na(lag(male)) &
             prev_partner_breeding == TRUE~ "4",

           ## 5: Focal captured, current partner not captured, previous partner not captured in current year or previous partner not known
           captured == "TRUE" &
             is.na(male) &
             (prev_partner_breeding == FALSE | is.na(prev_partner)) ~ "5",

           ## 6: Focal not captured, previous partner captured with another bird in current year
           captured == FALSE &
             prev_partner_breeding == TRUE &
             partner_known == TRUE ~ "6",

           ## Otherwise, 0
           TRUE ~ "0")) %>%

  filter(!is.na(female)) %>%
  arrange(pop, female, year, order) %>%

  ## Keep only the first record from each year
  distinct(pop, female, year, .keep_all = T) %>%
  select(pop, female, year, multi_partner, event_code)


## Pivot wider
pair_dat_test_wide <- pair_dat_complete_test %>%

  pivot_wider(id_cols = c(female, multi_partner),
              names_from = year,
              values_from = event_code)

## Save
write_csv(x = pair_dat_test_wide,
          file = paste0(here::here(),
                        paste("/R/data/test_data_CR_females.csv")))

############## Males

pair_dat_complete_test <- pair_test_dat %>%
  filter(pop %in% c("HOG")) %>%
  mutate(captured = ifelse(!is.na(male), TRUE, NA)) %>%
  group_by(pop) %>%
  tidyr::complete(male, year) %>%
  mutate(captured = ifelse(is.na(captured), FALSE, captured)) %>%

  ## Get previous partner for the male
  arrange(pop, male, year) %>%
  group_by(pop, male) %>%
  mutate(prev_partner = lag(female)) %>%

  ## In the current year, is the previous partner breeding?
  group_by(pop, year) %>%
  mutate(prev_partner_breeding = case_when((prev_partner %in% female ) & (!is.na(prev_partner)) ~ TRUE,
                                           !(prev_partner %in% female ) & (!is.na(prev_partner)) ~ FALSE)) %>%

  ## In the current year, is the previous partner breeding with a known individual?
  arrange(pop, prev_partner, year) %>%
  group_by(pop, prev_partner) %>%
  mutate(prev_partner_known_pair = !is.na(male)) %>%

  ## In the current year, is the previous partner breeding with a known individual?
  group_by(pop, year) %>%
  mutate(partner_known = case_when(prev_partner %in% female[captured == TRUE] & !is.na(prev_partner)  ~ TRUE,
                                   !(prev_partner %in% female[captured == TRUE] ) & !is.na(prev_partner) ~ FALSE)) %>%

  ## Was the female breeding in the previous year?
  group_by(pop, female) %>%
  arrange(pop, female) %>%
  mutate(breeding_previous_year = case_when(!is.na(lag(year)) & !is.na(female) ~ TRUE,
                                            is.na(lag(year)) & !is.na(female)~ FALSE,
                                            is.na(female) ~ NA)) %>%

  ## Does the female have the same partner as the previous year?
  group_by(pop, female, year) %>%
  mutate(same_female_partner = case_when(((male) == lag(male)) & !is.na(male) & !is.na(lag(female)) ~ TRUE,
                                         ((male) != lag(male)) & !is.na(male) & !is.na(lag(female))  ~ FALSE,
                                         TRUE ~ NA)) %>%

  arrange(pop, male, year) %>%
  group_by(pop, male) %>%
  mutate(caps = sum(captured, na.rm = T),
         # first_year = min(year[captured == TRUE], na.rm = T),
         # breeding_period = (max(year[captured == TRUE], na.rm = T) - min(year[captured == TRUE], na.rm = T)) + 1,
         # breaks = sum(captured, na.rm = T) < breeding_period,
         event_code = case_when(

           ## 0: not captured
           captured == FALSE &
             (is.na(lag(female)) | is.na(female)) &
             (isTRUE(prev_partner_breeding) | is.na(prev_partner_breeding)) ~ "0",

           ## 1: Focal captured, partner captured, same partner as previous partner
           captured == "TRUE" &
             (female == lag(female) & !is.na(female)) ~ "1",

           ## 2: Focal captured, current partner captured, different from previous partner
           captured == "TRUE" &
             !is.na(female) &
             (female != lag(female)) ~ "2",

           ## 3: Focal captured, partner captured, relationship with previous partner unknown
           captured == "TRUE" &
             is.na(lag(female)) &
             !is.na(female) ~ "3",

           ## 4: Focal captured, current partner not captured, previous partner breeding elsewhere in current year
           captured == "TRUE" &
             is.na(female) &
             !is.na(lag(female)) &
             prev_partner_breeding == TRUE~ "4",

           ## 5: Focal captured, current partner not captured, previous partner not captured in current year or previous partner not known
           captured == "TRUE" &
             is.na(female) &
             (prev_partner_breeding == FALSE | is.na(prev_partner)) ~ "5",

           ## 6: Focal not captured, previous partner captured with another bird
           captured == FALSE &
             prev_partner_breeding == TRUE  &
             partner_known == TRUE ~ "6",

           ## Otherwise, 0
           TRUE ~ "0")) %>%

  arrange(male, year) %>%
  filter(!is.na(male)) %>%
  select(pop, male, year, event_code)


## Pivot wider
pair_dat_test_wide <- pair_dat_complete_test %>%
  pivot_wider(id_cols = c(male),
              names_from = year,
              values_from = event_code)

## Save
write_csv(x = pair_dat_test_wide,
          file = paste0(here::here(),
                        paste("/R/data/test_data_CR_males.csv")))

