## Create CR for females in

library(tidyverse)
library(here)
library(glue)


## Read in pair data
pair_dat <- read_csv(paste0(here(),"/data/pop_pairs/all_pop_pairs.csv"))

## For each population, create CR for females, and save as a separate .csv
for(pop_code in unique(pair_dat$pop)){


  ## Get CR for each individual in selected pop
  pair_dat_complete <- pair_dat %>%
    filter(pop %in% pop_code) %>%
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
  pair_dat_wide <- pair_dat_complete %>%
    pivot_wider(id_cols = c(pop, female, multi_partner),
                names_from = year,
                values_from = event_code)

  # print(table(pair_dat_complete$event_code))

  ## Save
  write_csv(x = pair_dat_wide,
            file = paste0(here::here(),
                          paste("/data/event_codes/f/", pop_code, "_female_event_codes.csv", sep = "")))

}
