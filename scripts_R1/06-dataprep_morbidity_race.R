# GSL dust costs: morbidity incidence by race data prep
# albert.garcia@utah.edu
# created: 06/20/2025
# updated:

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Set up environment
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,
  stringr,
  readxl
)

options(scipen = 999)
options(dplyr.summarise.inform = FALSE)

`%ni%` <- Negate(`%in%`)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Census-tract population by age and race
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

agebyrace_sheet <- 11

ct_ACS_agebyrace_wide <- readxl::read_excel(
  "data/population/Demographic_Raw_Tables.xlsx",
  sheet = agebyrace_sheet
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Prepare race-by-age populations for morbidity
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# For morbidity, retain separate groups for ages 15--17,
# 18--19, and 20--24 because work-loss days apply only
# to adults ages 18--64.

ct_ACS_agebyrace_morbidity <- ct_ACS_agebyrace_wide %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = "agebyrace_group",
    values_to = "pop"
  ) %>%
  filter(!str_detect(agebyrace_group, "Total")) %>%
  mutate(
    agebyrace_group = str_remove_all(
      agebyrace_group,
      "[()]"
    ),
    
    race = case_when(
      str_detect(
        agebyrace_group,
        "American Indian and Alaska Native"
      ) ~ "American Indian and Alaska Native Alone",
      
      str_detect(
        agebyrace_group,
        "White"
      ) ~ "White NH",
      
      str_detect(
        agebyrace_group,
        "Hispanic"
      ) ~ "Hispanic or Latino",
      
      str_detect(
        agebyrace_group,
        "Black"
      ) ~ "Black or African American Alone",
      
      str_detect(
        agebyrace_group,
        "Asian"
      ) ~ "Asian alone",
      
      str_detect(
        agebyrace_group,
        "Hawaiian and PI"
      ) ~ "Hawaiian and PI Alone",
      
      str_detect(
        agebyrace_group,
        "Other Race"
      ) ~ "Other Race Alone",
      
      str_detect(
        agebyrace_group,
        "Multi Racial"
      ) ~ "Multi Racial",
      
      TRUE ~ NA_character_
    ),
    
    age_group = str_remove(
      agebyrace_group,
      paste0(race, " ")
    ),
    
    # Standardize this label to match the non-race table
    age_group = if_else(
      age_group == "18 and 19 years",
      "18 to 19 years",
      age_group
    ),
    
    # Combine these groups to match the non-race ACS table
    age_group = case_when(
      age_group %in% c(
        "25 to 29 years",
        "30 to 34 years"
      ) ~ "25 to 34 years",
      
      TRUE ~ age_group
    )
  ) %>%
  filter(!is.na(race)) %>%
  group_by(
    FIPS,
    County,
    age_group,
    race
  ) %>%
  summarise(
    pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lower_age = str_sub(age_group, 1, 2),
    
    lower_age = as.numeric(
      ifelse(
        lower_age == "Un",
        0,
        lower_age
      )
    ),
    
    age_group = str_remove(
      age_group,
      " years"
    ),
    
    upper_age = str_sub(age_group, -2, -1),
    
    upper_age = as.numeric(
      ifelse(
        upper_age == "er",
        99,
        upper_age
      )
    ),
    
    # "Under 5" corresponds to ages 0--4
    upper_age = ifelse(
      upper_age == 5,
      4,
      upper_age
    )
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Validate race population data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

missing_population_data <- ct_ACS_agebyrace_morbidity %>%
  filter(
    is.na(FIPS) |
      is.na(County) |
      is.na(race) |
      is.na(age_group) |
      is.na(lower_age) |
      is.na(upper_age) |
      is.na(pop)
  )

if (nrow(missing_population_data) > 0) {
  print(missing_population_data)
  
  stop(
    "Missing tract, county, race, age, or population data."
  )
}

duplicate_population_records <-
  ct_ACS_agebyrace_morbidity %>%
  count(
    FIPS,
    County,
    race,
    lower_age,
    upper_age
  ) %>%
  filter(n > 1)

if (nrow(duplicate_population_records) > 0) {
  print(duplicate_population_records)
  
  stop(
    "Duplicate race-specific population records found."
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Load non-race morbidity incidence and parameters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This file is produced by 05-dataprep_morbidity.R.
# It already contains the daily incidence rates and converted
# PM10 beta parameters.

ct_incidence_morbidity_base <- read.csv(
  "processed/ct_incidence_morbidity.csv",
  stringsAsFactors = FALSE
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Align non-race and race-specific age intervals
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The race-specific demographic table contains some broader
# age groups than the non-race table. For example, the race
# data contain ages 55--64, while the non-race data contain
# ages 55--59 and 60--64.

race_age_definitions <- ct_ACS_agebyrace_morbidity %>%
  distinct(
    FIPS,
    County,
    lower_age,
    upper_age,
    age_group
  ) %>%
  rename(
    race_lower_age = lower_age,
    race_upper_age = upper_age,
    race_age_group = age_group
  )


# Match each non-race age interval to the race-specific
# age interval that contains it.

base_to_race_age <- ct_incidence_morbidity_base %>%
  inner_join(
    race_age_definitions,
    by = c("FIPS", "County"),
    relationship = "many-to-many"
  ) %>%
  filter(
    lower_age >= race_lower_age,
    upper_age <= race_upper_age
  )


# Show any age intervals that required aggregation
base_to_race_age %>%
  filter(
    lower_age != race_lower_age |
      upper_age != race_upper_age
  ) %>%
  distinct(
    lower_age,
    upper_age,
    race_lower_age,
    race_upper_age
  ) %>%
  arrange(
    race_lower_age,
    lower_age
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Population-weighted mean helper
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weighted_mean_safe <- function(x, weights) {
  
  valid <- !is.na(x) & !is.na(weights)
  
  if (!any(valid)) {
    return(NA_real_)
  }
  
  if (sum(weights[valid]) == 0) {
    return(mean(x[valid]))
  }
  
  weighted.mean(
    x = x[valid],
    w = weights[valid]
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Aggregate incidence and beta parameters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Where a race-specific age interval contains multiple
# non-race intervals, calculate population-weighted incidence
# and beta parameters.

morbidity_parameters_by_race_age <-
  base_to_race_age %>%
  group_by(
    FIPS,
    County,
    race_lower_age,
    race_upper_age,
    race_age_group,
    endpoint
  ) %>%
  summarise(
    value = weighted_mean_safe(
      value,
      pop
    ),
    
    beta_pm10 = weighted_mean_safe(
      beta_pm10,
      pop
    ),
    
    beta_pm10_lower = weighted_mean_safe(
      beta_pm10_lower,
      pop
    ),
    
    beta_pm10_upper = weighted_mean_safe(
      beta_pm10_upper,
      pop
    ),
    
    .groups = "drop"
  ) %>%
  rename(
    lower_age = race_lower_age,
    upper_age = race_upper_age,
    age_group = race_age_group
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Join parameters to race-specific populations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct_incidence_morbidity_race <-
  ct_ACS_agebyrace_morbidity %>%
  inner_join(
    morbidity_parameters_by_race_age,
    by = c(
      "FIPS",
      "County",
      "lower_age",
      "upper_age",
      "age_group"
    ),
    relationship = "many-to-many"
  ) %>%
  select(
    FIPS,
    County,
    race,
    age_group,
    lower_age,
    upper_age,
    pop,
    endpoint,
    value,
    beta_pm10,
    beta_pm10_lower,
    beta_pm10_upper
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Validate final morbidity-by-race data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# There should be one observation per
# tract-race-age-endpoint combination.

duplicate_records <- ct_incidence_morbidity_race %>%
  count(
    FIPS,
    race,
    age_group,
    endpoint
  ) %>%
  filter(n > 1)

if (nrow(duplicate_records) > 0) {
  print(duplicate_records)
  
  stop(
    "Duplicate tract-race-age-endpoint records found."
  )
}


# Check for missing analysis inputs

missing_records <- ct_incidence_morbidity_race %>%
  filter(
    is.na(FIPS) |
      is.na(County) |
      is.na(race) |
      is.na(pop) |
      is.na(value) |
      is.na(beta_pm10) |
      is.na(beta_pm10_lower) |
      is.na(beta_pm10_upper)
  )

if (nrow(missing_records) > 0) {
  print(
    missing_records %>%
      count(
        endpoint,
        name = "n_missing"
      )
  )
  
  stop(
    "Missing race population, incidence, or beta values found."
  )
}


# Confirm that all intended endpoints are represented

expected_endpoints <- c(
  "HA, All Respiratory",
  "HA, All Cardiac Outcomes",
  "HA, Asthma",
  "ER visits, All Respiratory",
  "ER visits, All Cardiac Outcomes",
  "Work Loss Days"
)

missing_endpoints <- setdiff(
  expected_endpoints,
  unique(ct_incidence_morbidity_race$endpoint)
)

if (length(missing_endpoints) > 0) {
  stop(
    "Missing race morbidity endpoints: ",
    paste(
      missing_endpoints,
      collapse = ", "
    )
  )
}


# Review final record counts

ct_incidence_morbidity_race %>%
  count(
    endpoint,
    race,
    name = "n_records"
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Write output
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write.csv(
  ct_incidence_morbidity_race,
  file = "processed/ct_incidence_morbidity_race.csv",
  row.names = FALSE
)
