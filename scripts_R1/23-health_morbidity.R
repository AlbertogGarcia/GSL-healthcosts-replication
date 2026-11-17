# GSL dust costs: asthma impacts
# albert.garcia@utah.edu
# created: 05/28/2025
# updated: 

# Set up environment ########################################

# load or install necessary libraries. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(mapview,  # view spatial data in viewer
               measurements,  # convert units easily
               progress,  # progress bar
               scales, # add commas long numbers
               sf,  # shapefile read
               tidyverse, # tidyverse
               stringr,
               fuzzyjoin,
               ggplot2,
               readxl,
               cowplot
)

options(scipen=999)  # turn off sci notation
options(dplyr.summarise.inform = FALSE)  # turn off dplyr group by comments
options(java.parameters = "-Xmx8000m") 
`%ni%` <- Negate(`%in%`)  # "not in" function

# Color palette
palette <- list("white" = "#FAFAFA",
                "dark" = "#0c2230",
                "red" = "#d7191c",
                "blue" = "#2c7bb6",
                "orange" = "#fc8d62",
                "green" = "#66c2a5",
                "purple" = "#8da0cb",
                "bad" = "#d7191c",
                "current" = "#fdae61",
                "target" = "#abd9e9",
                "avg" = "#2c7bb6"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### set base parameters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_scenarios <- seq(4182, 4203, by = 1) 
current_scenario = 4192
relevant_scenarios <- c(4183, current_scenario, 4198, 4200) 

scenario_pal <- c(palette$bad, palette$current, palette$target, palette$avg)

n_years_storms = 6

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Load and merge processed data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#1 Emissions scenarios by water-level
scenario_pm_deltas <- read.csv("processed/scenario_pm_deltas_daily.csv", stringsAsFactors =  FALSE)%>%
  mutate(date = as.Date(date)) %>%
  filter(scenario %in% relevant_scenarios)

#2 Population and incidence
ct_incidence_morbidity <- read.csv("processed/ct_incidence_morbidity.csv", stringsAsFactors =  FALSE)


morbidity_valuations_2024 <- read.csv("data/health/morbidity_valuations_2024.csv",stringsAsFactors = FALSE) %>%
  select(Endpoint, COI_24)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Construct endpoint-specific exposure windows
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_morbidity_exposure <- function(
    daily_data,
    endpoints,
    outcome_day_offsets
) {
  
  window_length <- length(outcome_day_offsets)
  
  daily_data %>%
    select(
      scenario,
      date,
      centroid_name,
      FIPS,
      pm10_delta,
      current_pm10_delta
    ) %>%
    tidyr::crossing(
      endpoint = endpoints,
      outcome_day_offset = outcome_day_offsets
    ) %>%
    mutate(
      # Date on which the health outcome is realized
      outcome_date = date + outcome_day_offset,
      # Each exposure day contributes 1/k of its concentration
      # to a k-day moving average
      pm10_contribution = pm10_delta / window_length,
      current_pm10_contribution =
        current_pm10_delta / window_length
    ) %>%
    group_by(
      scenario,
      outcome_date,
      centroid_name,
      FIPS,
      endpoint
    ) %>%
    summarise(
      pm10_delta = sum(pm10_contribution),
      current_pm10_delta =
        sum(current_pm10_contribution),
      .groups = "drop"
    ) %>%
    mutate(
      relative_pm10_delta =
        pm10_delta - current_pm10_delta
    ) %>%
    rename(date = outcome_date)
}

# Same-day daily exposure
same_day_exposures <- make_morbidity_exposure(
  daily_data = scenario_pm_deltas,
  endpoints = c(
    "HA, All Respiratory",
    "HA, All Cardiac Outcomes",
    "HA, Asthma"
  ),
  outcome_day_offsets = 0
)

# Three-day moving average
ed_visit_exposures <- make_morbidity_exposure(
  daily_data = scenario_pm_deltas,
  endpoints = c(
    "ER visits, All Respiratory",
    "ER visits, All Cardiac Outcomes"
  ),
  outcome_day_offsets = 0:2
)

# Seven-day moving average
work_loss_exposures <- make_morbidity_exposure(
  daily_data = scenario_pm_deltas,
  endpoints = "Work Loss Days",
  outcome_day_offsets = 0:6
)

scenario_pm_morbidity <- bind_rows(
  same_day_exposures,
  ed_visit_exposures,
  work_loss_exposures
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Health-impact functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Total lakebed-attributable burden under scenario s,
# using current observed incidence as the baseline incidence
morbidity_hif_loglinear <- function(
    beta,
    scenario_exposure,
    current_exposure,
    baseline_cases
) {
  
  baseline_cases * (
    exp(beta * (scenario_exposure - current_exposure)) -
      exp(-beta * current_exposure)
  )
}

calculate_morbidity_impacts <- function(
    incidence_data,
    exposure_data
) {
  
  incidence_data %>%
    inner_join(
      exposure_data,
      by = c("FIPS", "endpoint"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      # Expected cases on one outcome day
      baseline_cases_daily = value * pop,
      
      # Log-linear primary estimate
      morbidity = morbidity_hif_loglinear(
        beta = beta_pm10,
        scenario_exposure = pm10_delta,
        current_exposure = current_pm10_delta,
        baseline_cases = baseline_cases_daily
      ) / n_years_storms,
      
      morbidity_at_beta_lower = morbidity_hif_loglinear(
        beta = beta_pm10_lower,
        scenario_exposure = pm10_delta,
        current_exposure = current_pm10_delta,
        baseline_cases = baseline_cases_daily
      ) / n_years_storms,
      
      morbidity_at_beta_upper = morbidity_hif_loglinear(
        beta = beta_pm10_upper,
        scenario_exposure = pm10_delta,
        current_exposure = current_pm10_delta,
        baseline_cases = baseline_cases_daily
      ) / n_years_storms,
      
      # The ordering can reverse for confidence intervals that cross zero
      morbidity_lower = pmin(
        morbidity,
        morbidity_at_beta_lower,
        morbidity_at_beta_upper
      ),
      
      morbidity_upper = pmax(
        morbidity,
        morbidity_at_beta_lower,
        morbidity_at_beta_upper
      ),
      pm_delta = pm10_delta
    )
}


ct_morbidity_age <- calculate_morbidity_impacts(
  incidence_data = ct_incidence_morbidity,
  exposure_data = scenario_pm_morbidity
) %>%
  select(
    FIPS,
    County,
    scenario,
    date,
    age_group,
    pop,
    pm_delta,
    endpoint,
    morbidity,
    morbidity_lower,
    morbidity_upper
  )



ct_morbidity <- ct_morbidity_age %>%
  group_by(FIPS, County, scenario, endpoint) %>%
  summarise(
    morbidity = sum(morbidity, na.rm = TRUE),
    morbidity_lower = sum(morbidity_lower, na.rm = TRUE),
    morbidity_upper = sum(morbidity_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    morbidity_valuations_2024,
    by = c("endpoint" = "Endpoint")
  ) %>%
  mutate(
    costs = morbidity * COI_24,
    costs_lower = morbidity_lower * COI_24,
    costs_upper = morbidity_upper * COI_24
  )

write.csv(ct_morbidity, file = "processed/ct_morbidity.csv", row.names = FALSE)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Totals morbidity
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_morbidity <- ct_morbidity_age %>%
  group_by(scenario, endpoint) %>%
  summarise(
    morbidity = sum(morbidity, na.rm = TRUE),
    morbidity_lower = sum(morbidity_lower, na.rm = TRUE),
    morbidity_upper = sum(morbidity_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    morbidity_valuations_2024,
    by = c("endpoint" = "Endpoint")
  ) %>%
  mutate(
    costs = morbidity * COI_24,
    costs_lower = morbidity_lower * COI_24,
    costs_upper = morbidity_upper * COI_24
  )

write.csv(total_morbidity, file = "processed/total_morbidity.csv", row.names = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# By race
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Keep only objects needed for the race-stratified analysis
objects_to_keep <- c(
  "scenario_pm_morbidity",
  "morbidity_valuations_2024",
  "n_years_storms",
  "morbidity_hif_loglinear",
  "calculate_morbidity_impacts",
  "%ni%"
)

rm(list = setdiff(ls(envir = .GlobalEnv), objects_to_keep))

gc()

#2 Population and incidence
ct_incidence_morbidity_race <- read.csv("processed/ct_incidence_morbidity_race.csv", stringsAsFactors =  FALSE)

ct_incidence_morbidity_race <- read.csv(
  "processed/ct_incidence_morbidity_race.csv",
  stringsAsFactors = FALSE
)

ct_morbidity_agebyrace <- calculate_morbidity_impacts(
  incidence_data = ct_incidence_morbidity_race,
  exposure_data = scenario_pm_morbidity
) %>%
  select(
    FIPS,
    County,
    scenario,
    date,
    race,
    age_group,
    pop,
    pm_delta,
    endpoint,
    morbidity,
    morbidity_lower,
    morbidity_upper
  )


ct_morbidity_race <- ct_morbidity_agebyrace %>%
  group_by(FIPS, County, scenario, race, endpoint) %>%
  summarise(morbidity = sum(morbidity, na.rm = T))%>%
  ungroup %>%
  left_join(morbidity_valuations_2024, by = c("endpoint" = "Endpoint")) %>%
  mutate(costs = morbidity*COI_24)

write.csv(ct_morbidity_race, file = "processed/ct_morbidity_race.csv", row.names = FALSE)


total_morbidity_race <- ct_morbidity_agebyrace %>%
  drop_na(race) %>%
  filter(race %ni% c("Other Race", "Multiple Race", "American Indian"))%>%
  group_by(scenario, endpoint, race) %>%
  summarise(morbidity = sum(morbidity, na.rm = T),
            morbidity_lower = sum(morbidity_lower, na.rm = T),
            morbidity_upper = sum(morbidity_upper, na.rm = T))%>%
  ungroup %>%
  left_join(morbidity_valuations_2024, by = c("endpoint" = "Endpoint")) %>%
  mutate(costs = morbidity*COI_24,
         costs_lower = morbidity_lower*COI_24,
         costs_upper = morbidity_upper*COI_24)
write.csv(total_morbidity_race, file = "processed/total_morbidity_race.csv", row.names = FALSE)

ct_morbidity_map <- ct_morbidity_agebyrace %>%
  group_by(FIPS, County, scenario, endpoint) %>%
  summarise(morbidity = sum(morbidity, na.rm = T))%>%
  ungroup %>%
  left_join(morbidity_valuations_2024, by = c("endpoint" = "Endpoint")) %>%
  mutate(costs = morbidity*COI_24)

write.csv(ct_morbidity_map, file = "processed/ct_morbidity_map.csv", row.names = FALSE)
