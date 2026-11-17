# GSL dust costs: School loss days impacts
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

all_scenarios <- seq(4182, 4202, by = 1) 
current_scenario = 4192
relevant_scenarios <- c(4183, current_scenario, 4198, 4200) 

scenario_pal <- c(palette$bad, palette$current, palette$target, palette$avg)

n_years_storms = 6

# cost of a school loss day
SLD_24 = 1673.504


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Load processed data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scenario_pm_deltas <- read.csv("processed/scenario_pm_deltas_daily.csv", stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  filter(scenario %in% relevant_scenarios)

ct_school <- read.csv("processed/ct_school.csv", stringsAsFactors = FALSE)

ct_school_race <- read.csv("processed/ct_school_race.csv", stringsAsFactors = FALSE)

# Read the valuation rather than hard-coding it
SLD_24 <- read.csv("data/health/morbidity_valuations_2024.csv", stringsAsFactors = FALSE) %>%
  filter(Endpoint == "School Loss Days") %>%
  pull(COI_24)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Construct five-day exposure: lags 0--4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

school_loss_exposures <- scenario_pm_deltas %>%
  select(
    scenario,
    date,
    FIPS,
    pm10_delta,
    current_pm10_delta
  ) %>%
  distinct() %>%
  tidyr::crossing(
    outcome_day_offset = 0:4
  ) %>%
  mutate(
    # An exposure on date t contributes to outcomes on dates t, t+1, ..., t+4.
    outcome_date = date + outcome_day_offset,
    # Each exposure day contributes one-fifth of its
    # concentration to the five-day moving average.
    pm10_contribution = pm10_delta / 5,
    current_pm10_contribution =
      current_pm10_delta / 5
  ) %>%
  group_by(
    scenario,
    outcome_date,
    FIPS
  ) %>%
  summarise(
    pm10_delta =
      sum(pm10_contribution, na.rm = TRUE),
    
    current_pm10_delta =
      sum(current_pm10_contribution, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  rename(date = outcome_date)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### School-loss concentration-response parameters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

schoolloss_parameters <- readxl::read_excel(
  "data/health/morbidity_parameters.xlsx"
) %>%
  filter(
    endpoint == "School Loss Days",
    pollutant == "pm10"
  ) %>%
  mutate(
    CI_lower = ifelse(
      is.na(CI_lower),
      parameter_value - 1.96 * se,
      CI_lower
    ),
    CI_upper = ifelse(
      is.na(CI_upper),
      parameter_value + 1.96 * se,
      CI_upper
    ),
    beta = case_when(
      parameter == "beta" ~ parameter_value,
      parameter == "RR" ~ log(parameter_value) / dose,
      TRUE ~ NA_real_
    ),
    beta_lower = case_when(
      parameter == "beta" ~ CI_lower,
      parameter == "RR" ~ log(CI_lower) / dose,
      TRUE ~ NA_real_
    ),
    beta_upper = case_when(
      parameter == "beta" ~ CI_upper,
      parameter == "RR" ~ log(CI_upper) / dose,
      TRUE ~ NA_real_
    )
  )

beta_pm10 <- schoolloss_parameters$beta
beta_pm10_lower <- schoolloss_parameters$beta_lower
beta_pm10_upper <- schoolloss_parameters$beta_upper


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Health-impact functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

schoolloss_hif_loglinear <- function(
    beta,
    scenario_exposure,
    current_exposure,
    baseline_cases
) {
  
  baseline_cases * (
    exp(
      beta *
        (scenario_exposure - current_exposure)
    ) -
      exp(-beta * current_exposure)
  )
}


schoolloss_hif_linear <- function(
    beta,
    scenario_exposure,
    baseline_cases
) {
  
  baseline_cases * beta * scenario_exposure
}


calculate_schoolloss_impacts <- function(
    incidence_data,
    exposure_data
) {
  
  incidence_data %>%
    inner_join(
      exposure_data,
      by = "FIPS",
      relationship = "many-to-many"
    ) %>%
    mutate(
      endpoint = "School Loss Days",
      # Expected school-loss days on one calendar day
      baseline_cases_daily =
        incidence_rate_daily * pop,
      # Primary log-linear specification
      SLD = schoolloss_hif_loglinear(
        beta = beta_pm10,
        scenario_exposure = pm10_delta,
        current_exposure = current_pm10_delta,
        baseline_cases = baseline_cases_daily
      ) / n_years_storms,
      SLD_lower =
        schoolloss_hif_loglinear(
          beta = beta_pm10_lower,
          scenario_exposure = pm10_delta,
          current_exposure = current_pm10_delta,
          baseline_cases = baseline_cases_daily
        ) / n_years_storms,
      SLD_upper =
        schoolloss_hif_loglinear(
          beta = beta_pm10_upper,
          scenario_exposure = pm10_delta,
          current_exposure = current_pm10_delta,
          baseline_cases = baseline_cases_daily
        ) / n_years_storms,
      # Linear sensitivity specification
      SLD_linear = schoolloss_hif_linear(
        beta = beta_pm10,
        scenario_exposure = pm10_delta,
        baseline_cases = baseline_cases_daily
      ) / n_years_storms,
      SLD_linear_lower =
        schoolloss_hif_linear(
          beta = beta_pm10_lower,
          scenario_exposure = pm10_delta,
          baseline_cases = baseline_cases_daily
        ) / n_years_storms,
      SLD_linear_upper =
        schoolloss_hif_linear(
          beta = beta_pm10_upper,
          scenario_exposure = pm10_delta,
          baseline_cases = baseline_cases_daily
        ) / n_years_storms,
      pm_delta = pm10_delta
    )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Calculate total and race-specific impacts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct_schoolloss <- calculate_schoolloss_impacts(
  incidence_data = ct_school,
  exposure_data = school_loss_exposures
)

ct_schoolloss_race <- calculate_schoolloss_impacts(
  incidence_data = ct_school_race,
  exposure_data = school_loss_exposures
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Population denominators
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

school_population <- ct_school %>%
  distinct(FIPS, age_group, pop) %>%
  summarise(student_pop = sum(pop, na.rm = TRUE)) %>%
  pull(student_pop)

reported_races <- c(
  "AfAm/Black",
  "Asian",
  "Hispanic/Latino",
  "Pacific Islander",
  "White"
)

school_population_race <- ct_school_race %>%
  filter(Race %in% reported_races) %>%
  distinct(FIPS, Race, age_group, pop) %>%
  group_by(Race) %>%
  summarise(
    student_pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Census-tract results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct_schoolloss_map <- ct_schoolloss %>%
  group_by(scenario, FIPS, County, endpoint) %>%
  summarise(
    SLD = sum(SLD, na.rm = TRUE),
    SLD_lower = sum(SLD_lower, na.rm = TRUE),
    SLD_upper = sum(SLD_upper, na.rm = TRUE),
    SLD_linear = sum(SLD_linear, na.rm = TRUE),
    SLD_linear_lower = sum(SLD_linear_lower, na.rm = TRUE),
    SLD_linear_upper = sum(SLD_linear_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    costs_SLD = SLD * SLD_24,
    costs_SLD_lower = SLD_lower * SLD_24,
    costs_SLD_upper = SLD_upper * SLD_24,
    costs_SLD_linear = SLD_linear * SLD_24,
    costs_SLD_linear_lower = SLD_linear_lower * SLD_24,
    costs_SLD_linear_upper = SLD_linear_upper * SLD_24
  )

write.csv(ct_schoolloss_map, "processed/ct_schoolloss_map.csv", row.names = FALSE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Total impacts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_schoolloss <- ct_schoolloss %>%
  group_by(scenario, endpoint) %>%
  summarise(
    SLD = sum(SLD, na.rm = TRUE),
    SLD_lower = sum(SLD_lower, na.rm = TRUE),
    SLD_upper = sum(SLD_upper, na.rm = TRUE),
    SLD_linear = sum(SLD_linear, na.rm = TRUE),
    SLD_linear_lower = sum(SLD_linear_lower, na.rm = TRUE),
    SLD_linear_upper = sum(SLD_linear_upper, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    student_pop = school_population,
    SLD_per_100k = 100000 * SLD / student_pop,
    costs_SLD = SLD * SLD_24,
    costs_SLD_lower = SLD_lower * SLD_24,
    costs_SLD_upper = SLD_upper * SLD_24,
    costs_SLD_linear = SLD_linear * SLD_24,
    costs_SLD_linear_lower = SLD_linear_lower * SLD_24,
    costs_SLD_linear_upper = SLD_linear_upper * SLD_24)

write.csv(total_schoolloss, "processed/total_schoolloss.csv", row.names = FALSE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Race-specific totals
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_schoolloss_race <- ct_schoolloss_race %>%
  filter(Race %in% reported_races) %>%
  group_by(scenario, endpoint, Race) %>%
  summarise(
    SLD = sum(SLD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    school_population_race,
    by = "Race"
  ) %>%
  mutate(
    SLD_per_100k =
      100000 * SLD / student_pop,
    costs_SLD = SLD * SLD_24)
write.csv(total_schoolloss_race, "processed/total_schoolloss_race.csv", row.names = FALSE)


schoolloss_race <- total_schoolloss_race %>%
  mutate(SLD_per_100k = SLD/student_pop*100000) %>%
  ggplot(aes(x=scenario, color = reorder(Race, -SLD_per_100k), y=SLD_per_100k))+
  geom_line(linewidth = 1)+
  geom_point(data = total_schoolloss_race %>% filter(scenario %in% relevant_scenarios)
             , size = 2.5)+
  # ggtitle("Distribution of mortality risk across race")+
  scale_y_continuous(name = "School Loss Days (per 100k students)")+
  #scale_x_reverse(
  scale_x_continuous(
    breaks = relevant_scenarios,
    name = "GSL water level (ftASL)")+
  scale_color_manual(values = c(palette$red, palette$green, palette$orange, palette$blue, palette$dark),
                     #name = "Race",
                     labels = c("Hawaiian/Pacific Islander", "Black/African American", "Hispanic/Latino", "White Non-hispanic", "Asian")
  )+
  ggtitle("School Loss Days by race")+
  theme_cowplot(16)+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.justification = "center"
  )
schoolloss_race
ggsave("figs_pdf/schoolloss_by_race.pdf", 
       width = 8.5, height = 5)
