# GSL dust costs: Morbidity impacts - projected
# albert.garcia@utah.edu
# created: 01/5/2026
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Load and merge processed data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#1 Emissions scenarios by water-level
scenario_pm_deltas <- read.csv("processed/scenario_pm_deltas_daily.csv", stringsAsFactors = FALSE) %>%
  filter(scenario %in% relevant_scenarios) %>%
  transmute(
    FIPS = as.character(FIPS),
    centroid_name = centroid_name,
    scenario,
    date = as.Date(date),
    event_days,
    pm10_delta,
    relative_pm10_delta
  )

current_exposure <- scenario_pm_deltas %>%
  filter(scenario == current_scenario) %>%
  select(centroid_name,date, FIPS,
    current_pm10_delta = pm10_delta)

scenario_exposure <- scenario_pm_deltas %>%
  left_join(current_exposure, by = c("centroid_name", "FIPS", "date")) %>%
  mutate(incremental_pm10_delta = if_else(scenario == current_scenario, 0, relative_pm10_delta)) %>%
  select(FIPS, scenario, date, event_days, incremental_pm10_delta, current_pm10_delta)

rm(scenario_pm_deltas, current_exposure)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Project morbidity one endpoint at a time
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

morbidity_valuations_projected <- read.csv("processed/morbidity_valuations_projected.csv", stringsAsFactors = FALSE) %>%
  filter(Year >= 2025)

endpoint_file_names <- read.csv("processed/morbidity_projection_endpoint_files.csv", stringsAsFactors = FALSE)

output_folder <- "processed/morbidity_projection_results"

dir.create(
  output_folder,
  showWarnings = FALSE,
  recursive = TRUE
)

endpoint_result_files <- character(
  nrow(endpoint_file_names)
)

for (i in seq_len(nrow(endpoint_file_names))) {
  
  e <- endpoint_file_names$endpoint[i]
  endpoint_file <- endpoint_file_names$endpoint_file[i]
  
  message(
    "Processing endpoint ", i,
    " of ",
    nrow(endpoint_file_names), ": ", e
  )
  
  #------------------------------------------------------------
  # Read only this endpoint's incidence and population data
  #------------------------------------------------------------
  
  endpoint_incidence <- readRDS(
    file.path(
      "processed/morbidity_projection_endpoints",
      paste0(endpoint_file, ".rds")
    )
  ) %>%
    filter(Year >= 2025) %>%
    mutate(
      FIPS = as.character(FIPS)
    ) %>%
    select(FIPS, County, Year, age_group, lower_age, upper_age, endpoint, pop, value, beta_pm10)
  
  if (nrow(endpoint_incidence) == 0) {
    stop("No incidence-projection observations found for endpoint: ", e)
  }
  
  missing_incidence_inputs <- endpoint_incidence %>%
    filter(
      is.na(FIPS) |
        is.na(Year) |
        is.na(pop) |
        is.na(value) |
        is.na(beta_pm10)
    )
  
  if (nrow(missing_incidence_inputs) > 0) {
    print(
      missing_incidence_inputs %>%
        count(
          is.na(pop),
          is.na(value),
          is.na(beta_pm10),
          name = "n"
        )
    )
    
    stop("Missing incidence inputs found for endpoint: ", e)
  }
  #------------------------------------------------------------
  # Calculate event-response totals for each distinct beta
  #
  # This collapses the date dimension before joining to every
  # age-by-year population observation.
  #------------------------------------------------------------
  
  endpoint_betas <- endpoint_incidence %>%
    distinct(beta_pm10) %>%
    filter(!is.na(beta_pm10))
  
  endpoint_exposure_response <- tidyr::crossing(
    scenario_exposure,
    endpoint_betas
  ) %>%
    mutate(
      current_response = 1 - exp(-beta_pm10 * current_pm10_delta),
      incremental_response = 1 - exp(-beta_pm10 * incremental_pm10_delta),
      combined_response = pmax(
        current_response + incremental_response,
        0
      ),
      event_response = combined_response * event_days
    ) %>%
    group_by(FIPS, scenario, beta_pm10) %>%
    summarise(
      annualized_response =
        sum(event_response, na.rm = TRUE) /
        n_years_storms,
      .groups = "drop"
    )
  
  #------------------------------------------------------------
  # Join the collapsed exposure response to projected baseline
  # cases and valuations
  #------------------------------------------------------------
  
  endpoint_results <- endpoint_incidence %>%
    inner_join(endpoint_exposure_response,
               by = c("FIPS", "beta_pm10"), relationship = "many-to-many"
    ) %>%
    inner_join(
      morbidity_valuations_projected %>%
        filter(Endpoint == e) %>%
        select(Endpoint, Year, COI_proj),
      by = c("endpoint" = "Endpoint",
             "Year"
             ),
      relationship = "many-to-one"
    ) %>%
    mutate(
      morbidity = annualized_response * value * pop,
      PV_costs_COI = COI_proj * morbidity / (1 + 0.03)^(Year - 2024) / 1e6
    ) %>%
    group_by(scenario, Year, endpoint) %>%
    summarise(
      morbidity = sum(morbidity, na.rm = TRUE),
      PV_costs_COI =
        sum(PV_costs_COI, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario, Year) %>%
    group_by(scenario, endpoint) %>%
    mutate(
      cum_morbidity = cumsum(morbidity),
      PV_cum_costs_COI = cumsum(PV_costs_COI)
    ) %>%
    ungroup()
  
  #------------------------------------------------------------
  # Save immediately and release memory
  #------------------------------------------------------------
  
  endpoint_result_file <- file.path(output_folder, paste0(endpoint_file, ".csv"))
  
  write.csv(endpoint_results, endpoint_result_file, row.names = FALSE)
  
  endpoint_result_files[i] <- endpoint_result_file
  
  rm(endpoint_incidence, endpoint_betas, endpoint_exposure_response, endpoint_results)
  
  gc()
}



total_morbidity_projections <- endpoint_result_files %>%
  purrr::map_dfr(~ read.csv(.x, stringsAsFactors = FALSE)
  ) %>%
  arrange(endpoint, scenario, Year)

write.csv(total_morbidity_projections, file = "processed/total_morbidity_projections.csv", row.names = FALSE)






# total_morbidity_projections <- data.frame()
# for(e in unique(morbidity_valuations_projected$Endpoint)){
#   
#   # Population and incidence projections
#   ct_morbidity_pollution <- read.csv("processed/ct_incidence_morbidity_projections.csv", stringsAsFactors =  FALSE) %>%
#     filter(endpoint == e) %>%
#     #Merge w/ pollution deltas
#     inner_join(scenario_pm_deltas, by = "FIPS") %>%
#     inner_join(morbidity_valuations_projected %>% filter(Endpoint == e)
#               , by = c("Year", "endpoint" = "Endpoint")) %>%
#     filter(Year >= 2025) 
#   
#   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   #### estimate morbidity impacts
#   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   
#   # test <- ct_morbidity_pollution %>% head(500)
#   
#   ct_morbidity_projections_temp <- ct_morbidity_pollution %>%
#     mutate(incidence_rate_event = value*event_days, # incidence rates are already daily for morbidity
#            pm10_delta = ifelse(scenario == current_scenario, pm10_delta, relative_pm10_delta),
#            morbidity_pm10 = ((1-(1/exp(beta_pm10*pm10_delta)))*incidence_rate_event*pop)/n_years_storms,
#            morbidity = morbidity_pm10,
#            pm_delta = pm10_delta
#     )%>%
#     drop_na(scenario)
#   
#   #### Get overall impacts, not just relative
#   ct_morbidity_projections_current <- ct_morbidity_projections_temp %>%
#     filter(scenario == current_scenario) %>%
#     rename(current_morbidity = morbidity,
#            current_pm_delta = pm_delta) %>%
#     select(FIPS, County, date, Year, age_group, lower_age, upper_age, endpoint, current_morbidity, current_pm_delta)
#   
#   ct_morbidity_projections <- ct_morbidity_projections_temp %>%
#     inner_join(ct_morbidity_projections_current, 
#                by = c("FIPS", "County", "date", "Year", "age_group", "lower_age", "upper_age", "endpoint"),
#                relationship = "many-to-many")%>%
#     mutate(relative_morbidity = ifelse(scenario == current_scenario, 0, morbidity),
#            morbidity = ifelse(relative_morbidity + current_morbidity >= 0, relative_morbidity + current_morbidity, 0),
#            relative_pm_delta = ifelse(scenario == current_scenario, 0, pm_delta),
#            pm_delta = relative_pm_delta + current_pm_delta,
#            PV_costs_COI = COI_proj*morbidity/(1+0.03)^(Year - 2024)
#     ) %>%
#     select(FIPS, County, scenario, date, Year, age_group, lower_age, upper_age, pop, pm_delta, endpoint, morbidity, PV_costs_COI)
#   
#   #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   ## Aggregation
#   
#   total_morbidity_projections <- ct_morbidity_projections %>%
#     group_by(scenario, Year, endpoint) %>%
#     summarise(morbidity = sum(morbidity, na.rm = T),
#               PV_costs_COI = sum(PV_costs_COI, na.rm = T)/1e6
#     )%>%
#     ungroup %>%
#     group_by(scenario) %>%
#     mutate(cum_morbidity = cumsum(morbidity),
#            PV_cum_costs_COI = cumsum(PV_costs_COI))%>%
#     ungroup %>%
#     rbind(total_morbidity_projections)
#   
#   
#   
#   rm(ct_morbidity_projections, ct_morbidity_projections_current, ct_morbidity_projections_temp, ct_morbidity_pollution)
#   print(e)
#   }
# 
# write.csv(total_morbidity_projections, file = "processed/total_morbidity_projections.csv", row.names = FALSE)
# 
