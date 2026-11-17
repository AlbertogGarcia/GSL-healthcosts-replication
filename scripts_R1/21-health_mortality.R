# GSL dust costs: Mortality impacts
# albert.garcia@utah.edu
# created: 05/28/2025
# updated: 09/26/2026

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

scenario_descrip <- c("No conservation",
                      "Current lake    ",
                      "Minimum healthy\nlake          ",
                      "Recent historical\naverage      "
                      )

n_years_storms <- 6

nonlinear_power <- 0.9
effect_increment <- 10

# Current scenario: attributable cases relative to the 4,203-ft counterfactual.
# Alternative scenarios: change in cases relative to the current scenario.
mortality_hif <- function(beta, pm_change, baseline_cases, is_current) {
  ifelse(
    is_current,
    (1 - exp(-beta * pm_change)) * baseline_cases,
    (exp(beta * pm_change) - 1) * baseline_cases
  )
}

mortality_hif_nonlinear <- function(
    beta,
    p,
    calibration_pm10,
    scenario_total_pm10,
    current_total_pm10,
    reference_pm10,
    baseline_cases,
    effect_increment = 10
) {
  
  # Match the published log-RR over a 10 µg/m3 increase
  # beginning at the representative calibration concentration.
  half_increment <- effect_increment / 2
  
  power_scale <- beta * effect_increment / (
    (calibration_pm10 + half_increment)^p -
      (calibration_pm10 - half_increment)^p
  )
  
  scenario_rr <- exp(
    power_scale * (
      scenario_total_pm10^p -
        current_total_pm10^p
    )
  )
  
  reference_rr <- exp(
    power_scale * (
      reference_pm10^p -
        current_total_pm10^p
    )
  )
  
  baseline_cases * (scenario_rr - reference_rr)
}

# main VSL and age-based VSL
VSL_24 = 12.57222

age_based_VSL_2024 <- read.csv("processed/age_based_VSL_2024.csv")%>%
  select(age, age_vsl_2024)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Load and merge processed data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#1 Emissions scenarios by water-level
scenario_pm_deltas <- read.csv("processed/scenario_pm_deltas_daily.csv", stringsAsFactors =  FALSE)%>%
  mutate(date = as.Date(date)) %>%
  filter(scenario %in% all_scenarios)

scenario_pm_deltas <- scenario_pm_deltas %>%
  mutate(
    reference_pm10 = non_gsl_dust_pm10,
    scenario_total_pm10 = reference_pm10 + pm10_delta,
    current_total_pm10 = reference_pm10 + current_pm10_delta
  )

nonlinear_calibration_pm10 <- scenario_pm_deltas %>%
  filter(scenario == current_scenario) %>%
  distinct(FIPS, date, reference_pm10) %>%
  summarise(
    calibration_pm10 = median(
      reference_pm10,
      na.rm = TRUE
    )
  ) %>%
  pull(calibration_pm10)


#2 Population and incidence
ct_incidence_mortality <- read.csv("processed/ct_incidence_mortality.csv", stringsAsFactors =  FALSE)

tract_population <- ct_incidence_mortality %>%
  filter(endpoint == "Mortality, All-cause") %>%
  group_by(FIPS) %>%
  summarise(
    tract_population = sum(pop, na.rm = TRUE),
    .groups = "drop"
  )

age_VSL_2024_grouped <- ct_incidence_mortality %>%
  select(lower_age, upper_age, age_group) %>%
  distinct()%>%
  tidyr::crossing(age_based_VSL_2024)%>%
  filter(age >= lower_age & age <= upper_age)%>%
  group_by(lower_age, upper_age, age_group) %>%
  summarise(age_vsl_2024 = mean(age_vsl_2024))

analysis_fips <- scenario_pm_deltas %>%
  distinct(FIPS)

study_population_by_county <- ct_incidence_mortality %>%
  filter(endpoint == "Mortality, All-cause") %>%
  semi_join(analysis_fips, by = "FIPS") %>%
  distinct(
    FIPS,
    County,
    lower_age,
    upper_age,
    pop
  ) %>%
  group_by(County) %>%
  summarise(
    population = sum(pop, na.rm = TRUE),
    .groups = "drop"
  )

study_population <- sum(
  study_population_by_county$population,
  na.rm = TRUE
)

study_population_by_county
study_population


#Merge w/ pollution deltas 
ct_mortality_pollution <- ct_incidence_mortality %>%
  left_join(scenario_pm_deltas, by = "FIPS")%>%
  left_join(age_VSL_2024_grouped, by = c("lower_age", "upper_age", "age_group"))

table(ct_mortality_pollution$County)

rm("ct_incidence_mortality")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Mortality impacts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Coefficients from Orellano et al., (2020)
###IMPORTANT: all are for a 10 microgram increase

#Mortality impact
ct_mortality_age_temp <- ct_mortality_pollution %>%
  mutate(
    beta_pm10 = case_when(
      endpoint == "Mortality, All-cause" ~ log(1.0041)/10,
      endpoint == "Mortality, Respiratory" ~ log(1.0091)/10,
      endpoint == "Mortality, Cardiovascular" ~ log(1.0060)/10
      ),
    beta_pm10_lower = case_when(
      endpoint == "Mortality, All-cause" ~ log(1.0034)/10,
      endpoint == "Mortality, Respiratory" ~ log(1.0063)/10,
      endpoint == "Mortality, Cardiovascular" ~ log(1.0044)/10
    ),
    beta_pm10_upper = case_when(
      endpoint == "Mortality, All-cause" ~ log(1.0049)/10,
      endpoint == "Mortality, Respiratory" ~ log(1.0119)/10,
      endpoint == "Mortality, Cardiovascular" ~ log(1.0077)/10
    ),
    pm10_change = ifelse(
      scenario == current_scenario,
      pm10_delta,
      relative_pm10_delta
    ),
    baseline_cases = incidence_rate_daily * pop / n_years_storms,
    mortality = mortality_hif(
      beta_pm10, pm10_change, baseline_cases,
      scenario == current_scenario
    ),
    mortality_pm10_lower = mortality_hif(
      beta_pm10_lower, pm10_change, baseline_cases,
      scenario == current_scenario
    ),
    mortality_pm10_upper = mortality_hif(
      beta_pm10_upper, pm10_change, baseline_cases,
      scenario == current_scenario
    ),
    mortality_lower = mortality_pm10_lower,
    mortality_upper = mortality_pm10_upper,
    # Sensitivity:
    mortality_linear = baseline_cases * beta_pm10 * pm10_change,
    mortality_nonlinear = mortality_hif_nonlinear(
      beta = beta_pm10,
      p = nonlinear_power,
      calibration_pm10 = nonlinear_calibration_pm10,
      scenario_total_pm10 = scenario_total_pm10,
      current_total_pm10 = current_total_pm10,
      reference_pm10 = reference_pm10,
      baseline_cases = baseline_cases,
      effect_increment = effect_increment
    ),
    pm_delta = pm10_change
    )%>%
  drop_na(scenario)

ct_mortality_age_current <- ct_mortality_age_temp %>%
  filter(scenario == current_scenario) %>%
  rename(current_mortality = mortality,
         current_mortality_lower = mortality_lower,
         current_mortality_upper = mortality_upper,
         current_mortality_linear = mortality_linear,
         current_pm_delta = pm_delta) %>%
  select(FIPS, centroid_name, County, date, age_group, lower_age, upper_age, endpoint, 
         current_mortality, current_mortality_linear,
         current_mortality_lower, current_mortality_upper, current_pm_delta)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Get overall mortality impacts
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm("ct_mortality_pollution")

ct_mortality_age <- ct_mortality_age_temp %>%
  left_join(ct_mortality_age_current, by = c("FIPS", "County", "date", "age_group", "lower_age", "upper_age", "endpoint"))%>%
  mutate(relative_mortality = ifelse(scenario == current_scenario, 0, mortality),
         mortality = relative_mortality + current_mortality,
         relative_mortality_linear = ifelse(scenario == current_scenario, 0, mortality_linear),
         mortality_linear = relative_mortality_linear + current_mortality_linear,
         relative_mortality_lower = ifelse(scenario == current_scenario, 0, mortality_lower),
         mortality_lower = relative_mortality_lower + current_mortality_lower,
         relative_mortality_upper = ifelse(scenario == current_scenario, 0, mortality_upper),
         mortality_upper = relative_mortality_upper + current_mortality_upper,
         relative_pm_delta = ifelse(scenario == current_scenario, 0, pm_delta),
         pm_delta = relative_pm_delta + current_pm_delta,
         costs_VSL = mortality*VSL_24,
         costs_VSL_lower = mortality_lower*VSL_24,
         costs_VSL_upper = mortality_upper*VSL_24,
         costs_VSL_linear = mortality_linear*VSL_24,
         costs_VSL_nonlinear = mortality_nonlinear*VSL_24,
         costs_age_VSL = mortality*age_vsl_2024) %>%
  select(FIPS, County, scenario, date, age_group, lower_age, upper_age, incidence_rate, pop, pm_delta, endpoint, 
         mortality, mortality_linear, mortality_nonlinear, mortality_lower, mortality_upper, costs_age_VSL, costs_VSL, costs_VSL_lower, costs_VSL_upper, costs_VSL_linear, costs_VSL_nonlinear)

mortality_shape_sensitivity <- ct_mortality_age %>%
  filter(endpoint == "Mortality, All-cause") %>%
  group_by(scenario) %>%
  summarise(
    mortality = sum(mortality, na.rm = TRUE),
    linear_mortality = sum(mortality_linear, na.rm = TRUE),
    nonlinear_mortality = sum(mortality_nonlinear, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    costs_VSL = mortality * VSL_24,
    linear_costs_VSL = linear_mortality * VSL_24,
    nonlinear_costs_VSL = nonlinear_mortality * VSL_24
  )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Mortality costs by concentration-response functional form
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mortality_shape_costs_plot_prep <-
  mortality_shape_sensitivity %>%
  filter(scenario %in% relevant_scenarios) %>%
  select(
    scenario,
    linear_costs_VSL,
    costs_VSL,
    nonlinear_costs_VSL
  ) %>%
  pivot_longer(
    cols = c(
      linear_costs_VSL,
      costs_VSL,
      nonlinear_costs_VSL
    ),
    names_to = "functional_form",
    values_to = "costs"
  ) %>%
  mutate(
    functional_form = case_when(
      functional_form == "linear_costs_VSL" ~
        "Strictly linear",
      
      functional_form == "costs_VSL" ~
        "Log-linear (primary)",
      
      functional_form == "nonlinear_costs_VSL" ~
        "Nonlinear"
    ),
    
    functional_form = factor(
      functional_form,
      levels = c(
        "Strictly linear",
        "Log-linear (primary)",
        "Nonlinear"
      )
    ),
    
    scenario = factor(
      scenario,
      levels = relevant_scenarios
    )
  )

shape_colors <- setNames(
  c(
    "grey45",
    palette$blue,
    palette$red
  ),
  c(
    "Strictly linear",
    "Log-linear (primary)",
    "Nonlinear"
  )
)

mortality_shape_costs_plot <-
  mortality_shape_costs_plot_prep %>%
  ggplot(
    aes(
      x = scenario,
      y = costs,
      fill = functional_form
    )
  ) +
  geom_col(
    position = position_dodge(width = 0.80),
    width = 0.70
  ) +
  geom_text(
    aes(
      label = paste0(
        "$",
        scales::number(
          costs,
          accuracy = 0.1
        ),
        "M"
      )
    ),
    position = position_dodge(width = 0.80),
    hjust = -0.10,
    size = 3.5
  ) +
  scale_x_discrete(
    labels = scenario_descrip
  ) +
  scale_y_continuous(
    name = "Mortality costs (USD millions)",
    breaks = scales::breaks_width(20),
    labels = scales::label_number(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.16)
    )
  ) +
  scale_fill_manual(
    values = shape_colors,
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "Expected annual mortality costs",
    title = paste0(
      "Mortality-cost sensitivity to ",
      "concentration–response form"
    )
  ) +
  coord_flip(
    clip = "off"
  ) +
  theme_cowplot(16) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    plot.title = element_text(
      hjust = 0.5
    ),
    plot.margin = margin(
      t = 5.5,
      r = 25,
      b = 5.5,
      l = 5.5
    )
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      reverse = FALSE
    )
  )

mortality_shape_costs_plot

ggsave("figs_pdf/mortality_costs_functional_form.pdf", plot = mortality_shape_costs_plot, width = 9, height = 7)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Create data for analyses
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# total mortality by census tract
ct_mortality <- ct_mortality_age %>%
  group_by(scenario, FIPS, County, endpoint) %>%
  summarise(pm_delta = mean(pm_delta, na.rm = T),
            mortality = sum(mortality, na.rm = T),
            population = sum(pop, na.rm = T)/n_distinct(date),
            costs_VSL = sum(costs_VSL, na.rm = T),
            costs_age_VSL = sum(costs_age_VSL, na.rm = T),
            ct_incidence_rate_annual = weighted.mean(incidence_rate, pop, na.rm = T)
            )%>%
  ungroup

write.csv(ct_mortality, file = "processed/ct_mortality.csv", row.names = FALSE)

county_mortality_4192 <- ct_mortality_age %>%
  group_by(County, scenario, endpoint) %>%
  summarise(mortality = sum(mortality, na.rm = T))%>%
  ungroup %>%
  filter(scenario == 4192)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Total overall mortality
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_mortality <- ct_mortality_age %>%
  group_by(scenario, endpoint) %>%
  summarise(mortality = sum(mortality, na.rm = T),
            mortality_lower = sum(mortality_lower, na.rm = T),
            mortality_upper = sum(mortality_upper, na.rm = T),
            mortality_linear = sum(mortality_linear, na.rm = T),
            mortality_nonlinear = sum(mortality_nonlinear, na.rm = T),
            costs_VSL = sum(costs_VSL, na.rm = T),
            costs_VSL_lower = sum(costs_VSL_lower, na.rm = T),
            costs_VSL_upper = sum(costs_VSL_upper, na.rm = T),
            costs_VSL_linear = sum(costs_VSL_linear, na.rm = T),
            costs_VSL_nonlinear = sum(costs_VSL_nonlinear, na.rm = T),
            costs_age_VSL = sum(costs_age_VSL, na.rm = T)
  )%>%
  ungroup

write.csv(total_mortality, file = "processed/total_mortality.csv", row.names = FALSE)


mortality_plot_prep <- total_mortality %>%
  select(scenario, mortality, endpoint) %>%
  filter(scenario %in% relevant_scenarios)%>%
  group_by(scenario)%>%
  mutate(mortality = ifelse(endpoint == "Mortality, All-cause", 2*mortality - sum(mortality), mortality),
         endpoint = case_when(endpoint == "Mortality, All-cause" ~ "All other causes",
                              endpoint == "Mortality, Respiratory" ~ "Respiratory",
                              endpoint == "Mortality, Cardiovascular" ~ "Cardiovascular"),
         total_mortality = sum(mortality),
         costs_VSL = mortality*VSL_24,
         total_costs_VSL = total_mortality*VSL_24
         ) %>%
  ungroup 

mortality_costs_plot <- mortality_plot_prep %>%
  ggplot(aes(x=reorder(scenario, scenario, order = T), y=mortality, fill = reorder(endpoint, - mortality))
  )+
  geom_bar(stat='identity')+
  scale_x_discrete(labels = scenario_descrip)+
  scale_y_continuous(
    "Expected premature mortalities", breaks = c(1,3,5),
    sec.axis = sec_axis(~ . * VSL_24, 
                        name = "Costs (millions USD)",
                        breaks = seq(0, 80, 20)
    )
  )+
  ggtitle("Current annual dust-induced mortality costs") + 
  scale_fill_manual(values = c(palette$blue, palette$red, palette$green))+
  theme_cowplot(16)+
  theme(legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)
  )+
  coord_flip()+
  guides(fill = guide_legend(title = NULL, reverse=T, hjust = 0.5))
mortality_costs_plot
ggsave("figs_pdf/mortality_costs.pdf", width = 9, height = 7)

#%%%%%%%%%%%%%%%%%%%%%%%%
## Plot by VSL type
#%%%%%%%%%%%%%%%%%%%%%%%

mortality_VSL_plot <- total_mortality %>%
  filter(
    scenario %in% relevant_scenarios,
    endpoint == "Mortality, All-cause"
  ) %>%
  pivot_longer(
    cols = c(
      costs_VSL,
      costs_age_VSL
    ),
    names_to = "VSL_type",
    values_to = "costs"
  ) %>%
  mutate(
    VSL_type = case_when(
      VSL_type == "costs_VSL" ~
        "EPA recommended",
      
      VSL_type == "costs_age_VSL" ~
        "Cohort-adjusted\n(Aldy & Viscusi 2008)"
    )
  ) %>%
  ggplot(
    aes(
      x = reorder(scenario, scenario, order = TRUE),
      y = costs,
      fill = reorder(VSL_type, costs)
    )
  ) +
  geom_col(
    width = 0.8,
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    aes(
      label = paste0(
        "$",
        scales::number(
          costs,
          accuracy = 0.1
        ),
        "M"
      )
    ),
    position = position_dodge(width = 0.9),
    hjust = -0.10,
    size = 3.5
  ) +
  scale_x_discrete(
    labels = scenario_descrip
  ) +
  scale_y_continuous(
    name = "Mortality costs (USD millions)",
    breaks = scales::breaks_width(20),
    labels = scales::label_number(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0.16)
    )
  ) +
  scale_fill_manual(
    values = c(
      palette$red,
      palette$blue
    ),
    name = "VSL"
  ) +
  labs(
    x = NULL,
    y = "Expected annual mortality costs"
  ) +
  theme_cowplot(16) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(
      t = 5.5,
      r = 25,
      b = 5.5,
      l = 5.5
    )
  ) +
  coord_flip(
    clip = "off"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      reverse = TRUE
    )
  )

mortality_VSL_plot

ggsave("figs_pdf/mortality_by_VSL.pdf", width = 9, height = 7)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Mortality by age group
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mortality_age <- ct_mortality_age %>%
  group_by(scenario, age_group, lower_age, upper_age, endpoint) %>%
  summarise(mortality = sum(mortality, na.rm = T),
            pop = sum(pop, na.rm = T)/n_distinct(date),
            costs_VSL = sum(costs_VSL, na.rm = T),
            costs_age_VSL = sum(costs_age_VSL, na.rm = T))%>%
  ungroup


mortality_age %>% 
  filter(scenario == current_scenario, 
         endpoint == "Mortality, All-cause"
         )%>%
  mutate(proportion_of_mortality = mortality/sum(mortality),
         age_group = ifelse(age_group == "85 and over", "Over 85", age_group))%>%
  group_by(age_group, lower_age, upper_age) %>% 
  summarise(proportion_of_mortality = mean(proportion_of_mortality))%>%
  ggplot(aes(x=reorder(age_group, lower_age), y=proportion_of_mortality))+
  geom_bar(stat='identity')+
  ggtitle("All-cause mortality by age")+
  scale_y_continuous(label = scales::percent, 
                    # breaks = seq(from = 0.01, to = .07, by = 0.01),
                     name = "Percent of total expected mortality")+
  theme_cowplot(16)+
  theme(axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 11)
  )
ggsave("figs_pdf/mortality_by_age.pdf",
       width = 9, height = 6)


mortality_county <- ct_mortality_age %>%
  group_by(scenario, County, endpoint) %>%
  summarise(mortality = sum(mortality, na.rm = T),
            pop = sum(pop, na.rm = T)/n_distinct(date),
            costs_VSL = sum(costs_VSL, na.rm = T),
            costs_age_VSL = sum(costs_age_VSL, na.rm = T))%>%
  ungroup
mortality_county


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Mortality by race
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rm(
  list = setdiff(
    ls(envir = .GlobalEnv),
    c(
      "mortality_hif",
      "current_scenario",
      "n_years_storms",
      "VSL_24",
      "palette",
      "relevant_scenarios",
      "age_based_VSL_2024",
      "scenario_pm_deltas"
    )
  ),
  envir = .GlobalEnv
)

# Same but disaggregated by race
ct_incidence_mortality_race <- read.csv("processed/ct_incidence_mortality_race.csv", stringsAsFactors =  FALSE) 

age_VSL_2024_grouped_race <- ct_incidence_mortality_race %>%
  select(lower_age, upper_age, age_group) %>%
  distinct()%>%
  tidyr::crossing(age_based_VSL_2024)%>%
  filter(age >= lower_age & age <= upper_age)%>%
  group_by(lower_age, upper_age, age_group) %>%
  summarise(age_vsl_2024 = mean(age_vsl_2024))

#Merge w/ pollution deltas 
ct_mortality_pollution_race <- ct_incidence_mortality_race %>%
  left_join(scenario_pm_deltas, by = "FIPS")%>%
  left_join(age_VSL_2024_grouped_race, by = c("lower_age", "upper_age", "age_group"))

rm("ct_incidence_mortality_race")

#All-cause coefficients
beta_pm25 <- log(1.0065)/10
beta_pm25_lower <- log(1.0044)/10
beta_pm25_upper <- log(1.0086)/10
beta_pm10 <- log(1.0041)/10
beta_pm10_lower <- log(1.0034)/10
beta_pm10_upper <- log(1.0049)/10
#Mortality impact
ct_mortality_agebyrace_temp <- ct_mortality_pollution_race %>%
  mutate(pm10_change = ifelse(
           scenario == current_scenario,
           pm10_delta,
           relative_pm10_delta
         ),
         baseline_cases = incidence_rate_daily * pop / n_years_storms,
         mortality = mortality_hif(
           beta_pm10, pm10_change, baseline_cases,
           scenario == current_scenario
         ),
         mortality_pm10_lower = mortality_hif(
           beta_pm10_lower, pm10_change, baseline_cases,
           scenario == current_scenario
         ),
         mortality_pm10_upper = mortality_hif(
           beta_pm10_upper, pm10_change, baseline_cases,
           scenario == current_scenario
         ),
         mortality_lower = mortality_pm10_lower,
         mortality_upper = mortality_pm10_upper,
         pm_delta = pm10_change
         )%>%
  drop_na(scenario)


#### Get overall impacts, not just relative
ct_mortality_agebyrace_current <- ct_mortality_agebyrace_temp %>%
  filter(scenario == current_scenario) %>%
  rename(current_mortality = mortality,
         current_mortality_lower = mortality_lower,
         current_mortality_upper = mortality_upper,
         current_pm_delta = pm_delta) %>%
  select(FIPS, County, date, race, age_group, lower_age, upper_age, endpoint, current_mortality, current_mortality_lower, current_mortality_upper, current_pm_delta)

ct_mortality_agebyrace <- ct_mortality_agebyrace_temp %>%
  left_join(ct_mortality_agebyrace_current, by = c("FIPS", "County", "date", "race", "age_group", "lower_age", "upper_age", "endpoint"))%>%
  mutate(relative_mortality = ifelse(scenario == current_scenario, 0, mortality),
         relative_mortality_lower = ifelse(scenario == current_scenario, 0, mortality_lower),
         relative_mortality_upper = ifelse(scenario == current_scenario, 0, mortality_upper),
         mortality = relative_mortality + current_mortality,
         mortality_lower = relative_mortality_lower + current_mortality_lower,
         mortality_upper = relative_mortality_upper + current_mortality_upper,
         relative_pm_delta = ifelse(scenario == current_scenario, 0, pm_delta),
         pm_delta = relative_pm_delta + current_pm_delta,
         costs_VSL = mortality*VSL_24,
         costs_age_VSL = mortality*age_vsl_2024) %>%
  select(FIPS, County, scenario, date, race, age_group, lower_age, upper_age, incidence_rate, pop, pm_delta, endpoint, mortality, mortality_lower, mortality_upper, costs_VSL, costs_age_VSL)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Aggregate to relevant scales
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct_mortality_race <- ct_mortality_agebyrace %>%
  group_by(scenario, FIPS, County, race, endpoint) %>%
  summarise(mortality = sum(mortality, na.rm = T),
            mortality_lower = sum(mortality_lower, na.rm = T),
            mortality_upper = sum(mortality_upper, na.rm = T),
            population = sum(pop, na.rm = T)/n_distinct(date),
            costs_VSL = sum(costs_VSL, na.rm = T),
            costs_age_VSL = sum(costs_age_VSL, na.rm = T))%>%
  ungroup

write.csv(ct_mortality_race, file = "processed/ct_mortality_race.csv", row.names = FALSE)

total_mortality_race <- ct_mortality_agebyrace %>%
  group_by(scenario, race) %>%
  summarise(mortality = sum(mortality, na.rm = T),
            population = sum(pop, na.rm = T)/n_distinct(date),
            pm_delta = weighted.mean(pm_delta, pop, na.rm = T),
            costs_VSL = sum(costs_VSL, na.rm = T),
            costs_age_VSL = sum(costs_age_VSL, na.rm = T)
            )%>%
  mutate(mortality_per_100k = mortality/population*100000,
         cost_per_capita = costs_VSL/population,
         age_cost_per_capita = costs_age_VSL/population
         )%>%
  ungroup %>%
  filter(race %in% c(
    "Asian alone", "White NH", "Hispanic or Latino", "Black or African American Alone", "Hawaiian and PI Alone"
  )
  )
write.csv(total_mortality_race, file = "processed/total_mortality_race.csv", row.names = FALSE)

mortality_race <- total_mortality_race %>%
ggplot(aes(x=scenario, color = reorder(race, -mortality_per_100k), y=mortality_per_100k))+
  geom_line(linewidth = 1)+
  geom_point(data = total_mortality_race %>% filter(scenario %in% relevant_scenarios)
             , size = 2.5)+
  # ggtitle("Distribution of mortality risk across race")+
  scale_y_continuous(name = paste("Premature mortality (per 100k)"),
                     limits = c(0, max(total_mortality_race$mortality_per_100k))
                     )+
  #scale_x_reverse(
  scale_x_continuous(
    breaks = relevant_scenarios,
                  name = "GSL water level (ftASL)")+
  scale_color_manual(values = c(palette$red, palette$blue, palette$dark, palette$orange, palette$green),
                     #name = "Race",
                     labels = c("Hawaiian/Pacific Islander", "White Non-hispanic", "Asian", "Hispanic/Latino", "Black/African American")
                     )+
  ggtitle("All-cause mortality by race")+
  theme_cowplot(16)+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.justification = "center"
        )
mortality_race

ggsave("figs_pdf/mortality_by_race.pdf", 
       width = 9, height = 5)


ct_mortality_map <- ct_mortality_agebyrace %>%
  group_by(scenario, FIPS, County, endpoint) %>%
  summarise(pm_delta = weighted.mean(pm_delta, pop, na.rm = T),
            mortality = sum(mortality, na.rm = T),
            population = sum(pop, na.rm = T)/n_distinct(date),
            costs_VSL = sum(costs_VSL, na.rm = T),
            ct_incidence_rate_annual = weighted.mean(incidence_rate,pop)
  )%>%
  ungroup
write.csv(ct_mortality_map, file = "processed/ct_mortality_map.csv", row.names = FALSE)

