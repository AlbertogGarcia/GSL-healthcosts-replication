# GSL dust costs: Lifetime costs of asthma onset
# albert.garcia@utah.edu
# created: 06/20/2025
# updated: 

# Set up environment ########################################

# load or install necessary libraries. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(measurements,  # convert units easily
               tidyverse, # tidyverse
               stringr,
               fuzzyjoin,
               ggplot2
)

options(scipen=999)  # turn off sci notation
options(dplyr.summarise.inform = FALSE)  # turn off dplyr group by comments
options(java.parameters = "-Xmx8000m") 
`%ni%` <- Negate(`%in%`)  # "not in" function




# CPI in 2020 and 2024
# (https://fred.stlouisfed.org/series/CPALTT01USA661S)
CPI_2024 <- 132.34870
CPI_2015 <- 100
CPI_1524 <- CPI_2024/CPI_2015

gdp_growth <- read.csv("data/health/realGDPpercapita_growth.csv", stringsAsFactors = F)%>%
  drop_na(rGDPpc)%>%
  mutate(growth_rate = rGDPpc/100,
         year = ifelse(str_sub(observation_date,-2,-1) > 25, paste0("19", str_sub(observation_date,-2,-1)), paste0("20", str_sub(observation_date,-2,-1)))
  )%>%
  select(year,growth_rate)

growth_rates <- gdp_growth %>%
  filter(year>2015)%>%
  mutate(cum_growth = cumprod(1 + growth_rate))

growth_GDPpc_1524 <- growth_rates %>% filter(year == 2024) %>% select(cum_growth) %>% pull()

income_elasticity <- 0.4

morbidity_valuations_2015 <- read.csv("data/health/morbidity_valuations_2015.csv", stringsAsFactors = F)

# EPA values used to calculate the wage-loss component of hospital admission COI estimates
median_weekly_wage_2015 <- 864
workdays_per_week <- 5
median_daily_wage_2015 <- median_weekly_wage_2015 / workdays_per_week

hospital_los <- tibble::tribble(
  ~Endpoint,  ~mean_los_days,
  "HA, All Cardiac Outcomes", 5.05,
  "HA, All Respiratory", 3.49,
  "HA, Asthma", 3.49,
  "HA, Stroke", 4.29
)

morbidity_valuations_2024 <- morbidity_valuations_2015 %>%
  left_join(hospital_los, by = "Endpoint") %>%
  mutate(lost_wages_2015 = if_else(
    str_starts(Endpoint, "HA,"),
    mean_los_days * median_daily_wage_2015,
    0),
    COI_2015_original = COI_2015,
    COI_2015 = COI_2015_original - lost_wages_2015,
    Qualifier = if_else(
      str_starts(Endpoint, "HA,"),
      "Medical costs",
      Qualifier
    ),
    COI_24 = COI_2015 *
      CPI_1524 *
      growth_GDPpc_1524^income_elasticity
  ) 

write.csv(morbidity_valuations_2024 %>%
            select(Endpoint, COI_24, Qualifier)
          
            , file = "data/health/morbidity_valuations_2024.csv", row.names = FALSE)