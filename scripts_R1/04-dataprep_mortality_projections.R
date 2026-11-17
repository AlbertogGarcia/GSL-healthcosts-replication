# GSL dust costs: Mortality and population projections
# albert.garcia@utah.edu
# created: 06/20/2025
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
               zoo # needed for interpolation of future adjusted mortality incidence
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
                "purple" = "#8da0cb"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### MORTALITY
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACS_year_final = 2023

# Projected county-level population growth by age

gardner_pop_projections <- readxl::read_excel("data/population/GardnerPI_projections.xlsx", sheet = 2)

gardner_pop_projections %>%
  filter(str_detect(Geography, " County$")) %>%
  mutate(
    "Under 18" =
      `Population Ages 0-4` + `Population Ages 5-17`,
    "18 to 64" =
      `Population Ages 18-24` + `Population Ages 25-64`,
    "Over 65" =
      `Population Ages 65-84` + `Population Ages 85+`
  ) %>%
  select(Year, "Under 18", "18 to 64", "Over 65") %>%
  pivot_longer(
    cols = 2:4,
    names_to = "age_group",
    values_to = "pop"
  ) %>%
  filter(Year >= ACS_year_final) %>%
  mutate(
    lower_age = case_when(
      age_group == "Under 18" ~ 0,
      age_group == "18 to 64" ~ 18,
      age_group == "Over 65" ~ 65
    )
  ) %>%
  group_by(age_group, lower_age, Year) %>%
  summarise(
    pop = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(age_group, Year) %>%
  group_by(age_group, lower_age) %>%
  mutate(
    baseline_pop = pop[Year == ACS_year_final][1],
    cum_growth = pop / baseline_pop
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = Year,
      y = cum_growth,
      color = reorder(age_group, lower_age)
    )
  ) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c(palette$green, palette$blue, palette$red)) +
  ylab("Population growth index") +
  ggtitle("Normalized Utah population projection by age group") +
  theme_classic(12) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(filename = "age_projections.pdf", path = "figs_pdf",
       height = 5, width = 6)

# Match projection age groups to mortality age groups
county_pop_projections <- gardner_pop_projections %>%
  filter(str_detect(Geography, " County$")) %>%
  mutate(
    `Population Ages 5-24` =
      `Population Ages 5-17` + `Population Ages 18-24`
  ) %>%
  select(Geography, Year,
    `Population Ages 0-4`,
    `Population Ages 5-24`,
    `Population Ages 25-64`,
    `Population Ages 65-84`,
    `Population Ages 85+`
  ) %>%
  pivot_longer(
    cols = `Population Ages 0-4`:`Population Ages 85+`,
    names_to = "age_group",
    values_to = "projected_pop"
  ) %>%
  rename(County = Geography) %>%
  mutate(
    age_group = str_remove(
      age_group,
      "Population Ages "
    ),
    age_group = str_replace(
      age_group,
      "-",
      " to "
    ),
    start.age = as.numeric(
      str_sub(age_group, 1, 2)
    ),
    end.age = ifelse(
      start.age == 85,
      100,
      as.numeric(str_sub(age_group, -2, -1))
    ),
    County = str_remove(
      County,
      " County"
    )
  ) %>%
  filter(Year >= ACS_year_final) %>%
  arrange(County, start.age, end.age, Year) %>%
  group_by(County, start.age, end.age) %>%
  mutate(
    baseline_projected_pop =
      projected_pop[Year == ACS_year_final][1],
    cum_growth =
      projected_pop / baseline_projected_pop
  ) %>%
  ungroup() %>%
  select(County, Year, start.age, end.age, cum_growth)


ct_incidence_mortality <- read.csv(
  "processed/ct_incidence_mortality.csv",
  stringsAsFactors = FALSE
) %>%
  filter(endpoint == "Mortality, All-cause")

ct_incidence_projections <- ct_incidence_mortality %>%
  fuzzyjoin::fuzzy_left_join(
    county_pop_projections,
    by = c(
      "County" = "County",
      "lower_age" = "start.age",
      "upper_age" = "end.age"
    ),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
  mutate(
    County = dplyr::coalesce(County.x, County.y),
    pop = pop * cum_growth) %>%
  select(-County.x, -County.y, -start.age, -end.age)


ct_incidence_projections_unadjusted <-
  ct_incidence_projections %>%
  select(-cum_growth)


# Were invalid incidence rates already present in the input?
ct_incidence_mortality %>%
  summarise(
    n = n(),
    n_missing = sum(is.na(incidence_rate)),
    n_nan = sum(is.nan(incidence_rate)),
    n_infinite = sum(is.infinite(incidence_rate))
  )

# Did the population projection join fail?
ct_incidence_projections_unadjusted %>%
  summarise(
    n = n(),
    n_missing_pop = sum(is.na(pop)),
    n_missing_year = sum(is.na(Year)),
    n_missing_incidence = sum(is.na(incidence_rate)),
    n_nan_incidence = sum(is.nan(incidence_rate))
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Projected mortality-rate adjustment factors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Rebase the EPA adjustment factors
mortality_adjustmentfactors <- readxl::read_excel(paste0("data/health/mortality/", "census_mortality_adjustmentfactors_2013relative.xlsx")) %>%
  group_by(start.age, end.age) %>%
  complete(Year = min(Year):max(Year)) %>%
  arrange(start.age, end.age, Year) %>%
  mutate(
    # Interpolate missing annual factors
    Factor = zoo::na.approx(
      Factor,
      x = Year,
      na.rm = FALSE
    ),
    # Normalize each EPA age group's series to base year
    baseline_factor = Factor[Year == 2025][1],
    Factor = Factor / baseline_factor
  ) %>%
  ungroup() %>%
  select(-baseline_factor)

mortality_age_groups <-
  ct_incidence_projections_unadjusted %>%
  distinct(
    lower_age,
    upper_age
  )

mortality_adjustment_by_age <-
  mortality_adjustmentfactors %>%
  tidyr::crossing(
    mortality_age_groups
  ) %>%
  mutate(
    overlap_ages = pmax(
      0,
      pmin(end.age, upper_age) -
        pmax(start.age, lower_age) + 1
    )
  ) %>%
  filter(overlap_ages > 0) %>%
  group_by(
    Year,
    lower_age,
    upper_age
  ) %>%
  summarise(
    Factor = weighted.mean(
      Factor,
      overlap_ages,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ct_incidence_projections_adj <-
  ct_incidence_projections_unadjusted %>%
  filter(Year >= 2025) %>%
  left_join(
    mortality_adjustment_by_age,
    by = c(
      "Year",
      "lower_age",
      "upper_age"
    ),
    relationship = "many-to-one"
  ) %>%
  mutate(
    incidence_rate = incidence_rate * Factor,
    incidence_rate_daily = incidence_rate_daily * Factor
  ) %>%
  select(all_of(names(ct_incidence_projections_unadjusted)))

write.csv(ct_incidence_projections_adj, file = "processed/ct_incidence_projections_adj.csv", row.names = FALSE)

