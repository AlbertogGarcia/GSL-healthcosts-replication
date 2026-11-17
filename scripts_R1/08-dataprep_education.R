# GSL dust costs: Mortality and asthma data processing
# albert.garcia@utah.edu
# created: 06/20/2025
# updated: 

# Set up environment ########################################

# packages=c("xlsx", "gdata", "dplyr","tidyr", "stringr", "fuzzyjoin", "stringr", 
#            "ggplot2", "stargazer", "plm", "cowplot", "sf", "lwgeom","data.table")

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
               purrr,
               ggplot2
)

options(scipen=999)  # turn off sci notation
options(dplyr.summarise.inform = FALSE)  # turn off dplyr group by comments
options(java.parameters = "-Xmx8000m") 
`%ni%` <- Negate(`%in%`)  # "not in" function

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### School-age populations: ages 5--14
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

age_sheet <- 3

ct_ACS_age_wide <- readxl::read_excel(
  "data/population/Demographic_Raw_Tables.xlsx",
  sheet = age_sheet
)

ct_ACS_5to14 <- ct_ACS_age_wide %>%
  select(-`Total population`) %>%
  pivot_longer(
    cols = `Under 5 years`:ncol(.),
    names_to = "age_group",
    values_to = "pop"
  ) %>%
  filter(
    age_group %in% c(
      "5 to 9 years",
      "10 to 14 years"
    )
  ) %>%
  mutate(
    age_group = str_remove(age_group, " years"),
    lower_age = as.numeric(str_sub(age_group, 1, 2)),
    upper_age = as.numeric(str_sub(age_group, -2, -1))
  )


agebyrace_sheet <- 11

ct_ACS_agebyrace_wide <- readxl::read_excel(
  "data/population/Demographic_Raw_Tables.xlsx",
  sheet = agebyrace_sheet
)

ct_ACS_race_5to14 <- ct_ACS_agebyrace_wide %>%
  pivot_longer(
    cols = 3:ncol(.),
    names_to = "agebyrace_group",
    values_to = "pop"
  ) %>%
  filter(!str_detect(agebyrace_group, "Total")) %>%
  mutate(
    agebyrace_group =
      str_remove_all(agebyrace_group, "[()]"),
    
    race = case_when(
      str_detect(
        agebyrace_group,
        "American Indian and Alaska Native"
      ) ~ "American Indian and Alaska Native Alone",
      
      str_detect(agebyrace_group, "Hispanic") ~
        "Hispanic or Latino",
      
      str_detect(agebyrace_group, "White") ~
        "White NH",
      
      str_detect(agebyrace_group, "Black") ~
        "Black or African American Alone",
      
      str_detect(agebyrace_group, "Asian") ~
        "Asian alone",
      
      str_detect(agebyrace_group, "Hawaiian and PI") ~
        "Hawaiian and PI Alone",
      
      str_detect(agebyrace_group, "Other Race") ~
        "Other Race Alone",
      
      str_detect(agebyrace_group, "Multi Racial") ~
        "Multi Racial"
    ),
    
    age_group = str_remove(
      agebyrace_group,
      paste0(race, " ")
    )
  ) %>%
  filter(
    age_group %in% c(
      "5 to 9 years",
      "10 to 14 years"
    )
  ) %>%
  mutate(
    age_group = str_remove(age_group, " years"),
    lower_age = as.numeric(str_sub(age_group, 1, 2)),
    upper_age = as.numeric(str_sub(age_group, -2, -1))
  ) %>%
  select(-agebyrace_group)






# Read in school district shapefiles
districts.shp <- read_sf("data/gis/UtahSchoolDist/SchoolDistricts.shp") %>% select(NAME, TOTENROLL) %>%
  mutate(NAME = str_remove(tolower(NAME), " county"))
# read in census tract shapefiles
tracts.shp <- read_sf("data/gis/tracts/CensusTracts2020.shp")%>%
  rename(FIPS = GEOID20)%>%select(FIPS)


# Read in absence rates and merge with cts
attendance <- readxl::read_excel(
  "data/schools/attendance23.xlsx",
  sheet = 5
) %>%
  filter(
    Measure == "Attendance Rate",
    `LEA Type` == "District",
    LEA_name != "Utah Schools for Deaf & Blind"
  ) %>%
  select(
    LEA_name,
    `K-12 Overall`,
    `AfAm/Black`,
    `American Indian`,
    `Asian`,
    `Hispanic/Latino`,
    `Multiple Race`,
    `Pacific Islander`,
    `White`
  ) %>%
  pivot_longer(
    cols = `K-12 Overall`:`White`,
    names_to = "Race",
    values_to = "attendance_rate"
  ) %>%
  filter(
    !is.na(attendance_rate),
    attendance_rate > 0,
    attendance_rate <= 1
  ) %>%
  mutate(
    # Probability of absence on a school day
    absence_probability_school_day =
      1 - attendance_rate,
    
    # Calendar-day probability used because modeled dust
    # events can occur on any calendar day
    incidence_rate_daily =
      absence_probability_school_day * 180 / 365,
    
    NAME = str_remove(
      tolower(LEA_name),
      " district"
    )
  ) %>%
  select(
    -attendance_rate,
    -absence_probability_school_day,
    -LEA_name
  )


jw <- comparator::JaroWinkler()

ct_attendance <- tracts.shp %>%
  st_centroid() %>%
  st_join(districts.shp) %>%
  fuzzyjoin::fuzzy_full_join(
    attendance, 
    by = "NAME",
    match_fun = function(x, y) { jw(x, y) > 0.1}
  )%>%
  mutate(jw = jw(NAME.x, NAME.y),
         FIPS = as.numeric(FIPS)) %>%
  group_by(NAME.x) %>%
  filter(jw == max(jw)) %>%
  ungroup %>%
  rename(District = NAME.y,
         District_enrollment = TOTENROLL) %>%
  select(FIPS, District, Race, District_enrollment, incidence_rate_daily) %>%
  st_drop_geometry() %>%
  mutate(race = case_when(
    str_detect(Race, "White") ~ "White NH",
    str_detect(Race, "Hispanic") ~ "Hispanic or Latino",
    str_detect(Race, "American Indian") ~ "American Indian and Alaska Native Alone",
    str_detect(Race, "Black") ~ "Black or African American Alone",
    str_detect(Race, "Asian") ~ "Asian alone",
    str_detect(Race, "Pacific Islander") ~ "Hawaiian and PI Alone",
    str_detect(Race, "Multi") ~ "Multi Racial"
  )
  )

ct_school <- ct_attendance %>%
  filter(Race == "K-12 Overall") %>%
  select(-Race, -race) %>%
  right_join(
    ct_ACS_5to14,
    by = "FIPS"
  )
sum(ct_school$pop)

# This will overcount students since the races overlap -- so we'll only use the race for disparities, not totals
ct_school_race <- ct_attendance %>%
  filter(!is.na(race)) %>%
  right_join(
    ct_ACS_race_5to14,
    by = c("FIPS", "race")
  ) %>%
  select(-race)

table(ct_school_race$Race)

write.csv(ct_school, file = "processed/ct_school.csv", row.names = FALSE)
write.csv(ct_school_race, file = "processed/ct_school_race.csv", row.names = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Baseline attendance rates by race for supplemental table
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Restrict the calculation to the census tracts included
# in the dust analysis
analysis_fips <- read.csv(
  "processed/scenario_pm_deltas_daily.csv",
  stringsAsFactors = FALSE
) %>%
  distinct(FIPS) %>%
  mutate(FIPS = as.character(FIPS))

attendance_baseline_by_race <- ct_school_race %>%
  mutate(FIPS = as.character(FIPS)) %>%
  semi_join(
    analysis_fips,
    by = "FIPS"
  ) %>%
  filter(
    Race %in% c("Pacific Islander", "Hispanic/Latino", "White", "AfAm/Black", "Asian")
  ) %>%
  distinct(
    FIPS,
    Race,
    age_group,
    pop,
    incidence_rate_daily
  ) %>%
  group_by(Race) %>%
  summarise(
    student_population =
      sum(pop, na.rm = TRUE),
    
    student_population_with_attendance =
      sum(
        pop[!is.na(incidence_rate_daily)],
        na.rm = TRUE
      ),
    
    baseline_incidence_daily =
      weighted.mean(
        incidence_rate_daily,
        pop,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Share of the race-specific student population for which
    # a district-race attendance rate was available
    attendance_data_coverage =
      student_population_with_attendance /
      student_population,
    
    # incidence_rate_daily was defined as annual expected
    # school-loss days divided by 365
    baseline_school_loss_days_annual =
      baseline_incidence_daily * 365,
    
    # Convert annual school-loss days back to the probability
    # of absence on a school day
    baseline_absence_rate =
      baseline_school_loss_days_annual / 180,
    
    baseline_attendance_rate =
      1 - baseline_absence_rate,
    
    baseline_attendance_percent =
      100 * baseline_attendance_rate
  ) %>%
  select(
    Race,
    attendance_data_coverage,
    baseline_attendance_rate,
    baseline_attendance_percent,
    baseline_school_loss_days_annual,
    baseline_incidence_daily
  )

print(attendance_baseline_by_race)
