# GSL dust costs: Mortality and asthma by race data processing
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
               cowplot,
               sjmisc
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
#### MORTALITY
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# census-tract population demographics by age AND RACE

agebyrace_sheet = 11
ct_ACS_agebyrace_wide <- readxl::read_excel("data/population/Demographic_Raw_Tables.xlsx", sheet = agebyrace_sheet)

ct_ACS_agebyrace <- ct_ACS_agebyrace_wide %>%
  pivot_longer(3:ncol(.), names_to = "agebyrace_group", values_to = "pop")%>%
  filter(!str_detect(agebyrace_group, "Total"))%>%
  mutate(agebyrace_group = str_remove_all(agebyrace_group, "[()]"),
         race = case_when(
           str_detect(agebyrace_group, "White") ~ "White NH",
           str_detect(agebyrace_group, "Hispanic") ~ "Hispanic or Latino",
           str_detect(agebyrace_group, "American Indian and Alaska Native") ~ "American Indian and Alaska Native Alone",
           str_detect(agebyrace_group, "Black") ~ "Black or African American Alone",
           str_detect(agebyrace_group, "Asian") ~ "Asian alone",
           str_detect(agebyrace_group, "Hawaiian and PI") ~ "Hawaiian and PI Alone",
           str_detect(agebyrace_group, "Other Race") ~ "Other Race Alone",
           str_detect(agebyrace_group, "Multi Racial") ~ "Multi Racial"
  ),
  age_group = str_remove(agebyrace_group, paste0(race, " ")),
  age_group = case_when(
    age_group %in% c("15 to 17 years", "18 and 19 years", "20 to 24 years") ~ "15 to 24 years",
    age_group %in% c("25 to 29 years", "30 to 34 years") ~ "25 to 34 years",
    age_group %ni% c("25 to 29 years", "30 to 34 years", "15 to 17 years", "18 and 19 years", "20 to 24 years") ~ age_group
  )
  ) %>%
  group_by(FIPS, County, age_group, race)%>%
  summarise(pop = sum(pop, na.rm = T))%>%
  ungroup %>%
  mutate(lower_age = str_sub(age_group,1,2),
  lower_age = as.numeric(ifelse(lower_age %in% "Un",0,lower_age)),
  age_group = str_remove(age_group," years"),
  upper_age = str_sub(age_group,-2,-1),
  upper_age = as.numeric(ifelse(upper_age %in% "er",100,upper_age)),
  upper_age = ifelse(upper_age==5,4,upper_age) ## What it actually should be (under 5 year old)
  )

# get county by age population totals
county_agebyrace_pop_ACS <- ct_ACS_agebyrace %>%
  group_by(County, lower_age, upper_age, age_group, race)%>%
  summarise(pop = sum(pop, na.rm = T))%>%
  ungroup

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot SL county age distributions
map_A <- c(
  "American Indian and Alaska Native Alone" = "American Indian",
  "Asian alone" = "Asian",
  "Black or African American Alone" = "Black/African American",
  "Hawaiian and PI Alone" = "Pacific Islander",
  "Hispanic or Latino" = "Hispanic/Latino",
  "White NH" = "White",
  "Multi Racial" = "Multiracial",
  "Other Race Alone" = "Other"
)

ct_ACS_agebyrace %>%
  filter(County %in% c("Salt Lake", 
                       "Weber", 
                       "Davis" 
                     #  "Utah" ,
                      # "Box Elder"
                       )) %>%
  group_by(County, lower_age, upper_age, age_group, race) %>%
  summarise(pop = sum(pop, na.rm = T)) %>%
  ungroup %>%
  group_by(race) %>%
  mutate(prop_age = pop/sum(pop, na.rm = T),
         race_std = recode(race, !!!map_A),
        age_group = ifelse(age_group == "85 and over", "85+", age_group)) %>%
  ungroup %>%
  filter(race_std %in% c("Pacific Islander", "Hispanic/Latino" , "White"
                    # "Asian",
                     #"Black/African American"
                    )) %>%
ggplot(aes(x = reorder(age_group, lower_age), y = prop_age, fill = race_std))+
  geom_bar(stat = 'identity', position="dodge")+
  scale_fill_manual(values = c(palette$orange, 
                                palette$red, 
                                palette$blue 
                               # palette$green, 
                              #  palette$dark
                               )
  )+
  xlab("Age group")+ylab("Proportion of population")+
  ggtitle("Age distributions across racial groups")+
  theme_cowplot()+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "top",
        legend.justification = "center")
ggsave("figs_pdf/race_ages.pdf", width = 9, height = 6.5)


# geom_histogram()# Get county level incidence data
county_incidence_agebyrace <- read.delim(
  "data/health/mortality/CDC_wonder/2018to23/county_race_age10yr.txt"
) %>%
  mutate(
    age_group = ifelse(
      Ten.Year.Age.Groups %in% c("< 1 year", "1-4 years"),
      "Under 5 years",
      Ten.Year.Age.Groups
    )
  ) %>%
  filter(
    age_group != "NS", 
    is.na(Notes) | Notes != "Total", 
    !is.na(County.Code)
  ) %>%
  # Construct the final mutually exclusive race/ethnicity categories
  mutate(
    race = case_when(
      Hispanic.Origin == "Hispanic or Latino" ~
        "Hispanic or Latino",
      Single.Race.6.Code == "NHOPI" ~
        "Hawaiian and PI Alone",
      Single.Race.6.Code == "A" ~
        "Asian alone",
      Single.Race.6.Code == "M" ~
        "Multi Racial",
      Single.Race.6.Code == "1002-5" ~
        "Black or African American Alone",
      Single.Race.6.Code == "2054-5" ~
        "American Indian and Alaska Native Alone",
      Single.Race.6.Code == "2106-3" ~
        "White NH",
      TRUE ~ "Other Race Alone"
    ),
    County = str_remove(County, " County, UT")
  ) %>%
  filter(County %in% county_agebyrace_pop_ACS$County) %>%
  # Combine the single-race rows that were recoded into the
  # common Hispanic or Latino category
  group_by(County,age_group,race) %>%
  summarise(
    Deaths = if (
      all(is.na(Deaths))
    ) {
      NA_real_
    } else {
      sum(Deaths, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  mutate(
    age_group = str_replace(age_group, "-", " to "),
    age_group = str_remove(age_group, " years"),
    lower_age = str_sub(age_group, 1, 2),
    lower_age = as.numeric(
      ifelse(lower_age == "Un", 0, lower_age)
    ),
    upper_age = str_sub(age_group, -2, -1),
    upper_age = as.numeric(ifelse(lower_age == 85, 100, upper_age)),
    upper_age = ifelse(upper_age == 5,4,upper_age)
  ) %>%
  left_join(
    county_agebyrace_pop_ACS,
    by = c("County", "race", "lower_age", "upper_age"),
    relationship = "many-to-one"
  ) %>%
  mutate(
    Deaths_annual = Deaths / 6,
    incidence_rate = Deaths_annual / pop
  ) %>%
  select(County, upper_age, lower_age, race, incidence_rate)

# match back to census tracts
ct_incidence_mortality <- ct_ACS_agebyrace %>%
  full_join(county_incidence_agebyrace, by = c("lower_age", "upper_age", "race", "County"))%>%
  mutate(incidence_rate = tidyr::replace_na(incidence_rate, 0),
         incidence_rate_daily = incidence_rate/365,
         endpoint = "Mortality, All Cause") %>%
  filter(pop>0)



write.csv(ct_incidence_mortality, file = "processed/ct_incidence_mortality_race.csv", row.names = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Baseline mortality rates by race for supplemental table
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseline_mortality_inputs <- ct_incidence_mortality %>%
  filter(endpoint == "Mortality, All Cause") %>%
  distinct(
    FIPS, County, race,
    age_group, lower_age, upper_age,
    pop, incidence_rate
  ) %>%
  filter(pop > 0)

baseline_mortality_crude <- baseline_mortality_inputs %>%
  group_by(race) %>%
  summarise(
    population = sum(pop, na.rm = TRUE),
    expected_baseline_deaths =
      sum(pop * incidence_rate, na.rm = TRUE),
    baseline_mortality_per_100k =
      100000 * expected_baseline_deaths / population,
    .groups = "drop"
  )

# Common age distribution for the entire study population
standard_age_weights <- baseline_mortality_inputs %>%
  group_by(lower_age, upper_age, age_group) %>%
  summarise(
    standard_population = sum(pop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    standard_age_weight =
      standard_population / sum(standard_population)
  )

# Population-weighted mortality rate within each race-age group
race_age_mortality_rates <- baseline_mortality_inputs %>%
  group_by(race, lower_age, upper_age, age_group) %>%
  summarise(
    population = sum(pop, na.rm = TRUE),
    expected_baseline_deaths =
      sum(pop * incidence_rate, na.rm = TRUE),
    mortality_rate =
      expected_baseline_deaths / population,
    .groups = "drop"
  )

baseline_mortality_standardized <-
  race_age_mortality_rates %>%
  left_join(
    standard_age_weights,
    by = c("lower_age", "upper_age", "age_group")
  ) %>%
  group_by(race) %>%
  summarise(
    age_standardized_mortality_per_100k =
      100000 * sum(
        standard_age_weight * mortality_rate,
        na.rm = TRUE
      ),
    age_weight_covered =
      sum(standard_age_weight[!is.na(mortality_rate)]),
    .groups = "drop"
  )
