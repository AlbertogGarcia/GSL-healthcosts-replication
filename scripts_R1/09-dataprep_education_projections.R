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
               ggplot2
)

options(scipen=999)  # turn off sci notation
options(dplyr.summarise.inform = FALSE)  # turn off dplyr group by comments
options(java.parameters = "-Xmx8000m") 
`%ni%` <- Negate(`%in%`)  # "not in" function

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### School Loss Days
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACS_year_final = 2023

# Implicitly going to apply 5-17 growth rate to the 5-14 population used in the school-loss analysis.

county_pop_projections <-
  readxl::read_excel(
    "data/population/GardnerPI_projections.xlsx",
    sheet = 2
  ) %>%
  filter(str_detect(Geography, " County$")) %>%
  select(Geography, Year, `Population Ages 5-17`) %>%
  transmute(
    County = str_remove(Geography, " County$"),
    Year = as.integer(Year),
    start.age = 5,
    end.age = 17,
    projected_pop = `Population Ages 5-17`
  ) %>%
  filter(Year >= ACS_year_final) %>%
  arrange(County, Year) %>%
  group_by(County, start.age, end.age) %>%
  mutate(
    baseline_projected_pop = projected_pop[Year == ACS_year_final][1],
    # Direct population index relative to 2023
    cum_growth = projected_pop / baseline_projected_pop
  ) %>%
  ungroup() %>%
  select(County, Year, start.age, end.age, cum_growth)




ct_school <- read.csv("processed/ct_school.csv", stringsAsFactors =  FALSE)

ct_school_projections <- ct_school %>%
  fuzzyjoin::fuzzy_left_join(county_pop_projections,
                             by = c("County" = "County",
                                    "lower_age" = "start.age" ,
                                    "upper_age" = "end.age"),
                             match_fun = list(`==`,`>=`, `<=`))%>%
  mutate(County = dplyr::coalesce(County.x,County.y),
         pop = pop*cum_growth)%>%
  select(-County.x,-County.y, -cum_growth, -start.age, -end.age)

write.csv(ct_school_projections, file = "processed/ct_school_projections.csv", row.names = FALSE)

