# GSL dust costs: Morbidity and population projections
# albert.garcia@utah.edu
# created: 01/07/2026
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
#### MORTALITY
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACS_year_final = 2023

ct_incidence_morbidity <- read.csv("processed/ct_incidence_morbidity.csv", stringsAsFactors =  FALSE)
table(ct_incidence_morbidity$age_group)

# Projected county-level population growth by age
gardner_pop_projections <- readxl::read_excel("data/population/GardnerPI_projections.xlsx", sheet = 2)
colnames(gardner_pop_projections)


# Match morbidity age groups to projected population age groups
county_pop_projections <- gardner_pop_projections %>%
  select(
    Geography,
    Year,
    `Population Ages 0-4`:`Population Ages 85+`
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
  arrange(
    County,
    start.age,
    end.age,
    Year
  ) %>%
  group_by(
    County,
    start.age,
    end.age
  ) %>%
  mutate(
    baseline_projected_pop =
      projected_pop[Year == ACS_year_final][1],
    cum_growth =
      projected_pop / baseline_projected_pop
  ) %>%
  ungroup() %>%
  select(
    County,
    Year,
    start.age,
    end.age,
    cum_growth
  )

ct_incidence_projections <- data.frame()

for (e in unique(ct_incidence_morbidity$endpoint)) {
  
  print(e)
  
  this_endpoint_projection <- ct_incidence_morbidity %>%
    filter(endpoint == e) %>%
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
      County = dplyr::coalesce(
        County.x,
        County.y
      ),
      pop = pop * cum_growth
    ) %>%
    select(
      -County.x,
      -County.y,
      -start.age,
      -end.age
    )
  
  ct_incidence_projections <- bind_rows(
    ct_incidence_projections,
    this_endpoint_projection
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Save separate projection file for each endpoint
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endpoint_file_names <- ct_incidence_projections %>%
  distinct(endpoint) %>%
  mutate(
    endpoint_file = endpoint %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", "_") %>%
      str_remove_all("^_|_$")
  )

dir.create(
  "processed/morbidity_projection_endpoints",
  showWarnings = FALSE,
  recursive = TRUE
)

for (i in seq_len(nrow(endpoint_file_names))) {
  
  this_endpoint <- endpoint_file_names$endpoint[i]
  this_file <- endpoint_file_names$endpoint_file[i]
  
  endpoint_data <- ct_incidence_projections %>%
    filter(endpoint == this_endpoint)
  
  saveRDS(
    endpoint_data,
    file = file.path(
      "processed/morbidity_projection_endpoints",
      paste0(this_file, ".rds")
    ),
    compress = TRUE
  )
  
  rm(endpoint_data)
  gc()
}

write.csv(
  endpoint_file_names,
  "processed/morbidity_projection_endpoint_files.csv",
  row.names = FALSE
)
