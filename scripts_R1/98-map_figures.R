# GSL dust costs: Figures using census tracts or other gis data
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
               tmap,
               basemaps,
               terra
)

options(scipen=999)  # turn off sci notation
options(dplyr.summarise.inform = FALSE)  # turn off dplyr group by comments
options(java.parameters = "-Xmx8000m") 
`%ni%` <- Negate(`%in%`)  # "not in" function

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Load and organize data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
current_scenario = 4192
map_baseline = 4203

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get relevant GSL lake elevation polygons
GSL_geotif <- terra::rast("data/gis/Great_Salt_Lake_TBDEM_100m_bilinear.tif")
## from-to-becomes
# classify the values into three groups 
# all values >= 0 and <= 0.25 become 1, etc.
m <- c(0, 4192, 4192,
       4192, map_baseline, map_baseline,
       map_baseline, Inf, NA
)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
GSL_rcl <- terra::as.polygons(classify(GSL_geotif, rclmat, include.lowest=FALSE))%>%
  st_as_sf %>%
  mutate(elevation = case_when(Band_1 == map_baseline ~ "Long-term avg.",
                               Band_1 == 4192 ~ "Current level")
  )

m <- c(0, 4198, 4198,
       4198, map_baseline, map_baseline,
       map_baseline, Inf, NA
)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
GSL_4198 <- terra::as.polygons(classify(GSL_geotif, rclmat, include.lowest=FALSE))%>%
  st_as_sf %>%
  mutate(elevation = case_when(Band_1 == map_baseline ~ "Long-term avg.",
                               Band_1 == 4198 ~ "Healthy lake: 4198 ftASL")
  )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cost results and integrate with gis files for tracts

ct_totals <- read.csv("processed/ct_totals.csv", stringsAsFactors =  FALSE) %>%
  mutate(costs_per_capita_pm = costs_per_capita/pm_delta)

scenarios_map <- read.csv("processed/scenarios_map.csv", stringsAsFactors =  FALSE) 

fips_county <- read.csv("data/gis/fips_codes/county_fips_master.csv")%>%
  filter(state_abbr %in% c("UT"))%>%
  mutate(County = str_remove(county_name, " County"),
         fips_county = as.character(county)
         )%>%
  select(County, fips_county)

tracts.shp <- read_sf("data/gis/tracts/CensusTracts2020.shp")%>%
  mutate(FIPS = as.double(GEOID20),
         fips_county = str_remove(COUNTYFP20, "^0+")
  )%>%
  select(FIPS, fips_county)

utah.shp <- tracts.shp %>% st_union() %>% st_as_sf()

ct_costs.shp <- tracts.shp %>%
  left_join(fips_county, by = "fips_county")%>%
  select(FIPS)%>%
  left_join(ct_totals, by = "FIPS")%>% 
  # filter(County %in% ct_totals$County)%>%
  mutate_at(vars(pm_delta:ncol(.)), ~tidyr::replace_na(., 0))%>%
  st_as_sf()

ct_pm10.shp <- tracts.shp %>%
  left_join(fips_county, by = "fips_county")%>%
  select(FIPS)%>%
  left_join(scenarios_map, by = "FIPS")%>% 
  mutate_at(vars(pm10_delta), ~tidyr::replace_na(., 0))%>%
  st_as_sf()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extra administrative boundaries

Utah_munic <- read_sf("data/gis/UtahMunic/Municipalities_Carto.shp")

SLC <- Utah_munic %>% # location of salt lake city
  filter(NAME %in% c("Salt Lake City"))

cities2 <- Utah_munic %>% # other relevant cities worth including on map
  filter(NAME %in% c("Logan", "Provo"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Backgrounds and extents
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

get_maptypes()

e_bg <- as(raster::extent(-12595000, -12370000, 4890000, 5125000), "SpatialPolygons")%>%
  st_as_sf()%>%
  st_set_crs(3857)
e_bg

bg <- basemaps::basemap_terra(ext= e_bg,
                                map_service = "carto", map_type = "light_no_labels"
)


bg_labels <- basemaps::basemap_terra(ext=ct_costs.shp, 
                                     map_service = "carto", map_type = "dark_only_labels"
)

ct_costs_current <- ct_costs.shp %>%
  filter(scenario == current_scenario | is.na(scenario)) %>%
  st_transform(terra::crs(bg)) %>% 
  st_as_sf()

ct_pm10_current <- ct_pm10.shp %>%
  filter(scenario == current_scenario | is.na(scenario)) %>%
  st_transform(terra::crs(bg)) %>% 
  st_as_sf()

hist(ct_costs_current$pm_delta)

ct_costs_4198 <- ct_costs.shp %>%
  filter(scenario == 4198 | is.na(scenario)) %>%
  st_transform(terra::crs(bg)) %>% 
  st_as_sf()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PM map
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct_pm10_current <- ct_pm10_current %>%
  mutate(pm10_plot = pmin(pm10_delta, 100))
ct_costs_current <- ct_costs_current %>%
  mutate(plot_costs_per_capita = pmin(costs_per_capita, 100))
# pm_palette <- c(
#   "#fdd49e",
#   "#fc8d59",
#   "#d7301f",
#   "#56000a"
# )

pm_breaks <- c(0, 20, 50, 100)

pm_labels <- c(
  "0",
  "",
  "",
  "100 or more"
)


pm_map <- tm_shape(e_bg) + tm_fill("#f5f9f9")+
  #tm_shape(bg) + tm_rgb(alpha = 0.9)+#tm_rgb()+
  tm_shape(ct_pm10_current) +
  tm_polygons(col = "pm10_delta",
          title = "Avg. event PM10 (\u00b5g/m\u00b3)",
          style = "cont",
          lwd = 0.2,
          palette = "YlOrRd",
          breaks = pm_breaks,
          labels = pm_labels
          )+
  tm_shape(GSL_rcl) + tm_fill("elevation",
                              title = "GSL elevation",
                              palette = c("#64a2b8", "#d5c5b5")
                              ) +
  tm_legend(frame = T,
            bg.color = "#f9f9f9")+
  tm_shape(SLC %>% st_centroid()) + tm_text("NAME", size = 1.2, col = "black", fontface = "bold")+
  tm_shape(cities2 %>% st_centroid()) + tm_text("NAME", size = 0.8, col = "black")+
  tm_layout(legend.outside = FALSE,
            legend.position = c("right", "top"),
            inner.margins = c(0,0,0,0),
            outer.margins = c(0,0,0,0))
pm_map
tmap_save(filename = "figs_pdf/pm_map_current.pdf", width = 6, height = 7)

tm_shape(GSL_rcl) + tm_fill("elevation",
                             title = "GSL elevation",
                             palette = c("#64a2b8", "#cdb9a6"))+
  tm_layout(legend.outside = FALSE,
            legend.position = c("right", "top"),
            inner.margins = c(0,0,0,0),
            outer.margins = c(0,0,0,0)) 
tmap_save(filename = "figs_pdf/lake_current.pdf", width = 5, height = 6)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Health Cost map
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cost_labels <- c(
  "0",
  "",
  "",
  "$100 or more"
)

CostPerCapita_map <- tm_shape(e_bg) + tm_fill("#f5f9f9")+
  #tm_shape(bg) + tm_rgb(alpha = 0.9)+#tm_rgb()+
  tm_shape(ct_costs_current) + 
  tm_polygons(col = "plot_costs_per_capita",
              title = "Dust-induced health\ncosts per capita",
              style = "cont",
              lwd = 0.2,
              pal = "YlOrRd",
              breaks = c(0, 20, 45, 100),
              labels = cost_labels
  )+
  tm_shape(GSL_rcl) + tm_fill("elevation",
                              palette = c("#64a2b8", "#cdb9a6"),
                              legend.show = FALSE
  ) + 
  tm_legend(frame = T,
            bg.color = "#f9f9f9")+
  # tm_shape(SLC %>% st_centroid()) + tm_text("NAME", size = 3/4, col = "white", fontface = "bold")+
  # tm_shape(cities2 %>% st_centroid()) + tm_text("NAME", size = 2/3, col = "black")+
  tm_layout(legend.outside = FALSE,
            legend.position = c("right", "top"),
            inner.margins = c(0,0,0,0),
            outer.margins = c(0,0,0,0)) +
  tm_scale_bar(breaks = c(0, 25, 50))+
  tm_compass(type = "4star", size = 1.75, position = c("right", "top"))
CostPerCapita_map
tmap_save(filename = "figs_pdf/cost_map_current.pdf", width = 6, height = 7)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Health costs per PM exposure
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary((ct_costs_current%>% filter(costs_per_capita_pm>0))$costs_per_capita_pm)

cost_by_delta_map <- tm_shape(e_bg) + tm_fill("#f5f9f9")+
  #tm_shape(bg) + tm_rgb(alpha = 0.9)+#tm_rgb()+
  tm_shape(ct_costs_current) + 
  tm_polygons(col = "costs_per_capita_pm",
              title = "Ratio of per capita costs\nto mg/m3 PM10 exposure",
              style = "cont",
              pal = "Purples",
              alpha = 0.75,
              lwd = 0.15,
              breaks = c(0, 1, 2)
              #labels = c(0, 50, 100, "$150 or more")
  )+
  tm_shape(GSL_rcl) + tm_fill("elevation",
                              palette = c("#64a2b8", "#cdb9a6"),
                              legend.show = FALSE
  ) + 
  tm_legend(frame = T,
            bg.color = "#f9f9f9")+
  tm_shape(SLC %>% st_centroid()) + tm_text("NAME", size = 1.2, col = "black", fontface = "bold")+
  tm_shape(cities2 %>% st_centroid()) + tm_text("NAME", size = 3/4, col = "black", fontface = "bold")+
  tm_layout(legend.outside = FALSE,
            legend.position = c("right", "top"),
            inner.margins = c(0,0,0,0),
            outer.margins = c(0,0,0,0)) +
  tm_scale_bar(breaks = c(0, 10, 25))+
  tm_compass(type = "4star", size = 1.25, position = c("right", "top"))
cost_by_delta_map
tmap_save(filename = "figs_pdf/cost_by_delta.pdf", width = 6, height = 7)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Inset and combined maps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# All basemaps packages return maps in EPSG:3857.
utah_3857 <- st_transform(utah.shp, 3857)

basemaps::get_maptypes()

# 2. Fetch the terrain basemap using free Esri services
bg_utah <- basemaps::basemap_terra(ext = utah_3857, 
                                   map_service = "esri", 
                                   map_type = "world_terrain_base"
                                   ) %>%
  terra::mask(utah_3857)

bg_ref <- basemaps::basemap_terra(ext = utah_3857,
                                   map_service = "esri",
                                   map_type = "world_reference_overlay"
                                    #"world_light_gray_reference"
) %>%
  terra::mask(utah_3857)


# e <- as(raster::extent(-12522831, -12340605, 4833572, 5076642), "SpatialPolygons")%>%
#   st_as_sf()%>%
#   st_set_crs(3857)
#proj4string(e) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

utah_inset  <- tm_shape(bg_utah) + tm_rgb()+
  tm_shape(bg_ref) + tm_rgb()+
  tm_shape(utah.shp) + tm_borders(col = "black")+
  tm_shape(e_bg) + tm_borders(lwd = 2, col = "red")+
  tm_layout(legend.outside = FALSE,
            inner.margins = c(0,0,0,0),
            outer.margins = c(0,0,0,0),
            frame = FALSE,
            bg.color = "transparent") 
utah_inset




#%%%%%%%%%%%%%%%%%%%%%%%%
#### Combine all maps into single pdf
#%%%%%%%%%%%%%%%%%%%%%%%%

library(grid)

pdf("figs_pdf/map.pdf", width = 13.25, height = 7, useDingbats = FALSE)

grid.newpage()

# Main maps: each gets half of the page, with a small gap
print(pm_map, vp = viewport(x = 0.25, y = 0.50, width = 0.495))

print(CostPerCapita_map, vp = viewport(x = 0.75, y = 0.50, width = 0.495))

print(utah_inset, vp = viewport(x = 0.06, y = 0.15, height = 0.25))

dev.off()

