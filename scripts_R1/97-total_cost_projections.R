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
               cowplot,
               ggpubr,
               purrr,
               ggplot2
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

scenario_descrip <- c("No conservation",
                      "Current lake    ",
                      "Minimum healthy\nlake          ",
                      "Recent historical\naverage      "
)

scenario_pal <- c(palette$bad, palette$current, palette$target, palette$avg)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Annual
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### Mortality
total_mortality <- read.csv("processed/total_mortality.csv", stringsAsFactors =  FALSE) %>%
  filter(endpoint == "Mortality, All-cause",
         scenario %in% relevant_scenarios)%>%
  mutate(costs = costs_VSL*1000000,
         costs_age = costs_age_VSL*1000000,
         costs_lower = costs_VSL_lower*1000000,
         costs_upper = costs_VSL_upper*1000000,
         costs_linear = costs_VSL_linear*1000000,
         costs_nonlinear = costs_VSL_nonlinear*1000000)%>%
  select(scenario, costs, costs_age, costs_linear, costs_nonlinear, costs_lower, costs_upper) %>%
  mutate(endpoint_category = "Mortality")

#### School Loss Days
total_schoolloss <- read.csv("processed/total_schoolloss.csv", stringsAsFactors =  FALSE) %>%
  mutate(costs = costs_SLD,
         costs_age = costs_SLD, # No age varying cost of school loss, so just giving standard
         costs_lower = costs_SLD_lower,
         costs_upper = costs_SLD_upper,
         costs_linear = costs_SLD,
         costs_nonlinear = costs_SLD)%>%
  select(scenario, costs, costs_age, costs_linear, costs_nonlinear, costs_lower, costs_upper)

#### All morbidity
total_morbidity <- read.csv("processed/total_morbidity.csv", stringsAsFactors =  FALSE) %>%
  group_by(scenario) %>%
  summarise(costs = sum(costs, na.rm = T),
            costs_age = sum(costs, na.rm = T), # no age varying COI
            costs_lower = sum(costs_lower, na.rm = T),
            costs_upper = sum(costs_upper, na.rm = T),
            costs_linear = sum(costs, na.rm = T),
            costs_nonlinear = sum(costs, na.rm = T))%>%
  ungroup %>%
  rbind(total_schoolloss) %>%
  group_by(scenario) %>%
  summarise(costs = sum(costs, na.rm = T),
            costs_age = sum(costs_age, na.rm = T),
            costs_lower = sum(costs_lower, na.rm = T),
            costs_upper = sum(costs_upper, na.rm = T),
            costs_linear = sum(costs_linear, na.rm = T),
            costs_nonlinear = sum(costs_nonlinear, na.rm = T))%>%
  ungroup %>%
  mutate(endpoint_category = "Morbidity")

totals_df <- total_morbidity %>%
  rbind(total_mortality) %>%
  mutate(costs_millions = costs / 1e6,
         costs_millions_lower = costs_lower / 1e6,
         costs_millions_upper = costs_upper / 1e6,
         costs_millions_age = costs_age / 1e6,
         costs_millions_linear = costs_linear / 1e6,
         costs_millions_nonlinear = costs_nonlinear / 1e6) %>%
  group_by(scenario) %>%
  summarise(total_costs = sum(costs_millions),
            total_costs_lower = sum(costs_millions_lower),
            total_costs_upper = sum(costs_millions_upper),
            total_costs_age = sum(costs_millions_age),
            total_costs_linear = sum(costs_millions_linear),
            total_costs_nonlinear = sum(costs_millions_nonlinear)
            ,.groups = "drop")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Manuscript annual totals
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annual_manuscript_totals <- totals_df %>%
  left_join(
    total_mortality %>%
      transmute(
        scenario,
        mortality_costs_m = costs / 1e6
      ),
    by = "scenario"
  ) %>%
  left_join(
    total_morbidity %>%
      transmute(
        scenario,
        morbidity_costs_m = costs / 1e6
      ),
    by = "scenario"
  ) %>%
  mutate(
    scenario_label = recode(
      as.character(scenario),
      "4183" = "No conservation",
      "4192" = "Current lake",
      "4198" = "Minimum healthy lake",
      "4200" = "Historical average"
    )
  ) %>%
  arrange(scenario)

current_annual_costs_m <- annual_manuscript_totals %>%
  filter(scenario == current_scenario) %>%
  pull(total_costs)

annual_manuscript_totals <- annual_manuscript_totals %>%
  mutate(
    change_from_current_m = total_costs - current_annual_costs_m,
    reduction_from_current_pct =
      100 * (current_annual_costs_m - total_costs) / current_annual_costs_m
  )

write.csv(
  annual_manuscript_totals,
  "processed/manuscript_annual_totals.csv",
  row.names = FALSE
)









total_annual_costs_plot <- total_morbidity %>%
  rbind(total_mortality) %>%
  mutate(costs_millions = costs/1000000)%>%
  ggplot(aes(x=reorder(scenario, scenario, order = T), y=costs_millions, fill = endpoint_category)
  )+
  geom_bar(stat='identity')+
  geom_text(
    data = totals_df,
    aes(
      x = reorder(scenario, scenario, order = TRUE),
      y = total_costs,
      label = round(total_costs, 2)
    ),
    inherit.aes = FALSE,
    hjust = -0.25,
    size = 4
  )+
  scale_y_continuous(limits = c(0, max(totals_df$total_costs)*1.07))+
  scale_x_discrete(labels = scenario_descrip)+
  geom_hline(yintercept = 0, linewidth = 0.25)+
  ggtitle("Current annual dust-induced health costs") + 
  xlab("GSL water level (ftASL)")+
  ylab("Costs (millions USD)") +
  scale_fill_manual(values = c(palette$blue, palette$red, palette$green))+
  theme_cowplot(14)+
  theme(legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #panel.grid.major.x = element_line(color = "gray80", size = 0.25),
        plot.title = element_text(hjust = 0.5)
  )+
  coord_flip()+
  guides(fill = guide_legend(title = NULL, reverse=T, hjust = 0.5))
total_annual_costs_plot
ggsave("figs_pdf/costs_mortality_vs_morbidity.pdf",
       width = 8, height = 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Sensitivity Analysis
# - age-based VSL
# - 5th and 95th pctl of conc-response function
# - linear CRF for mortality
# - nonlinear CRF for mortality
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annual_costs_sensitivity <- totals_df %>%
  select(
    scenario,
    `Main estimate` = total_costs,
    `Age-adjusted VSL` = total_costs_age,
    `All CRF 95% CI lower bounds` = total_costs_lower,
    `All CRF 95% CI upper bounds` = total_costs_upper,
    `Linear mortality CRF` = total_costs_linear,
    `Nonlinear mortality CRF` = total_costs_nonlinear
  ) %>%
  pivot_longer(
    cols = -scenario,
    names_to = "sensitivity",
    values_to = "annual_costs_millions"
  ) %>%
  mutate(
    sensitivity = factor(
      sensitivity,
      levels = c(
        "Main estimate",
        "Age-adjusted VSL",
        "All CRF 95% CI lower bounds",
        "All CRF 95% CI upper bounds",
        "Linear mortality CRF",
        "Nonlinear mortality CRF"
      )
    ),
    scenario = factor(
      scenario,
      levels = relevant_scenarios
    )
  )

costs_sensitivity_plot <- annual_costs_sensitivity %>%
  ggplot(
    aes(
      x = sensitivity,
      y = annual_costs_millions,
      fill = scenario
    )
  ) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = scenario_pal,
    name = "Lake elevation scenario",
    labels = scenario_descrip
  ) +
  scale_y_continuous(
    name = "Annual health costs (USD millions)",
    labels = scales::label_number(
      accuracy = 1,
      big.mark = ","
    ),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  labs(
    x = NULL,
    title = paste0(
      "Sensitivity of annual health-cost estimates\n",
      "to alternative parameterizations"
    )
  ) +
  coord_flip() +
  theme_cowplot(12) +
  theme(
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x =
      element_line(
        color = "gray70",
        linewidth = 0.25
      ),
    panel.grid.minor.x =
      element_line(
        color = "gray84",
        linewidth = 0.25
      )
  ) +
  guides(
    fill = guide_legend(
      reverse = TRUE,
      title.hjust = 0.5
    )
  )

costs_sensitivity_plot

ggsave("figs_pdf/costs_sensitivity.pdf",
       width = 8, height = 6)

sensitivity_lines <- totals_df %>%
  transmute(
    scenario = as.numeric(as.character(scenario)),
    `Main estimate` = total_costs,
    `Cohort-adjusted VSL` = total_costs_age,
    `Linear mortality CRF` = total_costs_linear,
    `Nonlinear mortality CRF` = total_costs_nonlinear
  ) %>%
  pivot_longer(
    cols = -scenario,
    names_to = "specification",
    values_to = "annual_costs"
  ) %>%
  mutate(
    specification = factor(
      specification,
      levels = c(
        "Main estimate",
        "Cohort-adjusted VSL",
        "Linear mortality CRF",
        "Nonlinear mortality CRF"
      )
    )
  )

sensitivity_bounds <- totals_df %>%
  transmute(
    scenario = as.numeric(as.character(scenario)),
    lower = total_costs_lower,
    upper = total_costs_upper
  )

costs_sensitivity_plot <- ggplot() +
  geom_ribbon(
    data = sensitivity_bounds,
    aes(
      x = scenario,
      ymin = lower,
      ymax = upper
    ),
    fill = "grey75",
    alpha = 0.55
  ) +
  geom_line(
    data = sensitivity_lines,
    aes(
      x = scenario,
      y = annual_costs,
      color = specification,
      linewidth = specification,
      linetype = specification
    )
  ) +
  geom_point(
    data = sensitivity_lines,
    aes(
      x = scenario,
      y = annual_costs,
      color = specification
    ),
    size = 2.5
  ) +
  scale_color_manual(
    values = c(
      "Main estimate" = palette$blue,
      "Cohort-adjusted VSL" = palette$green,
      "Linear mortality CRF" = "grey35",
      "Nonlinear mortality CRF" = palette$red
    ),
    name = NULL
  ) +
  scale_linewidth_manual(
    values = c(
      "Main estimate" = 0.8,
      "Cohort-adjusted VSL" = 0.8,
      "Linear mortality CRF" = 0.8,
      "Nonlinear mortality CRF" = 0.8
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "Main estimate" = "solid",
      "Cohort-adjusted VSL" = "dashed",
      "Linear mortality CRF" = "dashed",
      "Nonlinear mortality CRF" = "dashed"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = sort(unique(as.numeric(
      as.character(totals_df$scenario)
    )))
  ) +
  scale_y_continuous(
    labels = scales::label_dollar(
      prefix = "$",
      suffix = "M",
      accuracy = 1
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Great Salt Lake elevation (ft ASL)",
    y = "Annual health costs",
    title = "Sensitivity of annual health-cost estimates",
    subtitle = paste0(
      "Shaded region applies all lower or all upper ",
      "concentration–response bounds"
    )
  ) +
  theme_cowplot(13) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_legend(nrow = 2),
    linewidth = "none",
    linetype = "none"
  )

costs_sensitivity_plot


relative_sensitivity <- totals_df %>%
  transmute(
    scenario,
    `Cohort-adjusted VSL` =
      100 * (total_costs_age / total_costs - 1),
    `All CRF 95% CI lower bounds` =
      100 * (total_costs_lower / total_costs - 1),
    `All CRF 95% CI upper bounds` =
      100 * (total_costs_upper / total_costs - 1),
    `Linear mortality CRF` =
      100 * (total_costs_linear / total_costs - 1),
    `Nonlinear mortality CRF` =
      100 * (total_costs_nonlinear / total_costs - 1)
  ) %>%
  pivot_longer(
    cols = -scenario,
    names_to = "specification",
    values_to = "percent_change"
  )

relative_sensitivity_plot <- ggplot(
  relative_sensitivity,
  aes(
    x = percent_change,
    y = specification,
    color = factor(scenario)
  )
) +
  geom_vline(
    xintercept = 0,
    color = "grey55",
    linewidth = 0.5
  ) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 3
  ) +
  scale_x_continuous(
    labels = scales::label_number(
      accuracy = 1,
      suffix = "%"
    ),
    breaks = seq(-30, 30, 10)
  ) +
  scale_color_manual(
    values = scenario_pal,
    labels = scenario_descrip,
    name = "Lake elevation scenario"
  ) +
  labs(
    x = "Change relative to main estimate",
    y = NULL,
    title = "Sensitivity relative to the main estimate"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title.position = "top",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  )
relative_sensitivity_plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Separate plots by morbidity endpoint

# School Loss Days
total_schoolloss_endpoint <- read.csv("processed/total_schoolloss.csv", stringsAsFactors =  FALSE) %>%
  mutate(endpoint_subcategory = endpoint) %>%
  group_by(scenario, endpoint_subcategory) %>%
  summarise(costs = sum(costs_SLD, na.rm = T),
            morbidity = sum(SLD, na.rm = T),
            costs_lower = sum(costs_SLD_lower, na.rm = T),
            morbidity_lower = sum(SLD_lower, na.rm = T),
            costs_upper = sum(costs_SLD_upper, na.rm = T),
            morbidity_upper = sum(SLD_upper, na.rm = T))%>%
  ungroup

# All morbidity
total_morbidity_endpoint <- read.csv("processed/total_morbidity.csv", stringsAsFactors =  FALSE) %>%
  mutate(endpoint_subcategory = case_when(
    endpoint == "Work Loss Days" ~ endpoint,
    str_detect(endpoint, "HA") ~ "Hospital Admissions",
    str_detect(endpoint, "ER visits") ~ "Emergency Department Visits",
  )) %>%
  group_by(scenario, endpoint_subcategory) %>%
  summarise(costs = sum(costs, na.rm = T),
            morbidity = sum(morbidity, na.rm = T),
            costs_lower = sum(costs_lower, na.rm = T),
            morbidity_lower = sum(morbidity_lower, na.rm = T),
            costs_upper = sum(costs_upper, na.rm = T),
            morbidity_upper = sum(morbidity_upper, na.rm = T))%>%
  ungroup %>%
  rbind(total_schoolloss_endpoint)

endpoint_levels <- unique(total_morbidity_endpoint$endpoint_subcategory)
plot_list <- list()
i = 1

for(i in seq_along(endpoint_levels)){
  
  e <- endpoint_levels[i]
  
  this_endpoint <- total_morbidity_endpoint %>%
    filter(endpoint_subcategory == e) %>%
    mutate(cost_per_morbidity = costs / morbidity)
  
  c <- this_endpoint$cost_per_morbidity[1] / 1e6
  
  # Base plot
  p <- ggplot(
    this_endpoint,
    aes(
      x = reorder(scenario, scenario, order = TRUE),
      y = morbidity,
      fill = as.character(scenario)
    )
  ) +
    geom_bar(stat = "identity") +
    # geom_errorbar(
    #   aes(ymin = morbidity_lower,
    #       ymax = morbidity_upper),
    #   width = 0.25,
    #   linewidth = 0.5
    # ) +
    ggtitle(e) +
    scale_x_discrete(labels = scenario_descrip)+
    scale_y_continuous(
      "Dust-induced occurrences",
      labels = scales::label_comma(),
      sec.axis = sec_axis(~ . * c, name = "Costs (millions USD)")
    ) +
    scale_fill_manual(values = scenario_pal, 
                      name = "Lake elevation scenario",
                      labels = scenario_descrip)+
    theme_cowplot(14) +
    theme(
      plot.margin = margin(b = 20, r = 10),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 12)),
      legend.title.position = "top",
      legend.direction = "horizontal",
      legend.box.just = "center",  
      legend.justification = "center", 
      axis.title.x = element_text(size = 12),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    coord_flip() + 
    guides(fill = guide_legend(reverse = T, hjust = 0.5, title.hjust = 0.5))
  
  # Column logic (2 columns)
  is_left  <- i %% 2 == 1
  is_right <- i %% 2 == 0
  # Row logic (2 rows)
  is_top_row    <- i <= 2
  is_bottom_row <- i >  2
  
  # Remove bottom axis on top row
  
  # if (is_top_row) {
  #   p <- p + theme(
  #     axis.title.x.bottom = element_blank()
  #   )
  # }
  # 
  # # Remove top axis on bottom row
  # if (is_bottom_row) {
  #   p <- p + theme(
  #     axis.title.x.top = element_blank()
  #   )
  # }
  
  # Remove unwanted axes
  if (!is_left) {
    p <- p + theme(
      axis.text.y.left = element_blank()
    )
  }
  
  if (!is_right) {
    p <- p + theme(
      axis.text.y.right = element_blank()
    )
  }
  
  plot_list[[i]] <- p
}

ggarrange(
  plotlist = plot_list,
  ncol = 2, widths = c(1.35, 1.01),
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = TRUE,
  legend = "bottom"
)

ggsave("figs_pdf/morbidity_costs.pdf",
       width = 9, height = 10)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Projections through 2060
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### Mortality
total_mortality_projections <- read.csv("processed/total_mortality_projections.csv", stringsAsFactors =  FALSE) %>%
  filter(scenario %in% relevant_scenarios)%>%
  rename(PV_costs = PV_costs_VSL,
         PV_cum_costs = PV_cum_costs_VSL)%>%
  select(scenario, Year, 
         PV_costs, PV_cum_costs) %>%
  mutate(endpoint_category = "Mortality")


#### School Loss Days
total_schoolloss_projections <- read.csv("processed/total_schoolloss_projections.csv", stringsAsFactors =  FALSE) %>%
  rename(PV_costs = PV_costs_SLD,
         PV_cum_costs = PV_cum_costs_SLD)%>%
  select(scenario, Year, PV_costs, PV_cum_costs) 


#### Morbidity
total_morbidity_projections <- read.csv("processed/total_morbidity_projections.csv", stringsAsFactors =  FALSE) %>%
  group_by(scenario, Year) %>%
  summarise(PV_costs = sum(PV_costs_COI, na.rm = T),
            PV_cum_costs = sum(PV_cum_costs_COI, na.rm = T))%>%
  ungroup %>%
  rbind(total_schoolloss_projections) %>%
  group_by(scenario, Year) %>%
  summarise(PV_costs = sum(PV_costs, na.rm = T),
            PV_cum_costs = sum(PV_cum_costs, na.rm = T))%>%
  ungroup %>%
  mutate(endpoint_category = "Morbidity")


total_costs_proj <- total_morbidity_projections %>%
  rbind(total_mortality_projections) %>%
  group_by(scenario, Year) %>%
  summarise(PV_costs = sum(PV_costs),
            PV_cum_costs = sum(PV_cum_costs))

cum_costs_proj <- total_costs_proj %>%
  filter(Year == 2060)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Manuscript projected totals through 2060
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

projected_manuscript_totals <- cum_costs_proj %>%
  transmute(
    scenario,
    cumulative_pv_costs_m = PV_cum_costs,
    cumulative_pv_costs_b = PV_cum_costs / 1e3,
    scenario_label = recode(
      as.character(scenario),
      "4183" = "No conservation",
      "4192" = "Current lake",
      "4198" = "Minimum healthy lake",
      "4200" = "Historical average"
    )
  ) %>%
  arrange(scenario)

current_cumulative_pv_m <- projected_manuscript_totals %>%
  filter(scenario == current_scenario) %>%
  pull(cumulative_pv_costs_m)

projected_manuscript_totals <- projected_manuscript_totals %>%
  mutate(
    savings_vs_current_m =
      current_cumulative_pv_m - cumulative_pv_costs_m,
    savings_vs_current_b =
      savings_vs_current_m / 1e3,
    reduction_vs_current_pct =
      100 * savings_vs_current_m / current_cumulative_pv_m
  )

write.csv(
  projected_manuscript_totals,
  "processed/manuscript_projected_totals_2060.csv",
  row.names = FALSE
)

costs_proj <- ggplot(data = total_costs_proj,
       aes(x = Year, y = PV_costs, color = as.character(scenario))
)+
  geom_line(linewidth = 1)+
  geom_point(size = 2)+
  scale_y_continuous(name = "Present Value Costs (millions USD)") +
  ggtitle("Projected changes in annual costs")+
  scale_color_manual(values = scenario_pal)+
  scale_fill_manual(values = scenario_pal, guide="none")+
  theme_cowplot(12)+
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "gray80", size = 0.25),
        plot.title = element_text(hjust = 0.5)#, size=16)
  )
costs_proj

cum_costs_bar <- ggplot(data = cum_costs_proj,
       aes(x=reorder(scenario, scenario, order = T), y=PV_cum_costs/1e3, fill = as.character(scenario))
)+
  geom_bar(stat='identity')+
  ggtitle("Cumulative health costs through 2060") + 
  #xlab("GSL water level (ftASL)")+
  ylab("Present Value Costs (USD billions)") +
  scale_fill_manual(values = scenario_pal,
                    name = "Lake elevation scenario",
                     labels = scenario_descrip)+
  theme_cowplot(12)+
  theme(legend.title.position = "top",
    legend.direction = "horizontal",
    legend.box.just = "center",  
    legend.justification = "center", 
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.x = element_line(color = "gray80", size = 0.25)
  )+
  coord_flip()+
  guides(fill = guide_legend(reverse = T, title.hjust = 0.5, hjust = 0.5))
cum_costs_bar
legend <- get_legend(cum_costs_bar)

cum_costs_bar <- cum_costs_bar +
  theme(legend.position = "none")

p1 <- ggarrange(costs_proj, cum_costs_bar,
  ncol = 2,
  labels = c("A", "B")
)

ggarrange(p1, legend,
  ncol = 1,
  heights = c(1, 0.2)
)

ggsave("figs_pdf/total_costs_projected.pdf",
       width = 9, height = 5)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Racial disparities in costs per capita
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map_A <- c(
  "American Indian and Alaska Native Alone" = "American Indian",
  "Asian alone" = "Asian",
  "Black or African American Alone" = "Black/African American",
  "Hispanic or Latino" = "Hispanic/Latino",
  "Hawaiian and PI Alone" = "Pacific Islander",
  "White NH" = "White",
  "Multi Racial" = "Multiracial",
  "Other Race Alone" = "Other"
)
map_B <- c(
  "AfAm/Black" = "Black/African American",
  "Asian" = "Asian",
  "Hispanic/Latino" = "Hispanic/Latino",
  "Pacific Islander" = "Pacific Islander",
  "White" = "White"
)


total_mortality_race <- read.csv("processed/total_mortality_race.csv", stringsAsFactors =  FALSE) %>%
  filter(scenario %in% relevant_scenarios)%>%
  mutate(race_std = recode(race, !!!map_A),
         mortality_costs = costs_VSL*1000000)%>%
  select(scenario, race_std, mortality_costs, population, pm_delta)

total_schoolloss_race <- read.csv("processed/total_schoolloss_race.csv", stringsAsFactors =  FALSE) %>%
  mutate(race_std = recode(Race, !!!map_B)) %>%
  rename(costs = costs_SLD) %>%
  select(scenario, race_std, costs)

unique(total_mortality_race$race_std)
unique(total_schoolloss_race$race_std)

total_morbidity_race <- read.csv("processed/total_morbidity_race.csv", stringsAsFactors =  FALSE) %>%
  mutate(race_std = recode(race, !!!map_A)) %>%
  group_by(scenario, race_std) %>%
  summarise(costs = sum(costs, na.rm = T))%>%
  ungroup %>%
  rbind(total_schoolloss_race) %>%
  group_by(scenario, race_std) %>%
  summarise(morbidity_costs = sum(costs, na.rm = T))%>%
  ungroup 

total_costs_race <- total_mortality_race %>%
  full_join(total_morbidity_race, by = c("race_std", "scenario")) %>%
  mutate(costs = morbidity_costs + mortality_costs,
         costs_per_capita = costs/population,
         costs_per_delta = costs_per_capita/pm_delta) %>%
  drop_na(mortality_costs) %>%
  rename(Race = race_std)

total_costs_race %>%
  ggplot(aes(x=scenario, color = reorder(Race, -costs_per_capita), y=costs_per_capita))+
  geom_line(linewidth = 1)+
  geom_point(size = 2.5)+
  scale_y_continuous(name = "Health costs per capita",
                     labels = label_dollar())+
  scale_x_continuous(
    breaks = relevant_scenarios,
    name = "GSL water level (ftASL)")+
  scale_color_manual(values = c(palette$red, 
                                palette$blue, 
                                palette$orange, 
                                palette$green, 
                                palette$dark)
  )+
  ggtitle("Per-capita health costs by race")+
  theme_cowplot(14)+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.justification = "center"
  )
ggsave("figs_pdf/race_costs.pdf",
       width = 9, height = 6)

total_costs_long <- total_costs_race %>%
  pivot_longer(
    cols = c(morbidity_costs, mortality_costs),
    names_to = "endpoint_category",
    values_to = "cost"
  ) %>%
  mutate(
    damage_type = recode(
      endpoint_category,
      morbidity_costs = "Morbidity",
      mortality_costs = "Mortality"
    ),
    costs = cost/population
  )

total_costs_long %>%
  ggplot(aes(
    x = scenario,
    y = costs,
    color = reorder(Race, -costs_per_capita),
    group = Race
  )) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~ damage_type, nrow = 1) +
  scale_y_continuous(
    name = "Health costs per capita",
    labels = scales::label_dollar()
  ) +
  scale_x_continuous(
    breaks = relevant_scenarios,
    name = "GSL elevation (ftASL)"
  ) +
  scale_color_manual(values = c(
    palette$red,
    palette$blue,
    palette$orange,
    palette$green,
    palette$dark
  )) +
  ggtitle("Per-capita mortality and morbidity costs by race") +
  theme_bw(14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
ggsave("figs_pdf/race_costs_mortmorb.pdf",
       width = 9, height = 6)

cost_plot_current <- total_costs_race %>%
  filter(scenario == 4192) %>%
  ggplot(aes(
    x = reorder(Race, -costs_per_capita),
    y = costs_per_capita
  )) +
  geom_point(size = 5, color = palette$blue) +
  scale_y_continuous(
    name = "Annual per-capita health costs",
    labels = label_dollar()
  ) +
  scale_x_discrete(name = NULL) +
  ggtitle("Health costs by race") +
  theme_bw(13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
        axis.line = element_line(color = "grey50"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(
          angle = 40, hjust = 1, vjust = 1,
          size = 12
        )
  )
cost_plot_current

pm_plot_current <- total_costs_race %>%
  filter(scenario == 4192) %>%
  ggplot(aes(
    x = reorder(Race, -costs_per_capita),
    y = pm_delta
  )) +
  geom_point(size = 5, color = palette$blue) +
  scale_y_continuous(
    name = expression(paste("Pop. weighted PM10 exposure per event (", mu, "g/m"^3, ")"))
  ) +
  scale_x_discrete(name = NULL) +
  ggtitle("Dust exposure by race") +
  theme_bw(13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.25),
    axis.line = element_line(color = "grey50"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 40, hjust = 1, vjust = 1,
      size = 12
    )
  )
pm_plot_current

ggarrange(cost_plot_current, pm_plot_current,
          ncol = 2,
          labels = c("A", "B")
)

ggsave("figs_pdf/race_disparity_dotplot.pdf",
       width = 9, height = 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Output total costs per capita at census-tract for maps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### Mortality
ct_mortality <- read.csv("processed/ct_mortality_map.csv", stringsAsFactors =  FALSE) %>%
  filter(scenario %in% relevant_scenarios)%>%
  mutate(mortality_costs = costs_VSL*1000000)%>%
  select(FIPS, County, scenario, mortality_costs, population, pm_delta)

#### School Loss Days
ct_schoolloss <- read.csv("processed/ct_schoolloss_map.csv", stringsAsFactors =  FALSE) %>%
  rename(costs = costs_SLD)%>%
  select(FIPS, scenario, costs) 

#### All morbidity
ct_morbidity <- read.csv("processed/ct_morbidity_map.csv", stringsAsFactors =  FALSE) %>%
  group_by(FIPS, scenario) %>%
  summarise(costs = sum(costs, na.rm = T))%>%
  ungroup %>%
  rbind(ct_schoolloss) %>%
  group_by(FIPS, scenario) %>%
  summarise(morbidity_costs = sum(costs, na.rm = T))%>%
  ungroup 

ct_totals <- ct_mortality %>%
  full_join(ct_morbidity, by = c("FIPS", "scenario")) %>%
  mutate(costs = morbidity_costs + mortality_costs,
         costs_per_capita = costs/population) 

write.csv(ct_totals, file = "processed/ct_totals.csv", row.names = FALSE)
