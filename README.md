# Health damages and unequal economic burdens of dust events from the shrinking Great Salt Lake
This repository is for users to replicate the results of Garcia et al. 2027 (In XX): "Health damages and unequal economic burdens of dust events from the shrinking Great Salt Lake."

> **Authors**: Alberto Garcia (albert.garcia@utah.edu), Sara Grineski, Derek Mallia

> **Abstract:** Drying saline lakes are a major source of dust pollution around the world, yet the magnitude and distribution of the resulting health damages remain largely unquantified. Using an integrated damage-function framework, we combine atmospheric transport modeling with census-tract population data and epidemiological concentration–response functions to estimate the public health burden associated with dust emissions from the Great Salt Lake under several policy-relevant water-level scenarios. We then monetize these impacts using economic valuation methods used in regulatory benefit–cost analysis. At current lake levels, lakebed dust is associated with $87.64 million in annual health costs. Continued drying would exacerbate these costs, while attaining a minimum healthy lake level would avoid net present value damages of $1.45 billion through 2060 relative to current conditions. We further document pronounced distributional heterogeneity, as per-capita damages are highest among Pacific Islander communities and increase disproportionately as lake levels fall. Our results highlight the importance of accounting for air-quality externalities and recognizing that the distribution of pollution exposure may differ from that of the resulting economic damages when evaluating water management and ecosystem restoration policies.

## Repository structure  

In order to reproduce analyses, you must download the `data` folder at XX and add it to the repository.

```
GSL-healthcosts-replication
  |__ data (need to download)
  |__ scripts_R1
  |__ code
  |__ figs_pdf
  |__ results
```

## Software

This analysis was performed in R. The scripts for reproducing the paper analysis can be found in `scripts_R1`.

## License

The software code, data, and figures contained within this repository are made available under the [Creative Commons CC BY-NC-ND 2.0](https://creativecommons.org/licenses/by-nc-nd/2.0/) license.

## Code scripts
This section describes the code scripts developed to conduct this study. All code related to this study is available in the GitHub repository. The headings below represent folder names and paths.

### scripts_R1/
`02-dataprep_mortality.R`
prepares tract-by-age incidence data for mortality in the study region

`03-dataprep_mortality_race.R`
prepares tract-by-age-by-race incidence data for mortality in the study region
 
`04-dataprep_mortality_projections.R`
prepares mortality incidence data projected through 2060 in the study region

`05-dataprep_morbidity.R`
prepares tract-by-age-by-endpoint incidence data for morbidity endpoints excluding school absences in the study region

`06-dataprep_morbidity_race.R`
prepares tract-by-age-by-race-by-endpoint incidence data for morbidity endpoints excluding school absences in the study region
 
`07-dataprep_morbidity_projections.R`
prepares morbidity incidence data excluding school absences through 2060 in the study region

`08-dataprep_education.R`
prepares tract-by-age-by-race and tract-by-age incidence data for school absences in the study region

`09-dataprep_education_projections.R`
prepares school absence incidence data through 2060 in the study region

`11-VSL.R`
converts EPA VSL to 2024 dollars

`12-age_VSL.R`
converts Aldi and Viscusi cohort-adjusted VSL to 2024 dollars

`13-cost_morbidity.R`
converts mobidity cost of illness valuations to 2024 dollars

`14-costs_health_projected.R`
projects VSL and morbidity COIs through 2060

`21-health_mortality.R`
calculates annual mortality associated with dust events and the associated costs

`22-health_mortality_projections.R`
projects mortality associated with dust events and the associated costs through 2060

`23-health_morbidity.R`
calculates annual morbidity, excluding school absences, associated with dust events and the associated costs

`24-health_morbidity_projections.R`
projects morbidity, excluding school absences, associated with dust events and the associated costs through 2060

`25-health_schoolloss.R`
calculates annual school absence associated with dust events and the associated costs

`26-health_schoolloss_projections.R`
projects school absences associated with dust events and the associated costs through 2060

`97-total_cost_projections.R`
aggregates mortality, morbidity, and school absence costs, both annually and projected through 2060, for damage comparison. Most figures are generated in this script.

`98-map_figures.R`
creates map figures in the manuscript
