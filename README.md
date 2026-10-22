# Silvopasture and Bird Diversity in Colombia

This repository contains the code and data associated with the first chapter of my PhD dissertation, examining the impact of silvopasture on bird taxonomic and functional diversity in fragmented landscapes of Colombia. The project aims to explore how local landcover and land-use impact bird diversity, and how these relationships are moderated by landscape and environmental factors.

## Chapter Overview

As part of the Sustainable Cattle Ranching (SCR) project, this research uses data from 500+ unique point count sites collected over a decade to model multi-species occupancy and abundance. The study seeks to:

- Understand how silvopasture, a system combining productive trees and shrubs with pastoral land, influences bird taxonomic and functional(?) diversity through time.
- Generate recommendations regarding the optimal amount and configuration of on-farm tree cover to maintain key ecosystem services. 
- Explore key hypotheses such as the 'dominance of beta diversity' and 'intermediate landscape complexity' hypotheses in the context of fragmented landscapes in Colombia

The goal is to provide actionable insights for biodiversity conservation in the context of sustainable beef production and habitat restoration efforts in Colombia.

## Methodology

- **Data**: 500+ point count sites collected over a decade.
- **Analysis**: Bayesian multi-species occupancy / abundance models, alpha and beta diversity, community composition, functional diversity metrics

## Pipeline

Run from the repository root (open `Ch1-ssp-birds.Rproj`). Survey tables, point locations and site covariates are read from the sibling repository `../Ssp-bird-data-wrangling/`; scripts stop with a deliberate `stop()` before optional or export sections.

| Script | Purpose |
|----|----|
| `Scripts/00a_Extract_lcs.R` | Processes the digitized land-cover polygons for each period (middle, past, ubc) into snapped land-cover layers around each point count |
| `Scripts/00b_LSM.R` | Landscape metrics (e.g. `Dist_forest`, `ssp`, `te`) per period, then joins them to the survey table → `Derived/Excels/Event_covs_lsm.csv` |
| `Scripts/01_Wrangling_ms.R` | Builds the multi-species `spOccupancy` data list |
| `Scripts/Correlation_vis.R` | Screens candidate predictors for collinearity |
| `Scripts/spOccupancy_ms.R` | Fits the multi-species occupancy models and runs diagnostics |
| `Scripts/postHoc_ft.R`, `Scripts/NMDS.R`, `Scripts/Marginal_effects.R`, `Scripts/Plotting.R` | Trait analysis, community composition, marginal effects and the sampling map |
| `Scripts/spOcc_fns.R` | Shared functions, sourced by the scripts above |
| `Scripts/DAGs.R`, `Scripts/qmd/Math.qmd` | Causal diagrams and model notation |

`Scripts/_archive/` holds superseded scripts kept for reference. `Data/`, `Derived/`, `Figures/` and `Extra_scripts/` are local and not tracked.

## Usage

If you are interested in collaborating on this project please send me an email at skinnerayayron93 [at] gmail [dot] com. The data is not available publicly at present. For the time being the code can only be observed visually, although similar analyses could be performed using publicly available biodiversity data (e.g., from the vegan package). 

## License

Feel free to use, modify, and distribute this code however you wish. No restrictions apply.

## Acknowledgments

Thank you to the Sustainable Cattle Ranching project and to my advisors for providing data and helping me understand cattle ranching in Colombia. Thank you to the NGO SELVA for sponsoring a year long Fulbright grant in Colombia, where they provided logistical and conceptual guidance, as well as providing important training towards identification and manipulation of the Colombian avifauna. A special thanks to my advisors at UBC and TNC for providing funding, professional guidance, and both academic and personal support.

## Contact

For more information or questions, please contact skinnerayayron93 [at] gmail [dot] com
