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

Run from the repository root (open `Ch1-ssp-birds.Rproj`). Survey tables, point locations and site covariates are read from the sibling repository `../Ssp-bird-data-wrangling/`, and point-level landscape metrics (`Event_covs_lsm.csv`) from `../Farm_diversity/`, where the land-cover extraction and landscape-metric scripts now live; scripts stop with a deliberate `stop()` before optional or export sections.

| Script | Purpose |
|----|----|
| `Scripts/01_Wrangling_ms.R` | Builds the multi-species `spOccupancy` data list |
| `Scripts/Correlation_vis.R` | Screens candidate predictors for collinearity |
| `Scripts/spOccupancy_ms.R` | Fits the stacked multi-species occupancy models (one row per location × survey period) and runs diagnostics |
| `Scripts/spOccupancy_temporal.R` | Reshapes the stacked data list into the multi-season format (one site per location, year-group × season periods) and fits the multi-season spatial occupancy model (`svcTMsPGOcc` with a spatially varying intercept only, equivalent to `stMsPGOcc` but with range clipping) and its non-spatial counterpart (`tMsPGOcc`, location random effect); holds the primary model |
| `Scripts/postHoc_ft.R`, `Scripts/NMDS.R`, `Scripts/Marginal_effects.R`, `Scripts/Plotting.R` | Trait analysis (summary saved to `Derived/Excels/postHoc_<response>_<fit>.csv`), community composition (figures in `Figures/NMDS_meta/`; species counted inside each habitat ellipse saved to `Derived/Excels/nmds_ellipse_counts_<model_tag>.csv`), marginal effects and the sampling map |
| `Scripts/spOcc_fns.R` | Shared functions, sourced by the scripts above, including `fit_temporal()`, which fits and saves one multi-season model |
| `Scripts/DAGs.R`, `Scripts/qmd/Math.qmd` | Causal diagrams and model notation |

Data lists and fitted models are saved under `Derived/Rdata/`, split by data structure: `Model_input/Stacked/` and `Model_input/Temporal/` (data lists), `Model_output/Stacked/` and `Model_output/Temporal/` (fits).

`Scripts/_archive/` holds superseded scripts kept for reference. `Data/`, `Derived/`, `Figures/`, `Extra_scripts/` and `Docs/` (correspondence drafts) are local and not tracked.

## Usage

If you are interested in collaborating on this project please send me an email at skinnerayayron93 [at] gmail [dot] com. The data is not available publicly at present. For the time being the code can only be observed visually, although similar analyses could be performed using publicly available biodiversity data (e.g., from the vegan package). 

## License

Feel free to use, modify, and distribute this code however you wish. No restrictions apply.

## Acknowledgments

Thank you to the Sustainable Cattle Ranching project and to my advisors for providing data and helping me understand cattle ranching in Colombia. Thank you to the NGO SELVA for sponsoring a year long Fulbright grant in Colombia, where they provided logistical and conceptual guidance, as well as providing important training towards identification and manipulation of the Colombian avifauna. A special thanks to my advisors at UBC and TNC for providing funding, professional guidance, and both academic and personal support.

## Contact

For more information or questions, please contact skinnerayayron93 [at] gmail [dot] com
