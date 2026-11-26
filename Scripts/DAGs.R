## PhD birds in silvopastoral landscapes ##
# Directed acyclic graphs (DAG)
# This script specifies the causal assumptions, i.e. the mechanistic relationships between variables, for PhD dissertation. 

# Background --------------------------------------------------------------
# Goal: Estimate the direct effect of silvopasture on biodiversity

# When specifying the DAG, we focus on the mechanistic relationships between variables. Variables that are just statistically correlated (but no mechanistic link) do not end up in the DAG. Thus, causal modeling helps us determine the appropriate variables to include (or not include) to obtain unbiased estimates of the parameters of interest. But they are non-parametric: they do not make assumptions regarding the functional form or the distribution of the response variable

# Rationale ---------------------------------------------------------------

# Components of the detectability model that also influence the ecological process of interest can be included 
# Canopy height more relevant in the forest-only dataset. In the land-use gradient dataset / models, I already have binary forest (represented by lc_50m in DAG) in the detectability model 
# Time influences local landcover (quality) directly through tree characteristics, and indirectly through sampling
# Time influences biodiversity through forest growth extinction debts

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(dagitty)
library(ggdag)
# ggdag masks dplyr's filter(); the local tests filter ordinary data frames
filter <- dplyr::filter
# The local tests' Pillai-trace test for categorical nodes needs the CCP package installed (not attached)

# Simple DAG --------------------------------------------------------------
# Specify relationships for simplified dag
dag_simple <- dagify(
  Biodiversity ~ Local_LC + ssp_300m + Climate + Unknown, #Time_since_planting + Landscape_forest + Elevation 
  ssp_300m ~ Climate + Unknown,
  Unknown ~ Spatial_auto,
  exposure = c("Local_LC", "ssp_300m"), 
  outcome = "Biodiversity",
  latent = c("Unknown") 
)
# Adjustment set - Remove latent variables to see
adjustmentSets(dag_simple, type = "minimal") 

# Plot
tidy_dagitty(dag_simple, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")

# More complex DAG --------------------------------------------------------
## Specify relationships
dag_complex <- dagify(
  Biodiversity_obs ~ Biodiversity + Obs_skill + lc_50m, #+ Canopy_height
  Obs_skill ~ Sampling,
  Sampling ~ Time, 
  Biodiversity ~ lc_50m + Habitat_structure + ssp_300m + forest_300m + Climate + Landscape_forest + Dist_forest + Time + Unknown, #Time_since_planting  + Elevation
  lc_50m ~ Time + Sampling + Farmer_values,
  # Sampling protocols chose where points were surveyed, and a point's location fixes its climate, so the habitats sit at different places along the rainfall gradient (2026-11-25: pasture points wetter than live-fence points): a backdoor path from local land cover through sampling and climate, closed by conditioning on climate
  Climate ~ Sampling,
  # Local land cover is the habitat type; its structure (canopy density, tree height) grows with rainfall and time since planting, so climate shapes the habitat within a type rather than which type a point is; structure is a mediator, not adjusted for, and the type-specific climate response is the habitat-by-rainfall interaction
  Habitat_structure ~ lc_50m + Climate + Time,
  # The 50 m point-count disc lies inside the 300 m buffer, so the point's own land cover enters the buffer's shares (median forest within 300 m: 38% at forest points vs 11-16% at open points; silvopasture within 300 m: 22% at silvopasture points, 19% at pasture, 7% at forest); added 2026-11-25 after the local tests found both shares strongly dependent on sampling given climate, which these arrows explain (forest-heavy programs sampled forest points)
  # With these arrows the 300 m shares are partly mediators of the point's habitat, so conditioning on them targets the controlled direct effect of habitat (same landscape, different local habitat), identified given climate, landscape forest, time and the spatial unknowns (Unknown, the spatial random effect's stand-in), which also confound the shares and biodiversity
  ssp_300m ~ lc_50m + Climate + Farmer_values + Buffer_300m + Unknown,
  forest_300m ~ lc_50m + Farmer_values + Climate + Buffer_300m,
  Landscape_forest ~ Farmer_values + Sampling + Climate,
  # Distance to the forest edge depends on the habitat (a forest point is inside forest) and on where surveyors placed points relative to edges (silvopasture points sit closer to forest than pasture points): a mediator for the forest contrast, a confounder for silvopasture vs pasture, so adjusting for it is needed and habitat effects are read at a fixed distance
  # No pasture point lies closer than ~33 m to forest (a point at the edge is labelled forest), so forest vs pasture at the same distance is never observed; compare forest at its typical depth with the open habitats at a shared distance instead
  # More forest near a point puts the forest edge closer (local tests: partial correlation -0.58 between distance and forest within 300 m given sampling and habitat before this arrow was added, 2026-11-25)
  Dist_forest ~ lc_50m + Sampling + forest_300m,
  Unknown ~ Spatial_auto,
  exposure = c("lc_50m"), #, "ssp_300m"
  outcome = "Biodiversity", 
  latent = c("Buffer_300m", "Habitat_structure") #"Unknown", "Farmer_values"
)
# Adjustment set - Remove latent variables to see
adjustmentSets(dag_complex, type = "minimal") # all
# Silvopasture and forest within 300 m share one fixed buffer area, so one share cannot rise without another falling (digitization choices, e.g. drawing a stand as forest or as dispersed trees, add to this): a backdoor path through forest within 300 m, closed by conditioning on it; the buffer is marked latent because it cannot be adjusted for separately from its shares
adjustmentSets(dag_complex, exposure = "ssp_300m", type = "minimal")

# Plot
set.seed(123)
dag_complex_plot <- tidy_dagitty(dag_complex, layout = "fr") %>%
  ggdag_status(text_col = "black",
               text = TRUE,
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3,
               stylized = TRUE) +
  theme_dag() +
  guides(fill = "none", color = "none")
dir.create("Figures/DAGs", showWarnings = FALSE)
ggsave("Figures/DAGs/dag_complex.png", dag_complex_plot, width = 8, height = 6, bg = "white")
print(dag_complex_plot)

# >Test implied conditional independencies --------------------------------
### The DAG claims that some pairs of measured variables are independent given others; local tests check those claims against the surveyed points (Ankan, Wessels & Textor 2021, Current Protocols), and a clear dependence means a missing arrow or common cause
# Judged by effect size as well as the p-value: rows repeat locations across year-groups and points cluster on farms, so p-values are too small and confidence intervals too narrow
# Unit: one row per location and year-group (the land-cover shares and distances are updated per year-group; seasons within a year-group share them, except landscape canopy cover, averaged over seasons)
stacked_list_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm.rds"
lsm_path <- "../Farm_diversity/Derived/Excels/Event_covs_lsm.csv"
ring_covs_path <- "Derived/Excels/Scale_of_effect_fit_covs.csv"
landscape_radius_m <- 1500
# Partial correlation (or its categorical analogue) marking a negligible and a substantial dependence
ci_effect_small <- 0.1
ci_effect_large <- 0.3

## Measured nodes
# Sampling is coded as program-by-year-group, so it also fixes Time (the year-group), which is merged into it for testing; climate is annual rainfall, forest_300m and ssp_300m the forest and silvopasture shares within 300 m, Landscape_forest canopy cover within 1.5 km, Dist_forest the signed distance to the forest edge
occ_covs <- as_tibble(readRDS(stacked_list_path)$occ.covs) %>% mutate(across(c(Id_survey, Year_grp, Season), as.character))
forest_300 <- read_csv(lsm_path, show_col_types = FALSE) %>% distinct(Id_survey = as.character(Id_survey), Year_grp = as.character(Year_grp), forest)
landscape_cover <- read_csv(ring_covs_path, col_types = cols(Id_survey = "c", Year_grp = "c", Season = "c")) %>% filter(Radius_m == landscape_radius_m) %>% select(Id_survey, Year_grp, Season, Cover_pct)
dag_data <- occ_covs %>%
  left_join(forest_300, by = c("Id_survey", "Year_grp"), relationship = "many-to-one") %>%
  left_join(landscape_cover, by = c("Id_survey", "Year_grp", "Season"), relationship = "one-to-one") %>%
  group_by(Id_survey_no_dc, Year_grp) %>%
  summarise(Sampling = factor(paste(first(Uniq_db), first(Year_grp))), lc_50m = first(Habitat), Climate = first(Tot_prec), Dist_forest = first(Dist_forest), ssp_300m = first(ssp), forest_300m = first(forest), Landscape_forest = mean(Cover_pct), .groups = "drop") %>%
  select(-Id_survey_no_dc, -Year_grp) %>%
  as.data.frame()
stopifnot(!anyNA(dag_data))

## Testing DAG: Time merged into Sampling, every unmeasured node latent, so only claims among measured variables are listed
test_edges <- dagitty::edges(dag_complex) %>%
  mutate(across(c(v, w), \(x) if_else(x == "Time", "Sampling", x))) %>%
  filter(v != w) %>%
  distinct(v, w)
dag_test <- dagitty(paste0("dag{", paste(test_edges$v, "->", test_edges$w, collapse = "; "), "}"))
latents(dag_test) <- setdiff(names(dag_test), names(dag_data))
# Pillai-trace tests handle the categorical habitat and sampling nodes: for two continuous variables the estimate is their partial correlation with its 95% CI, and for a categorical node the root-mean-square canonical correlation (same 0-1 scale as an absolute correlation), for which dagitty gives no CI (it returns the estimate in both bounds, set to NA here)
ci_tests <- localTests(dag_test, dag_data, type = "cis.pillai", abbreviate.names = FALSE) %>%
  rownames_to_column("Claim") %>%
  rename(Lower = `2.5%`, Upper = `97.5%`, P_value = p.value, Estimate = estimate) %>%
  mutate(across(c(Lower, Upper), \(x) if_else(Lower == Upper, NA, x)),
         Dependence = case_when(abs(Estimate) >= ci_effect_large ~ "Substantial", abs(Estimate) >= ci_effect_small ~ "Weak", .default = "Negligible"),
         across(c(Estimate, Lower, Upper), \(x) round(x, 3)), P_value = signif(P_value, 3)) %>%
  arrange(desc(abs(Estimate)))
write_csv(ci_tests, "Derived/Excels/DAG_local_tests.csv")
print(ci_tests)

## Plot: one row per claim, the dependence with its 95% CI against the negligible and substantial bands
ci_tests_plot <- ggplot(ci_tests, aes(x = Estimate, y = fct_reorder(Claim, abs(Estimate)), colour = Dependence)) +
  annotate("rect", xmin = -ci_effect_small, xmax = ci_effect_small, ymin = -Inf, ymax = Inf, alpha = 0.15) +
  geom_vline(xintercept = c(-ci_effect_large, ci_effect_large), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper)) +
  scale_colour_manual(values = c(Negligible = "grey40", Weak = "#E69F00", Substantial = "#D55E00")) +
  labs(x = "Conditional dependence (partial correlation or Pillai analogue)", y = NULL)
ggsave("Figures/DAGs/dag_complex_local_tests.png", ci_tests_plot, width = 9, height = 0.3 * nrow(ci_tests) + 1.5, bg = "white")
print(ci_tests_plot)


# >No ssp_300m ----------------------------------------------------------
# Unclear whether we really have the data to include silvopasture in the matrix, so this DAG just focuses on the local landcover within 50m  

# Specify relationships
dag_LC <- dagify(
  Biodiversity ~ lc_50m + Climate + Landscape_forest + Tree_characteristics + Unknown, #Time_since_planting  + Elevation 
  lc_50m ~ Time + Climate + Farmer_values + Unknown, 
  Tree_characteristics ~ Time + Climate,
  Climate ~ Time,
  Landscape_forest ~ Farmer_values,
  Unknown ~ Spatial_auto,
  exposure = c("lc_50m"), 
  outcome = "Biodiversity" #,
  #latent = c("Unknown", "Farmer_values", "Tree_characteristics") 
)
# Adjustment set - Remove latent variables to see
adjustmentSets(dag_LC, type = "minimal") 

# Plot
tidy_dagitty(dag_LC, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")


# Interaction DAG ---------------------------------------------------------
# >Attia et al (2022) -----------------------------------------------------
# Specify relationships using Attia et al (2022) notation
idag <- dagitty('
dag {

  bb = "-.5,-.5,.5,.5"

  Time 
  lc_50m
  timeXlc_50m [exposure]
  Biodiversity      [outcome, pos="0,-0.45"]

  Time -> lc_50m
  Time -> timeXlc_50m
  lc_50m -> timeXlc_50m
  Time -> Biodiversity
  lc_50m -> Biodiversity
  timeXlc_50m -> Biodiversity
}
')
plot(idag, node.names = c("timeXlc_50m" = "Time x Local LC"), cex = 1.6)
adjustmentSets(idag)

# Specify relationships - Short term rainfall could be an important variable
# Sampling acts as a confounder on the interaction of timeXlc_50m, because it influences the chance of seeing pasture and SSP in later years (after 2019), and also indirectly influences biodiversity through the observation process
idag_complex <- dagify(
  Biodiversity_obs ~ Biodiversity + Obs_skill,
  Obs_skill ~ Sampling,
  Biodiversity ~ lc_50m + ssp_300m + Climate + Landscape_forest + Unknown + timeXlc_50m + Landscape_forestXlc_50m + climateXlc_50m + Time, #Time_since_planting  + Elevation
  timeXlc_50m ~ Time + lc_50m + Sampling,
  Landscape_forestXlc_50m ~ Landscape_forest + lc_50m, 
  climateXlc_50m ~ Climate + lc_50m,
  lc_50m ~ Time + Farmer_values + Sampling, 
  ssp_300m ~ Climate + Farmer_values + Unknown,
  Landscape_forest ~ Farmer_values + Time + Climate,
  Sampling ~ Time,
  Unknown ~ Spatial_auto,
  exposure = c("timeXlc_50m"), # "Landscape_forestXlc_50m"
  outcome = "Biodiversity" #,
  #latent = Sampling
)

# Adjustment set - Remove latent variables to see
adjustmentSets(idag_complex, type = "minimal") 

# Plot
tidy_dagitty(idag_complex, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")

# >Nilsson et al (2021) ---------------------------------------------------
# Specify relationships using Nilsson et al (2021) notation (delta Y)
idag <- dagify(
  Change_Biodiv_LC ~ Climate + Landscape_forest + Time_since_planting,  
  Landscape_forest ~ Farmer_values,
  outcome = "Change_Biodiv_LC",
  latent = c("Farmer_values") 
)
# Adjustment set - Remove latent variables to see
adjustmentSets(idag, type = "minimal") 

# Plot
tidy_dagitty(idag, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")

# All variables ------------------------------------------------------------
# Mathilde biases also include my biases, like that I know Meta better then other regions
# Climatic_extremes would be like strong El niño years, drought, etc. whereas climate refers more generally to fact that warmer & wetter tends to have higher species richness
# TIME - refer back to models in proposal
dag_phd <- dagify(
  Biodiv_Traits ~ Local_LC + Landscape_LC + Climate + Climatic_extremes + Elevation + Species_Biases + Planted_tree_species + Time_since_planting + Habitat_structure,
  Species_Biases ~ Detection + species_misID,
  Detection ~ Observer + Habitat + Weather,
  Local_LC ~~ GE_image, 
  Local_LC ~ Mathilde_biases, #Add + Planted_tree_species to illustrate a fork
  Landscape_LC ~ Brandt_raster, 
  GE_image ~ Clouds, 
  GE_image ~ Image_dates,
  labels = c("Biodiv_Traits" = "Biodiversity or \nTraits",
             "GE_image" = "Google \n Earth Image", 
             "species_misID" = "species \nmisidentification", 
             "Local_LC" = "Local \n landcover",
             "Landscape_LC" = "Landscape \n landcover"),
  exposure = "Local_LC", #
  outcome = "Biodiv_Traits",
  latent = c("species_misID", "Mathilde_biases", "Habitat_structure")
)
#impliedConditionalIndependencies(dag_phd)
adjustmentSets(dag_phd, type = "minimal", effect = "direct") #This means I don't need to include any of these in the model? 

tidy_dag_phd <- tidy_dagitty(dag_phd, layout = "fr") %>% 
  mutate(label = ifelse(is.na(label), str_replace_all(name, "_", "\n"), label), 
         #Manually create "status" column
         Status = case_when(
           name == "Local_LC" ~ "Exposure", #
           name == "Biodiv_Traits" ~ "Outcome",
           name %in% c("species_misID", "Mathilde_biases", "Habitat_structure") ~ "Latent"
         ))

# Using ggdag package
set.seed(123)
ggdag_status(tidy_dag_phd, text_col = "black",
             use_labels = "label",
             text = FALSE, 
             edge_type = "link_arc",
             node_size = 20,
             text_size = 3, 
             stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")

# Or for additional customizability 
set.seed(123)
ggplot(data = tidy_dag_phd, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_node(aes(color = Status)) +
  geom_dag_label_repel(aes(label = label), force = 10) + # Tried the padding arguments but didn't have luck 
  geom_dag_edges(curvature = 0) +
  theme_dag() +
  guides(fill = "none", color = "none")

