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

# Simple DAG --------------------------------------------------------------
# Specify relationships for simplified dag
dag_simple <- dagify(
  Biodiversity ~ Local_LC + Ssp_matrix + Climate + Unknown, #Time_since_planting + Landscape_forest + Elevation 
  Ssp_matrix ~ Climate + Unknown,
  Unknown ~ Spatial_auto,
  exposure = c("Local_LC", "Ssp_matrix"), 
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
  Biodiversity ~ lc_50m + Ssp_matrix + Climate + Landscape_forest + Dist_forest + Time + Unknown, #Time_since_planting  + Elevation 
  lc_50m ~ Time + Sampling + Farmer_values, 
  Ssp_matrix ~ Climate + Farmer_values + Unknown,
  #Tree_characteristics ~ Time,
  Landscape_forest ~ Farmer_values + Sampling,
  Dist_forest ~ lc_50m,
  Unknown ~ Spatial_auto,
  exposure = c("lc_50m"), #, "Ssp_matrix"
  outcome = "Biodiversity" #, 
  #latent = c("Unknown", "Farmer_values") #Tree_characteristics
)
# Adjustment set - Remove latent variables to see
adjustmentSets(dag_complex, type = "minimal") # all

# Plot
tidy_dagitty(dag_complex, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")


# >No Ssp_matrix ----------------------------------------------------------
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

