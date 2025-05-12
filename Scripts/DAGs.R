## PhD birds in silvopastoral landscapes ##
# Directed acyclic graphs (DAG)
# This script specifies the causal assumptions, i.e. the mechanistic relationships between variables, for PhD dissertation. 

# Background --------------------------------------------------------------
# Goal: Estimate the direct effect of silvopasture on biodiversity

# When specifying the DAG, we focus on the mechanistic relationships between variables. Variables that are just statistically correlated (but no mechanistic link) do not end up in the DAG. Thus, causal modeling helps us determine the appropriate variables to include (or not include) to obtain unbiased estimates of the parameters of interest. But it does not tell us how predictor variables may be related (i.e., interactions), or the functional form that variables should take.  

# Simple DAG --------------------------------------------------------------
# Specify relationships for simplified dag
dag_simple <- dagify(
  Biodiversity ~ Local_LC + Landscape_forest + Ssp_matrix + Precipitation + Elevation + Time_since_planting,
  Ssp_matrix ~ Precipitation,
  exposure = c("Local_LC", "Ssp_matrix", "Time_since_planting"), 
  outcome = "Biodiversity"
)
# Adjustment set 
adjustmentSets(dag_simple, type = "minimal", effect = "direct") 

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


# Complex DAG -------------------------------------------------------------

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

#Using ggdag package
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

#Or for additional customizability 
set.seed(123)
ggplot(data = tidy_dag_phd, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_node(aes(color = Status)) +
  geom_dag_label_repel(aes(label = label), force = 10) + #Tried the padding arguments but didn't have luck 
  geom_dag_edges(curvature = 0) +
  theme_dag() +
  guides(fill = "none", color = "none")
