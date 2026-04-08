## PhD birds in silvopastoral landscapes ##
# Post-hoc functional traits model 

## Objective: Understand if response to silvopasture is influenced by species' morphological and life-history functional traits 

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(posterior)
library(spOccupancy)
library(ggdist)

Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"
Ft <- read_csv(paste0(Wrangling_repo, Excels, "Traits/Functional_traits.csv"))

# Occupancy model output
#Occ_mod <- readRDS("Rdata/Model_output/Occ_mod_100spp_meta_restrict_phi2026-03-25.rds")

# Forest or non-forest? ---------------------------------------------------
Forest <- FALSE
if(Forest){
  variable_posthoc <- "Ano1"
} else{
  variable_posthoc <- "ssp"
  } 

# Format ------------------------------------------------------------------
# Rows correspond to the posterior MCMC samples and the columns correspond to species
y_draws <- as_draws_df(Occ_mod$beta.samples) %>% 
  select(starts_with(variable_posthoc)) %>% 
  rename_with(~str_remove(., paste0(variable_posthoc, "-")))

# >PCA --------------------------------------------------------------------
# PCA of bill morph
pca_bill <- prcomp(data = Ft, ~ Beak.Depth + Beak.Length_Nares + Beak.Width)
Ft$Beak_size <- pca_bill$x[,1]
Ft$Beak_shape <- pca_bill$x[,2]
biplot(pca_bill)

## If beak size axes are negative, make positive
if(all(pca_bill$rotation[,1] < 0)){
  Ft$Beak_size <- -Ft$Beak_size
} else if(all(pca_bill$rotation[,1] > 0)){
  Ft$Beak_size
} else{
  stop("Inconsistent axes directionality")
}

# PCA of body size
pca_size <- prcomp(data = Ft, ~ Mass + Tail.Length + Tarsus.Length)
Ft$Body_size <- pca_size$x[,1]
# Dominated by mass
biplot(pca_size)

## If body size axes are negative, make positive
if(all(pca_size$rotation[,1] < 0)){
  Ft$Body_size <- -Ft$Body_size
} else if(all(pca_size$rotation[,1] > 0)){
  Ft$Body_size
} else{
  stop("Inconsistent axes directionality")
}

# >General ----------------------------------------------------------------
## General format - Ensure species line up in Ft and y_draws
Spp <- str_replace(names(y_draws), "_", " ")
Ft_spp <- Ft %>%
  filter(Species_ayerbe %in% Spp) %>%
  rename(HWI = `Hand-Wing.Index`) %>%
  mutate(Species_ayerbe = factor(Species_ayerbe, levels = Spp),
         Trophic.Level = ifelse(
           Trophic.Level == "Scavenger", "Carnivore", Trophic.Level)
         ) %>%
  arrange(Species_ayerbe) %>%
  mutate(across(where(is.numeric), scale)) %>% 
  select(-Species_bl) %>% 
  distinct()

# >Morphology -------------------------------------------------------------
covariates_morph <- c("Body_size", "HWI", "Beak_size", "Beak_shape")
Ft_morph <- Ft_spp %>% 
  select(Species_ayerbe, all_of(covariates_morph)) %>%
  # postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(covariates_morph), ~ !is.na(.x)))

# Create list required by postHocLM
postHoc_l_morph <- list(y = y_draws, covs = Ft_morph)

# >Eye --------------------------------------------------------------------
Ft_eye <- Ft_spp %>% filter(!is.na(Eye_resid)) %>% 
  select(Species_ayerbe, Eye_resid)
Spp_eye <- Ft_eye %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_eye <- y_draws %>% select(all_of(Spp_eye)) 

# List required by postHocLM
postHoc_l_eye <- list(y = y_draws_eye, covs = Ft_eye)

# >Life-history -----------------------------------------------------------
covariates_lh <- c("Habitat.Density", "Migration", "Trophic.Level", "gen_length")
covariates_specialization <- c("Elev_range_final", "Range.Size", "Forest_bin", "Habitat_breadth", "Diet_breadth") # "Ecological_specialization"
Ft_lh <- Ft_spp %>% 
  select(Species_ayerbe, all_of(c(covariates_lh, covariates_specialization))) %>%
  # postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(c(covariates_lh, covariates_specialization)), ~ !is.na(.x)))
postHoc_l_lh <- list(y = y_draws, covs = Ft_lh)

# >Clutch -----------------------------------------------------------------
Ft_clutch <- Ft_spp %>% filter(!is.na(Clutch)) %>% 
  select(Species_ayerbe, Clutch)
Spp_clutch <- Ft_clutch %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_clutch <- y_draws %>% select(all_of(Spp_clutch))

# List required by postHocLM
postHoc_l_clutch <- list(y = y_draws_clutch, covs = Ft_clutch)

# >Nesting ----------------------------------------------------------------
# Nest location
covariates_nest_loc <- c("Nest_ground_bush", "N_nest_locs")
Ft_nesting <- Ft_spp %>% 
  select(Species_ayerbe, all_of(covariates_nest_loc)) %>%
  # postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(covariates_nest_loc), ~ !is.na(.x)))
Spp_nesting <- Ft_nesting %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_nesting <- y_draws %>% select(all_of(Spp_nesting))
postHoc_l_nesting <- list(y = y_draws_nesting, covs = Ft_nesting)

# Nest exposure
Ft_nest_exp <- Ft_spp %>% 
  filter(!is.na(Nest_exposure)) %>% 
  select(Species_ayerbe, Nest_exposure)
Spp_nest_exp <- Ft_nest_exp %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_nest_exp <- y_draws %>% select(all_of(Spp_nest_exp))
postHoc_l_nest_exp <- list(y = y_draws_nest_exp, covs = Ft_nest_exp)

# Examine correlations ----------------------------------------------------
# Morphology
Ft_spp %>% 
  select(all_of(covariates_morph), Mass, Tarsus.Length, Tail.Length) %>%
GGally::ggcorr(label = T, label_size = 2, label_round = 2, hjust = 0.75, size = 3, layout.exp = 1.01)

# Genearl life-history covariates (categorical)
if(FALSE){
  Ft_spp %>% 
    select(all_of(covariates_lh)) %>% 
    GGally::ggpairs()
}

# Specialization covariates
Ft_spp %>% 
  select(Species_ayerbe, all_of(covariates_specialization), Ecological_specialization) %>% 
  GGally::ggcorr(label = T, label_size = 2, label_round = 2, hjust = 0.75, size = 3, layout.exp = 1.01)

# Run model ---------------------------------------------------------------
## Morphology models
# Morph variables: Mass + HWI + Tarsus.Length + Tail.Length + Beak_size + Beak_shape
Morph_form <- as.formula(paste("~" , paste(covariates_morph, collapse = " + ")))
postHoc_morph <- postHocLM(formula = Morph_form , data = postHoc_l_morph)
# Residual eye size
postHoc_eye <- postHocLM(formula = ~Eye_resid, data = postHoc_l_eye)

## Life-history models
# Lh vars: Range.Size + Elev_range_final + Habitat.Density + Migration + Trophic.Level
lh_form <- as.formula(paste("~" , paste(covariates_lh, collapse = " + ")))
postHoc_lh <- postHocLM(formula = lh_form, data = postHoc_l_lh)

# Specialization
specialization_form <- as.formula(paste("~" , paste(covariates_specialization, collapse = " + ")))
postHoc_specialization <- postHocLM(formula = specialization_form, data = postHoc_l_lh)

# Clutch
postHoc_clutch <- postHocLM(formula = ~Clutch, data = postHoc_l_clutch)

## Nesting 
Nesting_form <- as.formula(paste("~" , paste(covariates_nest_loc, collapse = " + ")))
postHoc_nesting <- postHocLM(formula = Nesting_form , data = postHoc_l_nesting)
# Nest exposure
postHoc_nest_exp <- postHocLM(formula = ~Nest_exposure, data = postHoc_l_nest_exp)

# Plot --------------------------------------------------------------------
format.for.plotting <- function(samples){
  samples %>% as_draws_df() %>%
    #rename_with(~str_remove(., "Habitat")) %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    mutate(param = as.factor(param))
}

plot_half_eye <- function(df, y_var){
  df %>% 
    filter(!str_detect({{ y_var }}, "^\\.")) %>% # rm .iteration, .draw, and .chain
    ggplot(aes(x = value, y = {{ y_var }}, fill = {{ y_var }})) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(alpha = .8) +
    labs(x = NULL, y = NULL) +
    guides(fill = "none")
}

## Morphology
format.for.plotting(postHoc_morph$beta.samples) %>% plot_half_eye(y_var = param)
# Residual eye 
format.for.plotting(postHoc_eye$beta.samples) %>% plot_half_eye(y_var = param)

## Life-history
format.for.plotting(postHoc_lh$beta.samples) %>% plot_half_eye(y_var = param)
# Specialization
format.for.plotting(postHoc_specialization$beta.samples) %>% plot_half_eye(y_var = param)
# Nesting
format.for.plotting(postHoc_nesting$beta.samples) %>% plot_half_eye(y_var = param)
format.for.plotting(postHoc_nest_exp$beta.samples) %>% plot_half_eye(y_var = param)
# Clutch size
format.for.plotting(postHoc_clutch$beta.samples) %>% plot_half_eye(y_var = param)

# Inspect Threatened species ------------------------------------------------------
tibble(Species_ayerbe = str_replace(names(y_draws), "_", " "), 
       Response = colMeans(y_draws)) %>% 
  left_join(Ft[, c("Species_ayerbe", "gen_length")]) %>% 
  ggplot(aes(x = gen_length, y = Response)) + 
  geom_point() + 
  geom_smooth()

# Select T&E species
Spp_te <- Ft %>% 
  mutate(Species_ayerbe_ = str_replace(Species_ayerbe, " ", "_")) %>%
  filter(iucn_red_list != "LC") %>% 
  filter(Species_ayerbe_ %in% Occ_mod$sp.names) %>% 
  select(Species_ayerbe_, iucn_red_list)

# Plot
as_draws_df(Occ_mod$beta.samples) %>% 
  select(starts_with(variable_posthoc)) %>% 
  rename_with(~str_remove(., paste0(variable_posthoc, "-"))) %>%
  select(Spp_te$Species_ayerbe_) %>%
  pivot_longer(cols = everything(), 
               names_to = "Species_ayerbe", 
               values_to = "value") %>%
  plot_half_eye(y_var = Species_ayerbe)

# Examine sample sizes of T&E species
sort(obs_occ[Spp_te$Species_ayerbe_])
