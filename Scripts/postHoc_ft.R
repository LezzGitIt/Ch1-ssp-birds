## PhD birds in silvopastoral landscapes ##
# Post-hoc functional traits model 

## Objective: Understand if response to silvopasture is influenced by species' morphological and life-history functional traits 

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(posterior)
library(spOccupancy)
library(ggdist)
ggplot2::theme_set(cowplot::theme_cowplot())
# postHocLM samples its posterior, so a fixed seed keeps the saved summary reproducible
set.seed(20261105)

## Ch1 trait table built in the sibling Fd-pd repo (01_Ch1_traits.R): the deposit traits plus eye size, BirdLife clutch, nest and incubation traits, and Forest_specialization (1/rank-weighted wooded share of BIRDBASE habitat ranks, which replaced the AVONET-based Forest_bin)
Ft <- read_csv("../Fd-pd/Derived/Functional_traits_ch1.csv")

# Occupancy model output
Occ_mod_path <- "Derived/Rdata/Model_output/Temporal/Occ_mod_98spp_prec_hinge_f2_inits_2026-11-08.rds"
Occ_mod <- readRDS(Occ_mod_path)
# Stacked list the model's data were built from, for the distance-to-forest scaling and forest's reference depth
msOcc_l_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_98_hab_fix_jit1cm.rds"
source("Scripts/spOcc_fns.R")

# Contrasts ---------------------------------------------------------------
# Species' habitat effects regressed on traits: silvopasture vs pasture and silvopasture vs forest, each on the logit scale within each draw
contrasts <- c("Ssp_vs_pasture", "Ssp_vs_forest")
### Habitats compared at a distance to the forest edge each one's data support, as in the marginal-effects figures and the report: pasture and silvopasture at a shared distance outside forest, forest at its median depth inside forest
# Silvopasture vs pasture is the HabitatSsp coefficient, as both open habitats share the reference distance and distance enters additively; silvopasture vs forest is (HabitatSsp + Dist_open slope x reference distance) minus (HabitatBosque + Dist_inside slope x reference depth), since the raw HabitatSsp minus HabitatBosque compares them at the forest edge, where no pasture point lies
Ref_open_m <- 100

# Format ------------------------------------------------------------------
# Reference distances on the hinge terms' scale (metres divided by the distance-to-forest SD, without centring)
occ_covs <- readRDS(msOcc_l_path)$occ.covs
Dist_scaling <- read_csv(str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv"), show_col_types = FALSE) %>% filter(Process == "occupancy", Covariate == "Dist_forest")
Ref_depth_m <- forest_reference_depth(occ_covs$Dist_forest * Dist_scaling$SD + Dist_scaling$Mean, occ_covs$Habitat == "Bosque", occ_covs$Id_survey, occ_covs$Year_grp)
d_open <- Ref_open_m / Dist_scaling$SD
d_in <- Ref_depth_m / Dist_scaling$SD
# Rows correspond to the posterior MCMC samples and the columns correspond to species
beta_draws <- as.data.frame(as.matrix(Occ_mod$beta.samples))
species_draws <- function(term) beta_draws %>% select(starts_with(paste0(term, "-"))) %>% rename_with(~str_remove(., paste0(term, "-")))
y_draws_l <- list(
  Ssp_vs_pasture = species_draws("HabitatSsp"),
  Ssp_vs_forest = (species_draws("HabitatSsp") + species_draws("Dist_open") * d_open) - (species_draws("HabitatBosque") + species_draws("Dist_inside") * d_in)
)[contrasts]
# Species order (and the threatened-species section) from the silvopasture vs pasture draws
y_draws <- y_draws_l$Ssp_vs_pasture

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
covariates_lh <- c("Habitat.Density", "Migration_avo", "Trophic.Level", "Gen_length")
covariates_specialization <- c("Elev_range", "Range.Size", "Forest_specialization", "Habitat_breadth", "Diet_breadth") # "Ecological_specialization"
Ft_lh <- Ft_spp %>% 
  select(Species_ayerbe, all_of(c(covariates_lh, covariates_specialization))) %>%
  # Residents are the reference, so the migration coefficients compare long-distance and partial migrants with resident species
  mutate(Migration_avo = fct_relevel(Migration_avo, "Sedentary")) %>%
  # postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(c(covariates_lh, covariates_specialization)), ~ !is.na(.x)))
postHoc_l_lh <- list(y = y_draws, covs = Ft_lh)

# >Clutch -----------------------------------------------------------------
# Clutch_mean (95 of 99 species) rather than BirdLife's Clutch_bl (68 of 99; r = 0.85)
Ft_clutch <- Ft_spp %>% filter(!is.na(Clutch_mean)) %>%
  select(Species_ayerbe, Clutch_mean)
Spp_clutch <- Ft_clutch %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_clutch <- y_draws %>% select(all_of(Spp_clutch))

# List required by postHocLM
postHoc_l_clutch <- list(y = y_draws_clutch, covs = Ft_clutch)

# >Nesting ----------------------------------------------------------------
# Nest location
# N_nest_locs_sheard (Sheard et al. 2024 first; every location across three handbooks) or N_nest_locs_bb (BIRDBASE first; typical substrates); 89 of 99 species each, r = 0.78
covariates_nest_loc <- c("Nest_ground_bush", "N_nest_locs_sheard")
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
## One trait model per group of traits; each data list pairs a trait table with the occupancy-draw columns of its species
# Morph variables: Mass + HWI + Tarsus.Length + Tail.Length + Beak_size + Beak_shape
Morph_form <- as.formula(paste("~" , paste(covariates_morph, collapse = " + ")))
# Lh vars: Range.Size + Elev_range + Habitat.Density + Migration + Trophic.Level
lh_form <- as.formula(paste("~" , paste(covariates_lh, collapse = " + ")))
specialization_form <- as.formula(paste("~" , paste(covariates_specialization, collapse = " + ")))
Nesting_form <- as.formula(paste("~" , paste(covariates_nest_loc, collapse = " + ")))
postHoc_forms <- list(Morphology = Morph_form, Eye_size = ~Eye_resid, Life_history = lh_form, Specialization = specialization_form, Clutch_size = ~Clutch_mean, Nest_location = Nesting_form, Nest_exposure = ~Nest_exposure)
postHoc_data <- list(Morphology = postHoc_l_morph, Eye_size = postHoc_l_eye, Life_history = postHoc_l_lh, Specialization = postHoc_l_lh, Clutch_size = postHoc_l_clutch, Nest_location = postHoc_l_nesting, Nest_exposure = postHoc_l_nest_exp)
# postHocLM pairs trait rows with occupancy-draw columns by position, so the species must match in the same order
iwalk(postHoc_data, \(d, nm) if (!identical(str_replace(as.character(d$covs$Species_ayerbe), " ", "_"), colnames(d$y))) stop(nm, ": trait-table species do not match the occupancy-draw columns"))

## Every model for each contrast, swapping in that contrast's draws for the same species
postHoc_fits <- map(y_draws_l, \(y_contrast) imap(postHoc_data, \(d, nm) postHocLM(formula = postHoc_forms[[nm]], data = list(y = y_contrast[, colnames(d$y)], covs = d$covs))))

# Plot --------------------------------------------------------------------
# Formatting function
format.for.plotting <- function(samples){
  samples %>% as_draws_df() %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    mutate(param = as.factor(param))
}

# Plot posteriors function
plot_half_eye <- function(df, y_var, N_species = length(Occ_mod$sp.names), Intercept = FALSE){
  if(isFALSE(Intercept)) {df <-  df %>% filter(param != "(Intercept)") }
  df %>% filter(!str_detect({{ y_var }}, "^\\.")) %>% # rm .iteration, .draw, and .chain
    ggplot(aes(x = value, y = {{ y_var }}, fill = {{ y_var }})) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    stat_halfeye(alpha = .8) +
    labs(x = NULL, y = NULL, title = paste({{ N_species }}, "species")) +
    guides(fill = "none")
}

## Posterior of every coefficient, one plot per contrast x model
iwalk(postHoc_fits, \(fits, contrast) iwalk(fits, \(fit, nm) {
  print(format.for.plotting(fit$beta.samples) %>% plot_half_eye(y_var = param, N_species = nrow(postHoc_data[[nm]]$covs)) + labs(subtitle = paste(contrast, nm, sep = ": ")))
}))

# Save summary ------------------------------------------------------------
# One row per contrast x trait-model coefficient (posterior mean, 95% CrI, P(>0)) with the number of species and the median Bayesian R2, read by the progress report; the file name carries the occupancy fit's name so a summary is never attributed to another fit
postHoc_summary <- imap(postHoc_fits, \(fits, contrast) imap(fits, \(fit, nm) {
  b <- as.matrix(fit$beta.samples)
  # Trait = the trait-table column each coefficient comes from (a categorical trait gives one coefficient per non-reference level)
  traits <- setdiff(names(postHoc_data[[nm]]$covs), "Species_ayerbe")
  trait_of <- map_chr(colnames(b), \(term) { hits <- traits[startsWith(term, traits)]; if (length(hits) == 0) NA_character_ else hits[which.max(nchar(hits))] })
  tibble(Contrast = contrast, Model = nm, Trait = trait_of, Term = colnames(b), Mean = colMeans(b), Lower = apply(b, 2, quantile, 0.025), Upper = apply(b, 2, quantile, 0.975), P_positive = colMeans(b > 0), N_species = nrow(postHoc_data[[nm]]$covs), Bayes_R2 = median(as.numeric(fit$bayes.R2)))
}) %>% list_rbind()) %>%
  list_rbind() %>%
  mutate(across(c(Mean, Lower, Upper, P_positive, Bayes_R2), \(x) round(x, 3)))
write_csv(postHoc_summary, paste0("Derived/Excels/postHoc_traits_", sub("\\.rds$", "", basename(Occ_mod_path)), ".csv"))

# Inspect Threatened species ------------------------------------------------------
tibble(Species_ayerbe = str_replace(names(y_draws), "_", " "), 
       Response = colMeans(y_draws)) %>% 
  left_join(Ft[, c("Species_ayerbe", "Gen_length")]) %>%
  ggplot(aes(x = Gen_length, y = Response)) +
  geom_point() + 
  geom_smooth()

# Select T&E species
Spp_te <- Ft %>% 
  mutate(Species_ayerbe_ = str_replace(Species_ayerbe, " ", "_")) %>%
  filter(Iucn_red_list != "LC") %>% 
  filter(Species_ayerbe_ %in% Occ_mod$sp.names) %>% 
  select(Species_ayerbe_, Iucn_red_list)

# Sample sizes of T&E species: the number of point-count visits (replicates) with a detection, the same count the species-inclusion cutoff uses (Num_detect_cutoff in 01_Wrangling_ms.R)
N_det <- apply(Occ_mod$y, 1, sum, na.rm = TRUE)
N_det_te <- tibble(Species_ayerbe = Spp_te$Species_ayerbe_, N = N_det[Spp_te$Species_ayerbe_])
arrange(N_det_te, N)

# Plot of their silvopasture vs pasture effects, with each species' N to the right of the panel
y_draws %>%
  select(Spp_te$Species_ayerbe_) %>%
  pivot_longer(cols = everything(), 
               names_to = "Species_ayerbe", 
               values_to = "value") %>%
  plot_half_eye(y_var = Species_ayerbe, Intercept = NA) +
  geom_text(data = N_det_te, aes(x = Inf, y = Species_ayerbe, label = paste("N =", N)), inherit.aes = FALSE, hjust = -0.1, vjust = -0.5, size = 3.5) +
  coord_cartesian(clip = "off") +
  labs(caption = "N = point-count visits with a detection") +
  theme(plot.margin = margin(5.5, 50, 5.5, 5.5))