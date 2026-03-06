## PhD birds in silvopastoral landscapes ##
# Post-hoc functional traits model 

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(posterior)
library(spOccupancy)
library(ggdist)

Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"
Ft <- read_csv(paste0(Wrangling_repo, Excels, "Functional_traits.csv"))

# Occupancy model output
Occ_mod <- readRDS("Rdata/Model_output/Occ_mod_52spp_2026-03-05.rds")

# Forest or non-forest? ---------------------------------------------------
Forest <- TRUE

# Format ------------------------------------------------------------------
# Rows correspond to the posterior MCMC samples and the columns correspond to species
if(Forest){
  variable <- "Canopy_cover"
  y_draws <- as_draws_df(Occ_mod$beta.samples) %>% 
    select(starts_with(variable)) %>% 
    rename_with(~str_remove(., paste0(variable, "-")))
} else{
  variable <- "HabitatSsp"
  y_draws <- as_draws_df(Occ_mod$beta.samples) %>% 
    select(starts_with(variable)) %>% 
    rename_with(~str_remove(., paste0(variable, "-")))
}

# Subset functional traits
Spp <- str_replace(names(y_draws), "_", " ")
covariates <- c("Mass", "Elev_range_final", "Range.Size")

str_c(letters, collapse = ", ")
str_c(covariates, collapsae = "")
paste0(covariates, collapsae = " + ")

# PCA of bill morph
pca <- prcomp(data = Ft, ~ Beak.Depth + Beak.Length_Nares + Beak.Width)
Ft$Beak_size <- pca$x[,1]
Ft$Beak_shape <- pca$x[,2]

# Plot PCA


Ft_spp <- Ft %>%
  filter(Species_ayerbe %in% Spp) %>%
  rename(HWI = `Hand-Wing.Index`) %>%
  #select(Species_ayerbe, all_of(covariates)) %>%
  mutate(Species_ayerbe = factor(Species_ayerbe, levels = Spp)) %>%
  arrange(Species_ayerbe) %>%
# postHocLM cannot accept any missing values in covariates
  #filter(if_all(all_of(cqovariates), ~ !is.na(.x))) %>% 
  mutate(across(where(is.numeric), scale))

# Create list required by postHocLM
postHoc_l <- list(y = y_draws, covs = Ft_spp)

# >Eye --------------------------------------------------------------------
Ft_eye <- Ft_spp %>% filter(!is.na(Eye_resid)) %>% 
  select(-Clutch)
Spp_eye <- Ft_eye %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_eye <- y_draws %>% select(all_of(Spp_eye)) 

# List required by postHocLM
postHoc_l_eye <- list(y = y_draws_eye, covs = Ft_eye)

# >Clutch -----------------------------------------------------------------
Ft_clutch <- Ft_spp %>% filter(!is.na(Clutch)) %>% 
  select(-c(Eye_resid, Source_eye))
Spp_clutch <- Ft_clutch %>% pull(Species_ayerbe) %>% 
  str_replace(" ", "_")
y_draws_clutch <- y_draws %>% select(all_of(Spp_clutch))

# List required by postHocLM
postHoc_l_clutch <- list(y = y_draws_clutch, covs = Ft_clutch)

# Run model ---------------------------------------------------------------
## Morphology models
# Morph variables: Mass + HWI + Tarsus.Length + Tail.Length + Beak_size + Beak_shape
postHoc_morph <- postHocLM(
  formula = ~Mass , data = postHoc_l
)
# Residual eye size
postHoc_eye <- postHocLM(formula = ~Eye_resid, data = postHoc_l_eye)

## Life-history models
# Lh vars: Range.Size + Elev_range_final + Habitat.Density + Migration + Trophic.Level
postHoc_lh <- postHocLM(
  formula = ~Range.Size, data = postHoc_l
  )
# Clutch
postHoc_clutch <- postHocLM(formula = ~Clutch, data = postHoc_l_clutch, verbose = TRUE)

# Plot --------------------------------------------------------------------
format.for.plotting <- function(samples){
  samples %>% as_draws_df() %>%
    rename_with(~str_remove(., "Habitat")) %>%
    pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
    mutate(param = as.factor(param))
}

plot_half_eye <- function(df){
  df %>% 
    filter(!param %in% c(".iteration", ".draw", ".chain")) %>%
    ggplot(aes(x = value, y = fct_rev(param), fill = param)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::stat_halfeye(alpha = .8) +
    guides(fill = "none") +
    labs(x = NULL, y = NULL) +
    theme_minimal()
}

?postHocLM
## Morphology
format.for.plotting(postHoc_morph$beta.samples) %>% plot_half_eye()
# Residual eye 
format.for.plotting(postHoc_eye$beta.samples) %>% plot_half_eye()

## Life-history
format.for.plotting(postHoc_lh$beta.samples) %>% plot_half_eye()
# Clutch size
format.for.plotting(postHoc_clutch$beta.samples) %>% plot_half_eye()
 