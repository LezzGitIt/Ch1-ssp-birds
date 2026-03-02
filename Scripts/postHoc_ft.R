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
Occ_mod <- readRDS(paste0("Rdata/Model_output/Occ_mod", Sys.Date(), ".rds"))

# Format ------------------------------------------------------------------
# Rows correspond to the posterior MCMC samples and the columns correspond to species
y_draws <- as_draws_df(Occ_mod$beta.samples) %>% 
  select(starts_with("HabitatSsp")) %>% 
  rename_with(~str_remove(., "HabitatSsp-"))

# Subset functional traits
Spp <- str_replace(names(y_draws), "_", " ")
covariates <- c("Mass", "Habitat", "Range.Size")

Ft_spp <- Ft %>%
  filter(Species_ayerbe %in% Spp) %>%
  select(Species_ayerbe, all_of(covariates)) %>%
  mutate(Species_ayerbe = factor(Species_ayerbe, levels = Spp)) %>%
  arrange(Species_ayerbe) %>%
# postHocLM cannot accept any missing values in covariates
  filter(if_all(all_of(covariates), ~ !is.na(.x))) %>% 
  mutate(across(where(is.numeric), scale))

# Create list required 
postHoc_l <- list(y = y_draws, covs = Ft_spp)

# Run model ---------------------------------------------------------------
postHoc_mod <- postHocLM(formula = ~Mass + Habitat + Range.Size, data = postHoc_l)
summary(postHoc_mod)

# Plot --------------------------------------------------------------------
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

postHoc_mod$beta.samples %>% as_draws_df() %>%
  rename_with(~str_remove(., "Habitat")) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(param = as.factor(param)) %>%
  plot_half_eye()
