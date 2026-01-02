## PhD birds in silvopastoral landscapes##
# Analysis occupancy modeling -- Occupancy modeling using packages unmarked, ubms, and spOccupancy

# NOTES:: See Books/Applied hierarchical modeling notes Word doc for theory

# Contents
# 1)
# 2)
# 3)
# 4)
# 5)
# 6)

# Load libraries & data ---------------------------------------------------
library(tidyverse)
library(readxl)
library(unmarked)
library(spAbundance)
library(spOccupancy)
library(ubms) # function get_stancode() could be used to get code from model
library(flocker)
library(brms)
library(AICcmodavg)
library(furrr) # Use furrr::future_map()
library(coda)
library(tictoc)
library(ggpubr)
library(cowplot)
library(conflicted)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

load("Rdata/the_basics_05.10.25.Rdata")
load("Rdata/Taxonomy_11.14.24.Rdata")
load("Rdata/Occ_abu_inputs_05.16.25.Rdata")
load("Rdata/Occ_abu_models_01.09.25.Rdata")

source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")

#NOTE:: The unique row identifier  [e.g. ID x Year] is critical , this defines what a row is in your dataframes and how many rows each dataframe will have 
Row_identifier <- c("Id_muestreo", "Ano_grp")

# Distributions -----------------------------------------------------------
# Define function to create histogram for count data
hist_counts <- function(df, zeros = TRUE, title = title){
  if(zeros == FALSE){
    df <- df %>% filter(Count > 0)
  }
  histogram <- df %>% ggplot() +
    geom_histogram(aes(x = Count)) +
    labs(x = "Species count", y = "Frequency", title = title)
  if(zeros == FALSE){
    #breaks <- seq(1, max(df$Count, na.rm = TRUE), by = 3) # Manually set breaks
    histogram <- histogram + 
      scale_x_continuous(breaks = function(x) c(.5, pretty(x)), 
                         limits = c(1, NA))
  }
  return(histogram)
}

# Create histogram for each species
# Zeros included
Abund_zeros_hist <- imap(Abund_in_range, \(abund, name){
  hist_counts(df = abund, zeros = TRUE, title = name)
})

# Zeros NOT included
Abund_hist <- imap(Abund_in_range, \(abund, name){
  hist_counts(df = abund, zeros = FALSE, title = name)
})

## Examine mean relative to variance
# The Poisson distribution is controlled by a single variable lambda, which represents the mean (and should be equal to the variance parameter as well) of the Poisson distribution. 
# Examine if mean = variance, or degree of dispersion

# Define calculate mean & dispersion function
calc_mean_disp <- function(vec, name){
  tibble(
    Species_ayerbe_ = name,
    mean_count = mean(vec, na.rm = TRUE),
    dispersion = var(vec, na.rm = TRUE) / mean_count,
    zeros = any(vec == 0)
  )}

# Calculate mean & dispersion for dataset with 0s and without 0s. I would imagine that the ZIP model sees dispersions similar to the no_zeros dataset? 
# When zeros are present, this brings down the mean count but also brings down variance in many cases
Mean_var <- imap(Abund_in_range, \(abund, name) {
  No_zeros <- abund %>% filter(Count != 0)
  nz_mean_disp <- calc_mean_disp(No_zeros$Count, name)
  mean_disp <- calc_mean_disp(abund$Count, name)
  rbind(nz_mean_disp, mean_disp)
}) %>% list_rbind() %>% 
  arrange(dispersion)

# Plot histogram of dispersion values
Mean_var %>% arrange(zeros) %>%
  filter(zeros == TRUE) %>% 
  ggplot() + 
  geom_histogram(aes(x = dispersion))

# Simulate data assuming that the data were poisson distributed (as if variance = observed mean)
SimPois <- data.frame(rpois(233, lambda = .57))
hist(rpois(233, lambda = .68), breaks = 40, main = "Simulated Poisson") 

# Models ------------------------------------------------------------------
# >Unmarked: ML ---------------------------------------------------------------
# Set up multiple cores with future 
parallelly::availableCores()
future::plan("multisession", workers = 8)

## Maximum likelihood (ml) with unmarked
# Covs to play around with: 
# Detection: j_day + Pc_start + Uniq_db + Habitat_cons
# Abundance: Ano1 + Habitat_cons + Ecoregion + Elevation + Temperature + Brandt forest raster

# Specify mixture and formula 
mixture <- "ZIP"
formula_ml <- as.formula("~1 ~1 + (1 | Id_group)")

## Safely fit the N-mixture model
# Specify safe_pcount
safe_pcount <- safely(
  ~ pcount(formula_ml, data = .x, mixture = mixture)
)
# Track time
start_time <- Sys.time()
fit_ml_safe <- map(umf_abu_l, safe_pcount, .progress = TRUE) #future_map
run_time <- Sys.time() - start_time
units_rt <- units(run_time) # Units run time
run_time <- round(as.numeric(run_time), 1)
paste(run_time, units_rt)

fit_ml <- map(fit_ml_safe, "result")
fit_ml

# >ubms: Bayesian ---------------------------------------------------------
if(FALSE){
# Pull the locations of the most abundant species, good candidates to try
Top_abu_ord <- Top_abu_df2 %>% slice_max(order_by = Tot_count, n = 30) %>%
  arrange(Species_ayerbe) %>%
  arrange(desc(Tot_count)) %>% 
  #slice_head(n = 6) %>% 
  pull(Species_ayerbe_)

# Mixture = "P" is only thing allowed at present
formula_b <- as.formula("~1 ~1 + (1 | Id_group)")

## Safely fit the N-mixture model
safe_stan_pcount <- quietly(
  ~ ubms::stan_pcount(formula_b, data = .x, cores = 8, chains = 4, iter = 2000)
)


tic()
fit_stan[1] <- map(umf_abu_l[Top_abu_ord[1]], safe_stan_pcount, .progress = TRUE)
toc()

fit_stan_result <- map(fit_stan, "result")
fit_stan_result # Notice 30 min vs 1 min run times
}

# Bulk Effective Samples Size (ESS) is too low
# Some Pareto k diagnostic values are too high
map(fit_stan, "warnings")
help('pareto-k-diagnostic')

warn_function <- function(x) {
  if (x < 0) warning("Negative value!")
  x^2
}

# Wrap it with `quietly`
quiet_warn_function <- quietly(warn_function)

-2:2 %>%
  map(quiet_warn_function) %>%
  map_dfr(~ tibble(
    result = .x$result,
    warnings = paste(.x$warnings, collapse = "; ")
  ))

map(results, "warnings")
warnings


# Diagnostics -------------------------------------------------------------
# Low p-value indicats inadequate fit
gof_ppc <- gof(fit_stan_result[[1]])
plot(gof_ppc) 

# There should be a shaded gray area, & if model fits would expect 95% of data to be within shaded area
ubms::plot_residuals(fit_stan_result[[1]], submodel = "state")

# Chain convergence
traceplot(fit_stan_result[[1]], pars = c("beta_state", "beta_det"))

## Can also create your own 
# param = "y" draws observations from posterior, "z" draws latent abundance values (i.e. predicted abundances)
sim_y <- ubms::posterior_predict(fit_stan_result[[1]], param = "y", draws = 100)
dim(sim_y)

# calculate the proportion of zeros in each simulated dataset
prop0 <- apply(sim_y, 1, function(x) mean(x==0, na.rm=TRUE))

actual_prop0 <- mean(getY(fit_stan_result[[1]]) == 0, na.rm=TRUE)

#Compare
hist(prop0, col='gray', xlim = c(.7, actual_prop0 + .01))
abline(v=actual_prop0, col='red', lwd=2)


# Communicate with stan purists
sc <- get_stancode(fit_stan_result[[1]])
cat(sc)

# >Unmarked ML ------------------------------------------------------------
#parboot(), crossVal(), residuals(), vif(), functions in AICcmodavg package
# Parametric bootstrap
gof <- Nmix.gof.test(mod = fit_ml$Bubulcus_ibis, nsim = 10)
gof <- Nmix.chisq(fit_ml$Bubulcus_ibis, nsim = 50)
?Nmix.gof.test

safely(Nmix.gof.test)
# These should produce a c-hat estimate. If c-hat > 1 indicate overdispersion (variance > mean). Values much higher than 1 (i.e., > 4) probably indicate lack-of-fit. Values close to zero as well.


## Troubleshoot
# None of these functions work
fitted(fit_ml$Manacus_manacus, na.rm = TRUE)
simulate(fit_ml$Bubulcus_ibis, nsim = 1)
crossVal(fit_ml$Bubulcus_ibis)

class(fit_ml$Bubulcus_ibis)
length(fit_ml$Bubulcus_ibis@sitesRemoved) #52 

Nmix.gof.test(bi_fit)

Site_covs

table(umf_abu_l$Bubulcus_ibis@obsCovs$Uniq_db)


fitstats <- function(fm, na.rm=TRUE) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2, na.rm=na.rm)
  chisq <- sum((observed - expected)^2 / expected, na.rm=na.rm)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=na.rm)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

chisq <- function(fm) {
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

fm <- fit_ml$Bubulcus_ibis

(pb <- parboot(fit_ml$Bubulcus_ibis, statistic=chisq, nsim=100, parallel=FALSE))

pb <- parboot(fit_ml$Bubulcus_ibis, nsim=30)
plot(pb)

fit_ml$Bubulcus_ibis@


safe_gof <- safely(~Nmix.gof.test(.x, nsim = 10))
gof_metrics <- map(fit_ml, safe_gof )

pb_metrics <- map(fit_ml, \(fit){
  safe_pb(fit = fit)
})

safe_pb <- safely(~parboot(.x, fitstats, nsim=25))

pb_metrics <- map(fit_ml, safe_pb)

# >>Errors ------------------------------------------------------------------
## Create dataframe to keep track of errors

# Errors observed thus far: 
# <simpleError in optim(tmb_mod$par, fn = tmb_mod$fn, gr = tmb_mod$gr, method = method,     ...): initial value in 'vmmin' is not finite>
# <simpleError in solve.default(full): Lapack routine dgesv: system is exactly singular: U[3,3] = 0>
# <simpleError in solve.default(full): system is computationally singular: reciprocal condition number = 1.00645e-58>
#

##TO DO 
# Add K, other model diagnostics, structure of data (Unique identifier = Id_yr_season) -- Export to an Excel sheet where each sheet is a day / model run.

fit_errors <- map(fit_ml_safe, "error") 

# Visualize error messages
just_errors_l <- discard(fit_errors, is.null) 
just_errors_l

# Define date & time 
date <- format(Sys.Date(), "%m.%d.%y")
time <- format(Sys.time(), "%H:%M")

## Generate spp diagnostics df
# Define function
create_error_tbl <- function(fit_errors, just_errors_l) {
  # Create the diagnostics table
  error_tbl <- tibble(
    Date = date, 
    Time = time,
    Species_ayerbe_ = names(fit_errors),
    Error_msg = map_chr(fit_errors, \(error) {
      str_extract_all(as.character(error), "vmmin|singular") %>%
        unlist() %>%
        { if (length(.) > 0) .[1] else NA_character_ }
    })
  )
  # Check if all error messages are captured
  if (length(just_errors_l) != nrow(drop_na(error_tbl))) {
    print("Error messages not captured, check just_errors_l object")
  } else {
    print("All error messages captured")
  }
  return(error_tbl)
}

# Create error table, extract model AICc
Spp_mods_info <- create_error_tbl(fit_errors, just_errors_l) %>% 
  #Add model AICc
  mutate(Aic_c = map_dbl(fit_ml, \(mod) {
    if (length(mod) > 0) round(AICcmodavg::AICc(mod), 0) else NA 
  })) 

# Join with spp_info
Spp_diagnostics <- Top_abu_df2 %>% left_join(Spp_mods_info) %>% 
  select(-Species_ayerbe_) %>%
  relocate(c(Date, Time), .before = 1)

## Metadata df
Model_metadata <- tibble(
  Date = date, 
  Time = time,
  Description = "Null det and abu",
  Uniq_row_id = paste0(Row_identifier, collapse = " X "),
  Det = paste0(formula_ml)[2], 
  Abu = paste0(formula_ml)[3], 
  Mixture = mixture,
  Num_errors = length(just_errors_l),
  Run_time = run_time, 
  Units = units_rt
  #Comment = "Lots of NAs in Std errors"
)

## Export Excel files
# Put in list 
Diagnostics_names <- c("Species_diagnostics", "Model_metadata")
Diagnostics_names <- setNames(Diagnostics_names, Diagnostics_names)
Model_diagnostics_l <- list(Spp_diagnostics, Model_metadata)
names(Model_diagnostics_l) <- Diagnostics_names

# Define function to add data from new model to existing Excel sheet tabs
append_to_excel_xlsx <- function(file_path, sheet_name, new_data) {
  # Check if the file exists
  if (file.exists(file_path)) {
    # Load the workbook
    wb <- loadWorkbook(file_path)
    # Get all existing sheet names
    existing_sheets <- getSheets(wb)
    
    if (sheet_name %in% names(existing_sheets)) {
      # Read existing data from the target sheet
      existing_data <- read.xlsx(file_path, sheetName = sheet_name)
      # Combine existing data with new data
      combined_data <- bind_rows(existing_data, new_data)
    } else {
      # If the sheet doesn't exist, use only new data
      combined_data <- new_data
    }
  } else {
    # If the file doesn't exist, create a new workbook
    wb <- createWorkbook()
    combined_data <- new_data
  }
  
  # Remove the sheet if it already exists (to replace it with updated data)
  if (sheet_name %in% names(existing_sheets)) {
    removeSheet(wb, sheetName = sheet_name)
  }
  
  # Add updated data to the workbook
  sheet <- createSheet(wb, sheetName = sheet_name)
  addDataFrame(combined_data, sheet, row.names = FALSE)
  
  # Save the workbook
  saveWorkbook(wb, file_path)
}

# Export to two tabs in the Excel file
path <- "Derived/Excels/Model_diagnostics.xlsx"
imap(Model_diagnostics_l, \(df, name){
  append_to_excel_xlsx(
    file_path = path,
    sheet_name = name,
    new_data = df
  )
})

# Trends in model fitting -------------------------------------------------
## Visualize changes in model performance 
# Import the updated Excel sheets
Diagnostics_import <- map(Diagnostics_names, \(sheet_name){
  read_xlsx(path, sheet = sheet_name)
})

# Create Mod_round, where each combination of Date_time gets a unique number
Diagnostics_imp2 <- map(Diagnostics_import, \(df){
  df %>% mutate(Date_time = paste0(Date, Time)) %>%
    mutate(Mod_round = dense_rank(Date_time),
           Mod_round = fct(as.character(Mod_round)))
})

Spp_diagnostics <- Diagnostics_imp2$Species_diagnostics %>% 
  mutate(Species_ayerbe_ = str_replace_all(Species_ayerbe, " ", "_")) %>% 
  relocate(Species_ayerbe_, .after = Species_ayerbe)

# Examine trends with AICc as we iterate through models
Spp_diagnostics %>% 
  ggplot(aes(x = Mod_round, y = Aic_c, color = Species_ayerbe)) + 
  geom_point() + 
  geom_line(aes(group = Species_ayerbe)) +
  guides(color = "none")

# Examine metadata associated with model rounds
Diagnostics_imp2$Model_metadata %>% 
  distinct(Date, Description, Det, Abu, Mixture, Run_time, Mod_round) %>% 
  view()

## See if there are predictors of lack of model fit 
Mean_var_zeros <- Mean_var %>% filter(zeros == TRUE)
Spp_diagnostics2 <- Spp_diagnostics %>%
  mutate(N_err_msg = sum(!is.na(Error_msg)), .by = Species_ayerbe_) %>% 
  left_join(Mean_var_zeros) %>%
  arrange(Species_ayerbe_) 

Spp_diagnostics2 %>%
  select(N_err_msg, where(is.numeric)) %>%
  names(.) %>%
  setdiff(c("N_err_msg", "Aic_c")) %>%  # Exclude columns you don't want to plot
  map(~ {
    ggplot(Spp_diagnostics2, aes(x = N_err_msg, y = .data[[.x]])) +
      geom_point() +
      geom_smooth() +
      labs(
        title = .x,
        x = "N_err_msg",
        y = .x
      )
  })

# Multi-species model
#?occuMulti()


# ML selection ------------------------------------
## From Table 3 Kellner et al (2023)
# Model selection
#modSel() crossVal()

mods <- list(null=mod_null, hab=mod_hab, year=mod_year, habyear=mod_habyear, habxyear=mod_habxyear)
mods_l <- fitList(fits=mods)
mods_table <- modSel(mods)
mods_table <- show(mods_table) # get data frame 


# ML Model inference --------------------------------------------------------
#predict(), ranef(), posteriorSamples()
#lambda2psi(0:5)
# backTransform(linearComb(fm2, coefficients = c(1,0,0), type = 'det'))

# Can only visualize the effect of changing one thing at a time 
plotEffects(mod_habxyear, 'lambda', covariate='Year') 

#You can get the data used to make this plot (e.g. to customize the plot yourself) using `plotEffectsData`:
head(plotEffectsData(mod_habxyear, 'lambda', covariate='Year'))

# ranef() function generates posterior distributions of latent abundance or occurrence (abundance in this example) for each site-year. bup turns it into a point estimate
r <- ranef(mod_habxyear)
stopifnot(round(bup(r)[1],3) == 8.138)
head(bup(r))




## Diagnostics and model fit - Bayesian 
# See shinystan package too: Shiny web application works with the output of MCMC 
# Traceplots of all parameters 
imap(fit_stan_result,
     ~ .x %>% traceplot(alpha = 0.6) + 
       ggtitle(.y)
)

# Residuals against fitted values (default), or a covariate. Draws = number of posterior draws
imap(fit_stan_result, ~ .x %>% plot_residuals(submodel = "state", draws = 9) + 
       ggtitle(.y)) #covariate = "Elev"

# Summary estimates table
sum_tab <- map(fit_stan_result, \(fit_stan){
  sum_tab <- summary(fit_stan, "state")
  #coef(fit_stan)["state[Tot.prec]"] # Return a specific covariate value
  return(sum_tab)
})

# Extract residuals for custom plotting
stan_resid <- map(fit_stan_result,
     ~ .x %>% residuals("state")
)

# Extract draws from the posterior distributions for custom fit needs
Post_draws <- map(fit_stan_result, ~ .x %>% ubms::extract())

# Chi square test of goodness of fit 
# LOOK INTO THIS - NO PREDICTIVE POWER?
gof <- map(fit_stan_result,
    ~ .x %>% gof()
)

## Model selection
mods <- fitList(fit_null, fit_global)
waic()
modSel() # LOO cross-validation

## Model inference Bayesian 
#predict(), posterior_predict(), plot_marginal()

# From vignette - 
# Compare ML vs stan
cbind(unmarked=coef(mod_ml), stan=coef(mod_b))


# Vignette ----------------------------------------------------------------
library(MuMIn) 
##single season
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##run model
fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest,
                     mallardUMF, K=30)

##compute observed chi-square
obs <- Nmix.chisq(fm.mallard)

##round to 4 digits after decimal point
print(obs, digits.vals = 4)

##compute observed chi-square, assess significance, and estimate c-hat
obs.boot <- Nmix.gof.test(fm.mallard, nsim = 10)
##note that more bootstrap samples are recommended
##(e.g., 1000, 5000, or 10 000)
obs.boot

# Dredge functionality
mallard.y[mallard.y>1] <- 1
mallardUMF <- unmarkedFrameOccu(mallard.y, siteCovs = mallard.site,
                                obsCovs = mallard.obs)
(ufm.mallard.global <- occu(~ date + I(date^2) ~ forest + elev,mallardUMF))
# examine all combinations
(d1 <- dredge(ufm.mallard.global, rank = AIC))

# Fit model with psi held constant
(ufm.mallard <- occu(~ date + I(date^2) ~ 1, mallardUMF))
# Find the best detection model
(d2 <- dredge(ufm.mallard, rank = AIC))
# keep date and date^2 in comparisons (#2 model from above)
(d2 <- dredge(ufm.mallard.global, rank = AIC,subset=`p(date)`&`p(I(date^2))`))


# >spOccupancy ------------------------------------------------------------
# PROBLEMS: 
# 1) Not working because every site is not unique (gaica ubc vs gaica mbd). Would need to go back in and change subsetting process 
# 2) Does not allow NAs in occurrence (site) covs, so would also need to subset out all GAICA distancia data if we want to keep habitat as a covariate 
    
spOcc_mod <- map(spOcc_ll, \(spOcc){
  PGOcc(
    occ.formula = ~ Elev + Tot.prec + Habitat_cons + Ano1 + (1 | Id_group_fac),
    det.formula = ~ j_day + Pc_start,
    data = spOcc[1:3], n.burn = 2000, n.samples = 4000,
    n.thin = 4, n.chains = 3, verbose = FALSE, k.fold = 2
  )
}, .progress = TRUE)
summary(spOcc_mod[[5]])
map(spOcc_mod, summary)

# Multi-season models
?stPGOcc # Single species
?stMsPGOcc # Multi species
?svcTMsPGOcc # Svc multi species

## We want a model that is... In the following list, numbers 1 & 2 are our inferential goals, 3, 4, & 5 are not primary inference goals but are likely important to accurately estimate the effects we are interested in. 
# 1) svc: Allows for parameter estimates to vary through space, good as we want to nail the parameter estimates for silvopasture, but have 5 distinct (& heterogeneous) ecoregions. From vignette: 'The model reduces to a traditional single-species occupancy model when all covariate eﬀects are assumed constant across space and a spatial occupancy model (Johnson et al. 2013; Doser et al. 2022) when only the intercept is assumed to vary across space'. So this is essentially random slope & random intercept
# 2) multi-season: Important to understand how species occupancy varies through time in the different habitat types (forest, silvopasture, pasture)
# 3) spatial: Account for spatial autocorrelation. From vignette, this is essentially a random intercept
# 4) latent factor: Account for residual species correlations via a factor modeling approach, which is a dimension reduction technique that can account for correlations among a large number of species without a massive computational cost
# 5) integrated: We have different data collectors with different abilities to detect birds by sight & (especially) by ear

# NOTE there are several other relevant arguments if we add coords back into the mix (see below)

spOcc_ll$Amazona_ochrocephala$occ.covs %>% filter(is.na(Habitat_cons)) %>% 
  pull(Id_group)

tabyl(spOcc_ll$Amazona_ochrocephala$occ.covs$Habitat_cons)

# >spAbundance ------------------------------------------------------------
# Doser suggests that including spatial autocorrelation & negative-binomial distribution will help with zero-inflation

# Useful functions
?spAbundance::sfMsNMix
?NMix
postHocLM # Explore how species traits relate to species-specific parameter estimates from a multi-species occupancy model

# Save --------------------------------------------------------------------
# Save model output
#save(fit_ml, fit_ml_safe, fit_stan, Error_l, Error_meta, Ls_count, file = paste0("Rdata/Occ_abu_models_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))

# Data simulation ---------------------------------------------------------
N <- 100
Tree_cov <- runif(N, 0, 1)
alpha <- 0 
beta1 <- 2

lambda <- round(exp(alpha + beta1 * Tree_cov), 2)
lambda

True_abu <- rpois(N, lambda)

# Detectability
time <- rnorm(N, 0, 1)



data.frame(Tree_cov = Tree_cov, lambda = lambda, N = N) %>% 
  arrange(lambda)

# Model Manacus using Unmarked ----------------------------------------------
# >Formatting -------------------------------------------------------------
# At least one problem with unmarked is that it cannot account for spatial auto correlation
PC_date <- Bird_pcs %>%
  distinct(Id_muestreo, Nombre_institucion, Uniq_db, Departamento, Ano, Mes, Dia, Fecha) %>%
  mutate(surveyNum = 1:n()) # Add in time of day once standardized
PC_date %>% count(Uniq_db)

# Removed Habitat_ut b/c this was creating multiple rows for Distanciamiento points
uniqPCs <- Bird_pcs %>%
  distinct(Id_muestreo, Fecha, Ano, Mes, Departamento, elev_fnExtract, Latitud, Longitud_decimal) %>%
  group_by(Id_muestreo) %>%
  mutate(surveyNum = 1:n()) %>%
  arrange(Id_muestreo)
nrow(uniqPCs)

manacus <- Bird_pcs %>%
  filter(Especie == "Manacus manacus") %>%
  group_by(Id_muestreo) %>%
  mutate(count = sum(Count)) %>%
  distinct(Especie, Id_muestreo, count) %>%
  arrange(count)
nrow(manacus) # Not sure why 10 rows are added when using summarize () over mutate()
data.frame(manacus) # Confirmed for at least M9-04 on 3/15/17 that they did indeed report 22 manacus manacus
manacus2 <- merge(manacus, filter(uniqPCs, Departamento == "Meta"), by = "Id_muestreo", all.y = T)
manacus2 <- manacus2 %>%
  mutate(count = replace_na(count, 0)) %>%
  filter(year == 2013) # & Aaron_Hab != "SSP"
## Scale variables to aid in model convergence
manacus2[, c("eventDate", "elev")] <- manacus2 %>%
  dplyr::select(eventDate, elev) %>%
  lapply(scale)
table(manacus2$year, manacus2$month) # Assuming that there's no change

# Create df w/ counts per point count as well as site level covariates of interest
manacus.y <- manacus2 %>%
  dplyr::select(count, Id_muestreo, elev_fnExtract, surveyNum) %>%
  pivot_wider(names_from = surveyNum, values_from = count, values_fill = NA, values_fn = sum)

manacus.y$tot.abu <- manacus.y %>%
  dplyr::select(-c(parcelID, Aaron_Hab, elev)) %>%
  rowSums(., na.rm = T)
table(manacus.y$Aaron_Hab)


# >Create UnMarkedFrame ---------------------------------------------------
# Approaches: Hierarchical Modelling of Species Communities or unmarked
# What I learned from this initial round of playing around... See email to Cole Burton w/ subject line ("Issues w/ my analysis in Unmarked'). It seems the issues detailed there were using the full data set (across all years) and when setting mixture = "NB", the negative binomial distribution. The ZIP also seemed to have elevated abundance estimates but nothing like the NB. This did not happen w/ the Poisson distribution (estimates of 15.. still way too high , but not so bad). Notice that model selection changed drastically when you switched from AIC to AICc, and quick search shows that AICc is better for small sample sizes and approaches AIC at large sample sizes (thus it is likely better in all cases). The high counts (e.g. of 22 manacus in a single point count) seem to really destabilize the model.



gaica$uniqPC <- paste(gaica$parcelID, gaica$eventDate, sep = "_")


# Issue w/ multiple habitat types in single point is messing things up down the line
gaica$Aaron_Hab <- replace(gaica$Aaron_Hab, gaica$parcelID == "M4-06" & gaica$Aaron_Hab == "SSP", "Bosque ripario")
gaica$Aaron_Hab <- factor(gaica$Aaron_Hab, levels = c("Bosque", "Borde de bosque", "Bosque ripario", "SSP"))

gaica %>%
  group_by(uniqPC) %>%
  mutate(time.s = scale(eventTime)) %>%
  distinct(uniqPC, parcelID, eventDate, Department, Aaron_Hab, lat, long)

uniqPCs <- data.frame(gaica %>% distinct(uniqPC, parcelID, eventDate, year, month, Department, Aaron_Hab, lat, long) %>% group_by(parcelID) %>% mutate(surveyNum = 1:n()) %>% arrange(parcelID))
nrow(uniqPCs) # 1467
ColElev <- raster::getData("alt", country = "COL", mask = T)
uniqPCs$elev <- extract(x = ColElev, y = cbind(uniqPCs$long, uniqPCs$lat))
manacus <- gaica %>%
  filter(scientificName == "Manacus manacus") %>%
  group_by(uniqPC) %>%
  mutate(count = sum(individualCount)) %>%
  distinct(scientificName, uniqPC, count) %>%
  arrange(count)
nrow(manacus) # Not sure why 10 rows are added when using summarize () over mutate()
data.frame(manacus) # Confirmed for at least M9-04 on 3/15/17 that they did indeed report 22 manacus manacus
manacus2 <- merge(manacus, filter(uniqPCs, Department == "Meta"), by = "uniqPC", all.y = T)
manacus2 <- manacus2 %>%
  mutate(count = replace_na(count, 0)) %>%
  filter(year == 2013) # & Aaron_Hab != "SSP"
## Scale variables to aid in model convergence
manacus2[, c("eventDate", "elev")] <- manacus2 %>%
  dplyr::select(eventDate, elev) %>%
  lapply(scale)
table(manacus2$year, manacus2$month) # Assuming that there's no change


# Create df w/ counts per point count as well as site level covariates of interest
manacus.y <- manacus2 %>%
  dplyr::select(count, parcelID, Aaron_Hab, elev, surveyNum) %>%
  pivot_wider(names_from = surveyNum, values_from = count, values_fill = NA, values_fn = sum)

manacus.y$tot.abu <- manacus.y %>%
  dplyr::select(-c(parcelID, Aaron_Hab, elev)) %>%
  rowSums(., na.rm = T)
table(manacus.y$Aaron_Hab)


# Create dfs w/ detection covariates
man.obs <- list()
man.obs[[1]] <- manacus2 %>%
  dplyr::select(parcelID, surveyNum, eventDate) %>%
  pivot_wider(names_from = surveyNum, values_from = eventDate, values_fill = NA) %>%
  dplyr::select(-parcelID)
# man.obs[[1]] <- sapply(man.obs[[1]], scale) Unnecessary
names(man.obs) <- "date"
## STILL NEED TO SCALE TIME AND THEN CREATE DATA FRAME, BUT THERE ARE ISSUES W/ TIME IN THE DATA SET SO NO NEED TO DO THIS YET
# man.obs[[2]] <- gaica %>% dplyr::select(uniqPC, eventTime) %>% mutate(time.s = scale(eventTime)) %>% pivot_wider(names_from= parcelID, values_from = time.s, values_fill = NA) %>% dplyr::select(-parcelID)
# For now...
man.obs[[2]] <- data.frame(replicate(9, manacus.y$elev)) # 9 or 13
names(man.obs) <- c("date", "elev")

# Format final data frames
man.y <- as.matrix(manacus.y[, c(4:12)]) # 4:12 or 4:16
man.site <- manacus.y[, c(2, 3)]
manUMF <- unmarkedFramePCount(y = man.y, siteCovs = man.site, obsCovs = man.obs) # obsCovs = man.obs[[1]]
summary(manUMF)


# >Abund model selection & plotting ----------------------------------------------
# First determine the best distribution. For Poisson, mean should equal the variance
nrow(manacus2)
(mean(manacus2$count) + var(manacus2$count)) / 2
var(manacus2$count) / mean(manacus2$count) # This (happens to be?) almost the same as c-hat (1.36), the overdispersion parameter from goodness-of-fit test
table(manacus2$count)
SimPois <- data.frame(rpois(233, lambda = .57))
hist(rpois(233, lambda = .68), breaks = 40, main = "Simulated Poisson") # .68 = Avg of mean and variance
hist(manacus2$count, breaks = 40, main = "Manakins observed")
hist(rbinom(n = 233, size = 4, prob = .35)) # prob calculated as (524-319)/524
(233 - 151) / 233
ggplot(data = manacus2, aes(x = count)) +
  geom_histogram(bins = 60)


### Equivalent global models varying the assumed abundance distribution
GlPoi <- pcount(~ date + elev ~ Aaron_Hab + elev, data = manUMF, mixture = "P", K = 50)
GlNb <- pcount(~ date + elev ~ Aaron_Hab + elev, data = manUMF, mixture = "NB", K = 50)
GlZip <- pcount(~ date + elev ~ Aaron_Hab + elev, data = manUMF, mixture = "ZIP", K = 50) # Zero-inflated Poisson
summary(GlPoi)
summary(GlNb)
summary(GlZip)

# ZIP is best model w/ AIC, and Poisson best when AICc
aictab(c(GlPoi, GlNb, GlZip), modnames = c("GlPoi", "GlNb", "GlZip"))

### Next let's determine the detection covariates
DateElev <- pcount(~ date + elev ~ 1, data = manUMF, mixture = "P", K = 50)
Date <- pcount(~date ~ 1, data = manUMF, mixture = "P", K = 50)
Elev <- pcount(~elev ~ 1, data = manUMF, mixture = "P", K = 50)
Null <- pcount(~1 ~ 1, data = manUMF, mixture = "P", K = 50)
summary(DateElev)

# Date & Elevation is best predictor of detection
aictab(c(DateElev, Date, Elev, Null), modnames = c("DateElev", "Date", "Elev", "Null"))

### Finally, let's determine the best predictors of abundance
HabElevAbu <- pcount(~ date + elev ~ Aaron_Hab + elev, data = manUMF, mixture = "P", K = 50)
HabAbu <- pcount(~ date + elev ~ Aaron_Hab, data = manUMF, mixture = "P", K = 50)
ElevAbu <- pcount(~ date + elev ~ elev, data = manUMF, mixture = "P", K = 50)
NullAbu <- pcount(~ date + elev ~ 1, data = manUMF, mixture = "P", K = 50)
summary(HabAbu)

# Habitat alone is best model
aictab(c(HabElevAbu, HabAbu, ElevAbu, NullAbu), modnames = c("HabElev", "Hab", "Elev", "Null"))

# Let's confirm that Poisson is still best model
HabAbuZIP <- pcount(~ date + elev ~ Aaron_Hab, data = manUMF, mixture = "ZIP", K = 50)
aictab(c(HabAbuZIP, HabAbu), modnames = c("HabZIP", "HabPois"))


# Can't have detection covs if you're trying to backTransform detection. E.g.
backTransform(NullAbu, type = "det") # No go
# Thus, take the top model (regarding predictors of abundance) and set detection formula to ~1, then backTransform
Hab_det <- pcount(~1 ~ Aaron_Hab, data = manUMF, mixture = "P", K = 50)
# Notice detection estimates differ depending on covs present
backTransform(Null, type = "det")
backTransform(Hab_det, type = "det")

NullAbu <- pcount(~ date + elev ~ 1, data = manUMF, mixture = "P", K = 50)
backTransform(NullAbu, type = "state")



# How does our best model fit the data? Chi-square Goodness-of-fit Test
HabAbuGOF <- Nmix.gof.test(HabAbu, nsim = 50)
HabAbuGOF
warnings()
?Nmix.gof.test
?parboot
## Predict off of best model

newData <- data.frame(Aaron_Hab = unique(gaica$Aaron_Hab))
pred <- predict(HabAbu, type = "state", newdata = newData, appendData = TRUE)
pred

# Predict abundance by habitat type
ggplot(data = pred, aes(x = Aaron_Hab, y = Predicted)) +
  geom_bar(stat = "identity", fill = "skyblue", alpha = 0.7) +
  geom_errorbar(aes(x = Aaron_Hab, ymin = lower, ymax = upper), width = 0.4, colour = "orange", alpha = 0.9, size = 1.3) +
  ylab("Predicted abundance") +
  scale_x_discrete(labels = c("Forest", "Forest edge", "Riparian Forest", "Silvopasture")) +
  theme(axis.title.x = element_blank()) #+ ggpubr::stat_compare_means()

# Predict detection probability by date. Doesn't seem like the habitat type matters
newdata <- data.frame(Aaron_Hab = factor("Borde de bosque", levels = c("Bosque", "Borde de bosque", "Bosque ripario", "SSP")), date = seq(-2, 2, length = 50), elev = 0)
EpDate <- predict(HabAbu, type = "det", newdata = newdata, appendData = TRUE)
qplot(date, Predicted, data = EpDate, geom = "line", xlab = "Scaled date", ylab = "Probability of detection") + geom_ribbon(aes(x = date, ymin = lower, ymax = upper), alpha = 0.1) + theme_bw()

# Predict detection probability by date
newdata2 <- data.frame(Aaron_Hab = factor("Bosque", levels = c("Bosque", "Borde de bosque", "Bosque ripario", "SSP")), elev = seq(-2, 2, length = 50), date = 0)
EpElev <- predict(HabAbu, type = "det", newdata = newdata2, appendData = TRUE)
qplot(elev, Predicted, data = EpElev, geom = "line", xlab = "Scaled elevation", ylab = "Probability of detection") + geom_ribbon(aes(x = elev, ymin = lower, ymax = upper), alpha = 0.1) + theme_bw()


# >Occupancy modeling ------------------------------------------------------
manacus2occ <- manacus2
manacus2occ$count <- ifelse(manacus2$count > 0, 1, 0)
manacus.yocc <- manacus2occ %>%
  dplyr::select(count, parcelID, Aaron_Hab, elev, surveyNum) %>%
  pivot_wider(names_from = surveyNum, values_from = count, values_fill = NA, values_fn = sum)
man.yocc <- as.matrix(manacus.yocc[, c(4:12)])
manUMFocc <- unmarkedFrameOccu(y = man.yocc, siteCovs = man.site, obsCovs = man.obs)
summary(manUMFocc)

# First determine the best distribution. For Poisson, mean should equal the variance
nrow(manacus2)
(mean(manacus2$count) + var(manacus2$count)) / 2
table(manacus2$count)
SimPois <- data.frame(rpois(233, lambda = .57))
hist(rpois(233, lambda = .68), breaks = 40) # .68 = Avg of mean and variance
hist(manacus2$count, breaks = 40)
hist(rbinom(n = 233, size = 4, prob = .35)) # prob calculated as (524-319)/524
(233 - 151) / 233
ggplot(data = manacus2, aes(x = count)) +
  geom_histogram(bins = 60)

### Equivalent global models varying the assumed abundance distribution
Gl.logit <- occu(~ date + elev ~ Aaron_Hab + elev, data = manUMFocc, linkPsi = "logit")
Gl.cll <- occu(~ date + elev ~ Aaron_Hab + elev, data = manUMFocc, linkPsi = "cloglog")
summary(Gl.logit)
summary(Gl.cll)

# ZIP is best model
aictab(c(Gl.logit, Gl.cll), modnames = c("Gl.logit", "Gl.cll"))

### Next let's determine the detection covariates. Add 0 for occupancy
DateElevO <- occu(~ date + elev ~ 1, data = manUMFocc)
DateO <- occu(~date ~ 1, data = manUMFocc)
ElevO <- occu(~elev ~ 1, data = manUMFocc)
NullO <- occu(~1 ~ 1, data = manUMFocc)

# Date & Elev is best predictor of detection
aictab(c(DateElevO, DateO, ElevO, NullO), modnames = c("DateElev", "Date", "Elev", "Null"))

### Finally, let's determine the best predictors of abundance.
HabElevOcc <- occu(~ date + elev ~ Aaron_Hab + elev, data = manUMFocc)
HabOcc <- occu(~ date + elev ~ Aaron_Hab, data = manUMFocc)
ElevOcc <- occu(~ date + elev ~ elev, data = manUMFocc)
NullOcc <- occu(~ date + elev ~ 1, data = manUMFocc)
summary(DateElev)

# Null best model
aictab(c(HabElevOcc, HabOcc, ElevOcc, NullOcc), modnames = c("HabElev", "Hab", "Elev", "Null"))

backTransform(NullO, type = "det")
str(ElevAbu)

# How does our best model fit the data?
?mb.gof.test
HabOccGOF <- mb.gof.test(ElevOcc, nsim = 50)
HabOccGOF
warnings()

## Predict abundance from best model
preddata <- unmarked::predict(ElevOcc, type = "state", appendData = TRUE)
qplot(elev, Predicted, data = preddata, geom = "line", xlab = "Scaled elevation", ylab = "Probability of occupancy") + geom_ribbon(aes(x = elev, ymin = lower, ymax = upper), alpha = 0.1) + theme_bw()

# Predict detection probability by date
newDataOcc <- data.frame(Aaron_Hab = factor("Bosque", levels = c("Bosque", "Borde de bosque", "Bosque ripario", "SSP")), date = seq(-2, 2, length = 50), elev = 0)
EpDate <- predict(ElevOcc, type = "det", newdata = newDataOcc, appendData = TRUE)
qplot(date, Predicted, data = EpElev, geom = "line", xlab = "Scaled date", ylab = "Probability of detection") + geom_ribbon(aes(x = date, ymin = lower, ymax = upper), alpha = 0.1) + theme_bw()

# Predict detection probability by date
newDataOcc2 <- data.frame(Aaron_Hab = factor("Bosque", levels = c("Bosque", "Borde de bosque", "Bosque ripario", "SSP")), elev = seq(-2, 2, length = 50), date = 0)
EpElev <- predict(ElevOcc, type = "det", newdata = newDataOcc2, appendData = TRUE)
qplot(elev, Predicted, data = EpElev, geom = "line", xlab = "Scaled elevation", ylab = "Probability of detection") + geom_ribbon(aes(x = elev, ymin = lower, ymax = upper), alpha = 0.1) + theme_bw()


# STILL TO DO: SHOULD REALLY BE A DYNAMIC OCCUPANCY MODEL (COLEXT FUNCTION)? OR NEED TO RUN 2 SEPERATE SINGLE SEASON MODELS?
?colext
?pcountOpen

# spOccupancy vignette -------------------------------------------------------------
vignette(package = "spOccupancy") #For some reason the vignettes are not linked with R package
#Site: https://doserlab.com/files/spoccupancy-web/

# For MSOM - 12 species, 373 sites, 3 replicates 
data(hbef2015)
str(hbef2015)
?hbef2015
lapply(hbef2015[c(1, 2, 4)], head)
lapply(hbef2015[c(3)], head)

# Single species occupancy model for Ovenbird
sp.names <- dimnames(hbef2015$y)[[1]]
ovenHBEF <- hbef2015
ovenHBEF$y <- ovenHBEF$y[sp.names == "OVEN", , ]
nrow(ovenHBEF$y)
table(ovenHBEF$y)

# Single species occupancy model for Ovenbird
sp.names <- dimnames(hbef2015$y)[[1]]
btbwHBEF <- hbef2015
btbwHBEF$y <- btbwHBEF$y[sp.names == "BTBW", , ]

# Implement with own data 
lapply(btbwHBEF[c(1, 2, 4)], head)

# Specify model formulas
btbw.occ.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
btbw.det.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)

# Run the model
out <- spPGOcc(
  occ.formula = btbw.occ.formula,
  det.formula = btbw.det.formula,
  data = btbwHBEF, n.batch = 400, batch.length = 25,
  accept.rate = 0.43, cov.model = "exponential",
  NNGP = TRUE, n.neighbors = 5, n.burn = 2000,
  n.thin = 4, n.chains = 3, verbose = FALSE, k.fold = 2
)
summary(out)

ppc.out <- ppcOcc(out, fit.stat = "freeman-tukey", group = 1)
summary(ppc.out)

# Run a spatial multi-species occupancy model with real data. Need to structure the data in the same format as the Hubbard Brooke experimental forest data.
# Key functions
?spMsPGOcc()
?intMsPGOcc
