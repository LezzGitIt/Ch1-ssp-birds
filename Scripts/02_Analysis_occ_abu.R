## PhD birds in silvopastoral landscapes##
# Analysis occupancy modeling -- Occupancy modeling using packages unmarked, ubms, and spOccupancy

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
library(spOccupancy)
library(ubms) # function get_stancode() could be used to get code from model that could then be adapted
library(flocker)
library(brms)
library(ggpubr)
library(cowplot)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

load("Rdata/the_basics_12.24.24.Rdata")
load("Rdata/Taxonomy_11.14.24.Rdata")

source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")

# Data simulation ---------------------------------------------------------
# >Flocker -----------------------------------------------------------------
vignette("flocker_tutorial", package = "flocker")

# Try to simulate similar to AHM 10.4
simulate_flocker_data(
  n_rep = 3, n_pt = 100, n_sp = 1, n_season = 1, rep_constant = FALSE,
  covariates = ,
  params = list()
)

set.seed(1)
d <- simulate_flocker_data()
str(d)
names(d$params)
d$covariates
?simulate_flocker_data

d$event_covs

nrow(d$obs) # species-site combinations: 30 species * 50 sites
dim(d$obs) # 1500 x 4 repeat surveys

dim(d$unit_covs)
head(d$unit_covs)

## Fit repetition varying model
# Notice inclusion of event covars
fd_rep_varying <- make_flocker_data(
  obs = d$obs,
  unit_covs = d$unit_covs,
  event_covs = d$event_covs
)
?make_flocker_data

rep_varying <- flock(
  f_occ = ~ uc1 + (1 + uc1 | species),
  f_det = ~ uc1 + ec1 + (1 + uc1 + ec1 | species),
  flocker_data = fd_rep_varying,
  cores = 4,
  silent = 2,
  refresh = 0
)

## Fit repetition constant model (not event covars)
fd_rep_constant <- make_flocker_data(
  obs = d$obs,
  unit_covs = d$unit_covs
)


rep_constant <- flock(
  f_occ = ~ uc1 + (1 + uc1 | species),
  f_det = ~ uc1 + (1 + uc1 | species),
  flocker_data = fd_rep_constant,
  save_pars = save_pars(all = TRUE), # for loo with moment matching
  silent = 2,
  refresh = 0,
  cores = 4
)


# Modeling Manacus --------------------------------------------------------
# DELETE?
# At least one problem with unmarked is that it cannot account for spatial auto correlation
PC_date <- Bird_pcs %>%
  distinct(Id_muestreo, Nombre_institucion, Uniq_db, Departamento, Ano, Mes, Dia, Fecha) %>%
  mutate(surveyNum = 1:n()) # Add in time of day once standardized
PC_date %>% count(Uniq_db)

# Removed Habitat_homologado_ut b/c this was creating multiple rows for Distanciamiento points
uniqPCs <- Bird_pcs %>%
  distinct(Id_muestreo, Fecha, Ano, Mes, Departamento, elev_fnExtract, Latitud_decimal, Longitud_decimal) %>%
  group_by(Id_muestreo) %>%
  mutate(surveyNum = 1:n()) %>%
  arrange(Id_muestreo)
nrow(uniqPCs)

manacus <- Bird_pcs %>%
  filter(Especie == "Manacus manacus") %>%
  group_by(Id_muestreo) %>%
  mutate(count = sum(Numero_individuos)) %>%
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


# Create UnMarkedFrame ---------------------------------------------------
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


# Abund model selection & plotting ----------------------------------------------
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


# Occupancy modeling ------------------------------------------------------

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


# Habitat preferences -----------------------------------------------------
# Generate habitat preferences for 580 of 607 species observed in our point counts. Note some species are not able to be matched b/c of differences in taxonomy.. Would have to generate complete list of equivalents for Ayerbe -> Species_bl
library(traits) # Traits data
# traits::traitbank() #also see functions related to EOL

load("Hab_types.Rdata")

HBW <- data.frame(read_excel("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Analysis/Taxonomy/Handbook of the Birds of the World and BirdLife International Digital Checklist of the Birds of the World_Version_7.xlsx", sheet = "HBW-BirdLife v7 ", skip = 3))
HBWsp <- HBW %>% filter(`Subsp..Seq.` == 0) # Reduce file down just to recognized species
head(HBWsp)
# subset just relevant species of Colombia
HBWco <- distinct(Tax_df[c("Species_bl")]) %>%
  inner_join(HBWsp[c("Scientific name", "SISRecID")], join_by("Species_bl" == "Scientific name"))

HBW_hab_pref <- list()
for (i in 1:nrow(HBWco)) {
  print(i)
  HBW_hab_pref[[i]] <- birdlife_habitat(id = HBWco[i, 2])
}
Hab_types <- bind_rows(HBW_hab_pref) %>% inner_join(HBWco, join_by("id" == "SISRecID"))
head(Hab_types)
lapply(Hab_types[2:3], table)


# save(Hab_types, file = "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/R_Files/PhD/Hab_types.Rdata")

birdlife_threats(22689248)

# spOccupancy -------------------------------------------------------------
data(hbef2015)
?hbef2015
lapply(hbef2015[c(1, 2, 4)], head)
lapply(hbef2015[c(3)], head)

sp.names <- dimnames(hbef2015$y)[[1]]
btbwHBEF <- hbef2015
btbwHBEF$y <- btbwHBEF$y[sp.names == "BTBW", , ]

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
?spMsPGOcc()
