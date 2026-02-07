## Multi-species wrangling for spOccupancy

# Load libraries & data ---------------------------------------------------
## Load libraries
library(tidyverse)
library(hms)
library(janitor)
library(chron)
library(readxl)
library(sf)
library(cowplot)
library(conflicted)
library(ebirdst)
library(mgcv)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

## Load data
Wrangling_repo <- "../Ssp-bird-data-wrangling/"
Excels <- "Derived/Excels/"

## NOTE: Removing Hatico points until can digitize
Hatico_pc_nums <- str_pad(1:12, width = 2, pad = 0)

Pc_locs_dc_sf <- st_read(
  paste0(Wrangling_repo, "Derived_geospatial/shp/Pc_locs_dc.gpkg")
) %>% filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Bird_pcs_analysis <- read_csv(
  paste0(Wrangling_repo, Excels, "Bird_pcs/Bird_pcs_analysis.csv")
) %>% filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Event_covs <- read_csv(paste0(Wrangling_repo, Excels, "Event_covs.csv")) %>% 
  filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Site_covs <- read_csv(paste0(Wrangling_repo, Excels, "Site_covs.csv")) %>% 
  filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Ft <- read_csv(paste0(Wrangling_repo, Excels, "Functional_traits.csv"))

# Remove these points 
#aj_tbl <- Event_covs %>% filter(Uniq_db == "Gaica mbd" & Rep_season > 4)
#Event_covs <- Event_covs %>% anti_join(aj_tbl)
#Bird_pcs_analysis <- Bird_pcs_analysis %>% anti_join(aj_tbl)

# Breeding records matrix from Moreno-palacios
Breed_records <- read_csv("../Datasets_external/Palacios_breeding_database.csv") %>% 
  mutate(across(where(is.character), as.factor), 
         date = paste(year, month, day, sep = "-"),
         date = as.Date(date, format = "%Y-%m-%d"),
         doy = yday(date))

load(paste0(Wrangling_repo, "Rdata/NE_layers_Colombia.Rdata")) 
source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/Rcookbook/Themes_funs.R")

# Row identifier ----------------------------------------------------------
Event_covs %>% 
  select(Id_muestreo, Ano_grp, N_samp_periods, Season, 
         contains("rep", ignore.case = TRUE)) 

#NOTE:: The unique row identifier  [e.g. ID x Year] is critical , this defines what a row is in your dataframes and how many rows each dataframe will have 
Row_identifier <- c("Id_muestreo", "Ano_grp", "Season") #, "Season"

# NOTE:: This is a REMINDER, need to replace all 'Rep ano grp' with 'Rep season' and viceaversa depending on whether "Season" is included in the Row_identifier
Rep_identifier <- "Season" # Ano_grp
paste0("Rep_", Rep_identifier)

# Final formatting --------------------------------------------------------
Bird_pcs_analysis2 <- Bird_pcs_analysis %>%
  mutate(Species_ayerbe_ = str_replace_all(Species_ayerbe, " ", "_"))

## Add Rep_season to tbl 
# Pc_start & Rep_season are equivalent in terms of grouping (when combined with Id_muestreo & Ano_grp), but since we are using Rep_season for formatting of analysis dataframes it makes sense to include Rep_season as well
date_join_spp_obs <- Event_covs %>% 
  select(Id_muestreo, Ano_grp, Season, Fecha, Pc_start, Rep_season)
Bird_pcs_analysis3 <- Bird_pcs_analysis2 %>% 
  left_join(date_join_spp_obs)

# Species selection and ordering ----------------------------------------------
# From spOccupancy vignette: 'Place a common species first. The first species has all of its factor loadings set to fixed values, and so it can have a large influence on the resulting interpretation on the factor loadings and latent factors... For the remaining 𝑞−1 factors, place species that you believe will show different occurrence patterns than the first species'
# Personal - It seems that fitting widespread species results in better fit 

## Join with functional traits 
Ft_abu <- Bird_pcs_analysis3 %>% 
  count(Species_ayerbe, sort = T) %>% 
  left_join(Ft)

## Set minimum number of distinct point counts which in turn sets the number of species
Num_locs_cutoff <- 100
Ft_abu2 <- Ft_abu %>% filter(n > Num_locs_cutoff)
nrow(Ft_abu2) # 1 loc = 462 species, 10 = 251 species, 25 = 160, 40 = 117, 50 = 100 species, 100 = 39 species

## We want to find a few distinct species to place at the top of the species x site matrix. We will examine body size, habitat associations, trophic levels, and range size
Ft_abu3 <- Ft_abu2 %>% select(
  Species_ayerbe, n, Mass, contains(c("Habitat", "Trophic")), Range.Size
) %>% arrange(desc(n))
#Ft_abu3 %>% view()

# Habitat association - Forest, Shrubland, Wetland, grassland, human modified 
Ft_abu3 %>% count(Habitat, sort = T) %>% 
  filter(n > 1)
# Trophic levels - c("Carnivore", "Herbivore", "Omnivore") # Granivore, Nectarivore
Ft_abu3 %>% count(Trophic.Level, sort = T) 
Ft_abu3 %>% count(Trophic.Niche, sort = T) 

# Visualize body size
Ft_abu3 %>% ggplot() + 
  geom_histogram(aes(x = Mass))
# Select representative species - DONT CHANGE THIS ONE (WORKS)
Order_spp <- c("Tyrannus_melancholicus", "Synallaxis_azarae", "Milvago_chimachima", "Thraupis_episcopus", "Camptostoma_obsoletum", "Volatina_jacarina", "Troglodytes_aedon", "Pitangus_sulphuratus", "Bubulcus_ibis") # Other options could include "Milvago_chimachima", "Pteroglossus_castanotis", Ortalis_garrula, Manacus_manacus, Todirostrum_cinereum,  Columbina_talpacoti, Pitangus_sulphuratus
#Order_spp <- c("Tyrannus_melancholicus", "Synallaxis_azarae", "Bubulcus_ibis", "Troglodytes_aedon", "Volatina jacarina", "Troglodytes_aedon", "Camptostoma_obsoletum", "Columbina_talpacoti") 
Ft_abu3 %>% filter(Species_ayerbe %in% Order_spp)

# Remove species under the minimum number of distinct point counts 
Bird_pcs_analysis4 <- Bird_pcs_analysis3 %>% 
  filter(Species_ayerbe %in% Ft_abu2$Species_ayerbe)

# Site covariates ---------------------------------------------------------
## Covariates that are fixed for a given point count x Ano_grp combination
date_join_ano <- Event_covs %>% 
  distinct(Id_muestreo_no_dc, across(all_of(Row_identifier)), Ano, Uniq_db, Canopy_cover, Canopy_height_m, ssp, te)
Site_covs2 <- Site_covs %>%
  full_join(date_join_ano) %>% 
  # Change these two points to "Mosaic". This allows them to be included without biasing any of the landcovers we care more about (ssp, forest)
  mutate(
    Habitat = if_else(Habitat == "Cultivos", "Mosaic", Habitat)
  ) %>% 
  # Some GAICA points were sampled 2x in a single Ano_grp across different years. Take the measurement from 2016.
  slice_max(by = all_of(Row_identifier), order_by = Ano) %>%
  arrange(across(all_of(Row_identifier))) %>%
  mutate(Id_gcs = ifelse(Id_group_no_dc == "MB-R-OQ", "Oq_farm_id", Id_gcs)) %>%
  mutate(
    # Take the earlier of the two years so it is numeric
    Ano1 = as.numeric(str_split_i(Ano_grp, "-", 1)),
    across(where(is.character), as.factor),
    Habitat = fct_relevel(Habitat, c("Pastizales", "Mosaic", "Ssp", "Bosque")), #"Bosque ripario"
    Ecoregion = relevel(Ecoregion, ref = "Piedemonte")
  ) %>%
  select(all_of(Row_identifier), Ano1, Id_gcs, Id_group_no_dc, Uniq_db, Ecoregion, Elev, Tot_prec, Habitat, Uniq_db, Canopy_cover, Canopy_height_m, ssp, te)

# Scale numeric values, create categorical random effects
Site_covs3 <- Site_covs2 %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.))), 
# Categorical random effects have to be specified as numeric in spOccupancy
  Id_group_no_dc = as.numeric(Id_group_no_dc),
  Id_gcs = as.numeric(Id_gcs))

# >Rm habitat = cultivos?  -------------------------------------------------
# There are two points with Habitat == "Cultivos". These can be removed or can be changed to "Mosaic". Here I change them to "Mosaic", which allows them to be included without biasing any of the landcovers we care more about (ssp, forest)

# Adjust Row_identifier
#if(any(str_detect(Row_identifier, "Ano_grp"))){
#  Row_identifier <- str_replace_all(Row_identifier, "_grp", "1")
#}

#Rm_tbl <- Site_covs2 %>% filter(is.na(Habitat)) %>% 
  #distinct(across(all_of(Row_identifier)))

# Observation covariates ----------------------------------------------------
## Create and format key detection covariates

# >Breeding probability ---------------------------------------------------
# Moreno-palacios original model had guild , but remove here given that spOccupancy can't incorporate species level covariates
# To try - consider Ecoregion to region mapping; scale within each region instead of between regions; interaction between habitat type * breed_prob

if(FALSE){ # Slow
  GAM_region <- gam(
    breedingCor ~ s(doy, by=region, bs="cc", k=-1) + s(year, bs="re") + Elevation + latitude, # + guild
    family = "binomial", 
    data = Breed_records[Breed_records$status =="Native",], 
    method = "REML")
  saveRDS(GAM_region, file = "Rdata/Palacios_GAM_region.rds")
}
GAM_region <- readRDS("Rdata/Palacios_GAM_region.rds")

# Isolate the unique sites X julian day combinations along with relevant covariates
Sampled_sites <- Event_covs %>% 
  left_join(Site_covs) %>%
  distinct(Ecoregion, Julian_day, Elev, Lat)
# Rename to match with Moreno-Palacios model
newdat_sites <- Sampled_sites %>% 
  rename(doy = Julian_day,
         Elevation = Elev, 
         latitude = Lat) %>% 
  mutate(year = round(mean(Breed_records$year, na.rm = TRUE), 0))

# Generate table of equivalent names between SCR ecoregions and Moreno-palacios regions
Ecoregion_equiv <- tibble(
  Ecoregion = unique(Site_covs$Ecoregion), 
  region = c("Caribbean", "Andes", "Andes", "Caribbean", "Orinoquia")
)
newdat_sites2 <- newdat_sites %>% left_join(Ecoregion_equiv)

## Predict the Probability of breeding using Moreno-Palacios model
# The response here (Prob_breed) represents probability an individual (or observation) shows breeding activity on that day
Prob_breed_tbl <- newdat_sites2 %>% 
  mutate(Prob_breed = predict(GAM_region, newdat_sites2, type = "response"),
         Prob_breed = as.numeric(Prob_breed)) %>%
  # Change names back so that they match with Event_covs
  rename(Julian_day = doy, Elev = Elevation, Lat = latitude) %>%
  select(-c(year, region))

# Join to Event_covs
Event_covs2 <- Site_covs %>% 
  select(Id_muestreo_no_dc, Elev, Lat) %>%
  right_join(Event_covs) %>%
  left_join(Prob_breed_tbl) %>% 
  select(-c(Elev, Lat))

# >Binary forest ----------------------------------------------------------
# Create binary forest column for detection covariate
Site_covs_forest <- Site_covs %>% 
  mutate(Forest_bin = ifelse(Habitat == "Bosque", "Y", "N"))
Event_covs3 <- Event_covs2 %>% 
  left_join(Site_covs_forest[,c("Id_muestreo_no_dc", "Forest_bin")]) #%>%
  #mutate(Julian_day2 = poly(Julian_day, 2),
   #      Pc_start2 = poly(Pc_start, 2))

# >Format ------------------------------------------------------------------
# Generate Obs_covs_df where each row is a unique ID [Id_muestreo, Ano_grp], each column is a site visit, and the 'Variable' column identifies which Observation covariate the row corresponds to  
# Lengthen then widen dataframe
Obs_covs_df <- Event_covs3 %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(
    # Select detection covariates here
    cols = c(
      Nombre_institucion, Prob_breed, Pc_start, Pc_length, Sampling_day, Forest_bin #Julian_day2, Pc_start2
    ),      
    names_to = "Variable",          
    values_to = "Value"             
  ) %>%
  arrange(across(all_of(Row_identifier[-1])), Rep_season) %>%
  pivot_wider(
    id_cols = c(all_of(Row_identifier), Variable),
    names_from = c(Rep_season),
    names_glue = "Rep_season{Rep_season}",
    values_from = c(Value)
  ) 

# Split Obs_covs_df into a list, where each slot in the list is a covariate 
Obs_covs_grp <- Obs_covs_df %>%
  arrange(across(all_of(Row_identifier))) %>%
  group_by(Variable) 
Obs_covs_l <- Obs_covs_grp %>% group_split(.keep = FALSE)
keys <- Obs_covs_grp %>% 
  group_keys() %>% 
  pull(Variable)
names(Obs_covs_l) <- keys

# Define functions to convert type back to original type (date , time, character, etc.)
Type_tbl <- Event_covs3 %>% 
  select(Nombre_institucion, Prob_breed, Pc_start, Pc_length, Sampling_day, Forest_bin) %>% 
  map_dfr(., class) %>% 
  filter(Pc_start == "hms") %>% 
  pivot_longer(everything(), names_to = "Variable", values_to = "Type") %>% 
  arrange(Variable)
Type_vec <- Type_tbl %>% pull(Type)
convert_type_chr <- paste("as", Type_vec, sep = ".")
convert_type <- map(convert_type_chr, rlang::as_function)

# Apply convert_type to variables in list
Obs_covs_l2 <- map2(Obs_covs_l, convert_type, \(Obs_covs, type) {
  # Apply transformations on the first 3 elements
  Obs_covs %>% mutate(across(-all_of(Row_identifier), ~ type(.))) %>% 
    select(-all_of(Row_identifier)) %>%                    
    mutate(across(where(~ is.numeric(.) | inherits(., "hms")), ~ scale(.)))       
})

# Abundance formatting ----------------------------------------------------
## Join abundance data (just containing points where species was observed) with all the point counts that the species could have been observed at.

# Event_covs contains all points surveyed (including Spp_obs == 0)
date_join_rep <- Event_covs %>% 
  distinct(across(all_of(Row_identifier)), Rep_season)

# NOTE:: The right_join effectively adds NAs for Count where species could have been observed but weren't (lengthening the dataframe)
Abund_no_obs <- Bird_pcs_analysis4 %>%
  select(all_of(Row_identifier), Rep_season, Species_ayerbe_, Count) %>% 
  right_join(date_join_rep) %>% 
  # Change NAs to 0s
  mutate(Count = if_else(is.na(Count), 0, Count))

# NOTE:: There are only 56 point counts where Spp_obs == 0, but there are another 100+ points where all the species observed were >50m and thus were removed from the abundance data. 
Event_covs %>% filter(Spp_obs == 0)
Abund_no_obs %>% filter(is.na(Species_ayerbe_))

Abund_no_obs2 <- Abund_no_obs %>%
  summarize(
    Count = sum(Count),
    .by = c(Species_ayerbe_, all_of(Row_identifier), Rep_season)
  )

Abund_wide_l <- Abund_no_obs2 %>% 
  unite(Id_muestreo_ano, all_of(Row_identifier), sep = ".", remove = FALSE) %>%
  arrange(across(all_of(Row_identifier))) %>%
  select(-all_of(Row_identifier)) %>%
  pivot_wider(
    names_from = Id_muestreo_ano,
    values_from = Count,
    values_fill = 0
  ) %>%
  split(.$Rep_season)

# Generate a species list so that each replicate tbl (1-5) will have the full species list
Species_list <- Bird_pcs_analysis4 %>% distinct(Species_ayerbe_)
# Join with Species_list and convert NAs to 0s
Abund_rep_l <- map(Abund_wide_l, \(rep_df){
  rep_df %>% select(-Rep_season) %>% 
    right_join(Species_list) %>% 
    # Reorder species for more efficient estimation of latent factors
    arrange(match(Species_ayerbe_, Order_spp)) %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    column_to_rownames(var = "Species_ayerbe_")
  })

rownames(Abund_rep_l$`1`)[c(38, 14, 8, 43, 5, 6, 7)]

# Any of the observation covs (start time, pc length, or date) dataframes have NAs where a given point count wasn't surveyed at a given repetition. Use this data frame to assign NAs to the abundance dataframe
ncols_det <- ncol(Obs_covs_l$Pc_start)
start_det_cols <- length(Row_identifier) + 1

TF_mat <- !is.na(Obs_covs_l$Pc_start[,start_det_cols:ncols_det])  #659 rows

Abund_nas <- map(seq_along(Abund_rep_l), \(replicate){
  
  tbl <- Abund_rep_l[[replicate]]      # species × sites
  keep <- TF_mat[, replicate]     # length = sites
  
  tbl[, !keep] <- NA
  tbl
})

# Turn list (each slot is a replicate) into an array (each matrix slice is a replicate)
Abund_nas_ul <- unlist(Abund_nas)

# Ordering of unlist -> array requires filling the array by species, sites, and then replicates
sppXsite <- dim(Abund_rep_l[[1]])
n_replicates <- length(Abund_rep_l)

Abund_array <- array(Abund_nas_ul, dim = c(sppXsite[1], sppXsite[2], n_replicates))

# Understand Abundance array 
dim(Abund_array)
Abund_array[1,,1] # First species, all sites, first replicate
Abund_array[,1,1] # All species, first site, first replicate
Abund_array[,1,2] # Second replicate is all NAs because CIPAV only surveyed once

# Coordinates -------------------------------------------------------------

add_year_join <- Site_covs2 %>% distinct(across(all_of(Row_identifier)))

# Arrange coordinates 
Coords_arr <- Pc_locs_dc_sf %>% 
  distinct(Id_muestreo, geom) %>% 
  # Add year 
  left_join(add_year_join) %>%
  #anti_join(Rm_tbl) %>%
  arrange(across(all_of(Row_identifier)))
# Project coordinates for Colombia 
Coords <- Coords_arr %>% 
  st_transform(crs = st_crs("EPSG:32618")) %>% 
  st_jitter(amount = .001) %>% 
  st_coordinates()

# Biogeographic clipping --------------------------------------------------
## Idea from Socolar's (2022) paper: Biogeographic multi-species occupancy models for large-scale survey data. We only want to include point count locations that are within the species range (+ some buffer), differentiating a true zero (possible but not observed) vs points that are simply out of range (more like an NA).

# Create species list of those that have maps
No_maps <- c("Leptotila verreauxi", "Accipiter bicolor", "Thripadectes virgaticeps")
Species_list2 <- Species_list %>% 
  mutate(Species_ayerbe = str_replace(Species_ayerbe_, "_", " ")) 
Species_list_maps <- Species_list2 %>% 
  filter(!Species_ayerbe %in% No_maps) %>%
  pull(Species_ayerbe)

# Load in relevant shapefiles
safe_st_read <- safely(st_read)
Ayerbe_mod_spp_l <- map(Species_list2$Species_ayerbe, \(spp){
  safe_st_read(dsn = "../Geospatial_data/Ayerbe_shapefiles_1890spp", layer = spp)
})
Ayerbe_mod_spp_l2 <- map(Ayerbe_mod_spp_l, "result")
names(Ayerbe_mod_spp_l2) <- Species_list2$Species_ayerbe

## Leptotila verreauxi
# NOTE that Leptotila verreauxi is NULL because it is not in the Ayerbe shapefiles
map_vec(Ayerbe_mod_spp_l2, ~is.null(.x))

# Use ebird map for Leptotila verreauxi
whtdov_path <- ebirdst_download_status(
  species = "Leptotila verreauxi", download_ranges = TRUE, pattern = "range_raw_9km"
  ) 
whtdov_range <- load_ranges("whtdov", resolution = "9km", smoothed = F) %>% 
  select(scientific_name, geom) %>% 
  rename(Nombre = scientific_name, 
         geometry = geom)
# Crop to just Colombia (Not really necessary)
whtdov_range_col <- st_crop(whtdov_range, neCol)
Ayerbe_mod_spp_l2$`Leptotila verreauxi` <- whtdov_range_col

# Bind into single spatial object
Ayerbe_mod_spp <- do.call(rbind, Ayerbe_mod_spp_l2) %>% st_make_valid()

## Calculate the shortest distance (km) from each point count location to the species range
Dist_to_pcs <- Pc_locs_dc_sf %>%
  distinct(Id_muestreo, geom) %>%
  st_distance(Ayerbe_mod_spp) %>% 
  as_tibble() %>% 
  rename_with(~ Species_list2$Species_ayerbe) %>%
  Cap_snake %>%
  mutate(across(everything(), ~ (.x / 1000) %>% units::drop_units())) %>% 
  bind_cols(distinct(Pc_locs_dc_sf, Id_muestreo)) %>% 
  relocate(Id_muestreo, 1)

# Plot an example species to see distance (in km)
if(FALSE){
  # Add Pc_ids & their geometries back in
  Dist_pcs_sf <- Pc_locs_dc_sf %>% 
    distinct(Id_muestreo, geom) %>% 
    full_join(Dist_to_pcs)
  
  Dist_pcs_sf %>% select(Ramphocelus_carbo) %>%
    rename(Distance = Ramphocelus_carbo) %>%
    mutate(In_out = ifelse(Distance > 0, "Out", "In")) %>%
    ggplot() + 
    geom_sf(data = neCol, color = "green", fill = NA) +
    geom_sf(data = Ayerbe_mod_spp_l2$`Ramphocelus carbo`) +
    geom_sf(aes(color = In_out)) + # Distance
    scale_color_manual(values = c("In" = "black", "Out" = "red"))# +
  #geom_sf(data = Pc_locs_sf, color = "red") 
}

## Clipping
# If point is outside of range by more than the cutoff (in km) set to 0, otherwise set to 1. 
# TO DO - could consider using the actual observed distances in our dataset to inform the cutoff values for each species. At present, the code will eliminate some of our observations (anything further away than the cutoff distance) 
cutoff <- 15

# Format
Dist_to_pcs2 <- Dist_to_pcs %>%
  left_join(add_year_join) %>% # Add Ano_grp 
  arrange(across(all_of(Row_identifier))) %>%
  unite(Id_muestreo_ano, all_of(Row_identifier), sep = ".", remove = FALSE) %>%
  select(-all_of(Row_identifier))

# Clip using predefined cutoff distance
Biogeo_clip <- Dist_to_pcs2 %>%  
  mutate(across(where(is.numeric), ~ as.integer(.x < cutoff))) %>%
  column_to_rownames("Id_muestreo_ano") %>%
  t() %>% 
  as_tibble()

# Revert back to matrix 
Biogeo_clip_mat <- Biogeo_clip %>% as.matrix()

# spOccupancy frame --------------------------------------------------------
# Multi-species abundance list
msAbu_l <- list(y = Abund_array, abu.covs = Site_covs3, det.covs = Obs_covs_l2, coords = Coords, range.ind = Biogeo_clip_mat) 

# Convert abundance to occupancy
msOcc_l <- msAbu_l
msOcc_l$y <- ifelse(msAbu_l$y > 0, 1, 0)
names(msOcc_l)[2] <- "occ.covs"


# Save object -------------------------------------------------------------

saveRDS(msOcc_l, "Rdata/msOcc_l_season.rds")

# EXTRAS ------------------------------------------------------------------
stop()

# >Checks -----------------------------------------------------------------
dim(msOcc_l$y) # Matches expectations? 

# Ensure that row / column order is the same
check_l <- list(Site_covs3, Obs_covs_l[[1]], Coords_arr)
names(check_l) <- c("Site_covs3", "Obs_covs_l", "Coords_arr")

Order_id_ano <- map(check_l, \(df){
  df %>% 
    unite(Id_muestreo_ano, all_of(Row_identifier), sep = ".", remove = FALSE) %>% 
    pull(Id_muestreo_ano)
})

# All should be TRUE
table(Order_id_ano$Site_covs3 == Order_id_ano$Obs_covs_l)
table(Order_id_ano$Obs_covs_l == Order_id_ano$Coords_arr)
table(names(Abund_wide_l[[1]])[-c(1:2)] == Order_id_ano$Obs_covs_l)
table(names(Biogeo_clip) == Order_id_ano$Obs_covs_l)

# If any problems.. Diagnose
probs <- which(!Order_id_ano$Site_covs3 == Order_id_ano$Obs_covs_l)
Site_covs3 %>% slice(probs) 

# >Plot probability breeding ----------------------------------------------
## Plot in observed data
Prob_breed_tbl %>% 
  ggplot(aes(x = Julian_day, y = Prob_breed, color = Ecoregion)) +
  geom_point(alpha = .2) +
  geom_smooth()

## Plot smoothed 
# Average elevation and latitude by region
Avg_by_reg <- Breed_records %>% summarize(
  Elevation = mean(Elevation, na.rm = TRUE),
  latitude = mean(latitude, na.rm=TRUE),
  .by = region
)

# Generate smoothed data for plotting 
newdat_smooth <- expand.grid(
  doy = 1:365,
  region = levels(Breed_records$region),
  year = round(mean(Breed_records$year, na.rm = TRUE), 0) #,
  #guild = "Invertebrate" # Could select most common guild in our data
) %>% left_join(Avg_by_reg)

# The response here (prob_breed) represents probability an individual (or observation) shows signs of breeding on that day
newdat_smooth2 <- newdat_smooth %>%
  mutate(prob_breed = predict(GAM_region, newdat_smooth, type = "response"))

# Plot smoothed curves 
newdat_smooth2 %>% 
  select(region, doy, prob_breed) %>% 
  tibble() %>% 
  filter(!region %in% c("Unknown", "Amazonia", "Pacific")) %>% 
  ggplot(aes(x = doy, y = prob_breed, color = region)) +
  geom_line()
