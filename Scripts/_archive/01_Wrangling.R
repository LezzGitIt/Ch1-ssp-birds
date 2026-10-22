## PhD birds in silvopastoral landscapes ##
# Setting up occupancy & abundance data, the site (unit) covs, and observation (event) covariates

## This script generates the data structures that will be used as inputs for hierarchical modeling in future scripts 

# Script contents ---------------------------------------------------------
# 1) Subset: Filter to the 30 most abundant species
# 2) Biogeographic clipping: Remove point counts that are outside of the species distribution
# 3) Observation covariates: Create dataframes of covariates that vary by visit (time, date, etc.)
# 4) Site covariates: Create dataframe of covariates that are fixed for each unique row identifier (elevation, average temperature, landcover). 
# 4b) Merge the Lsm files, 
# 5) Abundance data: Create dataframe of observed abundances for each unique row identifier
# 6) Create unmarked::unmarkedFrame to store data in format understood by unmarked fitting function
# 7) spOccupancy: Use Obs, Site, & occupancy data to format for spOccupancy package. The main difference here is that spOccupancy also accepts coordinates of the sites 
# 8) Save & export important objects

# STILL TO DO: 
# Bring in the DW steps from the iNEXT script, look for overlap, create functions 

# Load libraries & data ---------------------------------------------------
## Load libraries
library(tidyverse)
library(hms)
library(janitor)
library(chron)
library(readxl)
library(sf)
library(unmarked)
library(conflicted)
library(cowplot)
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
Tax_df <- read_csv(paste0(Wrangling_repo, Excels, "Taxonomy/Taxonomy.csv"))
Bird_pcs_analysis <- read_csv(
  paste0(Wrangling_repo, Excels, "Bird_pcs/Bird_pcs_analysis.csv")
) %>% filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Event_covs <- read_csv(paste0(Wrangling_repo, Excels, "Event_covs.csv")) %>% 
  filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 
Site_covs <- read_csv(paste0(Wrangling_repo, Excels, "Site_covs.csv")) %>% 
  filter(!Id_muestreo_no_dc %in% paste0("MB-VC-EH_", Hatico_pc_nums)) 

load(paste0(Wrangling_repo, "Rdata/Lsm_l.Rdata"))
source("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Academia/Rcookbook/Themes_funs.R")
# vignette("unmarked") # Used vignette to help understand how data should be formatted

# Row identifier ----------------------------------------------------------
#NOTE:: The unique row identifier  [e.g. ID x Year] is critical , this defines what a row is in your dataframes and how many rows each dataframe will have 
Row_identifier <- c("Id_muestreo", "Ano_grp")

# Final formatting --------------------------------------------------------
Bird_pcs_analysis2 <- Bird_pcs_analysis %>%
  mutate(Species_ayerbe_ = str_replace_all(Species_ayerbe, " ", "_")) 

## Add Rep_ano_grp to tbl 
# Pc_start & Rep_ano_grp are equivalent in terms of grouping (when combined with Id_muestreo & Ano_grp), but since we are using Rep_ano_grp for formatting of analysis dataframes it makes sense to include Rep_ano_grp as well
date_join_spp_obs <- Event_covs %>% 
  filter(Spp_obs == 1) %>% # Only point counts with spp observed
  select(Id_muestreo, Ano_grp, Fecha, Pc_start, Rep_ano_grp)
Bird_pcs_analysis3 <- Bird_pcs_analysis2 %>% 
  left_join(date_join_spp_obs)

# Most abundant species ---------------------------------------------------
## Simplify this exercise by using only the most abundant species 
Spp_counts <- Bird_pcs_analysis3 %>%
  summarize(Count = sum(Count), .by = c(Species_ayerbe, Species_ayerbe_)) %>%
  arrange(desc(Count))

# Select the 30 most abundant species (that have maps)
No_maps <- c("Leptotila verreauxi", "Accipiter bicolor", "Sirystes sibilator", "Thripadectes virgaticeps")
Top_abu_df <- Spp_counts %>% 
  filter(!Species_ayerbe %in% No_maps) %>%
  slice_max(order_by = Count, n = 30) 
Top_abu <- Top_abu_df %>% 
  #arrange(Species_ayerbe) %>%
  pull(Species_ayerbe)

# Biogeographic_clipping --------------------------------------------------
## Idea from Socolar's (2022) paper: Biogeographic multi-species occupancy models for large-scale survey data. We only want to include point count locations that are within the species range (+ some buffer), differentiating a true zero (possible but not observed) vs points that are simply out of range (more like an NA).

# Load in relevant shapefiles 
Ayerbe_mod_spp_l <- map(Top_abu, \(spp){
  st_read(dsn = "../Geospatial_data/Ayerbe_shapefiles_1890spp", layer = spp) %>% st_make_valid()
})
names(Ayerbe_mod_spp_l) <- Top_abu
Ayerbe_mod_spp <- do.call(rbind, Ayerbe_mod_spp_l)

# Calculate the shortest distance (km) from each point count location to the species range
Dist_to_pcs <- Pc_locs_dc_sf %>%
  distinct(Id_muestreo, geom) %>%
  st_distance(Ayerbe_mod_spp) %>% 
  as_tibble() %>% 
  rename_with(~ Top_abu) %>%
  Cap_snake %>%
  mutate(across(everything(), ~ (.x / 1000) %>% units::drop_units())) %>% 
  bind_cols(distinct(Pc_locs_dc_sf, Id_muestreo)) %>% 
  relocate(Id_muestreo, 1)

# Add Pc_ids & their geometries back in
Dist_pcs_sf <- Pc_locs_dc_sf %>% 
  distinct(Id_muestreo, geom) %>% 
  full_join(Dist_to_pcs)

# Plot the rows with Pcs > 0 km
if(FALSE){
  Dist_pcs_sf %>% filter(if_any(where(is.numeric), ~ .x > 0)) %>%
    ggplot() + 
    geom_sf(data = Ayerbe_mod_spp_l$`Ixothraupis guttata`) +
    geom_sf(aes(color = Ixothraupis_guttata)) #+
  #geom_sf(data = Pc_locs_sf, color = "red") 
}

## Extract 'Id_muestreo' where the distance to the spp range == 0 for later subsetting
Pcs_in_range <- Dist_to_pcs %>%
  select(where(is.numeric)) %>% 
  map(~ Dist_to_pcs %>%
        filter(.x == 0) %>% 
        pull(Id_muestreo))

# Only retain species that have at least some observations within its range 
spp_rm <- Pcs_in_range %>%
  imap(~ if (length(.x) == 0) .y) %>% # Retain names of length-0 slots
  compact() 
if(length(spp_rm) > 0){
  spp_rm %>% list_simplify()
  Pcs_in_range <- Pcs_in_range %>% discard_at(spp_rm)
  
  spp_rm <- str_replace(spp_rm, "_", " ")
  Top_abu_df <- Top_abu_df %>% filter(!Species_ayerbe %in% spp_rm2)
  Top_abu <- Top_abu[!Top_abu %in% spp_rm2] 
}

# NOTE:: In the future will likely want to use some proportion of the total area via st_area()
st_area(Ayerbe_mod_spp) 

# Observation covariates ----------------------------------------------------
# Generate Obs_covs_df where each row is a unique ID [Id_muestreo, Ano_grp], each column is a site visit, and the 'Variable' column identifies which Observation covariate the row corresponds to  

# Lengthen then widen dataframe
Obs_covs_df <- Event_covs %>%
  mutate(across(everything(), as.character)) %>%
  # Select detection covariates here
  pivot_longer(
    cols = c(
      Nombre_institucion, Julian_day, Pc_start, Pc_length, Sampling_day
      ),      
    names_to = "Variable",          
    values_to = "Value"             
  ) %>%
  arrange(Ano_grp, Rep_ano_grp) %>%
  pivot_wider(
    id_cols = c(all_of(Row_identifier), Variable),
    names_from = c(Rep_ano_grp),
    names_glue = "Rep_ano_grp{Rep_ano_grp}",
    values_from = c(Value)
  ) 

# Split Obs_covs_df into a list, where each slot in the list is a covariate 
Obs_covs_spp <- map(Pcs_in_range, \(in_range){
  Obs_covs_grp <- Obs_covs_df %>%
    filter(Id_muestreo %in% in_range) %>% 
    arrange(across(all_of(Row_identifier))) %>%
    group_by(Variable) 
  Obs_covs <- Obs_covs_grp %>% group_split(.keep = FALSE)
  keys <- Obs_covs_grp %>% 
    group_keys() %>% 
    pull(Variable)
  names(Obs_covs) <- keys
  return(Obs_covs)
})

# Convert back to dates and times, respectively, and scale
# Define functions to convert type back to date & time
convert_type <- list(as.numeric, as.character, as_hms, as_hms, as.numeric)

Obs_covs <- map(Obs_covs_spp, \(Obs_covs) {
  # Apply transformations on the first 3 elements
  Obs_covs_hold <- map2(Obs_covs, convert_type, \(oc, type) {
    oc %>% mutate(across(-all_of(Row_identifier), ~ type(.))) %>% 
      select(-all_of(Row_identifier)) %>%                    
      mutate(across(where(~ is.numeric(.) | inherits(., "hms")), ~ scale(.)))             
  })
  return(Obs_covs_hold)
})

# Site covariates ---------------------------------------------------------
# Covariates that are fixed for a given point count x Ano_grp combination
date_join_ano <- Event_covs %>% 
  distinct(Id_muestreo_no_dc, across(all_of(Row_identifier)), Ano, Uniq_db, Canopy_cover, Canopy_height_m, ssp, te)
Site_covs_df <- Site_covs %>%
  full_join(date_join_ano) %>% 
  mutate(
    Habitat = if_else(Habitat == "Cultivos", NA, Habitat)
  ) %>% 
  # Some GAICA points were sampled 2x in a single Ano_grp across different years. Take the measurement from 2016.
  slice_max(by = all_of(Row_identifier), order_by = Ano) %>%
  mutate(
    # Take the earlier of the two years so it is numeric
    Ano1 = as.numeric(str_split_i(Ano_grp, "-", 1)),
    across(where(is.character), as.factor),
    Habitat = fct_relevel(Habitat, c("Pastizales", "Mosaic", "Ssp", "Bosque")), #"Bosque ripario"
    Ecoregion = relevel(Ecoregion, ref = "Piedemonte")
  ) %>% arrange(across(all_of(Row_identifier))) %>%
  select(Id_muestreo, Id_group_no_dc, Uniq_db, Ecoregion, Elev, Tot_prec, Habitat, Ano1, Uniq_db, Canopy_cover, Canopy_height_m, ssp, te)

# For each species, subset the relevant point counts
Site_covs_in_range <- map(Pcs_in_range, \(in_range){
  Site_covs_df %>% filter(Id_muestreo %in% in_range) 
}) 

# Extract Ids to match up with spOccupancy frame down the line
Ids_in_range_l <- map(Site_covs_in_range, \(sc_in_range){
  sc_in_range %>% pull(Id_muestreo)
})

# Remove Id_muestreo, scale numeric values
Site_covs <- map(Site_covs_in_range, \(sc_in_range){
  sc_in_range %>%
    select(-c(Id_muestreo)) %>% 
    mutate(across(where(is.numeric), ~ as.numeric(scale(.))))
}) 

# Abundance formatting  -----------------------------------------------------
# Include point counts with no spp observed 
add_pc_no_obs <- function(row_id, add_cols = NULL){ # additional columns
  Event_covs %>% 
    filter(Spp_obs == 0) %>% 
    distinct(across(all_of(c(row_id, add_cols)))) %>% 
    mutate(Species_ayerbe = NA, Count = NA)
}
No_obs_occ <- add_pc_no_obs(row_id = Row_identifier, add_cols = "Rep_ano_grp")

# Abundances -- Use rbind with No_obs_occ to add counts where no species were observed
Abund_no_obs <- Bird_pcs_analysis3 %>%
  distinct(across(all_of(Row_identifier)), Rep_ano_grp, Species_ayerbe, Count) %>% 
  rbind(No_obs_occ) 

# Create abundance list with just most abundant species 
Abund_l <- Abund_no_obs %>% 
  filter(Species_ayerbe %in% Top_abu) %>% 
  group_split(Species_ayerbe, .keep = FALSE)
names(Abund_l) <- names(Pcs_in_range)

# Join abundance data (just containing points where species was observed) with all the point counts that the species could have been observed with. 
date_join_rep <- Event_covs %>% 
  distinct(across(all_of(Row_identifier)), Rep_ano_grp)

# Add Rep_ano_grp & Ano_grp to the IDs in range to allow for a successful join with Abund_l
Abund_in_range <- map2(Pcs_in_range, Abund_l, \(in_range, abund){
  Date_in_range <- tibble(Id_muestreo = in_range) %>% 
    left_join(date_join_rep)
  #NOTE:: right_join effectively adds NAs for Count where species could have been observed but weren't (greatly lengthening the dataframe), and also does a clip, ensuring that there are not observations outside of buffer
  abund %>% right_join(Date_in_range) %>% 
    # Change NAs for 0s
    mutate(Count = if_else(is.na(Count), 0, Count))
})

# Pivot_wider and fill in with zeros 
# TO DO: Make a function here? Depending on iNEXT needs
Abund_zeros <- map(Abund_in_range, \(abund){
  abund %>% arrange(across(all_of(Row_identifier)), Rep_ano_grp) %>%
    pivot_wider(id_cols = all_of(Row_identifier),
                names_from = c(Rep_ano_grp),
                names_glue = "Rep_ano_grp{Rep_ano_grp}",
                values_from = Count, 
                values_fill = 0) %>% 
    arrange(across(all_of(Row_identifier))) %>% 
    select(-all_of(Row_identifier))
  #column_to_rownames(var = "Id_muestreo") 
})

# CHECK:: The number of point counts in range x Ano_grp 
nrow(Site_covs_df)
max(map_dbl(Abund_zeros, nrow))  # 664 should be max

# Use the Abund_zeros list to select only the relevant columns for the observation covariates
Obs_covs2 <- map2(Obs_covs, Abund_zeros, \(oc, abund){
  map(oc, \(oc){
    cols <- colnames(abund)
    oc %>% select(all_of(cols))
  })
})

# Any of the observation covs (start time, pc length, or date) dataframes have NAs where a given point count wasn't surveyed at a given repetition. Use this data frame to assign NAs to the abundance dataframe
Abund_nas <- map2(Abund_zeros, Obs_covs2, \(abund, oc){
  NAs <- oc$Pc_start
  abund %>% mutate(across(everything(), ~ if_else(is.na(NAs[[cur_column()]]), NA, .x)))
})

# Ensure the orders of Id_muestreo are the same in all input values to UMF
# Function to check if all values in a row are identical
is_same <- function(x) {
  all(x == x[1])
}

same_df <- tibble(
  az = Abund_nas$Bubulcus_ibis$Id_muestreo,
  sc = Site_covs$Bubulcus_ibis$Id_muestreo,
  oc = Obs_covs$Bubulcus_ibis$Fecha$Id_muestreo
) %>% rowwise() %>% 
  mutate(is_same = is_same(c_across(everything())))
table(same_df$is_same)

# unmarkedFrame -----------------------------------------------------------
## Create unmarkedFrames

# Multi-species
# NOTE:: Can do multi species occupancy in unmarked, but would only want to do this, for example, with UBC (& UniLlanos data).
#?unmarkedFrameOccuMulti # Set maxOrder = 1L to do zero higherorder interactions OR can manually set by putting covariates for the number of species present and 0s for everything else 
# Example with 3 species
stateformulas <- c("~Habitat", "~Habitat", "~Habitat", "0", "0", "0", "0")

# 'For temporary emigration fit Chandler (2011) model with gpcount(), Don't think this is necessary
# ?unmarked::gpcount()

# Abundance
umf_abu_l <- pmap(
  .l = list(Abund_nas, Site_covs, Obs_covs2), 
  \(abund, sc, oc){
    unmarkedFramePCount(y = abund, siteCovs = sc, obsCovs = oc)
  })

# Occupancy
Occ_nas <- map(Abund_nas, \(abund){
  abund %>% mutate(across(.cols = everything(), .fns = ~ ifelse(. > 0, 1, 0)))
})

umf_occ_l <- pmap(.l = list(Occ_nas, Site_covs, Obs_covs2), 
                  \(abund, sc, oc){
                    unmarkedFrameOccu(y = abund, siteCovs = sc, obsCovs = oc)
                  })

# spOccupancy frame --------------------------------------------------------
## Prep the data structure for spOccupancy package
# One benefit is that spOccupancy can account for spatial autocorrelation. 
# NOTE:: This is currently failing when fitting using spOccupancy because you cannot have multiple sites with the same coordinates. Quick solution would be to jitter coordinates a bit 

# Project coordinates for Colombia and filter by relevant PCs for each species
Coords <- map(Pcs_in_range,  \(in_range){
  Pc_locs_dc_sf %>% distinct(Id_muestreo, geom) %>% 
    filter(Id_muestreo %in% in_range) %>% 
    st_transform(crs = st_crs("EPSG:32618")) %>% 
    st_coordinates() # Remove as needed for plotting 
})

# Ensure adequate plotting 
if(FALSE){
  bi_sf <- Ayerbe_mod_spp %>% 
    filter(Nombre == "Bubulcus ibis") %>% 
    st_transform(crs = st_crs("EPSG:32618"))
  ggplot() + 
    geom_sf(data = bi_sf) +
    geom_sf(data = Coords[[1]])
}

## Combine inputs into the list structure needed by spOccupancy. 
# Adjust Row_identifier
if(any(str_detect(Row_identifier, "Ano_grp"))){
  Row_identifier <- str_replace_all(Row_identifier, "_grp", "1")
}

# Categorical random effects have to be specified as numeric in spOccupancy
Site_covs <- map(Site_covs, \(sc){
  sc %>% mutate(Id_group_no_dc_num = as.numeric(Id_group_no_dc))
})

# NOTE:: Different names required in the spOccupancy list, occ.covs = occurrence (or occupancy) covariates
spAbu_ll <- pmap(
  list(Abund_nas, Site_covs, Obs_covs2, Ids_in_range_l),
  \(abu, sc, oc, ids){
    # Create tbl with Ids & Ano1 to remove the habitats with NA
    rm_tbl <- sc %>% bind_cols(Id_muestreo = ids) %>% 
      filter(is.na(Habitat)) %>% 
      distinct(across(all_of(Row_identifier)))
    
    # Create list 
    spAbu_l <- list(y = abu, abu.covs = sc, det.covs = oc) 
    
    # Use rm_tbl to remove the necessary rows across all dataframes
    # modify_depth handles that det.covs is already a list 
    spAbu_l <- modify_depth(spAbu_l, 1, \(x){
      if (is.list(x) && all(map_lgl(x, is.data.frame))) {
        map(x, \(df){
          df %>% bind_cols(Id_muestreo = ids) %>% 
            anti_join(rm_tbl) %>% 
            select(-Id_muestreo) 
        })
      } else if (is.data.frame(x)) {
        x %>% bind_cols(Id_muestreo = ids) %>% 
          anti_join(rm_tbl) %>% 
          select(-Id_muestreo)
      } 
    })
    
    ## Add in coordinates
    add_year_join <- sc %>% bind_cols(Id_muestreo = ids) %>%
      distinct(Id_muestreo, Ano1)
    
    coords <- Pc_locs_dc_sf %>% 
      # Remove Ids out of range
      filter(Id_muestreo %in% unique(ids)) %>%
      distinct(Id_muestreo, geom) %>%
      # Add year to allow for anti_join
      left_join(add_year_join) %>%
      # Remove points with no LC
      anti_join(rm_tbl) %>% 
      # Order according to 'ids'
      arrange(factor(Id_muestreo, levels = tibble(ids))) %>%
      st_transform(crs = st_crs("EPSG:32618")) %>% 
      st_jitter(amount = .001) %>% 
      st_coordinates()
    # Add coords in 4th slot
    list_assign(spAbu_l, coords = coords)
  })

# Convert abundance to occupancy
spOcc_ll <- map(spAbu_ll, \(abu_l) {
  abu_l[[1]] <- abu_l[[1]] %>%
    mutate(across(.cols = everything(), .fns = ~ ifelse(. > 0, 1, 0)))
  names(abu_l)[2] <- "occ.covs"
  abu_l
})

# Top_abu_df --------------------------------------------------------------
Top_abu_df2 <- Top_abu_df %>% 
  mutate(Num_pcs = map_dbl(Pcs_in_range, length)) %>% 
  rename(Tot_count = Count)

# Save & Export -----------------------------------------------------------
stop()
# Export R data object for future analyses 
rm(list = ls()[!(ls() %in% c("umf_abu_l", "umf_occ_l", "spOcc_ll", "Top_abu_df2", "Abund_nas", "Abund_in_range", "Occ_nas"))]) #, "Site_covs", "Obs_covs2"
save.image(paste0("Rdata/Occ_abu_inputs_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))
#save.image("Rdata/Occ_abu_inputs_02.27.25.Rdata") # Manual

# Export Bird_pcs_analysis dataframe
if(FALSE){
  # Final data set I will analyze  w/ outliers reduced (max of 50 individuals) 
  Bird_pcs_analysis %>% 
    as.data.frame() %>%
    write.xlsx(file = "Derived/Excels/Bird_pcs_analysis.xlsx", 
               showNA = FALSE, row.names = FALSE)
  
  # Final dataset in csv format for data paper
  Birds_export <- Bird_pcs_analysis %>% 
    arrange(Id_muestreo, Id_muestreo_no_dc, Fecha, Pc_start) %>%
    select(-c(Species_ayerbe_))
  #select(Id_muestreo, Id_muestreo_no_dc, Fecha, Pc_start, Species_ayerbe, Count) 
  Birds_export %>% write_csv(file = "Derived/Excels/Bird_pcs_analysis.csv")
}

# Export R data object for BIOL314
#rm(list = ls()[!(ls() %in% c("Bird_pcs_all", "Bird_pcs_analysis", "Pc_hab", "Site_covs_df"))])
#save.image(paste0("Rdata/Biol314_inputs_", format(Sys.Date(), "%m.%d.%y"), ".Rdata"))

# Checks ----------------------------------------------------
## NOTE:: The UMF will throw an error if any of the dimensions don't match up, which is great. If the code runs this means that dimensions are correct 
# Error due to differing dimensions of SiteCovs
Bad_sites <- Site_covs$Amazona_ochrocephala %>% filter(Habitat != "Bosque")
unmarkedFramePCount(Abund_nas$Amazona_ochrocephala, Bad_sites, Obs_covs2$Amazona_ochrocephala)

# Error due to obsCov
Obs_covs3$Amazona_ochrocephala <- Obs_covs2$Amazona_ochrocephala
Obs_covs3$Amazona_ochrocephala[[3]] <- Obs_covs2$Amazona_ochrocephala[[3]] %>% select(-1)
unmarkedFramePCount(Abund_nas$Amazona_ochrocephala, Site_covs$Amazona_ochrocephala, Obs_covs3$Amazona_ochrocephala)


# There should be no NAs in the relevant columns at this point 
Event_covs %>% select(Uniq_db, Fecha, Pc_length, Pc_start) %>% 
  Na_rows_cols(distinct = TRUE)

# Visualize Pc_lengths by Uniq_db
Event_covs %>% distinct(Uniq_db, Pc_length) %>% 
  filter(!Uniq_db %in% c("Cipav mbd", "Gaica distancia"))

# Using the nested list Obs_covs 
map_depth(Obs_covs, 2, dim) # Troglodytes_aedon has all 551 rows possible, so we'll use it

# Add Id_muestreo back to the 3 dataframes
Example_l <- map(Obs_covs$Troglodytes_aedon, \(df){
  cbind(Id_muestreo = Site_covs_df$Id_muestreo, df) %>% 
    tibble()
})

# Could replace the nan for NA , but the coding is more complex than it should be 
Example_l$Pc_length[is.nan(Example_l$Pc_length)] <- NA

# This returns the proportion of values that are TRUE (NA in same loc) vs FALSE (NA in different loc)
Na_locs_equal(Example_l$Fecha, Example_l$Pc_start)
Na_locs_equal(Example_l$Fecha, Example_l$Pc_length)
Na_locs_equal(Example_l$Pc_start, Example_l$Pc_length)

# I think this must be causing the difference.. Likely due to the Pc_lengths that are 0. Fortunately I believe R interprets NA and NaN the same way in 99% of cases 
table(is.nan(Example_l$Pc_length[[2]]))
table(is.nan(Example_l$Pc_start[[2]]))

pcl <- Obs_covs_df %>% filter(Variable == "Pc_length")
pcs <- Obs_covs_df %>% filter(Variable == "Pc_start")
fecha <- Obs_covs_df %>% filter(Variable == "Fecha")

Na_locs_equal(pcl, pcs)
Na_locs_equal(fecha, pcs)
Na_locs_equal(fecha, pcl)

# DELETE 
TF_mat <- is.na(Example_l$Fecha) == is.na(Example_l$Pc_start) & is.na(Example_l$Fecha) == is.na(Example_l$Pc_length)
discrepancies <- which(!TF_mat, arr.ind = TRUE) %>% as_tibble()
discrepancies

ids_na_diff <- Example_l$Fecha[discrepancies$row, 1] %>% pull(Id_muestreo) %>% 
  unique()
Obs_covs_df %>% filter(Id_muestreo %in% ids_na_diff) %>% 
  select(Id_muestreo, unique(discrepancies$col)) %>%
  filter(if_any(everything(), is.na))

# Visualize example
df_meta %>% filter(Id_muestreo == "G-AD-M-LCA1_09") %>% 
  select(Id_muestreo, Fecha, Hora, Spp_obs)
