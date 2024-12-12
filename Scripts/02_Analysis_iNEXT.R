##PhD birds in silvopastoral landscapes##
#Analysis iNEXT -- Runs through iNEXT (4 steps and original) analysis
##This script creates various species accumulation curves to genearte comparable estimates of biodiversity across farms with different sampling efforts

#See associated documents in SNAPP -> Taxonomic diversity, as well as the associated Chao et al (2020) paper 

##NOTE: Some point counts are sampled multiple times while others sampled once. For example, a GAICA farm & CIPAV farm each may have 4 point counts, but GAICA could have sampled each farm 4x while CIPAV 1x##

#Contents
# 1) Prep data frame -- Sum the number of individuals by each species, and pivot data frame so spp is row & each column is a 1) point count or 2) farm X data base X year
# 2) iNEXT 4 steps -- Follow steps in Chao et al (2020) paper 
# 3) Preliminary analysis -- Goal: Understand how database (‘Uniq_db’), ecoregion, and the number of habitats surveyed (‘Num.hab’) influence taxonomic diversity & evenness for q = 0 (non-asymptotic at Cmax), as well as q = 1 & 2 (asymptotic). These are important variables to consider when interpreting taxonomic diversity results & thinking of potential biases.
# 4) Original iNEXT -- Before I discovered the iNEXT.4steps package I used the regular iNEXT package. I haven't revisited this code but it may not be necessary with improved iNEXT.4steps code
# 5) Vegan::specaccum() -- In Vegan we can create species accumulation curves based on the number of point counts on each farm. It seems like this is not necessary with the new iNEXT4steps framework, but Arturo recommended this at some point 

# Load libraries & data ---------------------------------------------------
library(iNEXT)
library(iNEXT.4steps)
library(vegan)
library(tidyverse)
library(janitor)
library(cowplot)
library(conflicted)
library(ggpubr)
ggplot2::theme_set(theme_cowplot())
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

load("Rdata/the_basics_09.25.24.Rdata")
load("Rdata/Taxonomy_09.26.24.Rdata")

#vignette(topic = "Introduction", package = "iNEXT")

# iNEXT curves -----------------------------------------------------------
# >Prep data frame ---------------------------------------------------------
#Remove species that weren't identified to species level
BirdPCs <- BirdPCs %>% filter(Nombre_ayerbe %in% TaxDf2$Species_ayerbe) %>% 
  #There are 6 rows in UBC data where # individuals wasn't recorded. We know there was at least 1 individual
  mutate(Numero_individuos = ifelse(is.na(Numero_individuos), 1, Numero_individuos))

#Sum number of individuals by each species, and pivot data frame so spp is row & each column is a.. 
#Point count
birds_wide_pc <- BirdPCs %>% #pc = point count
  group_by(ID_punto_muestreo_FINAL, Nombre_ayerbe) %>% 
  summarize(Numero_individuos = sum(Numero_individuos)) %>%
  #mutate(surveyNum = 1:n()) %>% 
  pivot_wider(names_from = ID_punto_muestreo_FINAL,
              values_from = Numero_individuos, 
              values_fn = mean,
              values_fill = 0)

#Each column is a farm X data base X year (ID_db_ano)
birds_wide_farm <- BirdPCs %>%
  mutate(ID_db_ano = paste(Id_gcs, Uniq_db, Ano_grp, sep = ".")) %>%
  group_by(ID_db_ano, Nombre_ayerbe) %>% 
  summarize(Numero_individuos = sum(Numero_individuos)) %>%
  #mutate(surveyNum = 1:n()) %>% 
  pivot_wider(names_from = ID_db_ano,
              values_from = Numero_individuos, 
              values_fn = mean,
              values_fill = 0)

#CHECK:: Should be 1 row per species
nrow(birds_wide_farm)
unique(BirdPCs$Nombre_ayerbe)

#Turn tibble to a dataframe w/ rownames
bw_df_pc <- birds_wide_pc %>% #birds wide data frame
  column_to_rownames(var = "Nombre_ayerbe")
bw_df_farm <- birds_wide_farm %>% #birds wide data frame
  column_to_rownames(var = "Nombre_ayerbe")

#Generate number of habitats per farm ID X Uniq_db X Ano_grp
Num.hab.df <- BirdPCs %>%
  mutate(ID_db_ano = paste(Id_gcs, Uniq_db, Ano_grp, sep = ".")) %>%
  #Adding Uniq_db & Ecoregion don't change number of rows in resulting df 
  group_by(ID_db_ano, Uniq_db) %>% #ADD ECOREGION HERE
  summarize(Num.hab = as.factor(length(unique(Habitat_homologado_ut)))) %>% 
  arrange(desc(Num.hab))

tabyl(Num.hab.df, Uniq_db, Num.hab)

# >iNEXT4steps ------------------------------------------------------------
i4steps <- iNEXT4steps(data = bw_df_farm, nboot = 0, details = TRUE) #TD = Taxonomic diversity

#Name "Assemblage' column, pivot_longer so Order.q is a single column. This permits joining with i4steps$summary[[2]] in next step
values_to <- c("SC_obs", "No_Asy_TD", "Evenness") #No_Asy_TD = Non asymptotic estimate of taxonomic div
dfs_long <- map2(i4steps$summary[c(1,3,4)], values_to, \(df, vt){
  df %>% rename_with(.cols = 1, ~"Assemblage") %>% 
    pivot_longer(cols = c(2:4), values_to = vt, names_to = "Order.q") %>% 
    mutate(Order.q = str_remove(Order.q, "q = "))
})

#Join the summary dataframes together with dfs_long data frames & Num.hab.df 
summary_df <- i4steps$summary[[2]]  %>% 
  mutate(Order.q = as.character(case_when(
    qTD == "Species richness" ~ 0,
    qTD == "Shannon diversity" ~ 1,
    qTD == "Simpson diversity" ~ 2
  ))) %>%
  left_join(dfs_long[[1]], by = c("Assemblage", "Order.q")) %>% 
  left_join(dfs_long[[2]], by = c("Assemblage", "Order.q")) %>% 
  full_join(dfs_long[[3]], by = c("Assemblage", "Order.q")) %>% 
  mutate(Evenness = ifelse(Order.q == "0", 1, Evenness)) %>% 
  left_join(Num.hab.df, by = join_by("Assemblage" == "ID_db_ano")) %>% 
  mutate(Id_gcs = str_split_fixed(Assemblage, pattern = "\\.", n = 3)[,1],
         Uniq_db = str_split_fixed(Assemblage, pattern = "\\.", n = 3)[,2],
         Ano = str_split_fixed(Assemblage, pattern = "\\.", n = 3)[,3]) %>% 
  relocate(c(Uniq_db, Id_gcs, Ano), .before = Assemblage) %>% 
  relocate(Order.q, .after = qTD)

#Export 
#write.xlsx(summary_df, file = paste0("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/SCR My products/Taxonomic diversity/Tax_div_Rd2_", format(Sys.Date(), "%m.%d.%y"), ".xlsx"), 
           #sheetName = "Diversity_estimates", row.names = F, showNA = F)

#Prepare for plotting#
i4_flat <- purrr::list_flatten(i4steps$details)

#Add metadata for plotting
i4_meta <- map(i4_flat, \(df){
  df %>% left_join(distinct(summary_df[,c("Assemblage", "Id_gcs", "Uniq_db", "Ano")])) %>% 
    left_join(Num.hab.df[,c("ID_db_ano", "Num.hab")], 
              by = join_by("Assemblage" == "ID_db_ano")) %>% 
    relocate(Assemblage, Id_gcs, Uniq_db, Ano)
})

##Plotting & summary stats##

#Sample Completeness#
#Plot sample completeness profile
facet_labels <- c("1" = "1 Habitat",
                  "2" = "2 Habitats",
                  "3" = "3 Habitats",
                  "4" = "4 Habitats")

ggplot(data = i4_meta$`Sample completeness`, 
                  aes(x = Order.q, y = Estimate.SC, color = Uniq_db, group = Assemblage)) + 
  geom_line(alpha = 0.2) + 
  facet_wrap(~Num.hab, labeller = labeller(Num.hab = facet_labels)) +
  labs(x = "Order q", y = "Sample coverage", color = "Database") + 
  theme(legend.position = "top")
#ggsave("Figures/Sample_completeness_prof.png", bg = "white", width = 9)

#Size based rarefaction#
#Create dfs of observed points & rarefied / extrapolated lines
size_plots <- list()
df.pt <- i4_meta$iNEXT_size_based %>% filter(Method == "Observed")
df.line <- i4_meta$iNEXT_size_based %>% filter(Method != "Observed") %>% 
  mutate(Method = factor(Method, levels = c("Rarefaction", "Extrapolation")))

#Option 1: Plot each panel separately #DELETE#
size_plots <- map(c(0,1,2), \(q){
  obs_rich <- df.pt %>% filter(Order.q == q)
  lines <- df.line %>% filter(Order.q == q)
  
  i4_meta$iNEXT_size_based %>% filter(Order.q == q) %>% 
ggplot(aes(x = m, y = qTD, color = Uniq_db,  group = interaction(Assemblage, Method))) + 
  geom_point(data = obs_rich) + 
  geom_line(data = lines, aes(linetype = Method), lwd=1.5, alpha = 0.2) 
})
ggarrange(size_plots[[1]], size_plots[[2]], size_plots[[3]], common.legend = TRUE, nrow = 1)

#Option 2: Plot using facets, this seems better
#Note the df.pt & df.line aren't from the iNEXT_coverage_based dataframe, but doesn't seem like it makes a difference 
facet_labels_q <- c("0" = "q = 0",
                  "1" = "q = 1",
                  "2" = "q = 2")
plot_sp.div <- function(iNEXT_df, x_var, xlab, coverage = TRUE){
  ggplot(data = iNEXT_df, aes(x = {{ x_var }}, y = qTD, 
                              color = Uniq_db, group = interaction(Assemblage, Method))) + 
    geom_point(data = df.pt, alpha = 0.8) + 
    geom_line(data = df.line, aes(linetype = Method), lwd=1.5, alpha = 0.3) +
    facet_wrap(~Order.q, labeller = labeller(Order.q = facet_labels_q)) + 
    labs(x = xlab, y = "Species diversity", color = "Database") + 
    guides(linetype = "none") + 
    theme(legend.position = "top") 
}

plot_sp.div(i4_meta$iNEXT_size_based, x_var = m, xlab = "Number of individuals", coverage = F)
#ggsave("Figures/Size_based_TD.png", bg = "white", width = 9)
plot_sp.div(i4_meta$iNEXT_coverage_based, x_var = SC, xlab = "Sample Coverage") 
#ggsave("Figures/Coverage_based_TD.png", bg = "white", width = 9)

#Asymptotic & empirical diversity profiles
i4_meta$`Observed and asymptotic diversity` %>% group_by(Uniq_db, Order.q, Method) %>% 
  summarize(mean.TD = mean(qTD)) %>% 
ggplot(aes(x = Order.q, y = mean.TD, color = Uniq_db, 
           group = interaction(Uniq_db, Method))) + 
  geom_line(aes(linetype = Method), lwd=1.5, alpha = 0.4) + 
  labs(x = "Order q", y = "Mean species diversity", color = "Database") + 
  guides(linetype = "none") + 
  theme(legend.position = "top")
#ggsave("Figures/Asy_emp_TD_prof.png", bg = "white", width = 9)

#Evenness profile
ggplot(data = i4_meta$Evenness_E3, 
       aes(x = Order.q, y = Evenness, color = Uniq_db, group = Assemblage)) +
  geom_line(lwd=1.5, alpha = 0.2) + 
  labs(x = "Order q", y = "Evenness", color = "Database") + 
  theme(legend.position = "top")
#ggsave("Figures/Evenness_prof.png", bg = "white", width = 9)

# >Preliminary analysis ---------------------------------------------------
#Goal: Understand how database (‘Uniq_db’), ecoregion, and the number of habitats surveyed (‘Num.hab’) influence taxonomic diversity & evenness for q = 0 (non-asymptotic at Cmax), as well as q = 1 & 2 (asymptotic). These are important variables to consider when interpreting taxonomic diversity results & thinking of potential biases.

df.a <- summary_df %>% filter(Id_gcs != 3133) %>% #Remove UBC farm with 2 forest points
  left_join(distinct(BirdPCs, Id_gcs, Ecoregion))
summary(lm(No_Asy_TD ~ Num.hab + Ecoregion, data = filter(df.a, Order.q == 0)))
summary(lm(TD_asy ~ Num.hab + Uniq_db + Ecoregion, data = filter(df.a, Order.q == 1)))
summary(lm(TD_asy ~ Num.hab + Uniq_db + Ecoregion, data = filter(df.a, Order.q == 2)))
summary(lm(Evenness ~ Num.hab + Uniq_db + Ecoregion, data = filter(df.a, Order.q == 1)))
summary(lm(Evenness ~ Num.hab + Uniq_db + Ecoregion, data = filter(df.a, Order.q == 2)))

# Original iNEXT --------------------------------------------------------
##Species diversity by farm##
start <- Sys.time() #Note this function takes time, can load time heavy objects above
#Main iNEXT function to generate Rarefaction & extrapolation sampling curves & sampling completeness curve
#In a loop could set endpoint = 3x sample size for each farm, and knots = 60
sp_div_farms <- iNEXT(bw_df_farm, q = c(0,1,2), datatype="abundance") 
Sys.time() - start

#Example farm. Note both report SC & m (abundance), & that they report the same estimates of species richness (b/c SC & m are linked)
sp_div_farms$iNextEst$size_based %>% filter(Assemblage == "CIPAV MBD_1262" & Order.q == 0)
sp_div_farms$iNextEst$coverage_based %>% filter(Assemblage == "CIPAV MBD_1262" & Order.q == 0)

#Fortify iNEXT object to create data frame for custom plotting
type <- c(1,2,3)
xlab <- c("Number of individuals", "Number of individuals", "Sample coverage")
ylab <- c("Species diversity", "Sample coverage", "Species diversity")
p_spp_div <- df.point <- df.line <- fort_sp_rich <- list()
for(i in 1:length(type)){
  print(i)
  fort_sp_rich[[i]] <- fortify(sp_div_farms, type = type[i])
  fort_sp_rich[[i]]$Uniq_db <- str_split_fixed(fort_sp_rich[[i]]$Assemblage, 
                                               pattern = "_", n = 2)[,1]
  fort_sp_rich[[i]]$Farm_Id_gcs <- str_split_fixed(fort_sp_rich[[i]]$Assemblage, 
                                              pattern = "_", n = 2)[,2]
  #Create point (observed data) & line (rarefied & extrapolated) data frames
  df.point[[i]] <- fort_sp_rich[[i]] %>% filter(Method == "Observed") #%>% filter(Farm_Id_gcs %in% c(2492, 1262))
  df.line[[i]] <- fort_sp_rich[[i]] %>% filter(Method != "Observed") #%>% filter(Farm_Id_gcs %in% c(2492, 1262))
  df.line[[i]]$Method <- factor(df.line[[i]]$Method,
                           c("Rarefaction", "Extrapolation"),
                           c("Rarefaction", "Extrapolation"))
  #Generate plots with facet_wrap for Hill numbers
p_spp_div[[i]] <- fort_sp_rich[[i]] %>% #filter(Farm_Id_gcs %in% c(2492, 1262)) %>% 
    ggplot(aes(x=x, y=y, group = interaction(Assemblage, Method))) +
    geom_point(data=df.point[[i]], aes(color = Uniq_db), size=5, alpha = 0.2) +
    geom_line(data=df.line[[i]], aes(color = Uniq_db, linetype = Method), lwd=1.5, alpha = 0.2) +
    {if(type[i] != 2)facet_wrap(~Order.q)} +
    #geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
    #                fill=Assemblage, colour=NULL), alpha=0.2) +
    labs(x = xlab[i], y= ylab[i]) + 
    scale_linetype_discrete(guide = "none") #+ 
  #ggrepel::geom_text_repel(data = df.point[[i]], aes(label = Farm_Id_gcs))
}

#Arrange and save plots
dir_figs_inext <- "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Figures/Proposal/iNEXT"

ggarrange(p_spp_div[[1]], p_spp_div[[3]], nrow = 2, common.legend = T, labels = "AUTO")
ggsave(paste0(dir_figs_inext, "/Sp_diversity_NI_SC.png"), width = 8, bg = "white") #Number individuals, sample coverage
ggarrange(p_spp_div[[2]] + theme(legend.position = "none"), labels = "C")
ggsave(paste0(dir_figs_inext, "/Sp_diversity_sampling_completeness.png"), bg = "white") #Number individuals, sample coverage

##Prepare Excel with relevant information plus accompanying figures##
#Add Uniq_db, FarmID & Order.q to the data frame that will be exported
AsyEst <- sp_div_farms$AsyEst %>% 
  mutate(Uniq_db = str_split_fixed(Assemblage, pattern = "_", n = 2)[,1],
         Farm_Id_gcs = str_split_fixed(Assemblage, pattern = "_", n = 2)[,2],
         Order.q = case_when(
           Diversity == "Species richness" ~ 0,
           Diversity == "Shannon diversity" ~ 1,
           Diversity == "Simpson diversity" ~ 2
         )) %>% 
  relocate(c(Uniq_db, Farm_Id_gcs), .before = Assemblage) %>% 
  relocate(Order.q, .after = Diversity) %>% 
  rename(AsyEst = Estimator)
head(AsyEst)

#Add SC & Number individuals (i.e., total abundance, column 'm') at the "Observed" point. I did check that 'coverage' and 'size' based dfs provide the same 'm' and 'SC' columns
SC <- sp_div_farms$iNextEst$coverage_based %>% filter(Method == "Observed")

AsyEst2 <- sp_div_farms$iNextEst$coverage_based %>% 
  group_by(Assemblage) %>% 
  slice_max(m) %>% 
  rename(Est2x = qD) %>% 
  select(Assemblage, Order.q, Est2x) %>% 
  left_join(AsyEst, by = c("Assemblage", "Order.q")) %>% 
  mutate(Dif_Asy_2x = AsyEst - Est2x) %>%
         #q0_Valid = ifelse(Dif_Asy_2x < 2, T, F)) %>% 
  relocate(Est2x, .before = AsyEst) %>% 
  relocate(Dif_Asy_2x, .after = AsyEst) %>% 
  ungroup()

#CHECK:: column 'm' is equal to the total abundance observed at a farm. 
bw_df_farm %>% summarize(across(
  .cols = everything(),
  .fn = ~sum(.x)))

AsyEst3 <- AsyEst2 %>% 
  left_join(SC[,c("Assemblage", "Order.q", "SC", "m")], by = c("Assemblage", "Order.q")) %>% 
  relocate(SC, m, .after = Farm_Id_gcs) %>% 
  rename(SC_obs = SC, Abundance_obs = m) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2)))

#Note low SC farm for UBC is Finca La Brisa
ggplot(AsyEst3, aes(x = SC_obs, y = after_stat(scaled))) + 
  geom_density(aes(color = Uniq_db)) + 
  labs(x = "Sample coverage", y = "Scaled frequency density")
ggsave(paste0(dir_figs_inext, "/SC_density.png"), bg = "white")

AsyEst3 %>% filter(Order.q == 0) %>% 
  ggplot(aes(x = Dif_Asy_2x)) + 
  geom_density() + 
  #geom_vline(xintercept = 2, linetype = "dashed") +
  labs(x = "Difference in sp. richness estimate between\n Asymptote and Estimate at 2x sample size")
ggsave(paste0(dir_figs_inext, "/Dif_Asy_2x.png"), bg = "white")

AsyEst3 %>% #pull(Farm_Id_gcs) %>% unique()
  arrange(desc(Abundance_obs))

#Plot Sample coverage & Dif in asymptote & the estimate at 2x the observed sample size
#This is the key thing I do not understand... I would expect that low sample coverage would lead to larger differences in the estimates at 2x the observed sample size
ggplot(AsyEst3, aes(x = SC_obs, y = Dif_Asy_2x)) + 
  geom_point(aes(color = Uniq_db)) + 
  geom_smooth() + 
  #xlim(.85, 1) +
  labs(x = "Sample coverage", 
       y = "Difference in sp. richness estimate between\n Asymptote and Estimate at 2x sample size")
ggsave(paste0(dir_figs_inext, "/SC_Dif_Asy_2x.png"), bg = "white")

#Create Sp_diversity_df with the Asymptotic estimates of species diversity (can change to 'Observed' or Est2x as well)
Sp_diversity_df <- AsyEst3 %>% select(Assemblage, Diversity, AsyEst) %>% 
  pivot_wider(names_from = Diversity, values_from = AsyEst)

Sp_diversity_df %>% ggplot(aes(x = `Species richness`, y = `Shannon diversity`)) + #Shannon
  geom_point() +
  geom_smooth() + 
  labs(title = "Asymptotic estimates of species diversity")
ggsave(paste0(dir_figs_inext, "/q0_q1_corr.png"), bg = "white") 

#Correlations between diversity metrics
Sp_diversity_df %>% select(-Assemblage) %>% 
  cor() %>% 
  data.frame() %>% 
  mutate(across(everything(), round, 2))

#Export AsyEst data frame
AsyEst3 %>% select(-Assemblage) %>% 
  relocate(Order.q, .after = Diversity) %>%
  rename(ObservedEst = Observed) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
  data.frame() #%>%
  #write.xlsx(file = "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/SCR My products/Tax_Div_Farms.xlsx", 
             #sheetName = "Tax_Div_Farm_Est", row.names = F)

#save(sp_div_farms, sp_rich, file = "/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/R_Files/PhD/Tax_div.Rdata")

#EXTRA
#Species richness by department plots
BirdPCs %>% filter(Nombre_institucion == "GAICA") %>% 
  group_by(Departamento, Ano) %>% 
  summarize(Richness = length(unique(Nombre_ayerbe))) %>% 
  ggplot(aes(x = Departamento, y = Richness)) +
  geom_bar()

BirdPCs %>%
  #filter(Nombre_institucion == "GAICA") %>%
  group_by(Departamento, Uniq_db) %>%
  summarize(Richness = length(unique(Nombre_ayerbe))) %>%
  ggplot(aes(x= fct_reorder(Departamento, Richness, .fun = max), y = Richness, 
             fill = fct_reorder(Uniq_db, Richness, .fun = max))) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Bird species richness",
       x = "Departament",
       y = "Richness",
       fill = "Data collector") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Users/aaronskinner/Library/CloudStorage/OneDrive-UBC/Grad_School/PhD/Figures/Proposal/Spp_rich_department.png", bg = "white")

##iNEXT by point count - Note there are so many point counts it's difficult to visualize anything. 
#Additionally, it may be necessary to do this process for each sampling institution or Uniq_db. 
start_letter <- c("C", "G", "U")
bw_df_dc <- sp_rich <- list() #Birds wide data frame for each data collector
for(i in 1:length(start_letter)){ 
  print(i)
  #Filter species (rows) to just those that were observed by each data collector to create more manageable data sets
  bw_df_dc[[i]] <- bw_df %>% select(starts_with(start_letter[i])) %>%
    filter(
      if_any(
        .cols = everything(),
        .fns  = ~ .x > 1,
      )
    )
  sp_rich[[i]] <- iNEXT(bw_df_dc[[i]], q = c(0), datatype="abundance") 
}


#Put this in Cookbook
data("bird")
bird #Data can be formatted like this but instead of North & South site it would be each point count (so hundreds of columns)
data("spider")
spider

#Uses bootstrapping to generate confidence intervals, so will want to set.seed() for reproducibility
?iNEXT
#Consider function div_profile in package hilldiv to create diversity profile curves (see Ausprey fig 1b)
iN_bird <- iNEXT(bird, q = c(0), datatype="abundance") #Change q number here to plot Hill numbers individually, or can set q = c(o,1,2) to plot them all using facet.var and col.var arugments

#Plot the iNEXT object as ggplot option, see vignette for examples of how we might want to customize 
ggiNEXT(iN_bird, type = 1, facet.var = "Assemblage") #Sample size
ggiNEXT(iN_bird, type = 3, facet.var="Order.q", color.var="Assemblage") 
ggiNEXT(iN_bird, type = 2) #Sample completeness curve
ggiNEXT(iN_bird, type = 3) #Sample coverage
?ggiNEXT

#Can manually set sample sizes
m <- c(1, 5, 20, 50, 100, 200, 400)
iNEXT(spider, q=0, datatype="abundance", size = m) #SC = sample coverage? 

#You can obtain estimates of species diversity at a specific sample size or coverage. Notice this is rarefaction for S site & extrapolation for N site
estimateD(bird, 
          datatype="abundance",
          base="size", #coverage
          level=225, #.985
          conf=0.95)


#Equation for q = 2, (i.e., the effective number of dominant species in the assemblage), why is the dominant species only worth 1.5 and non-dominant worth 400? 
(.8^2)^-1
(.05^2)^-1

# Vegan::specaccum curves by sampling effort ------------------------------
## In Vegan we can create species accumulation curves based on the number of point counts on each farm. It seems like this is not necessary with the new iNEXT4steps framework, but Arturo recommended this at some point 
#Create dataframe where Farm X PC X Datacollector are the rows#
#Add or remove Ano_grp depending on Christina's preference
srows <- BirdPCs %>% 
  group_by(Id_gcs, Uniq_db, ID_punto_muestreo_FINAL, Ano_grp, Nombre_ayerbe) %>% 
  summarize(Numero_individuos = sum(Numero_individuos)) %>%
  #mutate(surveyNum = 1:n()) %>% 
  pivot_wider(names_from = Nombre_ayerbe,
              values_from = Numero_individuos, 
              values_fn = mean,
              values_fill = 0)
#Create list 
srows_l <- srows %>% group_by(Id_gcs, Uniq_db, Ano_grp) %>%
  group_split()
#Name list
names(srows_l) <- map_chr(srows_l, \(df){
  unique(paste0(df$Id_gcs, ".", df$Uniq_db, ".", df$Ano_grp))
})

#Remove columns & create rownames
srows_l <- map(srows_l, \(df){
  df %>% select(-c(Id_gcs, Uniq_db, Ano_grp)) %>% 
    column_to_rownames(var = "ID_punto_muestreo_FINAL")
})

#Create species accumulation objects. Warning: the standard deviation is zero
sp_accum <- map(srows_l, \(df){
  specaccum(df)
})

#Create dataframes with key information for ggplotting 
sp_accum_df_l <- map(sp_accum, \(df){
  data.frame(Num_PCs = df$sites, Richness = df$richness, sd = df$sd)
})
names(sp_accum_df_l) <- names(srows_l)
#Combine dataframes 
sp_accum_df <- bind_rows(sp_accum_df_l, .id = "farm_dc_year") #dc = data collector

#Create function to decompose the ID x Uniq_db x Ano_grp column 
extract_metadata <- function(str, num){
  map_chr(str, \(chr){
    map_chr(strsplit(chr, "\\."), \(strgs){
      strgs[num]
    })
  })
}
IDs <- extract_metadata(sp_accum_df$farm_dc_year, 1)
Uniq_db <- extract_metadata(sp_accum_df$farm_dc_year, 2)
Ano_grp <- extract_metadata(sp_accum_df$farm_dc_year, 3)
sp_accum_df2 <- sp_accum_df %>% mutate(Id_gcs = IDs, Uniq_db = Uniq_db, Ano_grp = Ano_grp) 

#Plot 
#Create label data frame with a single row per farm X year X data collector
label_data <- sp_accum_df2 %>%
  group_by(Uniq_db, Id_gcs, Ano_grp) %>%
  slice_tail(n = 1)

ggplot(data = sp_accum_df2, aes(x = Num_PCs, y = Richness, group = interaction(Uniq_db, Id_gcs, Ano_grp), color = Uniq_db)) + 
  geom_line(alpha = .3) +
  labs(x = "Number of point counts", color = "Database") + 
  ggrepel::geom_text_repel(data = label_data, aes(label = Id_gcs)) + 
  theme(legend.position = "top")
#ggsave("Figures/Specaccum_curves_vegan.png", bg = "white", width = 9)

#No longer necessary to plot in base R# DELETE
plot(sp_accum[[1]], ci = 0, col = "red", ylim = c(0, max_rich))
imap(sp_accum[c(2:length(sp_accum))], \(accum_curve, name){
  col <- case_when(
    str_detect(name, "GAICA") ~ "red",
    str_detect(name, "CIPAV") ~ "blue",
    str_detect(name, "UNILLANOS") ~ "green",
    str_detect(name, "UBC") ~ "orange")
  plot(accum_curve, ci = 0, add = T, col = col, alpha = 0.1)
})
