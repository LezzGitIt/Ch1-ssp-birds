## PhD birds in silvopastoral landscapes ##
# NMDS of bird community composition by Habitat (Pastizales/Ssp/Bosque) in Meta, built from posterior draws of the fitted occupancy model's latent occupancy state (z.samples) rather than raw detections. A reference ordination is fit on the across-draw mean *dissimilarity matrix* (each draw's own Bray-Curtis distance matrix, averaged element-wise - not a distance matrix computed from the mean occupancy data, which isn't equivalent since Bray-Curtis is a nonlinear function of the input); each individual posterior draw is ordinated separately and Procrustes-rotated onto that reference, so the K per-habitat ellipses show how much the model's own estimation uncertainty shifts community composition, all on one comparable set of axes.

# User toggles ------------------------------------------------------------
K <- 50       # number of posterior draws to ordinate
nmds_k <- 3   # number of NMDS dimensions
seed <- 123   # for reproducible draw selection & metaMDS random starts
try_n <- 10   # metaMDS random starts (increase for a more stable final solution)
trymax_n <- 20

# TRUE = full species-name text cloud (every species at/above species_min_sites sites, the
# cluttered version); FALSE = drop that layer entirely and keep only the representative-species
# photos/labels, for a cleaner presentation figure. Also controls the export filename below.
show_species_labels <- TRUE

# TRUE = drop any individual posterior draw whose own metaMDS fit is degenerate (stress at/below
# stress_degenerate_threshold - an overfit, uninformative solution, not a genuinely good one) or
# exceeds stress_max_threshold (vegan's standard "not interpretable" ceiling), from every
# downstream calculation and plot layer (ellipses, species/site positions, stability diagnostics),
# not just the printed summary - FALSE keeps all K draws as fit, degenerate ones included.
exclude_bad_draws <- TRUE
# The observed degenerate draws in this model's K=50 fit sit at stress ~ 8e-05, orders of
# magnitude below the good cluster (~0.18-0.21) - 0.01 sits safely in the gap between them.
stress_degenerate_threshold <- 0.01
stress_max_threshold <- 0.2

if(nmds_k == 2){
  axes_names <- c("NMDS1", "NMDS2")
} 

if(nmds_k == 3){
  axes_names <- c("NMDS1", "NMDS2", "NMDS3")
}


# Validated 3-color categorical palette (dataviz skill, palette.md's first three slots -
# the only ones that clear all-pairs CVD/normal-vision floors together)
habitat_colors <- c(Pasture = "#eb6834", Silvopasture = "#2a78d6", Forest = "#1baf7a")

# Ellipse level - NULL = plain, unscaled SD ellipse (the standard "ellipses based on the SD of
# point scores" convention - a dispersion marker, not a coverage claim; points falling outside
# it are normal and expected). Set to e.g. 0.95 instead for a chi-square-scaled X% concentration
# ellipse - see conc_ellipse() below for why that's not the same as "X% SD" and why NULL is the
# better default for this particular plot.
ellipse_level <- NULL

# Drop species predicted present at fewer than this many sites (mean across the K draws) from
# the plotted labels, to declutter - set to 0 to keep everyone
species_min_sites <- 10

set.seed(seed)

# Libraries -----------------------------------------------------------------
library(tidyverse)
library(vegan)
library(ggrepel)
library(ggimage)
library(magick)

# Load data -----------------------------------------------------------------
msOcc_l <- readRDS("Derived/Rdata/Model_input/msOcc_l_meta_99_ref.rds")
#Occ_mod <- readRDS("Derived/Rdata/Model_output/Before_Oct_26/Occ_mod_99spp_1ecor_int_Dist_forest2026-05-06.rds")

# Build per-draw community matrices ------------------------------------------
# This data is already Meta-department-only (Department == "Meta" for every row); just exclude Mosaic
site_idx <- which(msOcc_l$occ.covs$Habitat != "Mosaic")
site_meta <- msOcc_l$occ.covs[site_idx, c("Id_survey", "Habitat")]
# Translate to English for the figure - recoded here (not just in the legend labels) so every
# downstream Habitat-grouped object (ellipses, species_candidates, etc.) already carries the
# English term, with nothing left keyed to the Spanish original.
site_meta$Habitat <- recode(site_meta$Habitat, Pastizales = "Pasture", Ssp = "Silvopasture", Bosque = "Forest")

# K draws evenly spaced across the posterior, rather than random, for even coverage
draw_idx <- round(seq(1, dim(Occ_mod$z.samples)[1], length.out = K))
z_draws <- Occ_mod$z.samples[draw_idx, , site_idx, drop = FALSE] # K x species x sites

# Every site has >= 1 species present in every individual draw, and every species is present
# at least once somewhere across the K draws (checked directly against this model/data) - so
# no rows or columns need to be dropped, and every draw's matrix has the same 460 x 99 shape,
# which is what makes the Procrustes alignment below valid (same points, same order, every time).
draw_comm <- function(i) {
  mat <- t(z_draws[i, , ])
  colnames(mat) <- Occ_mod$sp.names
  mat
}
comm_l <- map(seq_len(K), draw_comm)

# Reference dissimilarity matrix: average of the K draws' own Bray-Curtis distance matrices, not
# a distance matrix computed from the mean occupancy data (comm_mean) - averaging *after*
# computing dissimilarities avoids relying on distance-of-the-mean == mean-of-the-distances,
# which doesn't hold in general for a nonlinear dissimilarity like Bray-Curtis. Computed from all
# K draws as fit (not just the retained ones below) - a draw's fitting degeneracy is a property
# of how well metaMDS could resolve its distance matrix, not evidence that the distance matrix
# itself is invalid, so there's no principled reason to exclude it from the reference average.
dist_l <- map(comm_l, ~ vegdist(.x, method = "bray"))
dist_mean <- Reduce(`+`, dist_l) / length(dist_l)

# Fit one reference ordination + K per-draw ordinations ----------------------
# metaMDS's distance-matrix computation scales ~O(n^2), so K + 1 fits on ~460-row matrices are
# much cheaper in total than the single fit on a 4600-row stacked matrix this script used
# previously (~70 min - see Project_notes.md). Fitted independently, each draw's solution has
# an arbitrary axis rotation/reflection, so it's Procrustes-rotated onto the reference (below)
# before its ellipse is comparable to the others.
#
# fit_nmds() accepts either a raw community matrix (per-draw fits below, which still needs
# distance = "bray" to compute its own dissimilarity matrix internally) or an already-computed
# dissimilarity object (the reference fit, dist_mean, is passed as one) - metaMDS dispatches on
# the input's class either way and confirmed to give identical stress either path.
fit_nmds <- function(mat) {
  metaMDS(mat, k = nmds_k, distance = "bray", autotransform = FALSE, try = try_n, trymax = trymax_n, trace = 0)
}

nmds_cache <- paste0("Derived/Rdata/nmds/nmds_habitat_posterior_K", K, "_dim", nmds_k, "_distavg.rds")
if (file.exists(nmds_cache)) {
  cached <- readRDS(nmds_cache)
} else {
  cached <- list(nmds_ref = fit_nmds(dist_mean), nmds_draws = map(comm_l, fit_nmds))
  saveRDS(cached, nmds_cache)
}
nmds_ref <- cached$nmds_ref
nmds_draws <- cached$nmds_draws
nmds_ref # stress < 0.2 is interpretable, < 0.1 is ideal - bump nmds_k if too high

# Flag/exclude degenerate and high-stress individual draws ------------------
# A draw's own metaMDS fit can be degenerate (stress ~ 0) when that draw's community matrix
# happens to have many duplicate/near-identical site rows relative to nmds_k's dimensionality -
# vegan's own "stress is (nearly) zero: you may have insufficient data" warning at fit time is
# exactly this. A near-zero-stress solution isn't a good fit, it's an uninformative one (fit
# nearly perfectly because there was little real distance structure to resolve), so it shouldn't
# be averaged in alongside genuinely good draws. Flagged separately from the standard stress > 0.2
# "not interpretable" ceiling, since they're opposite failure modes, then both are excluded
# together if exclude_bad_draws is TRUE. The reference fit (nmds_ref, fit on dist_mean - the
# averaged dissimilarity matrix of all K original draws, cached) is left as-is either way - it's
# only the per-draw alignment target, and a single draw's fitting degeneracy doesn't make that
# draw's own distance matrix invalid.
draw_stress_all <- map_dbl(nmds_draws, "stress")
degenerate_draws <- which(draw_stress_all <= stress_degenerate_threshold)
high_stress_draws <- which(draw_stress_all > stress_max_threshold)
bad_draws <- union(degenerate_draws, high_stress_draws)

if (length(degenerate_draws) > 0) {
  warning(sprintf(
    "%d of %d posterior draws had a degenerate NMDS fit (stress <= %g): draw(s) %s.",
    length(degenerate_draws), K, stress_degenerate_threshold, paste(degenerate_draws, collapse = ", ")
  ), immediate. = TRUE, call. = FALSE)
}
if (length(high_stress_draws) > 0) {
  warning(sprintf(
    "%d of %d posterior draws exceeded the stress > %.2f 'not interpretable' threshold: draw(s) %s.",
    length(high_stress_draws), K, stress_max_threshold, paste(high_stress_draws, collapse = ", ")
  ), immediate. = TRUE, call. = FALSE)
}

if (exclude_bad_draws && length(bad_draws) > 0) {
  nmds_draws <- nmds_draws[-bad_draws]
  comm_l <- comm_l[-bad_draws]
  cat(sprintf(
    "Excluded %d of %d posterior draws (degenerate and/or stress > %.2f) from all downstream calculations and plots - %d draws remain.\n",
    length(bad_draws), K, stress_max_threshold, length(nmds_draws)
  ))
}

# Align each draw's site AND species scores onto the reference --------------
# An independently-fit NMDS solution's axes are only unique up to rotation/reflection, so each
# draw is Procrustes-rotated onto the reference before it's comparable to the others - the
# reference itself is only a stable, non-arbitrary alignment target now, not shown directly
# (fitting NMDS on the mean occupancy matrix does NOT give the average of the K draws' site
# positions, since NMDS is nonlinear - so "the mean" below is instead computed by averaging
# the K *aligned* draws' own positions). The same rotation fitted from a draw's site scores is
# applied to that draw's species scores via predict.procrustes(), since both live in the same
# pre-rotation coordinate space - this is what lets species positions be averaged the same way.
aligned <- map(nmds_draws, function(nmds_i) {
  pro <- procrustes(nmds_ref, nmds_i, symmetric = FALSE)
  list(sites = pro$Yrot, species = predict(pro, newdata = scores(nmds_i, display = "species")))
})

# Per-draw aligned site scores, for the K thin ellipses
draw_scores <- imap(aligned, function(a, i) {
  a$sites %>%
    as_tibble(.name_repair = ~ axes_names) %>%
    bind_cols(site_meta) %>%
    mutate(draw = i)
}) %>% list_rbind()

# Mean site/species positions, averaged across the K aligned draws (not read off the
# reference fit directly - see note above)
mean_site_scores <- map(aligned, "sites") %>%
  reduce(`+`) %>%
  (\(x) x / length(aligned))() %>%
  as_tibble(.name_repair = ~ axes_names ) %>%
  bind_cols(site_meta)

species_scores <- map(aligned, "species") %>%
  reduce(`+`) %>%
  (\(x) x / length(aligned))() %>%
  as_tibble(.name_repair = ~ axes_names) %>%
  mutate(Species = Occ_mod$sp.names)

# Mean number of sites each species is predicted present at, across the retained draws (fewer
# than K if exclude_bad_draws dropped any) - used to drop rarely-predicted species from the
# plotted labels below (species/site order is fixed and shared across comm_l, so this lines up
# positionally with species_scores$Species)
species_site_counts <- map(comm_l, colSums) %>% reduce(`+`) %>% (\(x) x / length(comm_l))()

# Representative species photos ----------------------------------------------
# Drop a PNG named "<Species_with_underscores>.png" into Figures/Species_images/ to illustrate it
# on the plot, then add its placement to `species_photo_offsets` below - an (x, y) nudge in NMDS
# units, added to that species' own true position, giving the offset position the photo (and its
# label) are drawn at. A thin line connects the photo back to the species' true point. Placement
# is manual on purpose: with only a handful of curated photos, picking where each one goes (e.g.
# alternating up/down within an ellipse, away from other species' labels) is a quick call to make
# by eye against the real plot, and not worth automating for this few. `species_true_coords`
# below prints each photographed species' real (unoffset) position every run, to read off when
# choosing/updating offsets. geom_image() preserves each photo's own aspect ratio automatically.
species_photo_dir <- "Figures/Species_images"
default_photo_size <- 0.15   # fraction of plot width - bump up/down to resize all photos at once
photo_label_gap <- 0.1        # NMDS units between a photo and its name label, drawn just below it

# photo_size defaults to default_photo_size for every row below; override individual rows (e.g.
# Tyrannus_melancholicus here) to size one photo differently from the rest.
species_photo_offsets <- tribble(
  ~Species,                 ~x_offset, ~y_offset, ~photo_size,
  "Sturnella_magna",        -0.15,     -0.6,      default_photo_size,
  "Vanellus_chilensis",     -0.2,      0.4,       default_photo_size,
  "Caracara_cheriway",       0,       0.8,     default_photo_size,
  "Tyrannus_melancholicus", -0.3,      0.55,      default_photo_size * 0.7,
  "Euphonia_laniirostris",   0.04,    -0.75,       default_photo_size,
  "Ramphocelus_carbo",       0.4,       0.65,    default_photo_size,
  "Thryophilus_rufalbus",    0.35,     -0.55,       default_photo_size,
  "Coragyps_atratus",       -0.26,     -0.9,      default_photo_size,
  "Cercomacroides_tyrannina", 0.4,      0.05,     default_photo_size,
)

species_photo_files <- list.files(species_photo_dir, pattern = "\\.png$", full.names = FALSE)
species_photo_names <- str_remove(species_photo_files, "\\.png$")

# Per-species spread across the K individual (pre-averaging) draws, plus how many sites each is
# predicted present at on average - a species whose position swings widely draw-to-draw is a
# signal its plotted position is more noise than compositional signal, most often because it's
# only rarely predicted present (few occupied sites means each draw's occupancy pattern for it
# can look quite different, and the mean position across draws inherits that instability).
species_score_spread <- imap(aligned, ~ as_tibble(.x$species, .name_repair = ~ axes_names) %>%
  mutate(Species = Occ_mod$sp.names, draw = .y)) %>%
  list_rbind() %>%
  summarise(NMDS1_sd = sd(NMDS1), NMDS2_sd = sd(NMDS2), .by = Species)

species_true_coords <- species_scores %>%
  filter(Species %in% intersect(species_photo_offsets$Species, species_photo_names)) %>%
  left_join(species_score_spread, by = "Species") %>%
  mutate(mean_sites_occupied = species_site_counts[match(Species, Occ_mod$sp.names)])
cat("Representative species - true NMDS position, draw-to-draw SD, and mean sites occupied (for setting/updating offsets, and for judging how stable a position is):\n")
print(species_true_coords, n = Inf)

species_photos <- species_photo_offsets %>%
  filter(Species %in% species_photo_names) %>%
  left_join(species_scores, by = "Species") %>%
  mutate(
    path = file.path(species_photo_dir, paste0(Species, ".png")),
    x_img = NMDS1 + x_offset,
    y_img = NMDS2 + y_offset,
    y_label = y_img - photo_label_gap,
    Species_label = str_replace_all(Species, "_", " ")
  )

# Species text labels - species with a photo get a fixed label anchored under their photo
# instead (below), so leave them out of the general repelled-label pool to avoid a second,
# competing label for the same species floating elsewhere on the plot.
species_scores_plot <- species_scores %>%
  filter(species_site_counts >= species_min_sites, !Species %in% species_photos$Species) %>%
  mutate(Species_label = str_replace_all(Species, "_", " "))

# Stress diagnostic ---------------------------------------------------------
# Per the vegan rule of thumb, stress < 0.2 is interpretable, < 0.1 is ideal - if the mean or
# several individual draws exceed 0.2, that's a signal nmds_k needs to go up (e.g. 2 -> 3).
# nmds_draws has already had any degenerate/high-stress draws removed above (if
# exclude_bad_draws), so this reports on the retained draws actually used below, not all K fit.
draw_stress <- map_dbl(nmds_draws, "stress")
cat(sprintf(
  "Stress across %d retained draws (of %d fit): mean = %.3f, range = [%.3f, %.3f], %d draws > 0.2\n",
  length(nmds_draws), K, mean(draw_stress), min(draw_stress), max(draw_stress), sum(draw_stress > 0.2)
))

# Ellipses: K per-draw (thin) + 1 averaged "mean" per habitat (bold) --------
# vegan::ordiellipse(kind = "sd") draws a *dispersion* ellipse of the site scores - NOT vegan's
# kind = "se", which is a much smaller ellipse for the uncertainty in the group *centroid's*
# location (shrinks as sample size grows) - a different quantity we don't want here.
#
# level = NULL (the default) reproduces the standard, widely-used convention (e.g. captioned as
# "ellipses based on the SD of point scores" in the literature): a plain, unscaled 1-SD ellipse
# - vegan's own default when `conf` isn't supplied to ordiellipse(). This is a purely descriptive
# dispersion marker, not a coverage claim, and it's normal/expected for many points to fall
# outside it.
# level = 0.95 (or any 0-1 value) instead scales by sqrt(qchisq(level, df = 2)) - the correct 2D
# concentration-ellipse radius (confirmed straight from vegan's source, `getAnywhere("ordiellipse")`
# - not literally "level SD"; 0.95 is ~2.45x SD, not 1.96, a 1D constant, not the 2D one). Tried
# on this project's actual data: because these ellipses are fit to the very same points they're
# plotted against (and posterior-draw-averaged NMDS scores are typically thinner-tailed than a
# true bivariate normal), a nominal 95% ellipse here enclosed ~97-99% of the points it was
# estimated from, not 95% - technically correct math, but not a very discriminating picture.
# Matching the field convention (unscaled SD) is the better default for this plot.
#
# ggordiplots::gg_ordiplot() (a thin wrapper around this same vegan machinery) was tried first
# but reproducibly errors even on this project's current, well-behaved data (traced into its
# own `with(df_ord[...], veganCovEllipse(...))` line - it expects a *live* ordination object to
# call vegan::scores()/vegan::ordiellipse() on directly, which we don't have here since these
# are Procrustes-aligned/posterior-averaged coordinates, not one fit's raw output). Calling
# vegan::ordiellipse() directly on the score matrix (draw = "none", so it just returns the
# per-group center/covariance/scale, doesn't plot) sidesteps that entirely and uses the same
# underlying vegan functions (ordiellipse + its internal veganCovEllipse polygon generator)
# gg_ordiplot itself would have used.
#
# ellipse_fit() stops at that per-group center/covariance/scale, without converting to a plotted
# polygon yet - so the bold "mean" ellipse below can be built by averaging the K per-draw fits
# themselves (center, covariance, and scale, each averaged element-wise across draws, per
# Habitat), not by re-fitting an ellipse to the K draws' already-averaged site positions
# (mean_site_scores). Averaging site positions before computing a dispersion ellipse mechanically
# shrinks the result - the same variance-reduction that makes a standard-error ellipse shrink
# with N - which is not what "the average of the 50 draws' own ellipses" means. ellipse_polygon()
# does the actual fit -> plotted-(x,y)-points conversion, shared by the thin per-draw ellipses
# and the bold averaged one.
ellipse_fit <- function(scores_mat, groups, level = NULL) {
  if (is.null(level)) {
    ordiellipse(scores_mat, groups = groups, kind = "sd", draw = "none")
  } else {
    ordiellipse(scores_mat, groups = groups, kind = "sd", conf = level, draw = "none")
  }
}

ellipse_polygon <- function(fit) {
  imap(fit, function(g, name) {
    vegan:::veganCovEllipse(g$cov, g$center, g$scale) %>%
      as_tibble(.name_repair = ~ c("x", "y")) %>%
      mutate(Group = name)
  }) %>% list_rbind()
}

draw_ellipse_fit <- ellipse_fit(
  as.matrix(draw_scores[, c("NMDS1", "NMDS2")]),
  groups = paste(draw_scores$Habitat, draw_scores$draw, sep = "_"),
  level = ellipse_level
)

draw_ellipses <- ellipse_polygon(draw_ellipse_fit) %>%
  mutate(Habitat = str_split_i(Group, "_", 1))

mean_ellipse_fit <- tibble(Group = names(draw_ellipse_fit)) %>%
  mutate(Habitat = str_split_i(Group, "_", 1)) %>%
  reframe(
    center = list(reduce(map(draw_ellipse_fit[Group], "center"), `+`) / n()),
    cov = list(reduce(map(draw_ellipse_fit[Group], "cov"), `+`) / n()),
    scale = mean(map_dbl(draw_ellipse_fit[Group], "scale")),
    .by = Habitat
  )

mean_ellipses <- pmap(mean_ellipse_fit, function(Habitat, center, cov, scale) {
  vegan:::veganCovEllipse(cov, center, scale) %>%
    as_tibble(.name_repair = ~ c("x", "y")) %>%
    mutate(Habitat = Habitat)
}) %>% list_rbind()

# Candidate species for photos -----------------------------------------------
# Ranks species by draw-to-draw position stability (lower combined SD = more consistent across
# the K draws - see species_score_spread above) within whichever Habitat's average ellipse
# center (mean_ellipse_fit$center, just computed) they sit closest to - a proxy for "near the
# middle of that habitat's ellipse" without needing to eyeball the plot. Restricted to species
# with at least candidate_min_sites occupied sites so rare/noisy species don't crowd out good
# candidates.
candidate_min_sites <- 20
habitat_centers <- mean_ellipse_fit %>%
  select(Habitat, center) %>%
  mutate(cx = map_dbl(center, 1), cy = map_dbl(center, 2))

species_candidates <- species_scores %>%
  left_join(species_score_spread, by = "Species") %>%
  mutate(mean_sites_occupied = species_site_counts[match(Species, Occ_mod$sp.names)]) %>%
  filter(mean_sites_occupied >= candidate_min_sites) %>%
  cross_join(habitat_centers) %>%
  mutate(dist_to_center = sqrt((NMDS1 - cx)^2 + (NMDS2 - cy)^2)) %>%
  slice_min(dist_to_center, by = Species, n = 1) %>%
  mutate(position_sd = sqrt(NMDS1_sd^2 + NMDS2_sd^2)) %>%
  select(Habitat, Species, mean_sites_occupied, position_sd, dist_to_center) %>%
  arrange(Habitat, position_sd)

cat(
  "\nCandidate species for photos - most stable position (lowest combined SD), nearest Habitat,",
  "min", candidate_min_sites, "sites occupied, top 5 per Habitat:\n"
)
species_candidates %>%
  slice_head(n = 5, by = Habitat) %>%
  print(n = Inf)

# Plot ------------------------------------------------------------------
# Layer order (each drawn on top of the last): faint site points, thin per-draw ellipses,
# bold mean ellipse, species labels on top of everything.
#
# Re-seeded here (same seed as the top of the script) because geom_text_repel's label layout
# draws its own random starting positions at render time - without re-seeding immediately before
# the plot is built, the labels' final arrangement (not the underlying species/site coordinates,
# which are already fixed by the cached NMDS fit) could still vary slightly run to run depending
# on how much of the seeded random stream earlier steps consumed (e.g. whether the NMDS cache was
# hit or refit).
set.seed(seed)

nmds_plot <- ggplot() +
  geom_point(
    data = mean_site_scores,
    aes(x = NMDS1, y = NMDS2, color = Habitat),
    alpha = 0.25, size = 1, show.legend = FALSE
  ) +
  geom_path(
    data = draw_ellipses,
    aes(x, y, group = Group, color = Habitat),
    alpha = 0.1, linewidth = 0.4
  ) +
  geom_path(
    data = mean_ellipses,
    aes(x, y, color = Habitat),
    linewidth = 1.4
  ) +
  # Anchored to the panel (data coords, x/y = -Inf/Inf), not the whole device like plot.tag was -
  # sits just above and right of the y-axis instead of the far top-left corner of the canvas.
  # Tracks nmds_draws (already filtered above if exclude_bad_draws is TRUE), so this always
  # reflects however many draws actually made it into the plot/calcs.
  annotate(
    "text", x = -Inf, y = Inf, hjust = -0.05, vjust = -0.8, size = 14 / .pt,
    label = sprintf("N = %d posterior draws", length(nmds_draws))
  ) +
  (if (show_species_labels) {
    geom_text_repel(
      data = species_scores_plot,
      aes(x = NMDS1, y = NMDS2, label = Species_label),
      size = 2.2, color = "grey20", max.overlaps = Inf, segment.size = 0.2
    )
  }) +
  geom_segment(
    data = species_photos,
    aes(x = NMDS1, y = NMDS2, xend = x_img, yend = y_img),
    color = "grey40", linewidth = 0.3
  ) +
  geom_image(
    data = species_photos,
    aes(x = x_img, y = y_img, image = path, size = I(photo_size))
  ) +
  geom_label(
    data = species_photos,
    aes(x = x_img, y = y_label, label = Species_label),
    size = 2.6, color = "grey20", fontface = "italic",
    fill = "white", label.size = 0, label.padding = unit(0.15, "lines")
  ) +
  scale_color_manual(values = habitat_colors) +
  # clip = "off" lets photos/labels spill into the plot margins rather than being cut off at the
  # panel edge - paired with the wider plot.margin below, which gives them literal canvas space
  # to render into
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    plot.margin = margin(t = 60, r = 60, b = 60, l = 60),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.caption = element_text(size = 14)
  ) +
  labs(
    x = "NMDS1", y = "NMDS2", color = "Habitat",
    # bottom-right by default (plot.caption's default hjust = 1) - a quick on-figure reference
    # for the stress value alongside the fuller diagnostic title above
    caption = sprintf("Average stress = %.3f", mean(draw_stress))
  )

# Not auto-printed here on purpose - printing a geom_image() plot to an interactive graphics
# device (e.g. RStudio's Plots pane) before that pane has a real, non-zero size throws "Error in
# grid.Call(C_convert, ...): Viewport has zero dimension(s)". ggsave() below is unaffected (it
# renders to its own device with explicit width/height), and the saved PNG is the actual output,
# so there's nothing lost by skipping the on-screen preview. Run `nmds_plot` yourself afterward,
# once your Plots pane is sized, if you want to see it inline.

# Export ------------------------------------------------------------------
# Larger canvas than the panel itself needs, so the margin-space photos (see coord_cartesian()
# above) have somewhere to actually land on the saved file instead of running off its edge.
# Filename tracks show_species_labels so the cluttered and clean versions don't overwrite each other.
# 16:9 aspect (widescreen PowerPoint slide) rather than the previous ~16:13 - wider canvas at the
# same height, so the figure fills more of the slide instead of leaving side margins once placed.
output_file <- if (show_species_labels) {
  "Figures/nmds_habitat_posterior.png"
} else {
  "Figures/nmds_habitat_posterior_no_text.png"
}
ggsave(output_file, nmds_plot, bg = "white", width = 20, height = 11.25)

# The 20x11.25 canvas above is sized generously so margin-space photos never get cut off, which
# leaves wasted white space wherever the actual content doesn't reach that far - trim it off the
# saved file directly rather than hand-tuning width/height to match the content's real extent.
magick::image_read(output_file) %>%
  magick::image_trim(fuzz = 2) %>%
  magick::image_write(output_file)
