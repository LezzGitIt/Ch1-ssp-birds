## PhD birds in silvopastoral landscapes ##
# Marginal-effects plots (predicted occupancy or detection along a covariate, optionally by a categorical covariate such as Habitat) for a fitted spOccupancy multi-species model.
# The plotted curve is the mean across species of each species' predicted probability, computed within each posterior draw (so the band is posterior uncertainty), with random effects (e.g. farm unit) and the spatial effect at 0; every non-focal covariate is held at its mean (numeric) or reference level (categorical).
# Same approach as Doser's stPGOcc-effects-plot-example.R (local copy in Extra_scripts/): a grid over the observed range, other covariates at 0 on the scaled axis, and a hand-built X.0 with the model's column names. Differences: multi-species (averaged across species per draw), the probabilities are computed directly from the species coefficients (the same fixed-effects prediction as predict() with ignore.RE = TRUE, but it also works for multi-season fits, whose X is a site x period x coefficient array), and the axis is back-transformed with the scaling saved by 01_Wrangling_ms.R.
# Works for single-season (stacked) and multi-season fits; for a multi-season fit, point the data list at the stacked list it was built from (same covariates, one row per surveyed location-period), which supplies the observed ranges and the scaling.

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(spOccupancy)
ggplot2::theme_set(cowplot::theme_cowplot())

# Configuration -----------------------------------------------------------
# Stacked list the multi-season fits were built from (spOccupancy_temporal.R)
msOcc_l_path <- "Derived/Rdata/Model_input/Stacked/msOcc_l_meta_278_hab_fix.rds"

## One figure per row: the fitted model, a short tag for the file name, and the focal variables. Either focal variable can be NA: no focal_cat gives one curve along focal_cont; no focal_cont gives a point + interval per level. Whether the two interact is decided by the loaded model's formula, not here. Every row must be fitted to the data list above
# 278-species primary model (multi-season spatial, no Farm_join; Meta_gradient row 78); the 98-species rows 86, 88, 89 used msOcc_l_meta_98_hab_fix_jit1cm.rds
plot_specs <- tribble(
  ~model_tag,          ~Occ_mod_path,                                    ~focal_cont,   ~focal_cat,
  "temporal_prec_278", "Derived/Rdata/Model_output/Temporal/278.rds", "Tot_prec",    "Habitat"
)

response <- "occupancy" # "occupancy" or "detection"
# Draw each categorical level in its own panel (with that level's species curves) rather than overlaid in one panel
facet_by_cat <- FALSE
# Faint per-species curves (or points); they can only be matched to a community estimate when there is one per panel, so by default (NULL) they are drawn only then
show_species <- NULL
# Predict focal_cont only across the values observed for each plotted level, and for the levels at which non-focal categorical covariates are held (e.g. no pasture points lie inside forest), rather than across the full range
restrict_to_observed <- TRUE
# In-figure caption giving the model tag and plotting conventions; FALSE for figures captioned elsewhere, e.g. in a report
show_caption <- TRUE

pred_n <- 200
ci_level <- 0.95
# Optional override of the default value of any non-focal covariate, on the model's (scaled) scale, e.g. list(Year1 = 0, Habitat = "Ssp")
fixed_vals <- list()

# Display names and colours (same palette as the NMDS figure)
level_labels <- c(Pastizales = "Pasture", Ssp = "Silvopasture", Bosque = "Forest", N = "Not in forest", Y = "In forest")
level_colors <- c(Pasture = "#eda100", Silvopasture = "#6250d6", Forest = "#008300")
axis_labels <- c(
  Dist_forest = "Distance to forest edge (m; negative = inside forest)", Tot_prec = "Annual precipitation (mm)",
  Year1 = "Year", ssp = "Silvopasture within 300 m (%)", Canopy_cover = "Forest cover within 1 km (%)",
  Habitat = "Habitat", Pc_start = "Start time (hours after midnight)", Prob_breed = "Breeding probability"
)
# Applied after back-transforming: offsets are added (Year1 is stored as the two-digit first year of the year group), divisors convert units (Pc_start is in seconds)
axis_offsets <- c(Year1 = 2000)
axis_divisors <- c(Pc_start = 3600)
# Axes labelled in whole units every 2 (Year1 would otherwise get ticks such as 2012.5)
integer_axes <- c("Year1")

# Functions ---------------------------------------------------------------
# Fallback for a missing (NA) lookup value
`%|%` <- function(a, b) if (is.na(a)) b else unname(a)

# Covariate names in a formula term or design-matrix column name (e.g. "I(Pc_start^2)" gives "Pc_start")
extract_vars <- function(term) {
  setdiff(str_extract_all(term, "[A-Za-z][A-Za-z0-9_.]*")[[1]], c("I", "scale", "log", "exp", "sqrt"))
}

# One figure: predictions from one fitted model along focal_cont (by focal_cat), axes in original units, saved under Figures/ and printed; covs and scaling come from the data list, settings from the configuration
make_marginal_plot <- function(Occ_mod_path, model_tag, focal_cont, focal_cat, covs, scaling, response, settings) {
  Occ_mod <- readRDS(Occ_mod_path)
  show_species <- settings$show_species %||% (is.null(focal_cat) || is.null(focal_cont) || settings$facet_by_cat)

  ## Model formula, rebuilt from the coefficient names so every fitted interaction is kept whatever the focal variables
  # Multi-season fits store X as a site x period x coefficient array
  X <- Occ_mod[[c(occupancy = "X", detection = "X.p")[[response]]]]
  model_cols <- if (length(dim(X)) == 3) dimnames(X)[[3]] else colnames(X)
  # Each side of a column name (split at ":") maps to its formula term: covariates and transformed terms stay intact, a dummy (e.g. "HabitatSsp") collapses to its covariate (the longest covariate name that prefixes it)
  model_terms <- map_chr(setdiff(model_cols, "(Intercept)"), \(col) {
    sides <- map_chr(strsplit(col, ":", fixed = TRUE)[[1]], \(side) {
      if (side %in% names(covs) || (str_detect(side, "\\(") && any(extract_vars(side) %in% names(covs)))) return(side)
      candidates <- names(covs)[startsWith(side, names(covs))]
      if (length(candidates) == 0) stop("Could not map model column '", side, "' to a covariate.")
      candidates[which.max(nchar(candidates))]
    })
    paste(sides, collapse = ":")
  }) %>% unique()
  form <- as.formula(paste("~", paste(model_terms, collapse = " + ")))
  model_vars <- unique(unlist(map(attr(terms(form), "term.labels"), extract_vars)))

  ## Prediction grid: focal_cont over its observed range (when restrict_to_observed, only among rows at each plotted level of focal_cat and at the held level of every other categorical covariate), crossed with the observed levels of focal_cat
  # Every other model variable is held at its fixed_vals value, else its mean (numeric), reference level (factor) or most common value
  held <- map(set_names(setdiff(model_vars, c(focal_cont, focal_cat))), \(v) {
    x <- covs[[v]]
    settings$fixed_vals[[v]] %||% if (is.numeric(x)) mean(x, na.rm = TRUE) else if (is.factor(x)) levels(x)[1] else names(which.max(table(x)))
  })
  held_cat <- keep(names(held), \(v) !is.numeric(covs[[v]]))
  # Rows observed at the held level of every non-focal categorical covariate
  in_held <- reduce(map(held_cat, \(v) as.character(covs[[v]]) == held[[v]]), `&`, .init = rep(TRUE, length(covs[[model_vars[1]]])))
  # Observed levels of focal_cat, in factor order
  cat_levels <- if (is.null(focal_cat)) NULL else intersect(levels(factor(covs[[focal_cat]])), unique(as.character(covs[[focal_cat]])))
  grid_over <- \(x) seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = settings$pred_n)
  pred_df <- if (is.null(focal_cont)) {
    tibble(!!focal_cat := cat_levels)
  } else if (is.null(focal_cat)) {
    tibble(!!focal_cont := grid_over(if (settings$restrict_to_observed) covs[[focal_cont]][in_held] else covs[[focal_cont]]))
  } else {
    map(cat_levels, \(k) {
      x <- if (settings$restrict_to_observed) covs[[focal_cont]][covs[[focal_cat]] == k & in_held] else covs[[focal_cont]]
      tibble(!!focal_cont := grid_over(x), !!focal_cat := k)
    }) %>% list_rbind()
  }
  # Categorical variables keep all their levels so model.matrix() builds their columns
  if (!is.null(focal_cat)) pred_df[[focal_cat]] <- factor(pred_df[[focal_cat]], levels = levels(factor(covs[[focal_cat]])))
  for (v in names(held)) pred_df[[v]] <- if (is.factor(covs[[v]])) factor(held[[v]], levels = levels(covs[[v]])) else held[[v]]

  ## Design matrix with exactly the model's columns
  # Stops if a level in the data has no column in the model (e.g. a Mosaic level left in a list whose model was fit without it), which would otherwise be predicted silently as the reference level
  X0 <- model.matrix(form, data = pred_df)
  missing_cols <- setdiff(model_cols, colnames(X0))
  extra_cols <- setdiff(colnames(X0), model_cols)
  if (length(missing_cols) > 0) stop("Prediction matrix lacks model columns: ", paste(missing_cols, collapse = ", "))
  if (length(extra_cols) > 0) stop("Prediction matrix has columns the model lacks (data and model do not match): ", paste(extra_cols, collapse = ", "))
  X0 <- X0[, model_cols, drop = FALSE]

  ## Predicted probabilities with random and spatial effects at 0, for every posterior draw and species (draws x species x predictions, the layout predict() returns); species coefficients are named "<coefficient>-<species>"
  coefs <- Occ_mod[[c(occupancy = "beta.samples", detection = "alpha.samples")[[response]]]]
  pred_samples <- vapply(Occ_mod$sp.names, \(sp) {
    plogis(as.matrix(coefs[, paste(colnames(X0), sp, sep = "-"), drop = FALSE]) %*% t(X0))
  }, matrix(0, nrow(coefs), nrow(X0))) %>%
    aperm(c(1, 3, 2))

  ## Per-species posterior means, and the community curve (mean across species within each draw, then quantiles)
  probs <- c((1 - settings$ci_level) / 2, 0.5, 1 - (1 - settings$ci_level) / 2)
  spp_mean <- t(apply(pred_samples, c(2, 3), mean))
  colnames(spp_mean) <- Occ_mod$sp.names
  species_df <- as_tibble(spp_mean) %>%
    mutate(.pred = row_number()) %>%
    pivot_longer(-.pred, names_to = "species", values_to = "mean") %>%
    left_join(mutate(pred_df, .pred = row_number()), by = ".pred")
  comm_quants <- apply(apply(pred_samples, c(1, 3), mean), 2, quantile, probs = probs)
  mean_df <- mutate(pred_df, .pred = row_number(), low = comm_quants[1, ], mean = comm_quants[2, ], high = comm_quants[3, ])

  ## Original units and display labels
  # Back-transform with the scaling saved next to the data list (left as is when focal_cont is not a scaled covariate); display names follow the order of the lookup, and levels without one keep their own
  focal_scaling <- if (is.null(focal_cont)) NULL else filter(scaling, Process == response, Covariate == focal_cont)
  to_display <- \(df) {
    if (!is.null(focal_cont) && nrow(focal_scaling) == 1) {
      df[[focal_cont]] <- (df[[focal_cont]] * focal_scaling$SD + focal_scaling$Mean) / (settings$axis_divisors[focal_cont] %|% 1) + (settings$axis_offsets[focal_cont] %|% 0)
    }
    if (!is.null(focal_cat)) {
      lab <- coalesce(unname(settings$level_labels[as.character(df[[focal_cat]])]), as.character(df[[focal_cat]]))
      df[[focal_cat]] <- factor(lab, levels = unique(c(intersect(settings$level_labels, lab), lab)))
    }
    df
  }
  species_df <- to_display(species_df)
  mean_df <- to_display(mean_df)
  # The configured palette when it covers every plotted level, otherwise ggplot's default hues
  plot_levels <- if (is.null(focal_cat)) NULL else levels(mean_df[[focal_cat]])
  plot_colors <- if (is.null(focal_cat)) NULL else if (all(plot_levels %in% names(settings$level_colors))) settings$level_colors[plot_levels] else setNames(scales::hue_pal()(length(plot_levels)), plot_levels)

  ## Plot
  x_var <- focal_cont %||% focal_cat
  if (is.null(focal_cont)) {
    # Categorical focal variable alone: species as grey jittered points, community mean + credible interval per level in the level colour
    marginal_plot <- ggplot(mean_df, aes(x = .data[[x_var]], y = mean, color = .data[[x_var]]))
    if (show_species) marginal_plot <- marginal_plot + geom_jitter(data = species_df, color = "grey40", alpha = 0.3, size = 0.8, width = 0.15, height = 0)
    marginal_plot <- marginal_plot +
      geom_pointrange(aes(ymin = low, ymax = high), linewidth = 1, size = 0.8) +
      scale_color_manual(values = plot_colors, guide = "none") +
      theme_bw()
  } else {
    # Lines along focal_cont: bold community curve + credible ribbon per level, optional faint species curves, optional one panel per level; with no categorical variable, a single black curve
    cat_var <- focal_cat %||% ".level"
    if (is.null(focal_cat)) {
      mean_df[[cat_var]] <- species_df[[cat_var]] <- "All"
      plot_colors <- c(All = "black")
    }
    marginal_plot <- ggplot(mean_df, aes(x = .data[[x_var]]))
    if (show_species) {
      marginal_plot <- marginal_plot + geom_line(data = species_df, aes(y = mean, group = interaction(species, .data[[cat_var]]), color = .data[[cat_var]]), alpha = 0.12, linewidth = 0.3)
    }
    marginal_plot <- marginal_plot +
      geom_ribbon(aes(ymin = low, ymax = high, fill = .data[[cat_var]]), alpha = 0.25) +
      geom_line(aes(y = mean, color = .data[[cat_var]]), linewidth = 1.2) +
      scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
      theme_bw()
    if (settings$facet_by_cat) marginal_plot <- marginal_plot + facet_wrap(vars(.data[[cat_var]])) + theme(legend.position = "none")
    if (is.null(focal_cat)) marginal_plot <- marginal_plot + guides(color = "none", fill = "none")
  }
  if (x_var %in% settings$integer_axes) marginal_plot <- marginal_plot + scale_x_continuous(breaks = \(lims) seq(ceiling(lims[1]), floor(lims[2]), by = 2))

  y_lab <- c(occupancy = "Mean occupancy across species", detection = "Mean detection across species")[[response]]
  has_re <- isTRUE(Occ_mod[[c(occupancy = "muRE", detection = "pRE")[[response]]]])
  re_label <- c(occupancy = "with no spatial effect", detection = "")[[response]]
  if (has_re) re_label <- paste(c(occupancy = "at an average farm unit", detection = "for an average collector team")[[response]], re_label)
  caption <- sprintf("%s. Bold: mean across %d species, %g%% CrI, %s; other covariates at their mean or reference level.%s", model_tag, length(Occ_mod$sp.names), 100 * settings$ci_level, str_trim(re_label), if (show_species) " Faint: posterior mean of each species." else "")
  marginal_plot <- marginal_plot +
    labs(x = if (is.null(focal_cont)) NULL else settings$axis_labels[x_var] %|% x_var, y = y_lab, color = NULL, fill = NULL, caption = if (settings$show_caption) str_wrap(caption, 110) else NULL) +
    theme(plot.caption = element_text(hjust = 0))

  out_file <- paste0("Figures/marginal_", paste(c(model_tag, tolower(c(focal_cat, focal_cont))), collapse = "_"), if (settings$facet_by_cat) "_facet" else "", ".png")
  ggsave(out_file, marginal_plot, width = 8, height = 5.5, bg = "white")
  print(marginal_plot)
  invisible(marginal_plot)
}

# Load data ---------------------------------------------------------------
msOcc_l <- readRDS(msOcc_l_path)
scaling <- read_csv(str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv"), show_col_types = FALSE)
# Covariates of the plotted process, one vector per covariate (detection covariates are site x replicate matrices, flattened column-wise so all covariates stay aligned)
covs <- list(occupancy = msOcc_l$occ.covs, detection = msOcc_l$det.covs)[[response]]
if (is.null(covs)) stop("response must be 'occupancy' or 'detection'.")
covs <- map(covs, \(x) if (is.data.frame(x) || is.matrix(x)) unlist(as.data.frame(x), use.names = FALSE) else x)

# Plot --------------------------------------------------------------------
plot_settings <- list(
  pred_n = pred_n, ci_level = ci_level, restrict_to_observed = restrict_to_observed, fixed_vals = fixed_vals,
  facet_by_cat = facet_by_cat, show_species = show_species, show_caption = show_caption, level_labels = level_labels, level_colors = level_colors,
  axis_labels = axis_labels, axis_offsets = axis_offsets, axis_divisors = axis_divisors, integer_axes = integer_axes
)
# NA in plot_specs means no focal variable of that kind
marginal_plots <- pmap(plot_specs, \(model_tag, Occ_mod_path, focal_cont, focal_cat) {
  make_marginal_plot(Occ_mod_path, model_tag, if (is.na(focal_cont)) NULL else focal_cont, if (is.na(focal_cat)) NULL else focal_cat, covs, scaling, response, plot_settings)
}) |> set_names(plot_specs$model_tag)
