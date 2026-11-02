## PhD birds in silvopastoral landscapes ##
# Marginal-effects plots (predicted occupancy or detection along a covariate, optionally by a categorical covariate such as Habitat) for a fitted spOccupancy multi-species model.
# The plotted curve is the mean across species of each species' predicted probability, computed within each posterior draw (so the band is posterior uncertainty), at an average farm unit (random effects set to 0, ignore.RE = TRUE) and with no spatial effect; every non-focal covariate is held at its mean (numeric) or reference level (categorical).
# Same approach as Doser's stPGOcc-effects-plot-example.R (local copy in Extra_scripts/): a grid over the observed range, other covariates at 0 on the scaled axis, a hand-built X.0 with the model's column names, and the model's class switched to its non-spatial counterpart so predict() needs no coordinates. Differences: multi-species (averaged across species per draw), ignore.RE = TRUE because this model has occupancy random effects, and the axis is back-transformed with the scaling saved by 01_Wrangling_ms.R.

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(spOccupancy)

# Configuration -----------------------------------------------------------
msOcc_l_path <- "Derived/Rdata/Model_input/msOcc_l_meta_99_bug_fixes.rds"
Occ_mod_path <- "Derived/Rdata/Model_output/Occ_mod_98spp_1ecor_int_prec_jit1mm_s42_2026-10-31.rds"

response <- "occupancy" # "occupancy" or "detection"
# Either focal variable can be NULL: no focal_cat gives one curve along focal_cont; no focal_cont gives a point + interval per level. Whether the two interact is decided by the loaded model's formula, not here.
focal_cont <- "Tot_prec"
focal_cat <- "Habitat"
# Draw each categorical level in its own panel (with that level's species curves) rather than overlaid in one panel
facet_by_cat <- FALSE
# Faint per-species curves (or points); they can only be matched to a community estimate when there is one per panel, so by default they are drawn only then
show_species <- is.null(focal_cat) || is.null(focal_cont) || facet_by_cat
# Predict focal_cont only across the values observed for each plotted level, and for the levels at which non-focal categorical covariates are held (e.g. no pasture points lie inside forest), rather than across the full range
restrict_to_observed <- TRUE

pred_n <- 200
ci_level <- 0.95
# Optional override of the default value of any non-focal covariate, on the model's (scaled) scale, e.g. list(Year1 = 0, Habitat = "Ssp")
fixed_vals <- list()

# Display names and colours (same palette as the NMDS figure)
level_labels <- c(Pastizales = "Pasture", Ssp = "Silvopasture", Bosque = "Forest", N = "Not in forest", Y = "In forest")
level_colors <- c(Pasture = "#eb6834", Silvopasture = "#2a78d6", Forest = "#1baf7a")
axis_labels <- c(
  Dist_forest = "Distance to forest edge (m; negative = inside forest)", Tot_prec = "Annual precipitation (mm)",
  Year1 = "Year", ssp = "Silvopasture within 300 m (%)", Canopy_cover = "Forest cover within 1 km (%)",
  Habitat = "Habitat", Pc_start = "Start time (hours after midnight)", Prob_breed = "Breeding probability"
)
# Applied after back-transforming: offsets are added (Year1 is stored as the two-digit first year of the year group), divisors convert units (Pc_start is in seconds)
axis_offsets <- c(Year1 = 2000)
axis_divisors <- c(Pc_start = 3600)

# Functions ---------------------------------------------------------------
# Covariates of one process as a list with one vector per covariate (detection covariates are site x replicate matrices, flattened column-wise so all covariates stay aligned)
get_covs <- function(data_list, response) {
  covs <- list(occupancy = data_list$occ.covs, detection = data_list$det.covs)[[response]]
  if (is.null(covs)) stop("response must be 'occupancy' or 'detection'.")
  map(covs, \(x) if (is.data.frame(x) || is.matrix(x)) unlist(as.data.frame(x), use.names = FALSE) else x)
}

# Default value for a non-focal covariate: mean for numeric, reference level for factors, otherwise the most common value
default_value <- function(x) {
  if (is.numeric(x)) return(mean(x, na.rm = TRUE))
  if (is.factor(x)) return(levels(x)[1])
  names(which.max(table(x)))
}

extract_vars <- function(term) {
  setdiff(str_extract_all(term, "[A-Za-z][A-Za-z0-9_.]*")[[1]], c("I", "scale", "log", "exp", "sqrt"))
}

# Map one side of a design-matrix column name (e.g. "HabitatSsp", "Dist_forest", "I(Pc_start^2)") to its formula term: transformed terms stay intact, a dummy collapses to its covariate (longest covariate name that prefixes it)
map_term_to_var <- function(term, covs) {
  if (str_detect(term, "\\(") && any(extract_vars(term) %in% names(covs))) return(term)
  if (term %in% names(covs)) return(term)
  candidates <- names(covs)[startsWith(term, names(covs))]
  if (length(candidates) == 0) stop("Could not map model column '", term, "' to a covariate.")
  candidates[which.max(nchar(candidates))]
}

# The formula the model was fit with, rebuilt from its design-matrix column names, so every fitted interaction is kept whatever the focal variables
derive_model_formula <- function(model_cols, covs) {
  terms <- unique(map_chr(setdiff(model_cols, "(Intercept)"), \(col) {
    paste(map_chr(strsplit(col, ":", fixed = TRUE)[[1]], map_term_to_var, covs = covs), collapse = ":")
  }))
  as.formula(paste("~", paste(terms, collapse = " + ")))
}

# Prediction grid: focal_cont over its observed range (when restrict_to_observed, only among rows at each plotted level of focal_cat and at the held level of every other categorical covariate), crossed with the observed levels of focal_cat; every other model variable at its default or fixed_vals value, with categorical variables keeping all their levels so model.matrix() builds their columns
build_pred_df <- function(covs, model_vars, focal_cont, focal_cat, pred_n, restrict_to_observed, fixed_vals) {
  held <- map(set_names(setdiff(model_vars, c(focal_cont, focal_cat))), \(v) fixed_vals[[v]] %||% default_value(covs[[v]]))
  held_cat <- keep(names(held), \(v) !is.numeric(covs[[v]]))
  # Rows observed at the held level of every non-focal categorical covariate
  in_held <- reduce(map(held_cat, \(v) as.character(covs[[v]]) == held[[v]]), `&`, .init = rep(TRUE, length(covs[[model_vars[1]]])))
  cat_levels <- NULL
  if (!is.null(focal_cat)) {
    cat_vec <- covs[[focal_cat]]
    cat_levels <- intersect(levels(factor(cat_vec)), unique(as.character(cat_vec))) # observed levels, in factor order
  }
  cont_grid <- function(x) seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = pred_n)
  pred_df <- if (is.null(focal_cont)) {
    tibble(!!focal_cat := cat_levels)
  } else if (is.null(focal_cat)) {
    tibble(!!focal_cont := cont_grid(if (restrict_to_observed) covs[[focal_cont]][in_held] else covs[[focal_cont]]))
  } else {
    map(cat_levels, \(k) {
      x <- if (restrict_to_observed) covs[[focal_cont]][covs[[focal_cat]] == k & in_held] else covs[[focal_cont]]
      tibble(!!focal_cont := cont_grid(x), !!focal_cat := k)
    }) %>% list_rbind()
  }
  if (!is.null(focal_cat)) pred_df[[focal_cat]] <- factor(pred_df[[focal_cat]], levels = levels(factor(covs[[focal_cat]])))
  for (v in names(held)) {
    pred_df[[v]] <- if (is.factor(covs[[v]])) factor(held[[v]], levels = levels(covs[[v]])) else held[[v]]
  }
  pred_df
}

# Design matrix with exactly the model's columns; stops if a level in the data has no column in the model (e.g. a Mosaic level left in a list whose model was fit without it), which would otherwise be predicted silently as the reference level
build_X0 <- function(pred_df, formula, model_cols) {
  X0 <- model.matrix(formula, data = pred_df)
  missing_cols <- setdiff(model_cols, colnames(X0))
  extra_cols <- setdiff(colnames(X0), model_cols)
  if (length(missing_cols) > 0) stop("Prediction matrix lacks model columns: ", paste(missing_cols, collapse = ", "))
  if (length(extra_cols) > 0) stop("Prediction matrix has columns the model lacks (data and model do not match): ", paste(extra_cols, collapse = ", "))
  X0[, model_cols, drop = FALSE]
}

# Per-species posterior means and the community curve (mean across species within each draw, then quantiles); predict.msPGOcc returns draws x species x predictions
summarize_predictions <- function(samples, pred_df, sp_names, ci_level) {
  probs <- c((1 - ci_level) / 2, 0.5, 1 - (1 - ci_level) / 2)
  spp_mean <- t(apply(samples, c(2, 3), mean))
  colnames(spp_mean) <- sp_names
  species_df <- as_tibble(spp_mean) %>%
    mutate(.pred = row_number()) %>%
    pivot_longer(-.pred, names_to = "species", values_to = "mean") %>%
    left_join(mutate(pred_df, .pred = row_number()), by = ".pred")
  comm_quants <- apply(apply(samples, c(1, 3), mean), 2, quantile, probs = probs)
  mean_df <- mutate(pred_df, .pred = row_number(), low = comm_quants[1, ], mean = comm_quants[2, ], high = comm_quants[3, ])
  list(species_df = species_df, mean_df = mean_df)
}

# Back-transform a standardized covariate to its original units with the scaling lookup saved next to the data list
unscale <- function(x, covariate, scaling, response, offsets, divisors) {
  s <- filter(scaling, Process == response, Covariate == covariate)
  if (nrow(s) != 1) return(x) # not a scaled covariate
  (x * s$SD + s$Mean) / (divisors[covariate] %|% 1) + (offsets[covariate] %|% 0)
}
# Fallback for a missing (NA) lookup value
`%|%` <- function(a, b) if (is.na(a)) b else unname(a)

# Display names for the levels of a categorical covariate, keeping the order of the lookup; levels without a display name keep their own
relabel_levels <- function(x, labels) {
  lab <- coalesce(unname(labels[as.character(x)]), as.character(x))
  factor(lab, levels = unique(c(intersect(labels, lab), lab)))
}

# Colours for the plotted levels: the configured palette when it covers every level, otherwise ggplot's default hues
level_palette <- function(levels, colors) {
  if (all(levels %in% names(colors))) colors[levels] else setNames(scales::hue_pal()(length(levels)), levels)
}

# Lines along the focal covariate: bold community curve + credible ribbon per level, optional faint species curves, optional one panel per level
plot_marginal_lines <- function(species_df, mean_df, x_var, cat_var, show_species, facet_by_cat, colors) {
  if (is.null(cat_var)) {
    cat_var <- ".level"
    mean_df[[cat_var]] <- species_df[[cat_var]] <- "All"
    colors <- c(All = "black")
  }
  p <- ggplot(mean_df, aes(x = .data[[x_var]]))
  if (show_species) {
    p <- p + geom_line(data = species_df, aes(y = mean, group = interaction(species, .data[[cat_var]]), color = .data[[cat_var]]), alpha = 0.12, linewidth = 0.3)
  }
  p <- p +
    geom_ribbon(aes(ymin = low, ymax = high, fill = .data[[cat_var]]), alpha = 0.25) +
    geom_line(aes(y = mean, color = .data[[cat_var]]), linewidth = 1.2) +
    scale_color_manual(values = colors, aesthetics = c("color", "fill")) +
    theme_bw()
  if (facet_by_cat) p <- p + facet_wrap(vars(.data[[cat_var]])) + theme(legend.position = "none")
  if (cat_var == ".level") p <- p + guides(color = "none", fill = "none")
  p
}

# Points per level for a categorical focal variable alone: species as jittered points, community mean + credible interval
plot_marginal_points <- function(species_df, mean_df, x_var, show_species, colors) {
  p <- ggplot(mean_df, aes(x = .data[[x_var]], y = mean, color = .data[[x_var]]))
  if (show_species) p <- p + geom_jitter(data = species_df, alpha = 0.3, size = 0.8, width = 0.15, height = 0)
  p +
    geom_pointrange(aes(ymin = low, ymax = high), linewidth = 1, size = 0.8, color = "black") +
    scale_color_manual(values = colors, guide = "none") +
    theme_bw()
}

# Load data ---------------------------------------------------------------
msOcc_l <- readRDS(msOcc_l_path)
Occ_mod <- readRDS(Occ_mod_path)
scaling <- read_csv(str_replace(msOcc_l_path, "\\.rds$", "_scaling.csv"), show_col_types = FALSE)

# Predict -----------------------------------------------------------------
covs <- get_covs(msOcc_l, response)
model_cols <- colnames(Occ_mod[[c(occupancy = "X", detection = "X.p")[[response]]]])
form <- derive_model_formula(model_cols, covs)
model_vars <- unique(unlist(map(attr(terms(form), "term.labels"), extract_vars)))

pred_df <- build_pred_df(covs, model_vars, focal_cont, focal_cat, pred_n, restrict_to_observed, fixed_vals)
X0 <- build_X0(pred_df, form, model_cols)

# Fixed effects only: the non-spatial class needs no coordinates, and ignore.RE = TRUE sets the random effects to 0
Occ_mod_ns <- Occ_mod
class(Occ_mod_ns) <- "msPGOcc"
out_pred <- predict(Occ_mod_ns, X.0 = X0, ignore.RE = TRUE, type = response)
pred_samples <- out_pred[[c(occupancy = "psi.0.samples", detection = "p.0.samples")[[response]]]]
plot_data <- summarize_predictions(pred_samples, pred_df, Occ_mod$sp.names, ci_level)

# Original units and display labels
plot_data <- map(plot_data, \(df) {
  if (!is.null(focal_cont)) df[[focal_cont]] <- unscale(df[[focal_cont]], focal_cont, scaling, response, axis_offsets, axis_divisors)
  if (!is.null(focal_cat)) df[[focal_cat]] <- relabel_levels(df[[focal_cat]], level_labels)
  df
})
plot_colors <- if (is.null(focal_cat)) NULL else level_palette(levels(plot_data$mean_df[[focal_cat]]), level_colors)

# Plot --------------------------------------------------------------------
x_var <- focal_cont %||% focal_cat
y_lab <- c(occupancy = "Occupancy probability (mean across species)", detection = "Detection probability (mean across species)")[[response]]
re_label <- c(occupancy = "at an average farm unit with no spatial effect", detection = "for an average collector team")[[response]]
caption <- sprintf("Bold: mean across %d species, %g%% CrI, %s; other covariates at their mean or reference level.%s", length(Occ_mod$sp.names), 100 * ci_level, re_label, if (show_species) " Faint: posterior mean of each species." else "")

marginal_plot <- if (is.null(focal_cont)) {
  plot_marginal_points(plot_data$species_df, plot_data$mean_df, x_var, show_species, plot_colors)
} else {
  plot_marginal_lines(plot_data$species_df, plot_data$mean_df, x_var, focal_cat, show_species, facet_by_cat, plot_colors)
}
marginal_plot <- marginal_plot +
  labs(x = axis_labels[x_var] %|% x_var, y = y_lab, color = NULL, fill = NULL, caption = str_wrap(caption, 110)) +
  theme(plot.caption = element_text(hjust = 0))

# Export ------------------------------------------------------------------
out_file <- paste0("Figures/marginal_", paste(tolower(c(focal_cat, focal_cont)), collapse = "_"), if (facet_by_cat) "_facet" else "", ".png")
ggsave(out_file, marginal_plot, width = 8, height = 5.5, bg = "white")
print(marginal_plot)
