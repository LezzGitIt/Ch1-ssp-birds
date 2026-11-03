## Define functions for working with multi-species spOccupancy products
# Extract parameter estimates and estimates of gof (ESS, rhat) from model list

extract_parms <- function(model, spatial = FALSE, ms = FALSE){
  
  beta <- tidyMCMC(model$beta.samples) %>% mutate(parameter = "Beta")
  alpha <- tidyMCMC(model$alpha.samples) %>% mutate(parameter = "Alpha")
  
  theta <- NULL
  if(spatial && !is.null(model$theta.samples)){
    theta <- tidyMCMC(model$theta.samples) %>% mutate(parameter = "Theta")
  }
  
  pars <- bind_rows(beta, alpha, theta)
  
  diag_tbl <- tibble(
    ess = c(model$ESS$beta, model$ESS$alpha, if(spatial) model$ESS$theta else NULL),
    rhat = c(model$rhat$beta, model$rhat$alpha, if(spatial) model$rhat$theta else NULL)
  )
  parms_df <- bind_cols(pars, diag_tbl)
  
  if(ms){
    parms_df <- parms_df %>% mutate(
      Species_num = str_split_i(term, "-", 2),
      Species_num = str_remove(Species_num, "sp"),
      term = str_split_i(term, "-", 1), 
      Species_num = ifelse(term == "phi", paste0("factor", Species_num), Species_num)
    ) %>% relocate(Species_num, .before = term)
  }
  return(parms_df)
}

# Rhat and effective sample size
plot_diagnostics <- function(parms_df){
  parms_df %>%
    ggplot(aes(ess, rhat, color = term, shape = parameter)) +
    geom_point() +
    geom_hline(yintercept = 1.1) +
    geom_vline(xintercept = 100) +
    guides(color = "none")
}

plot_species_diagnostics <- function(parms_df){
  parms_df %>%
    pivot_longer(cols = c(ess, rhat), names_to = "diagnostic") %>%
    ggplot(aes(species_ayerbe, value, color = term, shape = parameter)) +
    geom_point(position = "jitter") +
    facet_wrap(~diagnostic, scales = "free_y") +
    guides(color = "none") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

## GOF tests
# Run posterior predictive checks in spOccupancy
run_ppc <- function(mod){
  list(
    ppcOcc(mod, "freeman-tukey", 1),
    ppcOcc(mod, "freeman-tukey", 2) 
    #ppcOcc(mod, "chi-squared", 1),
    #ppcOcc(mod, "chi-squared", 2)
  )
}

## Bayesian p-value
# Function taking posterior predictive check object and calculating bayesian p-value
calc_bayes_p <- function(ppc_obj){
  fit_stat <- ppc_obj$fit.stat
  group <- ifelse(ppc_obj$group == 1, "Site", "Replicate")
  
  tot.post <- ppc_obj$n.chains * ppc_obj$n.post
  bayes_p <- sum(ppc_obj$fit.y.rep > ppc_obj$fit.y) / tot.post
  bayes_p <- round(bayes_p, 3)
  tibble(fit_stat, group, bayes_p)
}

# For multi_species objects
calc_bayes_p_ms <- function(ppc_obj){
  tot.post <- ppc_obj$n.chains * ppc_obj$n.post
  # Multi-species multi-season fits return the fit statistic as draws x species x period; sum it over periods so each species gets one Bayesian p-value, as for stacked fits (unsurveyed periods contribute 0 to both)
  if (ppc_obj$class %in% c("tMsPGOcc", "stMsPGOcc", "svcTMsPGOcc")) {
    ppc_obj$fit.y <- apply(ppc_obj$fit.y, c(1, 2), sum)
    ppc_obj$fit.y.rep <- apply(ppc_obj$fit.y.rep, c(1, 2), sum)
  }
  # Generate TF vector for each species (column)
  bayes_p <- colSums(ppc_obj$fit.y.rep > ppc_obj$fit.y) / tot.post
  bayes_p <- round(bayes_p, 3)
  # Create tibble
  fit_stat <- ppc_obj$fit.stat
  Species <- ppc_obj$sp.names
  group <- ifelse(ppc_obj$group == 1, "Site", "Replicate")
  tibble(fit_stat, group, Species, bayes_p)
} 

extract_gof <- function(ppc_list){
  imap(ppc_list, ~ map(.x, calc_bayes_p) %>%
         bind_rows() %>%
         mutate(species_ayerbe = .y)) %>%
    bind_rows()
}

# Visualize
plot_gof <- function(gof_tbl){
  gof_tbl %>%
    ggplot(aes(species_ayerbe, bayes_p, shape = group)) +
    geom_boxplot(position = position_dodge(.7)) +
    geom_point(aes(color = fit_stat), position = position_dodge(.7)) +
    geom_hline(yintercept = c(0.5, 0.1), linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "top")
}

## Fit summary
# Condense a fitted multi-species model (object or .rds path) into a small list holding everything the model log (spOcc_model_notes.xlsx) records, plus the community-level draws, so a fit made on another machine can be logged and its community effects explored without moving the multi-GB fit; saved as an .rds when save_path is given (default for a path input: the fit's name with a _summary suffix)
summarize_fit <- function(fit, save_path = NULL, gof = FALSE){
  if (is.character(fit)) {
    if (is.null(save_path)) save_path <- sub("\\.rds$", "_summary.rds", fit)
    fit_file <- basename(fit)
    fit <- readRDS(fit)
  } else {
    fit_file <- if (is.null(save_path)) NA_character_ else sub("_summary\\.rds$", ".rds", basename(save_path))
  }
  # Posterior median and 95% CrI of every column of a sample matrix
  post_quantiles <- function(samples) if (is.null(samples)) NULL else apply(as.matrix(samples), 2, stats::quantile, c(.5, .025, .975))
  # Posterior median within each chain, to spot a chain stuck in a different mode
  chain <- rep(seq_len(fit$n.chains), each = fit$n.post)
  chain_medians <- function(samples) if (is.null(samples)) NULL else t(sapply(seq_len(fit$n.chains), function(k) apply(as.matrix(samples)[chain == k, , drop = FALSE], 2, stats::median)))
  # Convergence: the share of species-level and spatial parameters with Rhat > 1.1 or ESS < 100 is the log's Percent 'problematic'; loadings (lambda) and spatial effects (w) are excluded from the minimum ESS
  ess_all <- unlist(fit$ESS[setdiff(names(fit$ESS), c("lambda", "w"))])
  rhat_param <- c(fit$rhat$beta, fit$rhat$alpha, fit$rhat$theta)
  ess_param <- c(fit$ESS$beta, fit$ESS$alpha, fit$ESS$theta)
  rhat_comm <- c(fit$rhat$beta.comm, fit$rhat$tau.sq.beta, fit$rhat$alpha.comm, fit$rhat$tau.sq.alpha)
  # Multi-season fits return WAIC per species; sum it to the community total
  waic <- spOccupancy::waicOcc(fit)
  if (!is.null(dim(waic))) waic <- colSums(waic[, c("elpd", "pD", "WAIC")])
  git_commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE), error = function(e) NA_character_, warning = function(w) NA_character_)
  out <- list(
    file = fit_file, class = class(fit), date = Sys.Date(),
    machine = paste(Sys.info()[c("sysname", "nodename")], collapse = " "), git_commit = git_commit, spOccupancy_version = as.character(utils::packageVersion("spOccupancy")),
    n_species = dim(fit$y)[1], n_sites = dim(fit$y)[2], n_periods = if (length(dim(fit$y)) == 4) dim(fit$y)[3] else NA,
    occ_terms = colnames(fit$beta.comm.samples), det_terms = colnames(fit$alpha.comm.samples),
    occ_re = colnames(fit$sigma.sq.psi.samples), det_re = colnames(fit$sigma.sq.p.samples),
    mcmc = c(n.samples = fit$n.samples, n.burn = fit$n.burn, n.thin = fit$n.thin, n.chains = fit$n.chains, n.post = fit$n.post),
    run_time_hrs = unname(fit$run.time["elapsed"]) / 3600,
    quantiles = list(beta.comm = post_quantiles(fit$beta.comm.samples), tau.sq.beta = post_quantiles(fit$tau.sq.beta.samples), alpha.comm = post_quantiles(fit$alpha.comm.samples), tau.sq.alpha = post_quantiles(fit$tau.sq.alpha.samples), sigma.sq.psi = post_quantiles(fit$sigma.sq.psi.samples), sigma.sq.p = post_quantiles(fit$sigma.sq.p.samples), theta = post_quantiles(fit$theta.samples)),
    prob_positive = list(beta.comm = colMeans(as.matrix(fit$beta.comm.samples) > 0), alpha.comm = colMeans(as.matrix(fit$alpha.comm.samples) > 0)),
    # Community-level draws (a few thousand numbers each) for contrasts, slopes within each habitat and marginal effects without the full fit
    draws = list(beta.comm = as.matrix(fit$beta.comm.samples), tau.sq.beta = as.matrix(fit$tau.sq.beta.samples), alpha.comm = as.matrix(fit$alpha.comm.samples), tau.sq.alpha = as.matrix(fit$tau.sq.alpha.samples), sigma.sq.psi = if (is.null(fit$sigma.sq.psi.samples)) NULL else as.matrix(fit$sigma.sq.psi.samples), sigma.sq.p = if (is.null(fit$sigma.sq.p.samples)) NULL else as.matrix(fit$sigma.sq.p.samples), theta = if (is.null(fit$theta.samples)) NULL else as.matrix(fit$theta.samples)),
    rhat = fit$rhat[intersect(c("beta.comm", "tau.sq.beta", "alpha.comm", "tau.sq.alpha", "sigma.sq.psi", "sigma.sq.p", "theta"), names(fit$rhat))],
    ess = fit$ESS[intersect(c("beta.comm", "tau.sq.beta", "alpha.comm", "tau.sq.alpha", "sigma.sq.psi", "sigma.sq.p", "theta"), names(fit$ESS))],
    min_ess = ess_all[which.min(ess_all)],
    loadings_median_ess = if (is.null(fit$ESS$lambda)) NA else stats::median(fit$ESS$lambda[fit$ESS$lambda > 0 & is.finite(fit$ESS$lambda)]),
    percent_problematic = 100 * mean(rhat_param > 1.1 | ess_param < 100),
    n_comm_rhat_gt_1.1 = sum(rhat_comm > 1.1), max_comm_rhat = max(rhat_comm),
    n_beta_rhat_gt_1.1 = sum(fit$rhat$beta > 1.1), n_beta = length(fit$rhat$beta), max_beta_rhat = max(fit$rhat$beta),
    mean_geweke = if (is.null(fit$lambda.samples)) NA else mean(abs(coda::geweke.diag(coda::as.mcmc(fit$lambda.samples))$z), na.rm = TRUE),
    waic = waic,
    chain_medians = list(beta.comm = chain_medians(fit$beta.comm.samples), tau.sq.beta = chain_medians(fit$tau.sq.beta.samples), sigma.sq.psi = chain_medians(fit$sigma.sq.psi.samples))
  )
  # Goodness of fit (slow: posterior predictive replicates of every species x site); the log's Site_bayes_p / Rep_bayes_p are the means of the species p-values
  if (gof) {
    gof_tbl <- dplyr::bind_rows(lapply(run_ppc(fit), calc_bayes_p_ms))
    out$gof_species <- gof_tbl
    out$gof_mean <- tapply(gof_tbl$bayes_p, gof_tbl$group, mean)
  }
  if (!is.null(save_path)) {
    saveRDS(out, save_path)
    cat("Saved fit summary", save_path, "\n")
  }
  out
}

plot_parm_estimates <- function(parms_df, ms = FALSE){
  parms_df <- parms_df %>%
    mutate(UCL = estimate + 1.96 * std.error,
           LCL = estimate - 1.96 * std.error)
  p <- ggplot(data = parms_df, aes(x = species_ayerbe, y = estimate, color = term))
  if(ms){ 
    p <- ggplot(data = parms_df, aes(x = Species, y = estimate, color = term))
  }
  p + geom_point(position = position_dodge(.7)) +
    geom_errorbar(position = position_dodge(.7), aes(ymin = LCL, ymax = UCL)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), 
          legend.position = "top") 
}

# Translate spOccupancy coefficient names (e.g. "HabitatSsp:Dist_forest") into readable labels, looking up each part of an interaction in a named vector of term labels; parts missing from the lookup keep their model name
label_params <- function(param, term_labels){
  purrr::map_chr(param, \(name){
    terms <- stringr::str_split_1(name, ":")
    labels <- dplyr::coalesce(unname(term_labels[terms]), terms)
    paste(labels, collapse = " × ")
  })
}

# Species-level occupancy coefficients against the community distribution they are drawn from: for each term, a histogram of the species' posterior medians (beta.samples) with the community normal N(mu, tau^2) overlaid, using the posterior medians of beta.comm and tau.sq.beta. The histogram is narrower than the curve because each species' estimate is shrunk toward mu (most for species with little data); tau^2 is the model's estimate of the true spread
plot_species_coefs <- function(model, terms = NULL, term_labels = NULL){
  spp_tbl <- tibble::tibble(
    name = colnames(model$beta.samples),
    estimate = apply(model$beta.samples, 2, stats::median)
  ) %>%
    dplyr::mutate(
      term = stringr::str_remove(name, "-[^-]+$"), # species names contain no "-", so the last "-" separates term and species
      Species = stringr::str_extract(name, "[^-]+$")
    )
  comm_tbl <- tibble::tibble(
    term = colnames(model$beta.comm.samples),
    mu = apply(model$beta.comm.samples, 2, stats::median),
    tau = sqrt(apply(model$tau.sq.beta.samples, 2, stats::median))
  )
  if(!is.null(terms)){
    spp_tbl <- dplyr::filter(spp_tbl, term %in% terms)
    comm_tbl <- dplyr::filter(comm_tbl, term %in% terms)
  }
  # Normal curve over mu +/- 3 tau for each term
  curve_tbl <- comm_tbl %>%
    dplyr::reframe(x = seq(mu - 3 * tau, mu + 3 * tau, length.out = 200), density = stats::dnorm(x, mu, tau), .by = term)
  # Facet labels carry the community mean and SD, so the plot can be read against summary(model)
  facet_labels <- comm_tbl %>%
    dplyr::mutate(
      label = if(is.null(term_labels)) term else label_params(term, term_labels),
      label = sprintf("%s\nmu = %.2f, tau = %.2f", label, mu, tau)
    ) %>%
    dplyr::select(term, label) %>%
    tibble::deframe()
  ggplot2::ggplot(spp_tbl, ggplot2::aes(x = estimate)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 25, fill = "grey75", color = "white") +
    ggplot2::geom_line(data = curve_tbl, ggplot2::aes(x = x, y = density), color = "#2a78d6", linewidth = 0.8) +
    ggplot2::geom_vline(data = comm_tbl, ggplot2::aes(xintercept = mu), color = "#2a78d6") +
    ggplot2::facet_wrap(~ term, scales = "free", labeller = ggplot2::as_labeller(facet_labels)) +
    ggplot2::labs(x = "Coefficient (logit scale)", y = "Density", caption = "Bars: species posterior medians. Blue: community distribution N(mu, tau^2) and its mean") +
    ggplot2::theme_bw()
}

# Posterior samples behind the random intercepts of one process ("occupancy" or "detection"): the species x level deviations (beta.star / alpha.star), the community-level coefficients (beta.comm / alpha.comm) and the random-effect variances (sigma.sq.psi / sigma.sq.p)
re_samples <- function(model, process = c("occupancy", "detection")){
  process <- match.arg(process)
  if(process == "occupancy"){
    list(star = model$beta.star.samples, comm = model$beta.comm.samples, sigma_sq = model$sigma.sq.psi.samples)
  } else {
    list(star = model$alpha.star.samples, comm = model$alpha.comm.samples, sigma_sq = model$sigma.sq.p.samples)
  }
}

# Readable labels for the collector-team random-effect levels, from the lookup that 01_Wrangling_ms.R writes next to each data list (<list>_collector_teams.csv): a named vector, names = Collector_team_num, values = team plus survey years, e.g. "Gaica distancia 1 (2019)", with the observers (Registered_by) on a second line when include_observers = TRUE
collector_team_labels <- function(lookup_path, include_observers = TRUE){
  lookup <- readr::read_csv(lookup_path, col_types = readr::cols(.default = "c"))
  labels <- sprintf("%s (20%s)", stringr::str_replace(lookup$Collector_team, "(\\D)(\\d+)$", "\\1 \\2"), lookup$Year_grp)
  if(include_observers){
    labels <- dplyr::if_else(is.na(lookup$Registered_by), labels, paste0(labels, "\n", lookup$Registered_by))
  }
  stats::setNames(labels, lookup$Collector_team_num)
}

# Species-specific random intercepts for one grouping variable: each grey point is one species' posterior median deviation for one group level, the black point is the mean across species for that level, and the dashed lines mark +/- 1.96 SD of the fitted random-effect distribution N(0, sigma^2). Each species' deviations are drawn independently, so the mean across species is pulled toward 0 by the model's structure (prior SD about sigma / sqrt(number of species)); a shift shared by all species at one level has no parameter of its own. Optional level_labels (named by level, e.g. from collector_team_labels()) replaces the level numbers on the axis
plot_re_species <- function(model, group, process = c("occupancy", "detection"), level_labels = NULL){
  process <- match.arg(process)
  samples <- re_samples(model, process)
  re_tbl <- tibble::tibble(
    name = colnames(samples$star),
    deviation = apply(samples$star, 2, stats::median)
  ) %>%
    tidyr::separate_wider_regex(name, c(Group = "[^-]+", "-", Level = "[^-]+", "-", Species = "[^-]+")) %>%
    dplyr::filter(Group == group) %>%
    dplyr::mutate(
      Level = if(is.null(level_labels)) Level else dplyr::coalesce(unname(level_labels[Level]), Level),
      Level = forcats::fct_reorder(Level, deviation, .fun = mean)
    )
  re_sd <- sqrt(stats::median(samples$sigma_sq[, group]))
  ggplot2::ggplot(re_tbl, ggplot2::aes(x = deviation, y = Level)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey60") +
    ggplot2::geom_vline(xintercept = c(-1.96, 1.96) * re_sd, linetype = "dashed", color = "grey60") +
    ggplot2::geom_jitter(height = 0.15, width = 0, alpha = 0.35, size = 1, color = "grey40") +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    ggplot2::labs(
      x = sprintf("Species-specific deviation from the species %s intercept (logit scale)", process), y = group,
      caption = sprintf("Grey: one point per species. Black: mean across species. Dashed: +/- 1.96 SD of N(0, sigma^2), sigma = %.2f", re_sd)
    ) +
    ggplot2::theme_bw()
}

# Community-level intercept for each level of one grouping variable, against the global community intercept: within each posterior draw, a level's intercept is the community mean intercept plus the mean of the species deviations for that level, so the intervals carry posterior uncertainty. Dotted lines mark the global mean +/- sigma, the SD of the species x level deviations, for scale. Optional group_sizes (named by level) adds the number of sites or surveys per level to the axis labels, and optional level_labels (named by level, e.g. from collector_team_labels()) replaces the level numbers
plot_re_level_means <- function(model, group, process = c("occupancy", "detection"), group_sizes = NULL, level_labels = NULL, ci_level = 0.95){
  process <- match.arg(process)
  samples <- re_samples(model, process)
  probs <- c((1 - ci_level) / 2, 0.5, 1 - (1 - ci_level) / 2)
  mu_draws <- samples$comm[, "(Intercept)"]
  star_names <- colnames(samples$star)
  star_group <- stringr::str_split_i(star_names, "-", 1)
  star_level <- stringr::str_split_i(star_names, "-", 2)
  levels_k <- unique(star_level[star_group == group])
  level_tbl <- purrr::map(levels_k, \(k){
    level_draws <- mu_draws + rowMeans(samples$star[, star_group == group & star_level == k, drop = FALSE])
    q <- stats::quantile(level_draws, probs)
    tibble::tibble(Level = k, low = q[[1]], median = q[[2]], high = q[[3]])
  }) %>%
    purrr::list_rbind()
  # Labels and sizes are both looked up by the raw level number
  level_tbl <- level_tbl %>%
    dplyr::mutate(
      Label = if(is.null(level_labels)) Level else dplyr::coalesce(unname(level_labels[Level]), Level),
      Label = if(is.null(group_sizes)) Label else sprintf("%s (n = %s)", Label, group_sizes[Level]),
      Level = forcats::fct_reorder(Label, median)
    )
  mu_q <- stats::quantile(mu_draws, probs)
  re_sd <- sqrt(stats::median(samples$sigma_sq[, group]))
  ggplot2::ggplot(level_tbl, ggplot2::aes(x = median, y = Level)) +
    ggplot2::annotate("rect", xmin = mu_q[[1]], xmax = mu_q[[3]], ymin = -Inf, ymax = Inf, fill = "#2a78d6", alpha = 0.15) +
    ggplot2::geom_vline(xintercept = mu_q[[2]], color = "#2a78d6") +
    ggplot2::geom_vline(xintercept = mu_q[[2]] + c(-1, 1) * re_sd, linetype = "dotted", color = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = low, xmax = high)) +
    ggplot2::labs(
      x = sprintf("Community %s intercept (logit scale)", process), y = group,
      caption = sprintf("Points: mean intercept across species per level, %g%% CrI. Blue: global mean and CrI. Dotted: global mean +/- sigma (%.2f)", 100 * ci_level, re_sd)
    ) +
    ggplot2::theme_bw()
}

# Landscape metrics (Dist_forest, ssp, te) per point and year, written by 00b_LSM.R in the sibling Farm_diversity repo
read_lsm <- function(lsm_path = "../Farm_diversity/Derived/Excels/Event_covs_lsm.csv"){
  readr::read_csv(lsm_path, show_col_types = FALSE) %>%
    dplyr::distinct(Id_survey, Year, Dist_forest, ssp, te)
}
