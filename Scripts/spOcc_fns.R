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