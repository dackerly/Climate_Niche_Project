
# sample a set of site climates from multivariate normal
generate_sites <- function(n_sites = 1000,
                           n_dims = 2,
                           eta = 1, # controls strength of correlations among variables
                           filter = NULL){
  rho <- trialr::rlkjcorr(1, n_dims, eta)
  g <- mvtnorm::rmvnorm(n_sites, sigma = rho)
  colnames(g) <- paste0("clim", 1:n_dims)
  g <- apply(g, 2, scale)
  if(!is.null(filter)) g <- filter(g)
  g
}

# nonlinear transformation function used to generate asymmetric niches
yeo_johnson <- function(y, l){
  # https://en.wikipedia.org/wiki/Power_transform
  dplyr::case_when(l != 0 & y >= 0 ~  ((y + 1) ^ l - 1) / l,
                   l == 0 & y >= 0 ~  log(y + 1),
                   l != 2 & y < 0  ~  -((1-y)^(2-l) - 1) / (2 - l),
                   l == 2 & y < 0  ~  -log(1 - y))
}

# randomly sample niche parameters for one species
generate_species <- function(id, sites, 
                             # all the following are hyperparameters controling distributions of niche params:
                             mu_mu = 0, mu_sigma = 2,
                             tau_shape = 3, tau_rate = 3,
                             rho_eta = 1,
                             lambda_mu = 0, lambda_sigma = .5,
                             alpha_shape1 = 1, alpha_shape2 = 20){
  n_dims <- ncol(sites)
  vars <- colnames(sites)
  mu <- rnorm(n_dims, mu_mu, mu_sigma) # niche means
  tau <- diag(rgamma(n_dims, tau_shape, tau_rate)) # variance
  rho <- trialr::rlkjcorr(1, n_dims, rho_eta) # correlation
  sigma <- tau %*% rho %*% tau # covariance matrix
  lambda <- rnorm(n_dims, lambda_mu, lambda_sigma) # niche asymmetry
  alpha <- rbeta(1, alpha_shape1, alpha_shape2) # max suitability at mu
  colnames(sigma) <- rownames(sigma) <- names(lambda) <- names(mu) <- vars
  list(id = id,
       mu = mu,
       sigma = sigma,
       lambda = lambda,
       alpha = alpha,
       dispersal = NULL, # placeholders
       competition = NULL)
}

# calculate occurrence probability as a funciton of site climate and one species' niche params
suitability <- function(spp, sites){
  centered <- sapply(1:length(spp$mu), 
                     function(i) sites[,i] - spp$mu[i])
  skewed <- sapply(1:length(spp$lambda), 
                   function(i) suppressWarnings(
                     yeo_johnson(centered[,i], spp$lambda[i])))
  d <- mvtnorm::dmvnorm(skewed, sigma = spp$sigma)
  d <- d / mvtnorm::dmvnorm(matrix(rep(0, length(spp$mu)), nrow = 1), sigma = spp$sigma)
  d * spp$alpha
}

# populate communities with species, generating a binary site by species matrix
# very simplistic for now, but could add spatial constraints, carrying capacities, abundances, etc.
assemble_community <- function(sites, species){
  suit <- sapply(species, suitability, sites = sites)
  occ <- apply(suit, 1, function(p) rbinom(length(p), 1, p))
  list(suit = suit,
       occ = t(occ))
}

# fit niche models to sample data
fit_sdms <- function(community, sites, method = "glm", filter = NULL){
  
  # build squared and interaction terms
  sites_df <- as.data.frame(sites)
  sites2 <- sites^2
  colnames(sites2) <- paste0(colnames(sites2), "sq")
  sites2 <- cbind(sites, sites2)
  # sitesx <- apply(combn(1:ncol(sites), 2), 2, 
  #       function(x) sites[,x[1]] * sites[,x[2]])
  # colnames(sitesx) <- paste0("x", 1:ncol(sitesx))
  # sites2 <- cbind(sites2, sitesx)
  sites2_df <- as.data.frame(sites2)
  
  # if requested, include only a subset of sites for model training
  if(!is.null(filter)){
    i <- filter(sites)
  } else {
    i <- 1:nrow(sites)
  }
  
  # fit and predict
  apply(community, 2, function(x){
    if(sum(x[i]) == 0) return(rep(NA, nrow(community)))
    if(method == "maxent"){
      fit <- dismo::maxent(sites_df[i,], x[i])
      return(dismo::predict(fit, sites_df))
    }
    if(method == "glm"){
      df <- as.data.frame(cbind(x = x[i], sites2[i,]))
      fit <- glm(x ~ ., data = df, family = binomial)
      return(predict(fit, sites2_df, type = "response"))
    }
  })
}

# calcualte niche summary stats for every species
summarize_niches <- function(sites, community, fits){
  
  # optimum across all sites
  opt_pa <- apply(sites, 2, function(x) apply(fits, 2, function(z){
    y <- x[z == max(z)]
    if(length(y) > 1) y <- sample(y, 1)
    y}))
  
  # mean of presences
  mean_p <- apply(sites, 2, function(x) apply(community$occ, 2, function(z) mean(x[z == 1])))
  
  # optimum of presences
  opt_p <- apply(sites, 2, function(x) apply(community$occ * fits, 2, function(z) mean(x[z == max(z)])))
  
  list(opt_pa = opt_pa,
       mean_p = mean_p,
       opt_p = opt_p)
}

# community-weighted means
summarize_communities <- function(community, niche_stats){
  lapply(niche_stats,
         function(x) t(apply(community$occ, 1, function(y){
           apply(x, 2, weighted.mean, w = y)
         })))
}

# calcualtes correlations and RMSE
evaluate <- function(a, b){
  funs <- list(r = function(x, y) cor(x, y, use = "pairwise.complete.obs"),
               rmse = function(x, y) sqrt(mean((x - y)^2, na.rm = T)))
  lapply(funs, function(f){
    y <- sapply(a, function(x){
      sapply(1:ncol(x), function(i) f(x[,i], b[,i]))
    })
    rownames(y) <- colnames(b)
    y
  })
}

# turn nested list into clean data frame
format_evals <- function(x){
  require(tidyverse)
  x %>%
    as.data.frame() %>%
    rownames_to_column("var") %>%
    gather(stat, value, -var) %>%
    separate(stat, c("eval", "eval_stat", "niche_metric"), sep = "\\.")
}

plot_evals <- function(x){
  require(tidyverse)
  if(!inherits(x, "data.frame")) x <- format_evals(x)
  if(! "sim" %in% names(x)) x$sim <- 1
  x %>%
    ggplot(aes(niche_metric, value, color = var, group = paste(var, sim))) +
    facet_grid(eval_stat ~ eval, scales = "free") +
    geom_line() + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  
}

# run a complete simulation and evaluation combining the above components
simulate <- function(n_sites = 10000, 
                     n_dims = 2,
                     n_species = 100,
                     eta = 1, 
                     filter = NULL,
                     training_filter = NULL,
                     sdm_method = "glm",
                     ... # args passed to generate_species
){
  # generate data
  sites <- generate_sites(n_sites, n_dims, eta, filter)
  species <- lapply(1:n_species, function(x) generate_species(x, sites, ...))
  community <- assemble_community(sites, species)
  
  # fit models for species/communities
  fits <- fit_sdms(community$occ, sites, method = sdm_method, filter = training_filter)
  niche_stats <- summarize_niches(sites, community, fits)
  cwm <- summarize_communities(community, niche_stats)
  
  # evaluate
  true_optima <- t(sapply(species, function(x) x$mu))
  list(niche = evaluate(niche_stats, true_optima),
       comm = evaluate(cwm, sites))
}

# run multiple simulations (in parallel, optionally)
rep_simulations <- function(n_sims = 10, 
                            n_cores = parallelly::availableCores() - 1,
                            seed = 1,
                            ...){
  require(furrr)
  if(n_cores > 1) plan(multisession, workers = n_cores)
  y <- future_map_dfr(1:n_sims, function(i){
    set.seed(i * seed)
    simulate(...) %>%
      format_evals() %>%
      mutate(sim = i) %>%
      suppressMessages()
  })
  if(n_cores > 1) plan(sequential)
  return(y)
}