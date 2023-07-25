
library(tidyverse)
source("scripts/simulation_functions.r")

# with defaults
# simulate() %>% plot_evals()

# with maxent (SLOW)
# simulate(sdm_method = "maxent") %>% plot_evals()

# with alternative params
# simulate(n_dims = 3, 
#          n_species = 500,
#          training_filter = function(x) which(x[,1] > 0)) %>%
#   plot_evals()

# multiple simulations (each with different climate covariance structures)
# rep_simulations(n_sims = 5, n_dims = 3) %>%
#   plot_evals()
# rep_simulations(n_sims = 8, n_dims = 3, sdm_method = "maxent") %>%
#   plot_evals()


# different range size distributions
d1 <- rep_simulations(n_sims = 100, alpha_shape1 = 1) # most species relatively rare
d2 <- rep_simulations(n_sims = 100, alpha_shape1 = 20) # relatively even range size distribution
p <- bind_rows(d1 %>% mutate(alpha_shape1 = 1),
          d2 %>% mutate(alpha_shape1 = 20)) %>%
  group_by(alpha_shape1, niche_metric, eval_stat, eval) %>%
  summarize(mean = mean(value),
            se = sd(value) / sqrt(length(value))) %>%
  ggplot(aes(niche_metric, mean, ymin = mean - se, ymax = mean + se,
             group = alpha_shape1, color = factor(alpha_shape1))) +
  geom_line() +
  geom_pointrange() +
  facet_wrap(~ eval + eval_stat, scales = "free")
ggsave("results/mk_simulations/alpha_shape.png", p)

d1 <- rep_simulations(n_sims = 100, tau_rate = 1) # broader niches
d2 <- rep_simulations(n_sims = 100, tau_rate = 10) # narrower niches
p <- bind_rows(d1 %>% mutate(tau_rate = 1),
               d2 %>% mutate(tau_rate = 10)) %>%
  group_by(tau_rate, niche_metric, eval_stat, eval) %>%
  summarize(mean = mean(value),
            se = sd(value) / sqrt(length(value))) %>%
  ggplot(aes(niche_metric, mean, ymin = mean - se, ymax = mean + se,
             group = tau_rate, color = factor(tau_rate))) +
  geom_line() +
  geom_pointrange() +
  facet_wrap(~ eval + eval_stat, scales = "free")
ggsave("results/mk_simulations/tau_rate.png", p)
