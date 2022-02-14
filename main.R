# packages ----------------------------------------------------------------
library(deSolve)
library(magrittr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jtools)
library(rsample)
library(tibble)
library(broom)
library(recipes)
library(parsnip)
library(workflows)
library(tune)
library(yardstick)
library(patchwork)
library(tidygraph)
library(ggraph)

theme_set(jtools::theme_nice())
library(InferInteractions)

get_classic_dynamics("chaos") # choose a dynamic

ts <- generate_time_series(eqns_per, time_range, state_initial, species_num) # simulate a dynamic
# Jacobian <- deqn_per(eqns_per, dataset)

plot_time_series(ts)

# fitted ------------------------------------------------------------------

reg_model <- choose_regression_model("linear")

fitted_models <- ts %>%
  preprocess_ts() %>%
  group_split(species) %>%
  map(fit_interaction_parameters) %>%
  bind_rows(.id = 'species') %>%
  mutate(species = paste0("x", species))

# fitted_models %>%
#   mutate(topology_label = as.factor(topology_label)) %>%
#   ggplot(aes(forcats::fct_reorder(topology_label, r2), r2, color = species)) +
#   geom_point()

topology_fitted <- fitted_models %>%
  filter(R2 > .8) %>%
  group_by(species) %>%
  sample_n(1) %>%
  # filter(R2 == max(R2)) %>%
  ungroup()

# topology_ground

ts_simu <- simualte_fitted_dynamics(topology_fitted)

plot_true_vs_simu(ts, ts_simu)
evaluate_fit(ts, ts_simu)

plot_interaction_topology(topology_ground)
plot_interaction_topology(topology_fitted)
