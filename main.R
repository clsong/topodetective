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
library(tidygraph)
library(pspline)

library(InferInteractions)

get_classic_dynamics("chaos") # choose a dynamic

ts <- generate_time_series(eqns_per, time_range, state_initial, species_num) # simulate a dynamic
# Jacobian <- deqn_per(eqns_per, dataset)

plot_time_series(ts)

# fitted ------------------------------------------------------------------

# reg_model <- choose_regression_model("linear")
reg_model <- choose_regression_model("linear")

topology_all <- rep(list(0:1), species_num) %>%
  expand.grid() %>%
  as_tibble() %>%
  mutate(topology_label = row_number()) %>%
  nest(topology = -topology_label)

fitted_models <- ts %>%
  differentiate_ts() %>%
  group_split(species) %>%
  map(~fit_interaction_parameters(., reg_model, topology_all)) %>%
  bind_rows(.id = "species") %>%
  mutate(species = paste0("x", species))

topology_fitted <- fitted_models %>%
  filter(R2 > .95) %>%
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

fitted_models %>%
  filter(
    (species == "x1" & topology_label == 8) |
    (species == "x2" & topology_label == 15) |
    (species == "x3" & topology_label == 14) |
    (species == "x4" & topology_label == 16)
  )
topology_ground

