# packages ----------------------------------------------------------------
library(deSolve)
library(magrittr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jtools)
library(rsample)
# library(here)
# library(Matrix)
# library(compiler)
# library(lars)
# library(elasticnet)
# library(tidyverse)
# library(reshape2)

# library(cowplot)
# library(rEDM)
# library(limSolve)
# library(Matrix)
# library(pspline)
# library(Deriv)
# library(deSolve)
# library(readxl)
# library(readr)
# library(ggrepel)
# library(corrplot)
# library(DiagrammeR)
# library(rlang)

source("R/equations.r")
source("R/solve_ODEs.R")
source("R/fit_trajectory.R")

get_classic_dynamics("chaos") # choose a dynamic

ts <- generate_time_series(eqns_per, time_range, state_initial, species_num) # simulate a dynamic
# Jacobian <- deqn_per(eqns_per, dataset)

plot_time_series(ts)

# fitted ------------------------------------------------------------------
ts_split <- initial_time_split(ts, prop = 3/4)

train_ts <- training(ts_split)
test_ts  <- testing(ts_split)

topology <- topology_ground
topology[topology != 0 ] <- 1

fitted <- fit_parameters(train_ts, topology)
simu <- fit_simulation(fitted)

evaluate_fit(simu, dataset)

plot_topology(fitted)
plot_fit_vs_simu(dataset, times, simu)


# exmaine all topologies-----------------------------------------------------------------
all_topologies <- rep(list(0:1), num^2 - num) %>%
  expand.grid() %>%
  as_tibble()

pb <- progress_estimated(nrow(all_topologies))
all_fitted <- 1:nrow(all_topologies) %>%
  map(~generate_topology(all_topologies, .)) %>%
  map(~fit_parameters(train, ., map = T))

pb <- progress_estimated(length(all_fitted))
all_fitted_evaluation <- all_fitted %>%
  map(~fit_simulation(., map = T)) %>%
  map(~evaluate_fit(., dataset))

topology_label_accurate <- all_fitted_evaluation %>%
  map_dfr(~summarise(.,
                 cor = mean(correlation),
                 NRMSE = mean(NRMSE))) %>%
  mutate(topology_label = row_number()) %>%
  arrange(NRMSE) %>%
  # filter(cor > .95) %>%
  filter(NRMSE < 0.005) %>%
  pull(topology_label)

topology_label_accurate %>%
  map(~plot_fit_vs_simu(dataset, times, fit_simulation(all_fitted[[.]]), save = T, topology_label = .))

topology_label_accurate %>%
  map(~plot_topology(all_fitted[[.]]))

topology_label_accurate %>%
  map(~plot_interactions(all_fitted[[.]], topology_ground))

# all_fitted_evaluation %>%
#   map_dfr(~summarise(.,
#                      cor = mean(correlation),
#                      NRMSE = mean(NRMSE))) %>%
#   mutate(topology_label = row_number()) %>%
#   arrange(NRMSE)
# dd

# # maynard -----------------------------------------------------------------
# adjacency_matrix <- matrix(c(-1, -1, -1, 0, 0, -1, -1,-1, -1, 0, -1, -1, -1, -1,-1,-1), byrow = T, ncol = 4)
#
# Jacobian_maynard <- get_Jacobian_maynard(dataset, times, adjacency_matrix)
#
# r_all <- Jacobian_maynard[[1]]
# J_all <- Jacobian_maynard[[2]]
#
# J_all %>%
#   ggplot(aes(x=time, y=J_maynard))+
#   geom_line()+
#   facet_grid(Var1~Var2)
#
# r_all %>%
#   mutate(time = row_number()) %>%
#   gather(key, value, -time) %>%
#   ggplot(aes(time, value, color = key))+
#   geom_line()
#
# Sigma <- J_all %>%
#   filter(time == times[length(times)]) %>%
#   select(-time) %>%
#   spread(Var2, J_maynard) %>%
#   select(-Var1) %>%
#   as.matrix()
#
# r <- r_all[nrow(r_all), ] %>%
#   unlist()
#
# simu <- SolutionOde(Sigma, r = r, N0 = state, MaxTime = 100) %>%
#   as_tibble()
#
# simu %>%
#   gather(variable, abundance, -time) %>%
#   ggplot(aes(time, abundance, group=variable, color= variable))+
#   geom_line()
#
#
#
