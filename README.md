
<!-- README.md is generated from README.Rmd. Please edit that file -->

# topodetective

<!-- badges: start -->
<!-- badges: end -->

The goal of topodetective is to infer topology of species interactions
from time series of species abundance. We use a regression-based method,
which can efficiently enumerate all possible topologies.

## Installation

You can install the development version of InferInteractions from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("clsong/InferInteractions")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tidyverse)
library(topodetective)

get_classic_dynamics("4_species_chaos") # choose a population dynamic
time_range <- seq(0, 50, by = .1)

set.seed(12345)
ts <- generate_time_series_LV(
  topology = topology_ground,
  state_initial = state_initial,
  time_range = time_range,
  noise = T,
  noise_level = .05
) # simulate a dynamic

plot_time_series(ts)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" /> We
then infer system parameters from time series alone. We choose all
possible topologies.

``` r
reg_model <- choose_regression_model("linear") # 'lasso' and 'ridge' are also available

fitted_models <- ts %>%
  differentiate_ts() %>%
  group_split(species) %>%
  map(fit_interaction_parameters) %>%
  bind_rows(.id = 'species') %>%
  mutate(species = paste0("x", species))
```

We then simulate the dynamics with the inferred parameters. We compare
the fit.

``` r
set.seed(123)
topology_fitted <- fitted_models %>%
  group_by(species) %>%
  top_n(3, R2) %>% 
  sample_n(1) %>%
  ungroup()

ts_simu <- simualte_fitted_dynamics(topology_fitted)
evaluate_fit(ts, ts_simu)
#> Warning in sim - obs: longer object length is not a multiple of shorter object
#> length

#> Warning in sim - obs: longer object length is not a multiple of shorter object
#> length

#> Warning in sim - obs: longer object length is not a multiple of shorter object
#> length

#> Warning in sim - obs: longer object length is not a multiple of shorter object
#> length
#> [1] 0.01567233

plot_true_vs_simu(ts, ts_simu)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

The fitted topology is different from the true topology.

``` r
 library(ggraph)

bind_rows(
  topology_ground %>% 
    mutate(facet = 'Groundtruth'),
  topology_fitted %>% 
    mutate(facet = 'Fitted'),
) %>% 
  plot_interaction_topology(facet = T) 
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
