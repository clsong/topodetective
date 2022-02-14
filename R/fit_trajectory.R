#' Get the deriative of a time series
#'
#' @return A tibble with original speecies abundance and its first derivative
#' @param ts Time series data
#' @export
differentiate_ts <- function(ts) {
  # log_diff <- ts %>%
  #   mutate_at(vars(matches("x")), log) %>%
  #   mutate_all(~ . - lag(.)) %>%
  #   drop_na() %>%
  #   mutate(row = row_number()) %>%
  #   gather(key, dlogN, -time, -row) %>%
  #   mutate(dlogN_dt = dlogN / time) %>%
  #   select(-time, -dlogN) %>%
  #   spread(key, dlogN_dt) %>%
  #   select(-row) %>%
  #   mutate(time = ts$time[-nrow(ts)])
  log_diff <-
    ts %>%
    mutate_at(vars(matches("x")), log) %>%
    gather(key, logN, -time) %>%
    group_by(key) %>%
    mutate(dlogN_dt = predict(sm.spline(time, logN), time, 1)) %>%
    ungroup() %>%
    select(-logN) %>%
    spread(key, dlogN_dt)

  log_diff %>%
    gather(species, log_change, -time) %>%
    left_join(ts, by = "time")
}

#' Find a regression method (linear, lasso, or ridge)
#'
#' @return A model specification
#' @param model_name Name of the model
#' @export
choose_regression_model <- function(model_name) {
  if (model_name == "linear") {
    model <- linear_reg() %>%
      set_engine("lm")
  } else if (model_name == "lasso") {
    model <- linear_reg(
      mode = "regression",
      penalty = 1,
      mixture = 1
    ) %>%
      set_engine("glmnet")
  } else if (model_name == "ridge") {
    model <- linear_reg(
      mode = "regression",
      penalty = 1,
      mixture = 0
    ) %>%
      set_engine("glmnet")
  }
  model
}

#' Fit interaction matrix with given topology and intrinsic growth rates from time series
#'
#' @return A tibble with fitted parameters
#' @param ts_species  Time series data of a single species
#' @param reg_model Which regression model to use
#' @param topology_all All possible topology that are used to fit
#' @export
fit_interaction_parameters <- function(ts_species,
                           reg_model = choose_regression_model("linear"),
                           topology_all) {
  species_num <- ncol(ts_species) - 3

  if (missing(topology_all)) {
    topology_all <- rep(list(0:1), species_num) %>%
      expand.grid() %>%
      as_tibble() %>%
      mutate(topology_label = row_number()) %>%
      nest(topology = -topology_label)
  }

  reg_recipe <- ts_species %>%
    select(-time, -species) %>%
    {
      recipe(log_change ~ ., data = .)
    }

  # df_split <- initial_time_split(ts_species, prop = 3 / 4)
  # df_train <- training(df_split)
  # df_test <- testing(df_split)

  fitted_models <- topology_all %>%
    mutate(workflow_fitted = map(topology, function(topology) {
      if (sum(topology == 0) > 0) {
        reg_recipe_local <- reg_recipe %>%
          step_rm(paste0("x", which(topology == 0)))
      } else{
        reg_recipe_local <- reg_recipe
      }

      workflow() %>%
        add_model(reg_model) %>%
        add_recipe(reg_recipe_local) %>%
        fit(data = ts_species)
    }))

  fitted_models %>%
    mutate(
      R2 = map_dbl(workflow_fitted, ~ glance(.)$r.squared),
      estimate = map(workflow_fitted, tidy)
    ) %>%
    unnest(estimate) %>%
    select(topology_label, topology, R2, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    rename(r = `(Intercept)`)
}

#' Simulate time series from fitted parameters
#'
#' @return A tibble with simulated time series of species abundances
#' @param topology_fitted The fitted parameter of the dynamics
#' @export
simualte_fitted_dynamics <- function(topology_fitted){
                           # state_initial = state_initial,
                           # time_range = time_range) {
  Sigma <- topology_fitted %>%
    select(starts_with("x")) %>%
    mutate(
      across(everything(), ~replace_na(.x, 0))
    ) %>%
    as.matrix()

  r <- topology_fitted %>%
    pull(r) %>%
    unlist()

  generate_time_series_LV(
    Sigma = Sigma,
    r = r,
    state_initial = state_initial,
    time_range = time_range
  )
}

#' Plot simulated time series with fitted parameters vs true time series
#'
#' @return A ggplot2 object
#' @param ts original time series
#' @param ts_simu simualted time series
#' @export
plot_true_vs_simu <- function(ts, ts_simu) {
  bind_rows(
    ts %>%
      mutate(type = 'true'),
    ts_simu %>%
      mutate(type = 'simulated')
  ) %>%
    gather(species, abundance, -time, -type) %>%
    ggplot(aes(time, abundance, group = type, color = type)) +
    geom_line()+
    facet_wrap(~species)+
    theme_nice() +
    theme(
      legend.position = 'top',
      legend.title = element_blank()
    )
}

#' Evaluate how close the simulated time series with fitted parameters is to true time series
#'
#' @return value of root-mean-square deviation
#' @param ts original time series
#' @param ts_simu simualted time series
#' @export
evaluate_fit <- function(ts, ts_simu) {
  calculate_NRMSE <- function(sim, obs) {
    sqrt(mean(sim - obs)^2) / (max(obs) - mean(obs))
  }

  2:ncol(ts) %>%
    map_dbl(~calculate_NRMSE(unlist(ts_simu[,.]), unlist(ts[,.]))) %>%
    mean()
}

#' Plot the interaction network and intrinsic growth rates
#'
#' @return A ggraph project
#' @param topology A tibble with species intrinsic growth rates and interaction strength
#' @export
plot_interaction_topology <- function(topology) {
  tbl_graph(
    nodes = topology %>%
      select(species, r),
    edges = topology %>%
      select(-r) %>%
      gather(from, strength, -species) %>%
      filter(from %in% topology$species) %>%
      unnest(strength) %>%
      filter(strength != 0) %>%
      filter(!is.na(strength)) %>%
      rename(to = species)
    ) %>%
    ggraph(layout = 'kk') +
    geom_edge_fan(aes(color = strength),
                  arrow = arrow(length = unit(4, 'mm')),
                  end_cap = circle(6, 'mm')) +
    geom_node_point(aes(size = r)) +
    geom_edge_loop(aes(colour = strength)) +
    scale_edge_colour_viridis() +
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white') +
    theme(aspect.ratio = 1)
}
