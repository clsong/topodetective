#' Transform a time series data into a tibble that is ready to fit the regression
#'
#' @return A tibble with log change of species abundance with original abundance
#' @param ts Time series data
#' @export
preprocess_ts <- function(ts) {
  log_diff <- ts %>%
    mutate_at(vars(matches("x")), log) %>%
    mutate_all(~ . - lag(.)) %>%
    drop_na() %>%
    mutate(row = row_number()) %>%
    gather(key, dlogN, -time, -row) %>%
    mutate(dlogN_dt = dlogN / time) %>%
    select(-time, -dlogN) %>%
    spread(key, dlogN_dt) %>%
    select(-row) %>%
    mutate(time = ts$time[-nrow(ts)])

  log_diff %>%
    gather(species, log_change, -time) %>%
    left_join(ts, by = "time")
  # log_diff %>%
  #   gather(species, log_change, -time) %>%
  #   ggplot(aes(time, log_change, color=species)) +
  #   geom_line()
}

#' Find a regression methods (linear, lasso, or ridge)
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
#' @return A tibble with simulated time series of species abundances
#' @param Sigma Interaction matrix
#' @param r intrinsic growth rates
#' @param state_initial
#' @param time_range
#' @export
fit_parameters <- function(ts_species,
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
      r2 = map_dbl(workflow_fitted, ~ glance(.)$r.squared),
      estimate = map(workflow_fitted, tidy)
    ) %>%
    unnest(estimate) %>%
    select(topology, r2, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    rename(r = `(Intercept)`)
}

fit_simulation <- function(fitted, map = F) {
  if (map == T) pb$tick()$print()

  Sigma <- fitted %>%
    select(-r) %>%
    as.matrix()

  r <- fitted[, 1] %>% unlist()

  generate_time_series_LV(Sigma, r = r, N0 = state_initial, MaxTime = max(time_range)) %>%
    as_tibble()
}

plot_interactions <- function(fitted, topology_ground) {
  bind_cols(
    fitted %>%
      select(-r) %>%
      as.matrix() %>%
      magrittr::set_colnames(1:4) %>%
      magrittr::set_rownames(1:4) %>%
      melt() %>%
      rename(fit = value),
    topology_ground %>%
      as.matrix() %>%
      magrittr::set_rownames(1:4) %>%
      melt() %>%
      rename(ground = value)
  ) %>%
    ggplot(aes(ground, fit)) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0) +
    theme_bw() +
    theme(aspect.ratio = 1)
}

plot_topology <- function(fitted) {
  num <- ncol(fitted) - 1

  n <- tibble(
    name = 1:num,
    width = .2
  ) %>%
    mutate(
      id_external = name
    )
  e <- fitted %>%
    select(-r) %>%
    as.matrix() %>%
    magrittr::set_colnames(1:num) %>%
    magrittr::set_rownames(1:num) %>%
    melt() %>%
    rename(
      from = Var1,
      to = Var2
    ) %>%
    mutate(color = if_else(value > 0, "dodgerblue", "firebrick2")) %>%
    mutate(value = abs(2 * value))

  create_graph() %>%
    add_nodes_from_table(
      table = n
    ) %>%
    add_edges_from_table(
      table = e,
      from_col = from,
      to_col = to,
      from_to_map = id_external
    ) %>%
    copy_edge_attrs(
      edge_attr_from = value,
      edge_attr_to = penwidth
    ) %>%
    render_graph(layout = "circular")
}

plot_fit_vs_simu <- function(dataset, times, simu, save = F, topology_label = F) {
  p <- bind_rows(
    dataset %>%
      as_tibble() %>%
      mutate(time = times) %>%
      gather(key, value, -time) %>%
      mutate(type = "ground truth"),
    simu %>%
      gather(key, value, -time) %>%
      mutate(type = "fitted")
  ) %>%
    ggplot(aes(time, value, group = type, color = type)) +
    geom_rect(aes(xmin = 0, xmax = times[split_end], ymin = -Inf, ymax = Inf),
      fill = "#F9F4FB", alpha = 1, linetype = 0
    ) +
    geom_rect(aes(xmin = times[split_end], xmax = max(times), ymin = -Inf, ymax = Inf),
      fill = "#FEF7F2", alpha = 1, linetype = 0
    ) +
    geom_line() +
    facet_wrap(~key) +
    scale_colour_manual(values = c("dodgerblue", "#FCBF4A")) +
    theme_classic() +
    theme(
      aspect.ratio = 1
    )

  if (save == T) {
    ggsave(paste0(topology_label, "-fit_vs_simu.pdf"), p)
  } else {
    p
  }
}

calculate_NRMSE <- function(sim, obs) {
  # (sim-obs)^2 %>%
  #   mean() %>%
  #   sqrt()/mean(obs)

  sqrt(mean(sim - obs)^2) / (max(obs) - mean(obs))
}

evaluate_fit <- function(simu, dataset, map = F) {
  if (map == T) pb$tick()$print()

  x <- simu %>%
    filter(time %in% time_range) %>%
    gather(key, value_simu, -time) %>%
    mutate_at(c("time", "value_simu"), as.numeric)

  y <- dataset %>%
    gather(key, value_ground, -time)

  left_join(x, y) %>%
    group_by(key) %>%
    summarise(
      NRMSE = calculate_NRMSE(value_simu, value_ground),
      correlation = cor(value_simu, value_ground)
    )
}

generate_topology <- function(all_topologies, j) {
  offdiag <- all_topologies[j, ] %>%
    unlist()

  topology <- matrix(NA, ncol = num, nrow = num)
  topology[lower.tri(topology)] <- offdiag[1:(length(offdiag) / 2)]
  topology[upper.tri(topology)] <- offdiag[-(1:(length(offdiag) / 2))]
  diag(topology) <- rep(1, num)

  topology
}
