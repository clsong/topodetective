#' Simulated time series of species abundances
#'
#' @return A tibble with simulated time series of species abundances
#' @param eqns_per Population dynamics
#' @param time_range Time range to run the simulation
#' @param state_initial Initial species abundances
#' @param species_num Number of all species
#' @export
generate_time_series <- function(eqns_per,
                                 time_range,
                                 state_initial,
                                 species_num) {
  ODE_system <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      expression <- paste0("x", 1:species_num) %>%
        paste(collapse = ",")

      res <- eval(parse(text = paste0("map_dbl(1:species_num, ~eqns_per[[.]](", expression, "))")))

      list(res * state)
    })
  }

  parameters <- c(a = 1)

  ode(
    y = state_initial, times = time_range,
    func = ODE_system, parms = parameters
  ) %>%
    as_tibble() %>%
    mutate_all(as.numeric)
}

#' Simulated time series of species abundances
#'
#' @return A tibble with simulated time series of species abundances
#' @param Sigma Interaction matrix
#' @param r intrinsic growth rates
#' @param state_initial Initial species abundances
#' @param time_range Time range to run the simulation
#' @export
generate_time_series_LV <- function(Sigma, r, state_initial, time_range) {
  alpha <- Sigma
  parms <- list(r = r, alpha = alpha)
  # delta_t <- 0.01 # time step
  # time_step <- seq(0, max(time_range), by = delta_t) # sequence of time
  model <- function(t, N, parms) {
    dN <- N * (parms$r + parms$alpha %*% N) + 1e-14
    list(dN)
  }
  ode(state_initial, time_range, model, parms, method = "ode45") %>%
    as_tibble() %>%
    mutate_all(as.numeric)
}


#' Plot time series of species abundance
#' @return A ggplot2 object
#' @param dataset
#' @export
plot_time_series <- function(dataset) {
  dataset %>%
    gather(species, abundance, -time) %>%
    ggplot(aes(time, abundance, group = species, color = species)) +
    geom_line() +
    theme_nice()
}

deqn_per <- function(eqns_per, dataset) {
  D <- function(i) {
    eqn <- eqns_per[[i]]

    get_D <- function(eqn, ...) {
      Deriv(eqn)(...) %>%
        as.matrix(nrow = 1) %>%
        t() %>%
        as_tibble()
    }

    expression <- paste0("..", 1:num) %>%
      paste(collapse = ",")
    eval(parse(text = paste0("pmap_dfr(as_tibble(dataset), ~get_D(eqn,", expression, "))"))) %>%
      mutate(time = row_number()) %>%
      gather(Var2, J_analytic, -time) %>%
      mutate(Var1 = paste0("x", i))
  }

  1:num %>%
    map_dfr(~ D(.x))
}

