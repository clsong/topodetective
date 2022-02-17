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
#' @param topology interaction strength and intrinsic growth rates
#' @param state_initial Initial species abundances
#' @param time_range Time range to run the simulation
#' @param noise Whether there should be noise
#' @export
generate_time_series_LV <- function(topology, state_initial, time_range,
                                    noise = F, noise_level = .01) {
  alpha <- topology %>%
    select(starts_with("x")) %>%
    mutate(
      across(everything(), ~replace_na(.x, 0))
    ) %>%
    as.matrix()

  r <- topology %>%
    pull(r) %>%
    unlist()

  parms <- list(r = r, alpha = alpha)

    model <- function(t, N, parms) {
      dN <- N * (parms$r + parms$alpha %*% N) + 1e-14
      list(dN)
    }

    simu <- ode(state_initial, time_range, model, parms, method = "ode45") %>%
      as_tibble() %>%
      mutate_all(as.numeric) %>%
      set_names(c("time", paste0("x", 1:length(state_initial))))

    if(noise){
      simu <- simu %>%
        mutate_at(vars(starts_with('x')), function(a){a * rnorm(length(a), 1, noise_level)})
    }

    simu
}


#' Plot time series of species abundance
#' @return A ggplot2 object
#' @param ts A tibble containing time series
#' @export
plot_time_series <- function(ts) {
  ts %>%
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

