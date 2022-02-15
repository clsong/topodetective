#' Get the parameters for some classic population dynamics
#'
#' @param equation_name Name of the classic population dynamics
#' @export
get_classic_dynamics <- function(equation_name,
                                 state_initial,
                                 time_range,
                                 species_num,
                                 topology_ground) {
  if (equation_name == "3_species_food_webs") {
    eqn1_per <- function(x1, x2, x3) (-x1 - 5 * x2 / (3 * x1 + 1) + 1)
    eqn2_per <- function(x1, x2, x3) (-0.1 * x3 / (2 * x2 + 1) + 5 * x1 / (3 * x1 + 1) - 0.4)
    eqn3_per <- function(x1, x2, x3) (.1 * x2 / (2 * x2 + 1) - .01)
    eqns_per <- list(eqn1_per, eqn2_per, eqn3_per)

    if (missing(state_initial)) {
      state_initial <- c(
        x1 = 0.003328141,
        x2 = 0.497268520,
        x3 = 0.868445870
      )
    }
    if (missing(time_range)) {
      time_range <- seq(0, 50, by = .01)
    }
    if (missing(species_num)) {
      species_num <- 3
    }
    if (missing(topology_ground)) {
      topology_ground <- matrix(c(
        -1, -1, 0,
        1, 0, -1,
        0, 1, 0
      ), byrow = T, ncol = 3)
    }
  }
  if (equation_name == "4_species_chaos") {
    eqn1_per <- function(x1, x2, x3, x4) (1 - x1 - 1.09 * x2 - 1.52 * x3)
    eqn2_per <- function(x1, x2, x3, x4) (1 - x2 - .44 * x3 - 1.36 * x4) * .72
    eqn3_per <- function(x1, x2, x3, x4) (1 - 2.33 * x1 - x3 - .47 * x4) * 1.53
    eqn4_per <- function(x1, x2, x3, x4) (1 - 1.21 * x1 - .51 * x2 - .35 * x3 - x4) * 1.27
    eqns_per <- list(eqn1_per, eqn2_per, eqn3_per, eqn4_per)
    # eqn1_per <- function(x1,x2,x3,x4) (1 -x1^1.2 - 1.09*x2 - 1.52*x3)
    # eqn2_per <- function(x1,x2,x3,x4) (1 -x2 - .44*x3 - 1.36 * x4 + x1*x3) * .72
    # eqn3_per <- function(x1,x2,x3,x4) (1- 2.33*x1/(x1+x2) - x3 - .47*x4) * 1.53
    # eqn4_per <- function(x1,x2,x3,x4) (1- 1.21*x1 - .51*x2 - .35*x3- x4) * 1.27

    if (missing(state_initial)) {
      state_initial <- c(x1 = .8, x2 = .4, x3 = .3, x4 = .7)
    }
    if (missing(time_range)) {
      time_range <- seq(0, 50, by = .01)
    }
    if (missing(species_num)) {
      species_num <- 4
    }
    if (missing(topology_ground)) {
      topology_ground <-
        matrix(c(
          -1, -1.09, -1.52, 0,
          0, -1 * .72, -0.44 * .72, -1.36 * .72,
          -2.33 * 1.53, 0, -1 * 1.53, -.47 * 1.53,
          -1.21 * 1.27, -.51 * 1.27, -.35 * 1.27, -1 * 1.27
        ), byrow = T, ncol = 4) %>%
        as_tibble() %>%
        set_names(paste0("x", 1:species_num)) %>%
        mutate(species = paste0("x", 1:species_num)) %>%
        mutate(r = c(1, .72, 1.53, 1.27))
    }
  }
  if (equation_name == "random_LV") {
    if (missing(state_initial)) {
      state_initial <- c(x1 = .8, x2 = .4, x3 = .3, x4 = .7)
    }
    if (missing(time_range)) {
      time_range <- seq(0, 50, by = .01)
    }
    if (missing(species_num)) {
      species_num <- 3
    }
    if (missing(topology_ground)) {
      Sigma <- runif(species_num^2, -1, 1) %>%
        matrix(nrow = species_num)
      diag(Sigma) <- 1.2*apply(Sigma, 1, function(x) sum(abs(x)))
      r_star <- Sigma %*% runif(species_num, 0, 1)
      topology_ground <- Sigma %>%
        as_tibble() %>%
        set_names(paste0("x", 1:species_num)) %>%
        mutate(species = paste0("x", 1:species_num)) %>%
        mutate(r = r_star[, 1])
    }
  }

  assign("state_initial", state_initial, envir = globalenv())
  assign("time_range", time_range, envir = globalenv())
  assign("species_num", species_num, envir = globalenv())
  if(equation_name != "random_LV"){
    assign("eqns_per", eqns_per, envir = globalenv())
  }
  assign("topology_ground", topology_ground, envir = globalenv())
}
