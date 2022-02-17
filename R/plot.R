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
    facet_wrap(~species, scales = 'free')+
    theme_nice() +
    theme(
      legend.position = 'top',
      legend.title = element_blank()
    )
}

#' Plot the interaction network and intrinsic growth rates
#' @return A ggraph project
#' @param topology A tibble with species intrinsic growth rates and interaction strength
#' @param title title of the ggplot
#' @export
plot_interaction_topology <- function(topology, title) {
  if(missing(title)) title <- ""
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
    ggtitle(title) +
    geom_edge_loop(aes(colour = strength)) +
    scale_edge_colour_viridis() +
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white') +
    theme(aspect.ratio = 1)
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
