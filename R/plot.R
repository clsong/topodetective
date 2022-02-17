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
#' @param facet Whether the graphs need to be faceted
#' @export
plot_interaction_topology <- function(topology, title, facet = F) {
  if(missing(title)) title <- ""
  if(facet){
    graph <- tbl_graph(
      nodes = topology %>%
        select(species) %>%
        distinct() %>%
        mutate(r = 1),
      edges = topology %>%
        select(-r) %>%
        gather(from, strength, -species, -facet) %>%
        filter(from %in% topology$species) %>%
        unnest(strength) %>%
        filter(strength != 0) %>%
        filter(!is.na(strength)) %>%
        rename(to = species) %>%
        left_join(
          topology %>%
            select(species, r, facet),
          by = c('from' = 'species', 'facet')
        ) %>%
        mutate(width = r)
    )

    p <-
      graph %>%
      ggraph(layout = 'kk') +
      geom_edge_fan(aes(color = strength, width = width),
                    arrow = arrow(length = unit(4, 'mm')),
                    start_cap = circle(2, 'mm'),
                    end_cap = circle(6, 'mm')) +
      geom_edge_loop(aes(colour = strength, width = width)) +
      scale_edge_width(range = c(.5, 2))+
      geom_node_point() +
      ggtitle(title) +
      facet_edges(~facet) +
      scale_edge_colour_viridis() +
      guides(
        edge_width = guide_legend(bquote(r[i])),
        edge_color = guide_legend(bquote(alpha[ij]))
      ) +
      theme_graph(foreground = 'steelblue', fg_text_colour = 'white') +
      theme(aspect.ratio = 1)
  } else{
    graph <- tbl_graph(
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
    )
    p <- graph %>%
      ggraph(layout = 'kk') +
      geom_edge_fan(aes(color = strength),
                    arrow = arrow(length = unit(4, 'mm')),
                    start_cap = circle(2, 'mm'),
                    end_cap = circle(6, 'mm')) +
      geom_node_point(aes(size = r)) +
      ggtitle(title) +
      geom_edge_loop(aes(colour = strength,)) +
      scale_edge_colour_viridis() +
      guides(
        node_size = guide_legend(bquote(r[i])),
        edge_color = guide_legend(bquote(alpha[ij]))
      ) +
      theme_graph(foreground = 'steelblue', fg_text_colour = 'white') +
      theme(aspect.ratio = 1)
  }
  p
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
