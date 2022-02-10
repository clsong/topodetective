fit_parameters <- function(data, topology, map = F){
  if(map == T) pb$tick()$print()

  species_num <- ncol(data) - 1

  log_diff <- data %>%
    mutate_at(vars(matches("x")), log) %>%
    mutate_all(~.-lag(.)) %>%
    drop_na() %>%
    mutate(row = row_number()) %>%
    gather(key, dlogN, -time, -row) %>%
    mutate(dlogN_dt = dlogN/time) %>%
    select(-time, -dlogN) %>%
    spread(key, dlogN_dt) %>%
    select(-row)

  abundance <- data[-nrow(data),-1]

  data_fit <- bind_cols(
    log_diff,
    abundance
  ) %>%
    set_colnames(c(paste0("dx", 1:num), paste0("x", 1:num)))

  get_coef <- function(variable, data_fit, topology){
    factors <- paste0("x", which(topology[variable, ] == 1))
    formula <- as.formula(paste(paste0("dx", variable), "~", paste(factors, collapse="+")))
    coef <- lm(formula, data= data_fit) %>%
      broom::tidy() %>%
      pull(estimate)
    names(coef) <- c('r', factors)
    if(length(coef) < num + 1){
      zero <- rep(0, num + 1 - length(coef))
      names(zero) <- setdiff(c('r', paste0("x", 1:num)), names(coef))
      coef <- c(coef, zero)
    }
    coef %>%
      enframe() %>%
      arrange(name) %>%
      select(-name)
  }

  1:species_num %>%
    map_dfc(~get_coef(variable = ., data_fit, topology)) %>%
    t() %>%
    as_tibble() %>%
    magrittr::set_colnames(c('r', paste0("a", 1:num)))
}

fit_simulation <- function(fitted, map = F){
  if(map == T) pb$tick()$print()

  Sigma <- fitted %>%
    select(-r) %>%
    as.matrix()

  r <- fitted[,1] %>% unlist()

  SolutionOde(Sigma, r = r, N0 = state, MaxTime = max(times)) %>%
    as_tibble()
}

plot_interactions <- function(fitted, topology_ground){
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
    geom_point(size =2) +
    geom_abline(slope=1,intercept = 0) +
    theme_bw()+
    theme(aspect.ratio = 1)
}

plot_topology <- function(fitted){
  num <- ncol(fitted)-1

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
    rename(from = Var1,
           to = Var2) %>%
    mutate(color = if_else(value > 0, 'dodgerblue', 'firebrick2')) %>%
    mutate(value = abs(2*value))

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

plot_fit_vs_simu <- function(dataset, times, simu, save = F, topology_label = F){

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

  if(save == T){
    ggsave(paste0(topology_label, '-fit_vs_simu.pdf'), p)
  } else{
    p
  }
}

calculate_NRMSE <- function(sim, obs){
  # (sim-obs)^2 %>%
  #   mean() %>%
  #   sqrt()/mean(obs)

  sqrt(mean(sim-obs)^2)/(max(obs)-mean(obs))
}

evaluate_fit <- function(simu, dataset, map = F){
  if(map == T) pb$tick()$print()

  x <- simu %>%
    filter(time %in% times) %>%
    gather(key, value_simu, -time)

  y <- dataset %>%
    as_tibble() %>%
    mutate(time = times) %>%
    gather(key, value_ground, -time)

  left_join(x,y) %>%
    group_by(key) %>%
    summarise(NRMSE = calculate_NRMSE(value_simu, value_ground),
              correlation = cor(value_simu, value_ground)
    )
}

generate_topology <- function(all_topologies, j){
  offdiag <- all_topologies[j,] %>%
    unlist()

  topology <- matrix(NA, ncol = num, nrow = num)
  topology[lower.tri(topology)] <- offdiag[1:(length(offdiag)/2)]
  topology[upper.tri(topology)] <- offdiag[-(1:(length(offdiag)/2))]
  diag(topology) <- rep(1, num)

  topology
}

