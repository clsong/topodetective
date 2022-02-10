get_fit <- function(A, r, x_obs, deria, adjacency_matrix){
  nspp <- length(x_obs)
  r <-  r %>% as.matrix() %>% as.vector()
  deria <-  deria %>% as.matrix() %>% as.vector()
  
  topology <- as.numeric(t(adjacency_matrix))
  topology[topology != 0] <- 1
  
  # implement quadratic programming for Lotka-Volterra, supplied with target r, A, and allowable tolerance
  fit_qp_LV <- function(A, r, x_obs, deria, A_tol= 1000, r_tol = 1000){
    nspp <- ncol(A)
    # convert A into a vector and create the target parameter vectors by appending r
    Avec <- as.numeric(t(A))
    Avec <- Avec * topology
    
    parvec <- c(Avec,r)
    # number of pars to fit
    npar <- length(parvec)
    # get the bounds
    Alow <- Avec-abs(Avec*A_tol) 
    Aupp <- Avec+abs(Avec*A_tol) 
    rlow <- r-abs(r*r_tol)
    rupp <- r+abs(r*r_tol)

    # inequality constraints
    h <- c(Alow,rlow,-Aupp,-rupp)
    G <- rbind(diag(npar),-diag(npar))
    # equality constraints, assuming at equilibirum
    E <- cbind(as.matrix(bdiag(replicate(nspp,matrix(x_obs,nrow=1),simplify = F))),diag(nspp))
    f <- deria
    # fit the qp model, returning a silent warning if it doesn't converge
    fit <- tryCatch(lsei(A=diag(npar),B=parvec,E=E,F=f,G=G,H=h),  error=function(e) e, warning=function(w) w)
    return(fit)
  }
  
  result_LV <- fit_qp_LV(A=A,r=r,x_obs=x_obs, deria)
  Afit <- t(matrix(result_LV$X[1:nspp^2], nspp,nspp))
  rfit <- result_LV$X[(nspp^2+1):(nspp^2+nspp)]
  return(list(Afit,rfit))
}

get_Jacobian_maynard <- function(dataset, times, adjacency_matrix){
  nspp <- ncol(dataset)
  
  deriative <- 1:ncol(dataset) %>% 
    map_dfc(~as_tibble(predict(sm.spline(times, dataset[,.x]), times, 1)/dataset[,.x])) %>% 
    magrittr::set_colnames(paste0('d', LETTERS[1:nspp], '/', LETTERS[1:nspp])) %>% 
    as.matrix()
  deriative[is.infinite(deriative) & deriative < 0] <- -1
  deriative[is.infinite(deriative) & deriative > 0] <- 1
  
  
  J <- list()
  r <- list()

  J[[1]] <- matrix(runif(nspp^2, 0, 1), nspp,nspp) * adjacency_matrix
  r[[1]] <- deriative[1,] - as.matrix(J[[1]]) %*% dataset[1,]
  for(i in 2:nrow(dataset)){
    res <- get_fit(A = J[[i-1]], r=r[[i-1]], x_obs = dataset[i,], deria=deriative[i,], adjacency_matrix) 
    J[[i]] <- res[[1]]
    r[[i]] <- res[[2]]
  }
  
  combine_J <- function(i){
    J[[i]] %>%
      magrittr::set_colnames(LETTERS[1:nspp]) %>% 
      magrittr::set_rownames(LETTERS[1:nspp]) %>% 
      melt() %>% 
      as_tibble() %>% 
      mutate(time = i) %>% 
      rename(J_maynard = value)
  }
  
  combined_r <- r %>% 
    map_dfc(~as_tibble(.x)) %>% 
    t() %>% 
    as_tibble() %>% 
    magrittr::set_colnames(LETTERS[1:nspp])
  
  combined_J <- 1:length(J) %>% 
    map_dfr(~combine_J(.x)) 
  
  list(combined_r, combined_J)
}

get_Jacobian_maynard_average <- function(dataset, times, num_trials){
  1:num_trials %>% 
    map_dfr(~mutate(get_Jacobian_maynard(dataset, times), trial = .x))
}