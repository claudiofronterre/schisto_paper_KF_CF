# Convert epsg to epsg KM ------------------------------------------------------
epsgKM <- function(x) paste(paste0("+init=epsg:", x), "+units=km")

floor_dec <- function(x, digits=1) round(x - 5*10^(-digits-1), digits)
ceiling_dec <- function(x, digits=1) round(x + 5*10^(-digits-1), digits)

# Envelope for variogram -------------------------------------------------------
variog_envelope <- function (geodata, coords = geodata$coords, 
                             data = geodata$data, 
                             obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  obj.variog$v <- NULL
  if ((is.matrix(data) | is.data.frame(data))) 
    if (ncol(data) > 1) 
      stop("envelops can be computed for only one data set at once")
  if (!is.null(obj.variog$estimator.type)) 
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) {
    return(data[sample(1:n.data)])
  }
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                       n.data = n.data)
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                               2])), as.double(as.vector(sim)), as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                         "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                 vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, 
             data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
             pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
             tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (save.sim == FALSE) 
    simula$data <- NULL
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                  v.upper = limits[2, ])
  if (save.sim) 
    res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

# Calculate and plot the variogram ---------------------------------------------
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = T) {
  require(ggplot2)
  require(geoR)
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))
  if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  envmc <- variog_envelope(coords = coords, data = data, 
                           obj.variog = empvario, nsim = nsim, messages = F)
  dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                        lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                        nbins = empvario$n)
  p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
    geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
    geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
    scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                       breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
    scale_y_continuous(name = "semivariance", 
                       #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                       limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
    ggtitle("Empirical semivariogram") +
    theme_classic()
  p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
  if(show_nbins) p2 else p1
}

# Empirical logit --------------------------------------------------------------
elogit <- function(num, den) {
  log((num + 0.5) / (den - num + 0.5))
}

# Check duplicated rows --------------------------------------------------------
check_dups <- function(cols, data) {
  data$row_id <- 1:nrow(data)
  dups <- data %>%
    group_by_at(cols) %>% 
    filter(n() > 1)  %>% 
    ungroup()
  if(nrow(dups) == 0) {
    dups <- NULL
  } else {
    dups$dup_id <- unclass(factor(apply(dups[, cols], 1, paste0, collapse = "")))
    dups <- dups[order(dups$dup_id), ]
  }
  return(dups)
}


# Calculate and extract UTM km coordinates from data ---------------------------

extract_coords <- function(x) {
  x[c("long", "lat")] <- st_coordinates(x)
  x[c("utm_x", "utm_y")] <- st_coordinates(st_transform(x, crs = epsgKM(unique(x$crs))))
  x$geometry <- NULL
  x$crs <- NULL
  as.data.frame(x)
}

# Fit binomial geostatistical models to list of datasets -----------------------

fitALL <- function(data, path, message = T) {
  
  all_countries <- unique(data$country)
  data_list <- split(data, f = data$country)
  
  # Fit a linear geostatistical model to the empirical logit to get starting
  # value for the variance covariance parameters
  
  map(data_list, safely(function(x) {
    x <- as.data.frame(x)
    country <- unique(x$country)
    cat("Analysing", nrow(x), "survey data from", country, "(",
        match(country, all_countries), "out of", length(all_countries), " countries)\n\n") 
    phi_start <- median(dist(x[, c("utm_x", "utm_y")])) / 3
    f <- elogit ~ 1
    cat("Fitting Linear Geostatistical Model to", country, "data\n")
    geo_linear <- linear.model.MLE(formula = f,
                                   coords =  ~ utm_x + utm_y,
                                   data = x, start.cov.pars = c(phi_start, 1),
                                   kappa = .5, messages = message, method = "BFGS")
    
    par0 <- as.numeric(coef(geo_linear))
    cat("Estimated parameters", par0, "\n\n")
    
    # Fitting binomial geostatistical model
    cmcmc <- list()
    cmcmc[[1]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
    cmcmc[[2]] <- control.mcmc.MCML(n.sim = 10000, burnin = 2000, thin = 8)
    # cmcmc[[3]] <- control.mcmc.MCML(n.sim = 65000, burnin = 5000, thin = 6)
    
    par0 <- as.numeric(coef(geo_linear))
    p <- length(par0) - 3
    f <- positive ~ 1
    niter <- length(cmcmc)
    for(i in 1:niter) {
      cat("\nFitting Geostatistical Binomial Model to", country, "data:",
          "Iteration", i, "of", niter, "\n")
      init_cov_pars <- c(par0[p + 2], par0[p + 3] / par0[p + 1])
      geo_binomial <- binomial.logistic.MCML(formula = f, units.m = ~ examined,
                                             coords = ~ utm_x + utm_y,
                                             data = x,
                                             par0 = par0,
                                             start.cov.pars = init_cov_pars,
                                             control.mcmc = cmcmc[[i]], 
                                             kappa = 0.5,
                                             method = "nlminb", messages = message, 
                                             plot.correlogram = message)
      par0 <- as.numeric(coef(geo_binomial))
      llik <- geo_binomial$log.lik
      cat("Estimated parameters", par0, "Log-likelihood:", llik, "\n\n")
    }
    cat("Saving final results for", country, "\n\n")
    fpath <- paste0(path, "geo_binomial_", gsub(" ", "_", country), ".rds")
    geo_binomial$country <- country
    saveRDS(geo_binomial, fpath)
  }))
}


# Calculate distances ----------------------------------------------------------

ggdist <- function(coords) {
  require(ggplot2)
  dd <- sort(as.numeric(dd <- dist(x[, c("utm_x", "utm_y")], diag = F, upper = F)))
  df <- data.frame(prop = 1:length(dd) / length(dd), dists = dd)
  ggplot(df, aes(x = dists, y = prop)) +
    geom_line() +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Distance between pairs of points (km)", 
         y = "Proportions of pairs of points within distance x")
    
}

# Extract parameters -----------------------------------------------------------
extract_params <- function(x) {
  fit <- readRDS(x)
  params <- as.numeric(PrevMap:::coef.PrevMap(fit))
  n <- length(fit$y)
  return(c(n, params))
}

# Generate scenarios -----------------------------------------------------------

generate_scenarios <- function(params, mu, quantiles) {
  q <- quantiles
  q_sigma2 <- q_phi <- q_nugget <- q
  
  if(is.list(q)) {
    q_mu <- q[[1]]
    q_sigma2 <- q[[2]]
    q_phi <- q[[3]]
    q_nugget <- q[[4]]
  }
  
  quant <- function(x, q) as.numeric(quantile(x, p = q))
  
  # Generate all possible scenarios
  scenarios <- expand.grid(mu = mu,
                           sigma2 = quant(params$`sigma^2`, q_sigma2),
                           phi = quant(params$phi, q_phi),
                           nugget = quant(params$`tau^2`, q_nugget)) %>% 
    as_tibble()
  return(scenarios)
}

# Simulate spatial data --------------------------------------------------------
simulate_spdata <- function(coords, nsims, mu, cov_pars, nugget, cov_model = "matern", kappa = 0.5) {
  require(geoR, quietly = T)
  coords <- as.matrix(coords)
  npts <- nrow(coords)
  dist_mat <- as.matrix(dist(coords))
  Sigma <- diag(nugget, nrow = npts) + cov.spatial(dist_mat, cov.pars = cov_pars, cov.model = cov_model, kappa = kappa)
  p <- length(mu)
  if(any(is.na(match(dim(Sigma), p))))
    stop("Dimension problem!")
  D <- chol(Sigma)
  sim <- t(matrix(rnorm(npts * nsims), ncol = p) %*% D + rep(mu, rep(nsims, p)))
  rownames(sim) <- NULL
  return(sim)
}

# Run scenarios ----------------------------------------------------------------

run_scenarios <- function(scenarios, nsims, schools, n_district, n_subd) {
  
  # Create progress bar
  require(progress, quietly = T)
  
  pb <- progress_bar$new(
    format = "simulating [:bar] :percent eta: :eta",
    total = nrow(scenarios), clear = FALSE, width= 60)
  
  
  # Metrics to be calculated
  metrics <-  c("correct_d", "correct_subd", 
                "correct_school_d", "correct_school_subd",
                "overtreated_d", "overtreated_subd",
                "undertreated_d", "undertreated_subd",
                "se_d", "se_subd", "sp_d", "sp_subd",
                "ppv_d", "ppv_subd", "npv_d", "npv_subd") 
  
  # Output for internal simulations
  output <- matrix(NA, nrow = nsims, ncol = length(metrics))
  colnames(output) <- metrics
  
  # Summary results
  results <- matrix(NA, nrow = nrow(scenarios), ncol = ncol(output) + 2)
  colnames(results) <- c(metrics, "prev_sim", "prevl10")
  
  
  
  for (i in 1:nrow(scenarios)) {
    # Simulate fake data
    mu <- scenarios$mu[i]
    sigma2 <- scenarios$sigma2[i]
    phi <- scenarios$phi[i] 
    nugget <- scenarios$nugget[i]
    fake_data <- simulate_spdata(coords = schools[, c("x", "y")], 
                                 nsims = nsims, 
                                 mu = rep(0, nrow(schools)), 
                                 cov_pars = c(sigma2, phi), 
                                 nugget = nugget)  
    fake_prev <- plogis(mu + fake_data)
    prev_sim <- mean(colMeans(fake_prev))
    prevl10 <- mean(apply(fake_prev, 2, function(x) mean(x < 0.1)))
    for (j in 1:nsims) {
      
      # True underlying prevalence for each school
      schools$true_prev <- fake_prev[, j]
      schools$true_class <- schools$true_prev >= 0.1
      
      # Calculate true prevalence at both district and sub-district level
      districts_true <- schools %>% 
        group_by(ID_district) %>% 
        summarise(d_true_prev = mean(true_prev),
                  d_true_class = d_true_prev >= 0.1)
      
      subd_true <- schools %>% 
        group_by(ID_subd) %>% 
        summarise(subd_true_prev = mean(true_prev),
                  subd_true_class = subd_true_prev >= 0.1)
      
      # Simulate prevalence during survey (test 50 children)
      schools$sample_prev <- rbinom(n = nrow(schools), size = 50, prob = schools$true_prev) / 50
      
      # Select 5 schools per district and 5 per sub-district and calcualte the 
      # observed/empirical prevalence
      districts_sample <- schools %>% 
        group_by(ID_district) %>% 
        sample_n(n_district) %>% 
        summarise(d_sample_prev = mean(sample_prev),  
                  d_sample_class = d_sample_prev >= 0.1)
      
      subd_sample <- schools %>% 
        group_by(ID_subd) %>% 
        sample_n(n_subd) %>% 
        summarise(subd_sample_prev = mean(sample_prev),  
                  subd_sample_class = subd_sample_prev >= 0.1)
      
      districts <- districts_true %>% 
        inner_join(districts_sample, by = "ID_district") %>% 
        mutate(correct_d = d_true_class == d_sample_class)
      
      subd <- subd_true %>% 
        inner_join(subd_sample, by = "ID_subd") %>% 
        mutate(correct_subd = subd_true_class == subd_sample_class)
      
      schools_final <- schools %>% 
        inner_join(districts, by = "ID_district") %>% 
        inner_join(subd, by = "ID_subd") %>% 
        mutate(correct_school_d = d_sample_class == true_class,
               correct_school_subd = subd_sample_class == true_class,
               overtreated_d = d_sample_class > true_class,
               undertreated_d = d_sample_class < true_class,
               overtreated_subd = subd_sample_class > true_class,
               undertreated_subd = subd_sample_class < true_class)
      
      pred <- factor(as.numeric(schools_final$d_sample_class), levels = c(0, 1))
      true <- factor(as.numeric(schools_final$true_class), levels = c(0, 1))
      schools_final$se_d <- caret::sensitivity(pred, true, positive = levels(pred)[2])
      schools_final$sp_d <- caret::specificity(pred, true, negative = levels(pred)[1])
      schools_final$npv_d <- caret::negPredValue(pred, true, negative = levels(pred)[1])
      schools_final$ppv_d <- caret::posPredValue(pred, true, positive = levels(pred)[2])
      
      pred <- factor(as.numeric(schools_final$subd_sample_class), levels = c(0, 1))
      true <- factor(as.numeric(schools_final$true_class), levels = c(0, 1))
      schools_final$se_subd <- caret::sensitivity(pred, true, positive = levels(pred)[2])
      schools_final$sp_subd <- caret::specificity(pred, true, negative = levels(pred)[1])
      schools_final$npv_subd <- caret::negPredValue(pred, true, negative = levels(pred)[1])
      schools_final$ppv_subd <- caret::posPredValue(pred, true, positive = levels(pred)[2])
      
      
      output[j, ] <- apply(schools_final[, metrics], 2, mean)
      
    }
    results[i, ] <- c(apply(output, 2, mean), prev_sim, prevl10)
    pb$tick()
  }
  return(as_tibble(cbind(scenarios, results)))
}



plot_results <- function(data, x, y, col, xlab, ylab, lx, ly, lcol) {
  x_vals <- data[[x]]
  y_vals <- data[[y]]
  col_vals <- unique(data[[col]])
  
  breaks_x <- create_breaks(x_vals, lx)
  breaks_y <- create_breaks(y_vals, ly)
  breaks_col <- create_breaks(col_vals, lcol)
  
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(aes_string(col = col), alpha = 1, shape = 19, size = .85) +
    scale_color_viridis_c("Prevalence", direction = -1, 
                          limits = c(breaks_col[1], breaks_col[length(breaks_col)]), 
                          breaks = breaks_col) +
    labs(x = xlab, y = ylab) +
    # scale_x_continuous(breaks = breaks_x,
    #                    limits = c(breaks_x[1], breaks_x[length(breaks_x)])) +
    # scale_y_continuous(breaks = breaks_y,
    #                    limits = c(breaks_y[1], breaks_y[length(breaks_y)])) +
    theme_ipsum_rc(axis_title_size = 14, base_size = 14)
}

create_breaks <- function(x, l) {
  breaks <- seq(min(x), max(x), l = l) 
  n <- length(breaks)
  breaks[2:(n - 1)] <- round(breaks[2:(n - 1)])
  breaks[1] <- floor(breaks[1])
  breaks[n] <- ceiling(breaks[n])
  breaks
}


plot_results2 <- function(data, x, y, col, shape, xlab, ylab, lx, ly, lcol, exp) {
  x_vals <- data[[x]]
  y_vals <- data[[y]]
  col_vals <- unique(data[[col]])
  
  breaks_x <- create_breaks(x_vals, lx)
  breaks_y <- create_breaks(y_vals, ly)
  breaks_col <- create_breaks(col_vals, lcol)
  
  data$shape <- as.factor(data[[shape]])
  
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(aes_string(col = col, shape = "shape"), alpha = 1) +
    scale_color_viridis_c("Prevalence", direction = -1, 
                          limits = c(breaks_col[1], breaks_col[length(breaks_col)]), 
                          breaks = breaks_col) +
    labs(x = xlab, y = ylab) +
    scale_x_continuous(breaks = breaks_x, 
                       limits = c(breaks_x[1], breaks_x[length(breaks_x)])) +
    scale_y_continuous(breaks = breaks_y, 
                       limits = c(breaks_y[1], breaks_y[length(breaks_y)])) +
    scale_shape_manual(exp, values = c(1, 2, 3), 
                       labels = round(unique(data[[shape]]), 2)) +
    theme_ipsum_rc(axis_title_size = 14, base_size = 14)  
}

create_labels <- function(x, greater = F, smaller = F, text_sep = " - ") {
  n <- length(x)
  x <- gsub(" ", "", format(x))
  labs <- paste(x[1:(n - 1)], x[2:(n)], sep = text_sep)
  if (greater) {
    labs[length(labs)] <- paste("\u2265", x[n - 1])
  }
  if (smaller) {
    labs[1] <- paste("<", x[2])
  }
  
  return(labs)
}

floor_dec <- function(x, digits=1) round(x - 5*10^(-digits-1), digits)
ceiling_dec <- function(x, digits=1) round(x + 5*10^(-digits-1), digits)
