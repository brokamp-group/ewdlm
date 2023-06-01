#### Exposure Window Simulation ####
# define the beginning and end of a window of specific type
window_bounds <- function(type, num_lags) {
  if (type == "null" || type == 0) {
    start <- NA
    end <- NA
    
  } else if (type == "early" || type == 1) {
    start <- 1
    end <- floor(num_lags * 0.2)
    
  } else if (type == "middle" || type == 2) {
    start <- floor(num_lags * 0.4)
    end <- floor(num_lags * 0.6)
    
  } else if (type == "late" || type == 3) {
    start <- floor(num_lags * 0.8)
    end <- num_lags
    
  } else if (type == "full" || type == 4) {
    start <- 1
    end <- num_lags
    
  } else if (type == "dual" || type == 5) {
    start <- c(floor(num_lags * 0.2), floor(num_lags * 0.6))
    end <- c(floor(num_lags * 0.4), floor(num_lags * 0.8))
    
  } else stop("Invalid window type")
  
  bounds <- tibble(start = start,
                   end = end)
  return(bounds)
}


# create a triangular peak of lag weights/coefficients matching the length of the window
create_peak <- function(start, end, peak_effect) {
  if (!is.na(start) && !is.na(end)) {
    mid <- floor((end-start)/2) + start
    lag_coef <- c(seq(0, peak_effect, length.out = mid - start + 1),
                  seq(peak_effect, 0, length.out = end - mid))
    return(lag_coef)
  } else {
    return(NA)
  }
}


# generate lag weights based on window type and peak effect
lag_effects <- function(type, num_lags, peak_effect) {
  
  if (type == "null" || type == 0) {
    effects <- rep(0, num_lags)
  } else if (type == "full" || type == 4) {
      effects <- rep(peak_effect, num_lags)
  } else {
    effects <- rep(0, num_lags)
    bounds <- window_bounds(type, num_lags)
    effects[bounds$start[1]:bounds$end[1]] <- create_peak(bounds$start[1], bounds$end[1], peak_effect)
    
    if (type == "dual" || type == 5) {
      effects[bounds$start[2]:bounds$end[2]] <- create_peak(bounds$start[2], bounds$end[2], peak_effect * -1)
    }
  }
  
  return(effects)
}


#### Methylome Simulation ####
# for one site, simulate M values that are associated with exposure in a specific window
sim.m <- function(m_mean, m_sd, exposure, lag_coefs, prop_noise) {
  num_obs <- dim(exposure)[1]
  
  if (sum(lag_coefs) != 0) {
    effect_raw <- exposure %*% lag_coefs
    effect_raw_mean <- mean(effect_raw)
    effect_raw_sd <- sd(effect_raw)
    
    effect_std <- (effect_raw - effect_raw_mean)/effect_raw_sd
    
    R2 <- 1 - prop_noise
    ssr <- sum(effect_std^2)
    e <- resid(lm(rnorm(dim(pm_wide)[1]) ~ effect_std))
    noise <- e * sqrt((1-R2)/R2 * ssr/(sum(e^2)))
    
  } else {
    effect_std <- rep(0, num_obs)
    noise <- rnorm(num_obs)
  }
  
  m_sim_raw <- effect_std + noise
  m_sim_std <- m_sim_raw / sd(m_sim_raw)
  m_sim <- (m_sim_std * m_sd) + m_mean
  
  return(m_sim)
}


# use lag_effects() and sim.m() to simulate epigenome-wide methylation for each subject
# allows for random assignment of lag_form and coef_max with "rand" input
sim.methylome <- function(m_summary, # requires columns named "site", "mean", and "sd"
                          exposure, # requires lags as columns and observations as rows
                          prop_effects, 
                          coef_max,
                          lag_form,
                          noise) {
  
  if (typeof(lag_form) == "character") {
    lag_form_num <- switch(lag_form,
                           "early" = 1,
                           "middle" = 2,
                           "late" = 3,
                           "full" = 4,
                           "dual" = 5,
                           "rand" = 99)
  } else {
    lag_form_num <- lag_form
  }
  
  num_sites <- dim(m_summary)[1]
  num_lags <- dim(exposure)[2]
  
  assoc_exists <- rbinom(num_sites, 1, prop_effects)
  
  if (lag_form_num == 99) {
    window_type <- assoc_exists * sample(1:5, num_sites, replace = T)
  } else {
    window_type <- assoc_exists * rep(lag_form_num, num_sites)
  }
  
  if (coef_max == "rand") {
    window_peak <- assoc_exists * sample(seq(0.2, 0.6, 0.1), num_sites, replace = T)
  } else {
    window_peak <- assoc_exists * as.numeric(coef_max)
  }
  
  lag_coefs <- map2(window_type, 
                    window_peak, 
                    ~lag_effects(.x, num_lags, .y),
                    .progress = T)
  
  m_value <- pmap(list(m_summary$mean,
                       m_summary$sd,
                       lag_coefs,
                       noise),
                  ~sim.m(..1, ..2, exposure, ..3, ..4),
                  .progress = T)
  
  m_sims <- tibble(site = m_summary$site,
                   assoc_exists = assoc_exists,
                   window_type = window_type,
                   window_peak = window_peak,
                   lag_coefs = lag_coefs,
                   m_value = m_value)
  
  return(m_sims)
}