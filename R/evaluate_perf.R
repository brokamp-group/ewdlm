#### Evaluate EWDLM Performance ####
# helper function to get performance metrics from a comparison of predictions and truth 
get_perf <- function(true, pred, cen_good = NULL, output = c("stats", "counts")) {
  tp <- true & pred
  fn <- true & !pred
  tn <- !true & !pred
  fp <- !true & pred
  
  if (!is_empty(cen_good)) {
    orig_tp <- tp
    orig_fp <- fp
    tp <- orig_tp & cen_good
    fp <- orig_fp | (orig_tp & !cen_good)
  }  
  
  if (output == "stats") {
    out <- tibble(ppv = sum(tp) / sum(pred),
                  npv = sum(tn) / sum(!pred), 
                  sens = sum(tp) / sum(true),
                  spec = sum(tn) / sum(!true))
    
  } else if (output == "counts") {
    out <- tibble(num_tp = sum(tp),
                  num_fn = sum(fn),
                  num_tn = sum(tn),
                  num_fp = sum(fp))
  }
  
  return(out)
}


# helper function to identify lags which are acceptable window centers based on window type
window_range <- function(type, num_lags, coverage) {
  lags <- 0:(num_lags - 1)
  
  get_win_range <- function(start, end) {
    cen <- median(lags[start:end])
    len <- length(lags[start:end])
    wing <- (len * coverage) / 2
    win_range <- lags[floor(cen - wing):ceiling(cen + wing)]
    return(win_range)
  }
  
  if (type == "null" || type == 0) {
    win_range <- NA
  } else {
    if (type == "full" || type == 4) {
      win_range <- lags
    } else {
      bounds <- window_bounds(type, num_lags)
      win_range <- get_win_range(bounds$start[1], bounds$end[1])
      
      if (type == "dual" || type == 5) {
        win_range <- c(win_range, get_win_range(bounds$start[2], bounds$end[2]))
      }
    }
  }
  
  return(win_range)
}


# assess performance of the analysis 
# outputs performance metrics, confusion matrix values, and a summary of windows
ewdlm_perf <- function(output_path, # path to directory where run_EWDLM output
                       simulation_path, # path to directory where simulated M values were saved
                       empirical_cutoff = F, # option to use user-defined alpha cutoffs after FDR correction
                       alpha = 0.05) # IF empirical_cutoff == T, alpha should be a named numerical vector with cutoffs for "ewas", "dlm_all", and "dlm_abs"
  {
  if (empirical_cutoff) {
    alpha_ewas <- alpha["ewas"]
    alpha_dlm_all <- alpha["dlm_all"]
    alpha_dlm_abs <- alpha["dlm_abs"]
  } else {
    alpha_ewas <- alpha_dlm <- alpha_dlm_abs <- alpha
  }
  
  print(paste0("Evaluating output at ", output_path, " (", Sys.time(), ")"))
  cumul_fits <- open_dataset(output_path, partitioning = c("batch")) |> 
    select(-starts_with("mat"), -batch) |>
    collect() |> 
    mutate(ewas_p_adj = p.adjust(ewas_p_raw, method = "fdr"),
           sig_ewas = (ewas_p_adj < alpha_ewas),
           dlm_all_p_raw = pnorm(abs(all_coef/all_se), lower.tail = F),
           dlm_all_p_adj = p.adjust(dlm_all_p_raw, method = "fdr"),
           sig_dlm_all = (dlm_all_p_adj < alpha_dlm_all),
           dlm_abs_p_raw = pnorm(abs(abs_coef/abs_se), lower.tail = F),
           dlm_abs_p_adj = p.adjust(dlm_abs_p_raw, method = "fdr"),
           sig_dlm_abs = (dlm_abs_p_adj < alpha_dlm_abs))
  
  sites_sig_dlm_all <- cumul_fits$site[cumul_fits$sig_dlm_all]
  sites_sig_dlm_abs <- cumul_fits$site[cumul_fits$sig_dlm_abs]
  sites_combined <- union(sites_sig_dlm_all, sites_sig_dlm_abs)
  
  true_fits <- read_parquet(paste0(simulation_path, "m_sims.parquet"),
                            col_select = c("site", "window_type"))
  
  mat_fits <- open_dataset(output_path, partitioning = c("batch")) |> 
    select(site, mat_coef, mat_ci_lower, mat_ci_upper) |> 
    dplyr::filter(site %in% sites_combined) |> 
    collect()
  
  num_lags <- length(mat_fits$mat_coef[[1]])
  
  mat_fits_long <- mat_fits |> 
    mutate(lag = list(0:(num_lags - 1))) |> 
    unnest(cols = c(mat_coef, mat_ci_lower, mat_ci_upper, lag)) |> 
    mutate(signSum = factor(as.character(sign(mat_ci_lower) + sign(mat_ci_upper)),
                            levels = c("-2", "0", "2")))
  
  windows <- mat_fits_long |> 
    dplyr::filter(signSum != "0") |> 
    group_by(site) |> 
    nest() |> 
    mutate(window_id = map(data, ~cumsum(c(1, abs(.$lag[-length(.$lag)] - .$lag[-1]) > 1)))) |> 
    unnest(cols = c(data, window_id))
  
  window_summary <- windows |> 
    group_by(site, window_id) |> 
    summarize(start = min(lag),
              end = max(lag),
              length = n(),
              midpoint = start + floor(length/2),
              max_effect_coef = mat_coef[abs(mat_coef) == max(abs(mat_coef))],
              max_effect_lag = lag[mat_coef == max_effect_coef],
              .groups = "drop") |> 
    ungroup() |> 
    group_by(site) |> 
    nest() |> 
    left_join(true_fits, by = "site") |> 
    mutate(acceptable_cen = map(window_type, ~window_range(.x, num_lags, 0.5))) |> 
    unnest(cols = "data") |> 
    mutate(in_range = map2_lgl(midpoint, acceptable_cen, ~ .x %in% .y)) |> 
    ungroup() |> 
    select(-acceptable_cen)
  
  window_matches <- window_summary |> 
    group_by(site) |> 
    summarize(cen_good = sum(as.logical(in_range)) > 0)
  
  comparison_tibble <- true_fits |> 
    mutate(assoc_exists = (window_type != 0)) |> 
    left_join(select(cumul_fits, site, starts_with("sig")),
              by = "site") |> 
    left_join(select(window_matches, site, cen_good),
              by = "site")
  
  perf_stats <- comparison_tibble |> 
    summarize(ewas = get_perf(assoc_exists, sig_ewas, output = "stats"),
              dlm_all = get_perf(assoc_exists, sig_dlm_all, output = "stats"),
              dlm_all_loc = get_perf(assoc_exists, sig_dlm_all, cen_good, output = "stats"),
              dlm_abs = get_perf(assoc_exists, sig_dlm_abs, output = "stats"),
              dlm_abs_loc = get_perf(assoc_exists, sig_dlm_abs, cen_good, output = "stats"))
  
  perf_by_window <- comparison_tibble |> 
    group_by(window_type) |>
    summarize(ewas = get_perf(assoc_exists, sig_ewas, output = "counts"),
              dlm_all = get_perf(assoc_exists, sig_dlm_all, output = "counts"),
              dlm_all_loc = get_perf(assoc_exists, sig_dlm_all, cen_good, output = "counts"),
              dlm_abs = get_perf(assoc_exists, sig_dlm_abs, output = "counts"),
              dlm_abs_loc = get_perf(assoc_exists, sig_dlm_abs, cen_good, output = "counts"))
  
  return(list(perf_stats, perf_by_window, window_summary))
}