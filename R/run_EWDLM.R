#### Epigenome-Wide Distributed Lag Model ####
# create distributed lag models for each site in batches and save relevant outputs from crosspred_abs
# runs EWAS, allfit EWDLM, and abs(allfit) EWDLM
run_EWDLM <- function(mvalues, # tibble with columns "site" and "m_value"
                      mdl_dlm, # model object for individual DLMs 
                      mdl_lm, # model object for individual lms creating EWAS
                      x.cb, # output of exposure data fed into crossbasis function from dlnm
                      x.pen, # output of x.cb fed into cbPen function from dlnm
                      x.means, # mean of exposure across lag period for use in EWAS
                      covariates, # tibble of covariates needed in models
                      batches, # tibble with columns start and end
                      dir_path) # path to directory in which M value database should be stored
  { 
  pb <- progress_bar$new(
    format = "  running [:bar] (:current/:total) | Time Elapsed: :elapsed | Estimated Time to Completion: :eta",
    total = dim(batches)[1])
  
  for (i in 1:dim(batches)[1]) {
    dlm_output <- mvalues[batches$start[i]:batches$end[i],] |> 
      mutate(model_dlm = map(m_value, ~mdl_dlm(.x, x.cb, x.pen, covariates)),
             pred = map(model_dlm, ~crosspred_abs(x.cb, .x, at = 1, abs_sum = T)),
             model_lm = map(m_value, ~mdl_lm(.x, x.means, covariates))) |> 
      select(-m_value) |> 
      mutate(mat_edf = map_dbl(model_dlm, ~sum(.x$edf))) |> 
      select(-model_dlm) |> 
      mutate(ewas_p_raw = map_dbl(model_lm, ~summary(.x)$coefficients[2,4])) |> 
      select(-model_lm) |> 
      mutate(all_coef = map_dbl(pred, ~c(.x$allfit)),
             all_se = map_dbl(pred, ~c(.x$allse)),
             all_ci_lower = map_dbl(pred, ~c(.x$alllow)),
             all_ci_upper = map_dbl(pred, ~c(.x$allhigh)),
             abs_coef = map_dbl(pred, ~c(.x$absfit)),
             abs_se = map_dbl(pred, ~c(.x$absse)),
             abs_ci_lower = map_dbl(pred, ~c(.x$abslow)),
             abs_ci_upper = map_dbl(pred, ~c(.x$abshigh)),
             mat_coef = map(pred, ~c(.x$matfit)),
             mat_ci_lower = map(pred, ~c(.x$matlow)),
             mat_ci_upper = map(pred, ~c(.x$mathigh))) |>
      select(-pred)
    
    dir.create(paste0(dir_path, i, "/"), recursive = T)
    write_parquet(dlm_output, paste0(dir_path, i, "/data.parquet"))
    
    rm(dlm_output)
    pb$tick()
  }
}