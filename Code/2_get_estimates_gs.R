get_estimates <- function(tr_input, stan_file, thin = 2, iter = 2500) {
  
  key <- tr_input %$% unique(key)
  
  triangle <- tr_input %>%
    select(AccidentYear, DevelopmentLag, inc_paid) %>%
    rename(x1 = AccidentYear,
           x2 = DevelopmentLag,
           y = inc_paid)
  
  mean_sd <- triangle %>%
    mutate(y = ifelse(x1 + x2  > 1996, NA, y)) %$%
    {
      c(attr(scale(y), "scaled:center"),
        attr(scale(y), "scaled:scale"))
    }
  
  triangle_stan <-
    triangle %>%
    filter(x1 + x2 < 1999) %>% 
    mutate(y = ifelse(x1 + x2  > 1996, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(
      x1 = rescale(x1, to = c(0, 1)),
      x2 = rescale(x2, to = c(0, 1)),
      y = (y - mean_sd[1]) / mean_sd[2]
    )
  
  obs_res <- triangle %>%
    mutate(
      lower = ifelse(x1 + x2 == 1997 | x1 + x2 == 1998, y, NA)
    ) %$%
    {
      sum(lower, na.rm = TRUE)
    }
  
  data_stan_fit <- list(
    N = nrow(triangle_stan),
    N1 = triangle_stan %>% filter(!is.na(y)) %>% nrow(),
    N2 = triangle_stan %>% filter(is.na(y)) %>% nrow(),
    x1 = triangle_stan %>% pull(x1),
    x2 = triangle_stan %>% pull(x2),
    y_input = triangle_stan %>% filter(!is.na(y)) %>% pull(y)
  )
  
 # fit <-
 #   stan(
 #     file = stan_file,
 #     data = data_stan_fit,
 #     chains = 4,
 #     thin = thin,
 #     iter = iter
 #   )
  
 
 model <- stan_model(stan_file)

 fit <- sampling(model, 
		 data = data_stan_fit,
 		 chains = 4,
 		 thin = thin,
 		 iter = iter)
  
 
  distr_res <-
    as.array(fit, pars = paste0("y[",  37:55, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    rowSums()
  
  emp_cdf <- ecdf(distr_res)
  
  perc_res <- emp_cdf(obs_res)
  
  result <- tibble(
    key = key,
    obs_res = obs_res,
    lower_25 = quantile(distr_res, .025),
    est_res = mean(distr_res),
    upper_25 = quantile(distr_res, .975),
    perc = perc_res,
    prior_fun = str_extract(stan_file, "(?<=Functions/).+(?=\\.stan)")
  )

  dso_filename <- model@dso@dso_filename
  loaded_dlls <- getLoadedDLLs()
   if (dso_filename %in% names(loaded_dlls)) {
       message("Unloading DLL for model dso ", dso_filename)
       model.dll <- loaded_dlls[[dso_filename]][['path']]
       dyn.unload(model.dll)
       } else {
       message("No loaded DLL for model dso ", dso_filename)
       }

     loaded_dlls <- getLoadedDLLs()
     loaded_dlls <- loaded_dlls[str_detect(names(loaded_dlls), '^file')]
       if (length(loaded_dlls) > 10) {
           for (dll in head(loaded_dlls, -10)) {
	        message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
                dyn.unload(dll[['path']])
	   }
       }
     message("DLL Count = ", length(getLoadedDLLs()))

  return(result)
  
}
