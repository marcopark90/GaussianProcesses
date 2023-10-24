get_estimates_lr <- function(tr_input, stan_file, thin = 2, iter = 2500, type = "df") {
  
  key <- tr_input %$% unique(key)
  
  premiums <- tr_input %>% 
    rename(x1 = AccidentYear) %>% 
    mutate(x1 = rescale(x1, to = c(0, 1))) %>% 
    select(x1, earned_prem) %>% 
    distinct(x1, earned_prem)

  triangle <- tr_input %>% 
    mutate(inc_LR = inc_paid / earned_prem) %>% 
    select(AccidentYear, DevelopmentLag, inc_LR) %>%
    rename(x1 = AccidentYear,
           x2 = DevelopmentLag,
           y = inc_LR)
  
  mean_sd <- triangle %>%
    mutate(y = ifelse(x1 + x2  > 1998, NA, y)) %$%
    {
      c(attr(scale(y), "scaled:center"),
        attr(scale(y), "scaled:scale"))
    }
  
  triangle_stan <-
    triangle %>%
    mutate(y = ifelse(x1 + x2  > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(
      x1 = rescale(x1, to = c(0, 1)),
      x2 = rescale(x2, to = c(0, 1)),
      y = (y - mean_sd[1]) / mean_sd[2]
    )
  
  obs_res <- triangle %>%
    mutate(lower = ifelse(x1 + x2 > 1998, NA, y),
           x1 = rescale(x1, to = c(0, 1)),
           x2 = rescale(x2, to = c(0, 1))
    ) %>%
      arrange(desc(!is.na(lower)), x1, x2) %>%
      left_join(premiums, by = "x1") %>% 
      mutate(y = y * earned_prem) %>%
      filter(is.na(lower)) %>% 
      pull(y)

  data_stan_fit <- list(
    N = nrow(triangle_stan),
    N1 = triangle_stan %>% filter(!is.na(y)) %>% nrow(),
    N2 = triangle_stan %>% filter(is.na(y)) %>% nrow(),
    x1 = triangle_stan %>% pull(x1),
    x2 = triangle_stan %>% pull(x2),
    y_input = triangle_stan %>% filter(!is.na(y)) %>% pull(y)
  )
  
  fit <-
    stan(
      file = stan_file,
      data = data_stan_fit,
      chains = 4,
      thin = thin,
      iter = iter
    )
  
  
  prem_56_100 <- triangle_stan %>% 
                  left_join(premiums, by = "x1") %>%
                  filter(is.na(y)) %>% 
                  pull(earned_prem)
  
  distr_res <-
    as.array(fit, pars = paste0("y[",  56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(1, function(x) x * prem_56_100) %>% 
    t()
  
  emp_cdf <- ecdf(rowSums(distr_res))
  
  perc_res <- emp_cdf(sum(obs_res))
  
  avg_pred_tr <-
    as.array(fit, pars = paste0("y[",  56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(1, function(x) x * prem_56_100) %>% 
    t() %>%
    colMeans()
  
  lower_pred_tr <-
    as.array(fit, pars = paste0("y[",  56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>% 
    apply(1, function(x) x * prem_56_100) %>% 
    t() %>% 
    apply(2, function(x) quantile(x, .025))
  
  upper_pred_tr <-
    as.array(fit, pars = paste0("y[",  56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>% 
    apply(1, function(x) x * prem_56_100) %>% 
    t() %>% 
    apply(2, function(x) quantile(x, .975))
  
  
  tr_join <- tr_input %>% 
    select(AccidentYear, DevelopmentLag, earned_prem) %>% 
    rename(x1 = AccidentYear, 
           x2 = DevelopmentLag)
  
  avg_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2  > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], avg_pred_tr)) %>% 
    left_join(tr_join, by = c("x1", "x2")) %>%
    mutate(y = ifelse(x1 + x2  <= 1998, y * earned_prem, y)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>% 
    incr2cum() %>% 
    ata() %>% 
    attr("vwtd") %>% 
    unname()
  
  lower_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2  > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], lower_pred_tr)) %>% 
    left_join(tr_join, by = c("x1", "x2")) %>%
    mutate(y = ifelse(x1 + x2  <= 1998, y * earned_prem, y)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>% 
    incr2cum() %>% 
    ata() %>% 
    attr("vwtd") %>% 
    unname()
  
  upper_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2  > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], upper_pred_tr)) %>%
    left_join(tr_join, by = c("x1", "x2")) %>%
    mutate(y = ifelse(x1 + x2  <= 1998, y * earned_prem, y)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>% 
    incr2cum() %>% 
    ata() %>% 
    attr("vwtd") %>% 
    unname()
  
    if (type == "list") {
    result <- list(
      stan_fit = fit,
      obs_res = sum(obs_res),
      distr_res = rowSums(distr_res),
      res_lower_25 = quantile(rowSums(distr_res), .025),
      est_mean = mean(rowSums(distr_res)),
      res_upper_975 = quantile(rowSums(distr_res), .975),
      perc = perc_res,
      rmse = rmse(obs_res, colMeans(distr_res)),
      lower_ldf = lower_ldf,
      avg_ldf = avg_ldf,
      upper_ldf = upper_ldf)
  
  } else if (type == "df") {
  
    result <- tibble(
      key = key,
      obs_res = sum(obs_res),
      lower_25 = quantile(rowSums(distr_res), .025),
      est_res = mean(rowSums(distr_res)),
      upper_25 = quantile(rowSums(distr_res), .975),
      perc = perc_res,
      rmse = rmse(obs_res, colMeans(distr_res)),
      prior_fun = str_extract(stan_file, "(?<=Functions/).+(?=\\.stan)")
  )
  }
  
  return(result)
  
}