generate_init_list <- function(N) {
  A <- matrix(rnorm(N^2), N, N)
  KW_init <- crossprod(A)

  KW_init <- KW_init / max(abs(KW_init))

  list(
    KW_init = KW_init,
    alpha = 1.0,
    rho = 1.0,
    a1 = 2.0,
    b1 = 2.0,
    a2 = 2.0,
    b2 = 2.0,
    psi1 = 1.0,
    psi2 = 1.0,
    sigma = 0.1
  )
}

get_estimates_lr <- function(
  tr_input,
  stan_file,
  thin = 2,
  iter = 2500,
  type = "df"
) {
  key <- tr_input %$% unique(key)

  premiums <- tr_input %>%
    pull(earned_prem)

  triangle_lr <- tr_input %>%
    mutate(inc_paid_lr = inc_paid / earned_prem) %>%
    select(AccidentYear, DevelopmentLag, inc_paid_lr, earned_prem) %>%
    rename(x1 = AccidentYear, x2 = DevelopmentLag, y = inc_paid_lr)

  mean_sd <- triangle_lr %$%
    {
      c(attr(scale(y), "scaled:center"), attr(scale(y), "scaled:scale"))
    }

  triangle_stan_lr <-
    triangle_lr %>%
    mutate(type = ifelse(x1 + x2 > 1998, "lower", "upper")) %>%
    arrange(desc(type), x1, x2) %>%
    mutate(
      x1_rescale = rescale(x1, to = c(0, 1)),
      x2_rescale = rescale(x2, to = c(0, 1)),
      y = (y - mean_sd[1]) / mean_sd[2]
    )

  obs_res <- tr_input %>%
    mutate(
      lower = ifelse(AccidentYear + DevelopmentLag > 1998, NA, inc_paid)
    ) %>%
    filter(is.na(lower)) %>%
    pull(inc_paid) %>%
    sum()

  data_stan_fit_lr <- list(
    N = nrow(triangle_stan_lr),
    N_upper = triangle_stan_lr %>% filter(type == "upper") %>% nrow(),
    N_lower = triangle_stan_lr %>% filter(type == "lower") %>% nrow(),
    x1 = triangle_stan_lr %>% pull(x1_rescale),
    x2 = triangle_stan_lr %>% pull(x2_rescale),
    y_input = triangle_stan_lr %>% pull(y),
    y_input_upper = triangle_stan_lr %>% filter(type == "upper") %>% pull(y)
  )

  fit <-
    stan(
      file = stan_file,
      data = data_stan_fit_lr,
      init = function() generate_init_list(data_stan_fit_lr$N),
      chains = 4,
      thin = thin,
      iter = iter,
      control = list(max_treedepth = 15)
    )

  prem_56_100 <- triangle %>%
    mutate(type = ifelse(x1 + x2 > 1998, "lower", "upper")) %>%
    arrange(desc(type), x1, x2) %>%
    filter(type == "lower") %>%
    pull(earned_prem)

  distr_res <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(1, function(x) x * prem_56_100) %>%
    t()

  emp_cdf <- ecdf(rowSums(distr_res))

  perc_res <- emp_cdf(sum(obs_res))

  avg_pred_tr <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(1, function(x) x * prem_56_100) %>%
    t() %>%
    colMeans()

  lower_pred_tr <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(1, function(x) x * prem_56_100) %>%
    t() %>%
    apply(2, function(x) quantile(x, .025))

  upper_pred_tr <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(1, function(x) x * prem_56_100) %>%
    t() %>%
    apply(2, function(x) quantile(x, .975))

  tr_join <- tr_input %>%
    select(AccidentYear, DevelopmentLag, earned_prem) %>%
    rename(x1 = AccidentYear, x2 = DevelopmentLag)

  avg_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], avg_pred_tr)) %>%
    left_join(tr_join, by = c("x1", "x2")) %>%
    mutate(y = ifelse(x1 + x2 <= 1998, y * earned_prem, y)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>%
    incr2cum() %>%
    ata() %>%
    attr("vwtd") %>%
    unname()

  lower_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], lower_pred_tr)) %>%
    left_join(tr_join, by = c("x1", "x2")) %>%
    mutate(y = ifelse(x1 + x2 <= 1998, y * earned_prem, y)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>%
    incr2cum() %>%
    ata() %>%
    attr("vwtd") %>%
    unname()

  upper_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], upper_pred_tr)) %>%
    left_join(tr_join, by = c("x1", "x2")) %>%
    mutate(y = ifelse(x1 + x2 <= 1998, y * earned_prem, y)) %>%
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
      upper_ldf = upper_ldf
    )
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
