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

get_estimates <- function(
  tr_input,
  stan_file,
  thin = 2,
  iter = 2500,
  type = "df"
) {
  key <- tr_input %$% unique(key)

  triangle <- tr_input %>%
    select(AccidentYear, DevelopmentLag, inc_paid) %>%
    rename(x1 = AccidentYear, x2 = DevelopmentLag, y = inc_paid)

  mean_sd <- triangle %$%
    {
      c(attr(scale(y), "scaled:center"), attr(scale(y), "scaled:scale"))
    }

  triangle_stan <-
    triangle %>%
    mutate(type = ifelse(x1 + x2 > 1998, "lower", "upper")) %>%
    arrange(desc(type), x1, x2) %>%
    mutate(
      x1_rescale = rescale(x1, to = c(0, 1)),
      x2_rescale = rescale(x2, to = c(0, 1)),
      y = (y - mean_sd[1]) / mean_sd[2]
    )

  obs_res <- triangle %>%
    mutate(lower = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(lower)), x1, x2) %>%
    filter(is.na(lower)) %>%
    pull(y)

  data_stan_fit <- list(
    N = nrow(triangle_stan),
    N_upper = triangle_stan %>% filter(type == "upper") %>% nrow(),
    N_lower = triangle_stan %>% filter(type == "lower") %>% nrow(),
    x1 = triangle_stan %>% pull(x1_rescale),
    x2 = triangle_stan %>% pull(x2_rescale),
    y_input = triangle_stan %>% pull(y),
    y_input_upper = triangle_stan %>% filter(type == "upper") %>% pull(y)
  )

  fit <-
    stan(
      file = stan_file,
      data = data_stan_fit,
      init = function() generate_init_list(data_stan_fit$N),
      chains = 4,
      thin = thin,
      iter = iter,
      control = list(max_treedepth = 15)
    )

  samples_fit <- rstan::extract(fit)

  distr_res <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1])

  emp_cdf <- ecdf(rowSums(distr_res))

  perc_res <- emp_cdf(sum(obs_res))

  avg_pred_tr <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    colMeans()

  lower_pred_tr <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(2, function(x) quantile(x, .025))

  upper_pred_tr <-
    as.array(fit, pars = paste0("y_new[", 56:100, "]")) %>%
    apply(3, rowMeans) %>%
    multiply_by(mean_sd[2]) %>%
    add(mean_sd[1]) %>%
    apply(2, function(x) quantile(x, .975))

  avg_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], avg_pred_tr)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>%
    incr2cum() %>%
    ata() %>%
    attr("vwtd") %>%
    unname()

  lower_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], lower_pred_tr)) %>%
    as.triangle(origin = "x1", dev = "x2", value = "y") %>%
    incr2cum() %>%
    ata() %>%
    attr("vwtd") %>%
    unname()

  upper_ldf <- triangle %>%
    mutate(y = ifelse(x1 + x2 > 1998, NA, y)) %>%
    arrange(desc(!is.na(y)), x1, x2) %>%
    mutate(y = c(y[!is.na(y)], upper_pred_tr)) %>%
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
