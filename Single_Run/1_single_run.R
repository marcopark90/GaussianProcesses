tr_input <- filter(full_data, key == "620_ppauto")

stan_file = "./Stan/Output/matern_5_2_sq_exp_prod_v2.stan" # stan file with best porameters
stan_file_niwp = # stan file with best porameters
  results_am <- get_estimates(
  triangle,
  stan_file = stan_file,
  iter = 5000,
  type = "list"
)

results_lr <- get_estimates_lr(
  triangle,
  stan_file = stan_file,
  iter = 5000,
  type = "list"
)

results_am_niw <- get_estimates(
  triangle,
  stan_file = stan_file_niwp,
  iter = 5000,
  type = "list"
)

results_lr_niw <- get_estimates_lr(
  triangle,
  stan_file = stan_file_niwp,
  iter = 5000,
  type = "list"
)

write_rds(results_am, "./Single_Run/results_am.rds")
write_rds(results_lr, "./Single_Run/results_lr.rds")
write_rds(results_am_niw, "./Single_Run/results_am_niwp.rds")
write_rds(results_lr_niw, "./Single_Run/results_lr_niwp.rds")
