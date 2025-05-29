tr_list <- map(unique(full_data$key), ~ filter(full_data, key == .x)) %>% 
            as_tibble_col(column_name = "triangle_data") %>% 
            mutate(map_df(triangle_data, ~unique(.x %>% magrittr::extract("key"))))

df_input <- read_csv("./Grid_Search/final_results.csv") %>%
              mutate(cov_function = paste0(file.path("./Stan/Output", prior_fun), ".stan")) %>%
              left_join(tr_list)

df_input_niwp <- read_csv("./Grid_Search/final_results.csv") %>% 
              mutate(cov_function = paste0(file.path("./Stan/Output", prior_fun), "_niwp.stan")) %>% 
              left_join(tr_list)

results_am <- map2(df_input$triangle_data,
                     df_input$cov_function,
                     ~get_estimates(.x, stan_file = .y, thin = 2, iter = 3000)) %>%
              bind_rows()

results_lr <- map2(df_input$triangle_data,
                   df_input$cov_function,
                   ~get_estimates_lr(.x, stan_file = .y, thin = 2, iter = 3000)) %>%
              bind_rows()

results_am_niwp <- map2(df_input_niwp$triangle_data,
                   df_input_niwp$cov_function,
                   ~get_estimates(.x, stan_file = .y, thin = 2, iter = 3000)) %>%
              bind_rows()

results_lr_niwp <- map2(df_input_niwp$triangle_data,
                   df_input_niwp$cov_function,
                   ~get_estimates_lr(.x, stan_file = .y, thin = 2, iter = 3000)) %>%
              bind_rows()

# rmse lower is bettter
# ks_pvalue higher is better

write_csv(results_am, "./Full_Run/full_results_am.csv")
write_csv(results_lr, "./Full_Run/full_results_lr.csv")
write_csv(results_am_niwp, "./Full_Run/full_results_am_niwp.csv")
write_csv(results_lr_niwp, "./Full_Run/full_results_lr_niwp.csv")
