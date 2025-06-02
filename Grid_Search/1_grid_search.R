tr_list <- map(unique(full_data$key), ~ filter(full_data, key == .x))

prior_scripts <- file.path("./Stan/Output") %>% 
  list.files(full.names = T, pattern = "*.stan")

gs_input <- expand_grid(functions = prior_scripts,
			triangles = tr_list)

grid_results <- map2(gs_input$triangles,
                     gs_input$functions,
                     ~get_estimates(.x, stan_file = .y, thin = 2, iter = 2500)) %>%
                     bind_rows()

final_results <- grid_results %>% 
                 mutate(rmse = sqrt((obs_res - est_res)^2)) %>% 
                 group_by(key) %>% 
                 arrange(rmse) %>% 
                 filter(row_number() == 1) %>% 
                 select(key, prior_fun)

# rmse lower is bettter
# ks_pvalue higher is better

write_csv(grid_results, "./Grid_Search/full_results.csv")
write_csv(final_results, "./Grid_Search/final_results.csv")
