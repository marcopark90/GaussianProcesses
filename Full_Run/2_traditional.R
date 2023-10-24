tr_list <- map(unique(full_data$key), ~ filter(full_data, key == .x)) %>% 
  as_tibble_col(column_name = "triangle_data") %>% 
  mutate(map_df(triangle_data, ~unique(.x %>% magrittr::extract("key"))))

trad_res <- tr_list %>% 
  mutate(obs_res = map(triangle_data, ~filter(., .$AccidentYear + .$DevelopmentLag > 1998) %>% 
                         pull(inc_paid)),
         boot = map(triangle_data, ~mutate(., cum_paid = if_else(.$AccidentYear + .$DevelopmentLag > 1998, NA, .$cum_paid)) %>% 
                      as.triangle(origin = "AccidentYear", dev = "DevelopmentLag", value = "cum_paid") %>% 
                      BootChainLadder(R = 1000)),
         boot_res_distr = map(boot, ~extract2(., "IBNR.Totals")),
         boot_inc_pmts = map(boot, ~apply(.$IBNR.Triangles, c(1,2), function(x) mean(x)) %>% 
                               as_tibble() %>% 
                               pivot_longer(cols = everything()) %>% 
                               filter(value != 0)  %>% 
                               pull(value))) %>%
  mutate(results = pmap(list(obs_res, boot_res_distr, boot_inc_pmts),
                        ~tibble("obs_res" = sum(..1),
                                "lower25" = quantile(..2, .025),
                                "est_res" = mean(..2),
                                "upper25" = quantile(..2, .975),
                                "perc" = ecdf(..2)(sum(..1)),
                                "rmse" = rmse(..1, ..3)))) %>% 
  select(key, results) %>% 
  unnest(results)

write_csv(trad_res, "./Full_Run/trad_res.csv")


