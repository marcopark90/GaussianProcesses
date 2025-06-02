library(tidyverse)

funs <- file.path("./Stan/Functions") %>% 
          list.files(full.names = T, pattern = "*.stan") %>% 
          map(~tibble(fun = read_file(.),
                    name = str_extract(., "(?<=Functions/).+(?=\\.stan)"))) %>% 
          list_rbind()

pars <- file.path("./Stan/Parameters") %>% 
          list.files(full.names = T) %>% 
          map(~tibble(pars = read_file(.) %>% str_trim(),
                    V = str_extract(., "_v[:digit:]+"))) %>% 
          list_rbind()

df <- expand_grid(funs, pars) %>% 
        mutate(fun_1 = str_split_i(fun, "// Parameters", i = 1),
               fun_2 = str_split_i(fun, "// Parameters", i = 2),
               cov_fun = str_c(fun_1, pars, fun_2),
               file_name = str_c(name, V)) %>% 
        select(cov_fun, file_name)

map2(df$cov_fun, df$file_name, ~write_file(.x, file = paste0("./Stan/Output/", .y,".stan")))