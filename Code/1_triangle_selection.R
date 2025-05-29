# Download data -----------------------------------------------------------

download_triangle <- function(link, type = str_extract(link, "[^/_]+(?=_)")){
  read_csv(link) %>% 
    mutate(lob = type) %>% 
    dplyr::select(AccidentYear,
           DevelopmentLag,
           starts_with("CumPaidLoss_"),
           starts_with("EarnedPremNet_"),
           GRCODE,
           lob)  %>% 
    rename(cum_paid = starts_with("CumPaidLoss_"),
           earned_prem = starts_with("EarnedPremNet_")) %>% 
    group_by(GRCODE, lob, AccidentYear) %>% 
    mutate(inc_paid = cum_paid - lag(cum_paid, 1, 0)) %>%
    ungroup()
  
}

pp_auto_link <- "https://www.casact.org/sites/default/files/2021-04/ppauto_pos.csv"
work_comp_link <- "https://www.casact.org/sites/default/files/2021-04/wkcomp_pos.csv"
comm_auto_link <- "https://www.casact.org/sites/default/files/2021-04/comauto_pos.csv"
med_mal_link <- "https://www.casact.org/sites/default/files/2021-04/medmal_pos.csv"
prod_liab_link <- "https://www.casact.org/sites/default/files/2021-04/prodliab_pos.csv"
oth_liab_link <- "https://www.casact.org/sites/default/files/2021-04/othliab_pos.csv"

pp_auto <- download_triangle(pp_auto_link)
work_comp <- download_triangle(work_comp_link)
comm_auto <- download_triangle(comm_auto_link)
med_mal <- download_triangle(med_mal_link)
prod_liab <- download_triangle(prod_liab_link)
oth_liab <- download_triangle(oth_liab_link)

# Data Manipulation -------------------------------------------------------

full_data <- bind_rows(list(pp_auto,
                            work_comp,
                            comm_auto,
                            med_mal,
                            prod_liab,
                            oth_liab))

lob_codes <- full_data %>% 
              distinct(GRCODE, lob)

keep_df <- tibble(lob = character(),
                  GRCODE = numeric())

for (i in 1:nrow(lob_codes)) {
  if (full_data %>%
      filter(lob == lob_codes$lob[i],
             GRCODE == lob_codes$GRCODE[i]) %>%
      group_by(AccidentYear) %>%
      summarise(n_cells = sum(inc_paid > 0)) %$% sum(n_cells) == 100) {
    
    keep_df <- keep_df %>%
      bind_rows(tibble(lob = lob_codes$lob[i],
                       GRCODE = lob_codes$GRCODE[i]))
  }
}

keep_df <- keep_df %>% 
  unite("key", GRCODE:lob)

full_data <- full_data %>% 
             unite("key", GRCODE:lob, remove = TRUE) %>% 
             filter(key %in% keep_df$key) %>% 
             select(AccidentYear, DevelopmentLag, cum_paid, inc_paid, earned_prem, key)

rm(list = c("comm_auto",
            "comm_auto_link",
            "download_triangle",
            "i",
            "keep_df",
            "lob_codes",
            "med_mal",
            "med_mal_link",
            "oth_liab",
            "oth_liab_link",
            "pp_auto",
            "pp_auto_link",
            "prod_liab"        ,
            "prod_liab_link",
            "work_comp",
            "work_comp_link" ))