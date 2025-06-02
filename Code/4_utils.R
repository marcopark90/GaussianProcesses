plot_res <- function(data_list, title = "Predictive Reserve Distribution"){
  data_list %>% 
  pluck("distr_res") %>% 
  as_tibble() %>% 
  ggplot() +
  geom_histogram(aes(x = value, 
                     y = after_stat(density)), 
                 col = "blue", 
                 fill = "deepskyblue4", 
                 alpha = .6,
                 bins = 50) +
  scale_x_continuous('Reserve values') +
  scale_y_continuous('Frequency') + 
  theme(plot.title = element_text(face = 'bold', size = 10, hjust = .5),
        axis.title.x = element_text(face = "bold", colour = "black", size = 10),
        axis.title.y = element_text(face = "bold", colour = "black", size = 10)) +
  geom_vline(aes(xintercept = data_list$est_mean, color = 'Estimated'), linetype = 1, linewidth = 2) +
  geom_vline(aes(xintercept = data_list$obs_res, color = 'Observed'), linetype = 8, linewidth = 2) +  
  scale_color_manual(name = "", values = c(Estimated = "blue", Observed = "red")) +
  ggtitle(title)
}
