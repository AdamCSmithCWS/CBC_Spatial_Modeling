plot_variance_spatial_pattern <- function(indices = inds,
                          selected_draws = right_draws,
                          end_year = 2021,
                          start_year = 1966,
                          base_map_blank = base_map,
                          single_map = TRUE,
                          title = ""){
  
  max_year_num <- 1+(end_year-indices$meta_data$start_year)
  min_year_num <- 1+(start_year-indices$meta_data$start_year)
  nyears = max_year_num - min_year_num
  
  trends_plot <- purrr::map_dfr(indices$sample, 
                                .f = ~ 100*(((.x[selected_draws,max_year_num] / .x[selected_draws,min_year_num])^(1/nyears))-1),
                                .id = "strata_name") %>% 
    pivot_longer(.,cols = -strata_name,
                 names_to = ".draw") %>% 
    mutate(strata_name = gsub("stratum_","",
                              strata_name))
  
  if(single_map){
    trends_plot <- trends_plot %>% 
      group_by(strata_name) %>% 
      summarise(value = mean(value)) %>% 
      mutate(t_plot = cut(value, breaks = c(-Inf, breaks, Inf),
                          labels = labls))
  }else{
    trends_plot <- trends_plot %>% 
      mutate(t_plot = cut(value, breaks = c(-Inf, breaks, Inf),
                          labels = labls))
  }
  
  pal <- setNames(
    c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
      "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"),
    levels(trends_plot$t_plot))
  
  
  trend_map <- base_map_blank %>% 
    right_join(.,trends_plot,
               by = "strata_name",
               multiple = "all")
  
  plot_out <- ggplot()+
    geom_sf(data = trend_map,
            aes(fill = t_plot))+
    scale_fill_manual(values = pal)+
    labs(title = title)+
    theme_minimal()
  if(!single_map){
    plot_out <- plot_out +
      facet_wrap(vars(.draw))
  }
  
  return(plot_out)
  
}
