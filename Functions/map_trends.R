map_trends <- function(trends = trends,
                       variable = "trend",
                       plot_trend = TRUE,
                          base_map_blank = base_map,
                          title = "",
                       legend_title = NULL,
                       zoom_out = 0.05,
                       add_base = TRUE,
                       region_name = "region",
                       facetgrid = FALSE,
                       facet_1 = NULL,
                       facet_2 = NULL,
                       ...){
  
  zoom_out <- zoom_out+1
  
  if(!is_tibble(trends)){
    trends <- trends[["trends"]] %>% 
      filter(region_type == "stratum") 
  }

  if(facetgrid){
    if(!is.null(facet_2) &
       !is.null(facet_1)){
    trends <- trends %>% 
      rename_with(.,~gsub(facet_1,
                          "fac1",
                          .x)) %>% 
      rename_with(.,~gsub(facet_2,
                          "fac2",
                          .x))
    }else{
    if(!is.null(facet_1)){
      trends <- trends %>% 
        rename_with(.,~gsub(facet_1,
                            "fac1",
                            .x))  
    }
      if(!is.null(facet_2)){
        trends <- trends %>% 
          rename_with(.,~gsub(facet_2,
                              "fac2",
                              .x))  
      }
    }
  }
 
  
if(plot_trend){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls <- c(paste0("< ", breaks[1]),
             paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),
             paste0("> ",breaks[length(breaks)]))
  labls <- paste0(labls, " %")
  
  trends_plot <- trends %>% 
    rename_with(.,~gsub(variable,"value",.x)) %>% 
    mutate(t_plot = cut(value, breaks = c(-Inf, breaks, Inf),
                        labels = labls,
                        ordered_result = TRUE))
  
  pal <- setNames(
    c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
      "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"),
    levels(trends_plot$t_plot))
  
  start_year <- min(trends_plot$start_year)
  end_year <- max(trends_plot$end_year)
  
  
  if(!is.null(legend_title)){
    l_title <- paste(legend_title,"\n",start_year,"-",end_year)
    
  }else{
  l_title <- paste("Trend \n",start_year,"-",end_year)
  }
  if(length(unique(trends_plot$start_year)) > 1){
    l_title <- "Trend"
  }
  }else{
    

    trends_plot <- trends %>% 
      rename_with(.,~gsub(variable,"t_plot",.x)) 
    
    pal <- scale_colour_viridis_c()
    start_year <- min(trends_plot$start_year)
    end_year <- min(trends_plot$end_year)
    
    if(!is.null(legend_title)){
      l_title <- paste(legend_title,"\n",start_year,"-",end_year)
      
    }else{
    l_title <- paste(variable,"\n",start_year,"-",end_year)
    }
    if(length(unique(trends_plot$start_year)) > 1){
      if(!is.null(legend_title)){
        l_title <- paste(legend_title)
        
      }else{
        l_title <- paste(variable)
      }
    }
  }
  
  
  
  

  
  trend_map <- base_map_blank %>% 
    right_join(.,trends_plot,
               by = c("strata_name" = region_name),
               multiple = "all")
  
  map_ext <- sf::st_bbox(trend_map)
  
  
  plot_out <- ggplot()+
    geom_sf(data = trend_map,
            aes(fill = t_plot))+
    labs(title = title)+
    theme_minimal()
  
  if(plot_trend){
    plot_out <- plot_out +
    scale_fill_manual(values = pal,
                      guide = guide_legend(title = l_title, reverse = TRUE))
  }else{
    plot_out <- plot_out +
      scale_fill_viridis_c(guide = guide_colourbar(title = l_title))
  }
  if(add_base){
    plot_out <-  plot_out+
      geom_sf(data = base_map_blank,
              alpha = 0)
  }
  
  if(facetgrid){
    if(is.null(facet_2)){
  plot_out<- plot_out +
    facet_grid(rows = vars(fac1))
    }
    if(is.null(facet_1)){
      plot_out<- plot_out +
        facet_grid(cols = vars(fac2))
    }
    if(!is.null(facet_2) &
       !is.null(facet_1)){
      plot_out<- plot_out +
        facet_grid(rows = vars(fac1),
                   cols = vars(fac2))
    }
  }
  
  plot_out<- plot_out +
    coord_sf(xlim = map_ext[c("xmin","xmax")]*zoom_out,
             ylim = map_ext[c("ymin","ymax")]*zoom_out)+
    theme(plot.margin = unit(rep(1,4),"mm"),
          panel.spacing = unit(2,"mm"),
          axis.text = element_text(size = 6),
          title = element_text(size = 8))
  
  return(plot_out)
  
}
