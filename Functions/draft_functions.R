


# index function ---------------------------------------------------------------
index_function <- function(fit = stanfit,
                           parameter = "n",
                           strat = NULL, # or strat = "Stratum_Factored",
                           year = "Year",
                           first_dim = "s",
                           quant = 0.95,
                           weights_df = NULL, # or weights_df = cntry_wts, if supplied a composite summary estimate is generated
                           area = NULL, # or "Area",
                           summary_regions = NULL, # or summary_regions = "country", optional column to provide summary regions, if not supplied summarises across all strata if weights_df is upplied
                           year_1 = 1966,
                           to_summarise = FALSE # or to_summarise = TRUE, if TRUE provides summary indices across regions, otherwise jsut strata estimates
){
  
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)
  year_1 <- year_1
  if(is.null(strat)){
    dims <- year
    
  }else{
    if(!is.null(weights_df)){
      if(is.null(strat)){
        stop("Weights were supplied but no stratum dimension")
        
        return(NULL)
      }
      if(!any(grepl(pattern = strat,x = names(weights_df)))){
        stop("weights_df must include a column with a name that matches strat")
        return(NULL)
      }
      if(!to_summarise){
        to_summarise <- TRUE
        warning("weights were supplied by to_summarise is FALSE, changing to_summaries to TRUE")
      }
    }
    if(to_summarise){
      if(is.null(weights_df) | is.null(strat)){
        stop("to_summarise is TRUE but no weights were supplied or strat is NULL")
        return(NULL)
      }
    }
    if(tolower(substr(first_dim,1,1)) == "s"){
      dims <- c(strat,year)
    }
    if(tolower(substr(first_dim,1,1)) == "y"){
      dims <- c(year,strat)
    }
  }
  
  smpls <- posterior_samples(fit = fit,
                             parm = parameter,
                             dims = dims)
  
  if(!is.null(weights_df)){
    nstrata_fit <- length(unique(smpls[[strat]]))
    nstrata_w <- nrow(weights_df)
    if(nstrata_fit != nstrata_w){
      stop("Lengths of strata and weights are different")
      return(NULL)
    }
    
    if(!is.null(summary_regions)){
      
      weights_df <- weights_df %>% 
        rename_with(.,~gsub(pattern = area,
                            replacement = "stratum_area",
                            x = .x,
                            fixed = TRUE)) %>% 
        rename_with(., ~gsub(pattern = summary_regions,
                             replacement = "summary_region",.x,
                             fixed = TRUE))
      tmp <- weights_df %>% 
        group_by(summary_region) %>% 
        summarise(regional_area_sum = sum(stratum_area),.groups = "keep")
      
      weights_df <- weights_df %>% 
        left_join(.,tmp,by = "summary_region") %>% 
        mutate(stratum_weight = stratum_area/regional_area_sum)
      
    }else{
      weights_df <- weights_df %>% 
        rename_with(.,~gsub(pattern = area,
                            replacement = "stratum_area",
                            x = .x,
                            fixed = TRUE)) %>% 
        mutate(stratum_weight = stratum_area/sum(stratum_area),
               summary_region = "Survey_wide") 
      summary_regions <- "Survey_wide"
      
    }
    
    smpls <- smpls %>% 
      left_join(.,weights_df,
                by = strat)
    
    # if(summary_regions != "Survey_wide"){
    inds <- smpls %>% 
      mutate(.value = .value*stratum_weight) %>% 
      rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                           fixed = TRUE)) %>% 
      group_by(yyy,summary_region,.draw) %>% 
      summarise(.vsum = sum(.value),
                .groups = "keep") %>% 
      group_by(yyy,summary_region) %>% 
      summarise(mean = mean(.vsum),
                median = median(.vsum),
                lci = quantile(.vsum,lu),
                uci = quantile(.vsum,uu),
                .groups = "keep") %>% 
      mutate(true_year = yyy+year_1-1) %>% 
      rename_with(., ~gsub(replacement = summary_regions,pattern = "summary_region",.x,
                           fixed = TRUE)) %>% 
      rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                           fixed = TRUE))
    
    smpls <- smpls %>% 
      mutate(.value = .value*stratum_weight) %>% 
      rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                           fixed = TRUE)) %>% 
      group_by(yyy,summary_region,.draw) %>% 
      summarise(.value = sum(.value),
                .groups = "keep") %>%
      rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                           fixed = TRUE)) %>% 
      rename_with(., ~gsub(replacement = summary_regions,
                           pattern = "summary_region",.x,
                           fixed = TRUE)) %>% 
      mutate(true_year = yyy+year_1-1) %>% 
      rename_with(., ~gsub(pattern = "yyy",replacement = year,.x,
                           fixed = TRUE))
    
    weights_df <- weights_df %>% 
      rename_with(., ~gsub(replacement = area,
                           pattern = "stratum_area",.x,
                           fixed = TRUE))%>% 
      rename_with(., ~gsub(replacement = summary_regions,
                           pattern = "summary_region",.x,
                           fixed = TRUE))
  }else{ ## everything up to here requires a weights_df data frame
    if(!is.null(strat)){
      inds <- smpls %>% 
        rename_with(., ~gsub(pattern = strat,replacement = "sss",.x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>% 
        group_by(sss,yyy) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  lci = quantile(.value,lu),
                  uci = quantile(.value,uu),
                  .groups = "keep") %>%
        mutate(true_year = yyy+year_1-1) %>% 
        rename_with(., ~gsub(replacement = strat,pattern = "sss",.x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                             fixed = TRUE))
      summary_regions <- strat
      
      smpls <- smpls %>% 
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>% 
        rename_with(., ~gsub(replacement = strat,pattern = "sss",.x,
                             fixed = TRUE)) %>% 
        mutate(true_year = yyy+year_1-1) %>% 
        rename_with(., ~gsub(pattern = "yyy",replacement = year,.x,
                             fixed = TRUE))
      
    }else{
      summary_regions <- "Complete_pooling"
      
      inds <- smpls %>%
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>%
        group_by(yyy) %>%
        summarise(mean = mean(.value),
                  median = median(.value),
                  lci = quantile(.value,lu),
                  uci = quantile(.value,uu),
                  .groups = "keep") %>%
        mutate(true_year = yyy+year_1,
               Survey_wide = "Survey_wide") %>%
        rename_with(., ~gsub(replacement = year,pattern = "yyy",.x,
                             fixed = TRUE))
      
      
      smpls <- smpls %>%
        rename_with(., ~gsub(pattern = year,replacement = "yyy",.x,
                             fixed = TRUE)) %>%
        mutate(true_year = yyy+year_1-1,
               Survey_wide = "Survey_wide") %>%
        rename_with(., ~gsub(pattern = "yyy",replacement = year,.x,
                             fixed = TRUE))
      
    }
  }
  
  inds <- inds %>% 
    mutate(region_type = summary_regions)

  return(list(indices = inds,
              samples = smpls,
              parameter = parameter,
              strat = strat,
              year = year,
              dims = dims,
              quant = quant,
              weights_df = weights_df,
              area = area,
              summary_regions = summary_regions,
              to_summarise = to_summarise,
              year_1 = year_1
  ))
}
# ------------------------------------------------------------------------------


# add smooths to indices -------------------------------------------------------
get_smooth_ns <- function(idx_list){
  
  # establish quantiles
  lu <- ((1-(as.numeric(idx_list$quant)))/2)
  uu <- 1-((1-(as.numeric(idx_list$quant)))/2)
  
  # get n samples and smooth n
  dat_i <- idx_list$samples %>% ungroup() %>% 
    group_by(across(c(2, 3))) %>% 
    mutate(.value = smooth_sample(x=Year, y=.value))

  # make summaries of smooth n
  sum_smooth <- dat_i %>% ungroup() %>% 
    group_by(across(c(1, 2, 5))) %>% 
    summarise(mean = mean(.value),
              median = median(.value),
              lci = quantile(.value, lu),
              uci = quantile(.value, uu)) %>% 
    select(1, 2, 4, 5, 6, 7, 3) %>% 
    mutate(region_type=idx_list$summary_regions)
  
  # return list with two dataframes
  return(list(samples_smooth=dat_i, indices_smooth=sum_smooth))
}
# ------------------------------------------------------------------------------





# new map function -------------------------------------------------------------
map_function <- function(trds = bcr_10yr_trds, spunit = "bcr"){
  breaks1 <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls1 <- c(paste0("< ", breaks1[1]),
              paste0(breaks1[-c(length(breaks1))],":", breaks1[-c(1)]),
              paste0("> ",breaks1[length(breaks1)]))
  labls1 <- paste0(labls1, "%")
  cols1 <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", 
             "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", 
             "#313695")
  if(spunit == "bcr"){
    map1 <- bbsBayes2::load_map(stratify_by = spunit, type = "strata") %>% 
      mutate(stratum = as.numeric(str_extract(strata_name, pattern = "(\\d)+")))
  }
  if(spunit == "prov_state"){
    map1 <- bbsBayes2::load_map(stratify_by = spunit, type = "strata") %>% 
      mutate(stratum = strata_name)
  }
  map_ext1 <- sf::st_bbox(right_join(map1, trds)) * 1.2
  trds_map1 <- left_join(map1, trds, by = "stratum") %>% 
    mutate(t_plot = cut(trend, breaks = c(-Inf, breaks1, Inf),
                        labels = labls1, ordered_result = TRUE)) %>% 
    mutate(t_plot = factor(t_plot, levels = labls1)) %>% 
    mutate(sig = as.character(sig)) %>% 
    mutate(sig = ifelse(is.na(sig), "0", sig))
  pal1 <- setNames(cols1, levels(trds_map1$t_plot))
  title1 <- paste0(species, ": ", 
                   unique(na.omit(trds_map1$start_year)), " through ", 
                   unique(na.omit(trds_map1$end_year)))
  tp <- ggplot() + 
    geom_sf(data = map1, fill = "grey80") + 
    geom_sf(data = trds_map1 %>% filter(!is.na(t_plot)), 
            aes(fill = t_plot), show.legend = TRUE) + 
    geom_sf_pattern(data = trds_map1,
                    aes(fill = t_plot, pattern_type = sig),
                    pattern = 'magick',
                    pattern_scale = 0.5,
                    pattern_fill = "black",
                    show.legend = F) +
    scale_pattern_type_discrete(choices = c("gray100", "left45")) +
    scale_fill_manual(values = pal1,
                      drop = F,
                      na.translate = F,
                      na.value = "grey80",
                      guide = guide_legend(title = "Trend (%/yr)", 
                                           reverse = TRUE)) +
    coord_sf(xlim = map_ext1[c("xmin", "xmax")],
             ylim = map_ext1[c("ymin", "ymax")]) +
    labs(title = title1)+
    theme_bw() +
    guides(pattern_type = "none") + 
    theme(plot.margin = unit(rep(1, 4),"mm"),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 8),
          title = element_text(size = 11),
          plot.title = element_text(hjust = 0.5))
  return(tp)
}
# ------------------------------------------------------------------------------





# NEW new trends function based on gam endpoints -----------------------------------
test <- trends_function(ind_list = ps_idx_lst,
                            start_year = 1966,
                            end_year = 2019,
                            quant = 0.95)

trends_function <- function(ind_list = ps_idx_lst,
                            start_year = 1966,
                            end_year = 2019,
                            quant = 0.95){
  
  samples <- ind_list$samples
  samples_smooth <- ind_list$samples_smooth
  #strat <- ind_list$strat
  year <- ind_list$year
  #dims <- ind_list$dims
  #weights_df <- ind_list$weights_df
  #area <- ind_list$area
  summary_regions <- ind_list$summary_regions
  #to_summarise <- ind_list$to_summarise
  
  if(is.null(end_year)){
    end_year <- max(samples$true_year)
  }
  if(is.null(start_year)){
    start_year <- min(samples$true_year)
  }
  
  nyrs <- end_year-start_year
  lu <- ((1-(quant))/2)
  uu <- 1-((1-(quant))/2)
  
  # get samples
  indt <- samples %>% 
    ungroup() %>% 
    filter(true_year %in% c(start_year, end_year)) %>% 
    # ungroup() %>% 
    select(-matches(match = year, ignore.case = FALSE)) %>% 
    pivot_wider(names_from = true_year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",start_year),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",end_year),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "stratum_trend",
                         pattern = summary_regions,.x,
                         fixed = TRUE))
  
  # compute endpoint trends
  tt <- indt %>% 
    group_by(.draw, stratum_trend) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end/start,ny = nyrs),
              ch = chng(end/start),
              .groups = "keep") %>% 
    group_by(stratum_trend) %>% 
    summarise(trend = mean(t),
              lci = quantile(t,lu,names = FALSE),
              uci = quantile(t,uu,names = FALSE),
              percent_change = median(ch),
              p_ch_lci = quantile(ch,lu,names = FALSE),
              p_ch_uci = quantile(ch,uu,names = FALSE),
              prob_decline = prob_dec(ch,0),
              prob_decline_GT30 = prob_dec(ch,-30),
              prob_decline_GT50 = prob_dec(ch,-50),
              prob_decline_GT70 = prob_dec(ch,-70))%>% 
    rename_with(., ~gsub(replacement = summary_regions,
                         pattern = "stratum_trend",.x,
                         fixed = TRUE)) %>% 
    mutate(region_type = summary_regions)
  
  # save
  out_lst <- list(endpoint_trends=tt)
  
  # get smooth samples
  indt <- samples_smooth %>%
    ungroup() %>% 
    filter(true_year %in% c(start_year, end_year)) %>% 
    select(-matches(match = year, ignore.case = FALSE)) %>% 
    pivot_wider(names_from = true_year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",start_year),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",end_year),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "stratum_trend",
                         pattern = summary_regions,.x,
                         fixed = TRUE))
  
  # compute smooth endpoint trends
  tt <- indt %>% 
    group_by(.draw, stratum_trend) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end/start,ny = nyrs),
              ch = chng(end/start),
              .groups = "keep") %>% 
    group_by(stratum_trend) %>% 
    summarise(trend = mean(t),
              lci = quantile(t,lu,names = FALSE),
              uci = quantile(t,uu,names = FALSE),
              percent_change = median(ch),
              p_ch_lci = quantile(ch,lu,names = FALSE),
              p_ch_uci = quantile(ch,uu,names = FALSE),
              prob_decline = prob_dec(ch,0),
              prob_decline_GT30 = prob_dec(ch,-30),
              prob_decline_GT50 = prob_dec(ch,-50),
              prob_decline_GT70 = prob_dec(ch,-70))%>% 
    rename_with(., ~gsub(replacement = summary_regions,
                         pattern = "stratum_trend",.x,
                         fixed = TRUE)) %>% 
    mutate(region_type = summary_regions)
  
  # save
  out_lst$smooth_endpoint_trends <- tt
  
  # return
  return(out_lst)
}
# ------------------------------------------------------------------------------









# # new trends function based on gam endpoints -----------------------------------
# trends_function <- function(ind_list,
#                             start_year,
#                             end_year,
#                             quant = 0.95){
#   # get info
#   indices <- ind_list$indices
#   samples <- ind_list$samples
#   parameter <- ind_list$parameter
#   strat <- ind_list$strat
#   year <- ind_list$year
#   dims <- ind_list$dims
#   weights_df <- ind_list$weights_df
#   area <- ind_list$area
#   summary_regions <- ind_list$summary_regions
#   to_summarise <- ind_list$to_summarise
#   
#   # fill in star and end years
#   if(is.null(end_year)){
#     end_year <- max(samples$true_year)
#   }
#   if(is.null(start_year)){
#     start_year <- min(samples$true_year)
#   }
#   
#   # calc quantiles
#   lu <- ((1 - (quant)) / 2)
#   uu <- 1 - ((1 - (quant)) / 2)
#   nyrs <- end_year - start_year + 1 # have to add 1 right?????????????
# 
#   # get samples
#   indt <- samples %>% 
#     ungroup() %>% 
#     rename_with(., ~gsub(replacement = "stratum_trend",
#                        pattern = summary_regions,.x,
#                        fixed = TRUE))
#   
#   # test data
#   # indt <- indt %>% filter(.draw<100)
#   
#   # generate trends
#   tt <- indt %>% 
#     group_by(.draw, stratum_trend) %>% 
#     summarise(
#       t = tgam(y = .value, x = true_year, t1 = start_year, t2 = end_year),
#       ch = ((((t / 100) + 1)^nyrs) - 1) * 100,
#       .groups = "keep") %>% 
#     group_by(stratum_trend) %>% 
#     summarise(trend = mean(t),
#             lci = quantile(t, lu, names = FALSE),
#             uci = quantile(t, uu, names = FALSE),
#             percent_change = median(ch),
#             p_ch_lci = quantile(ch, lu, names = FALSE),
#             p_ch_uci = quantile(ch, uu, names = FALSE),
#             prob_decline = prob_dec(ch, 0),
#             prob_decline_GT30 = prob_dec(ch, -30),
#             prob_decline_GT50 = prob_dec(ch, -50),
#             prob_decline_GT70 = prob_dec(ch, -70)) %>% 
#     rename_with(., ~gsub(replacement = summary_regions,
#                        pattern = "stratum_trend", .x,
#                        fixed = TRUE))
#   
#   # add metadata
#   tt <- tt %>% 
#     mutate(region_type = summary_regions,
#            start_year = start_year, end_year = end_year) 
#   return(tt)
# }
# # ------------------------------------------------------------------------------




# helper functions --------------------------------------------------------
smooth_sample <- function(y, x){
  require(mgcv)
  ny <- length(x) # full time series
  kts <- as.numeric(round(ny/10, 0)) # knot every five years
  dat <- data.frame(logy = log(y), x = x)
  yhat <- exp(fitted(gam(logy ~ 1 + s(x, k=kts, bs="ps"),
                         family = "gaussian", data = dat)))
  return(yhat)
}


# gam trend function. p-spline with nyears/4 knots. gaussian on log Y
# this function is used in the trends_function()
# tgam <- function(y, x, t1, t2){
#   require(mgcv)
#   ny <- length(x) # full time series
#   kts <- as.numeric(round(ny/4, 0)) # knot every four years
#   df <- data.frame(y = y, x = x)
#   # predict from p-spline with ny/4 knots, gaussian on log Y
#   yhat <- exp(fitted(gam(y ~ 1 + s(x),
#                          family = "gaussian", data = df)))
#   nt1 <- yhat[which(x == t1)] # pick start year
#   nt2 <- yhat[which(x == t2)] # pick end year
#   trd <- 100 * (((nt2 / nt1)^(1 / (t2 - t1))) - 1) # calc annl percent chg
#   return(trd)
# }
#

p_lt <- function(x,th){
  length(which(x < th))/length(x)
}

p_neg <- function(x){
  length(which(x < 0))/length(x)
}

texp <- function(x,ny = 2019-1974){
  (x^(1/ny)-1)*100
}

chng <- function(x){
  (x-1)*100
}

# this function is used in the trends_function()
prob_dec <- function(ch,thresh){
  length(which(ch < thresh))/length(ch)
}

### function to extract the dimension values from an bayesian model fit
### works within the gather_samples function
# this function is used in posterior_samples()
dim_ext <- function(dim = 1,
                    var = "",
                    cl = "Parameter",
                    dat = NULL){
  ##3 function to extract the indicator values from cmdstanr output
  require(stringr)
  
  pat = paste0("(?<=",var,"\\[")
  
  if(dim > 1){
    for(j in 1:(dim-1)){
      
      pat2 = paste0(pat,")[:digit:]+")
      cl2 = str_extract(unlist(dat[,cl]),pattern = pat2)
      
      d = max(nchar(cl2))
      
      pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
    }
  }
  
  
  pat = paste0(pat,")[:digit:]+")
  dds = as.integer(str_extract(unlist(dat[,cl]),pattern = pat))
  return(dds)
  
}

### function to generate the same tidy output as gather-draws in tidybayes package
## dims should be a character vector defining the dimensions of the parameter
## e.g., parm = "nsmooth", dims = c("stratum","year"),
## function works with cmdstanr output by default and rstan fits
## with is_rstan == TRUE
# this function is used in index_function()
posterior_samples <- function(fit = cmdstanfit,
                              parm = "nsmooth",
                              dims = NULL,
                              is_rstan = FALSE,
                              is_mcmc = FALSE){
  require(posterior)
  require(tidyverse)
  
  if(length(dims) > 0){
    parm_ex <- paste0(parm,"[")
  }else{
    parm_ex <- parm
  }
  if(class(fit)[1] == "stanfit"){is_rstan <- TRUE}
  
  if(class(fit)[1] == "mcmc"){is_mcmc <- TRUE}
  
  if(is_rstan | is_mcmc){
    samples <- as_draws_df(as.array(fit)) %>% 
      dplyr::select(starts_with(parm_ex,ignore.case = FALSE),
                    .chain,
                    .iteration,
                    .draw)#,pars = c(parm)))
  }else{
    
    samples <- as_draws_df(fit$draws(variables = c(parm)))
    
  }

  plong <- suppressWarnings(samples %>% pivot_longer(
    cols = starts_with(parm_ex,ignore.case = FALSE),
    names_to = c(".variable"),
    values_to = ".value",
    values_drop_na = TRUE
  )) 
  
  for(dn in 1:length(dims)){
    dd = dims[dn]
    plong[,dd] = dim_ext(dim = dn,
                         var = parm,
                         cl = ".variable",
                         dat = plong)
    
  }
  
  plong <- plong %>% mutate(.variable = parm)
  return(plong)
  
}

# posterior summaries function
posterior_sums <- function(samples = n_samples,
                           quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                           ci = 0.95,
                           dims = NULL){
  
  cis = c((1-ci)/2,1-((1-ci)/2))
  
  if(!is.null(dims)){
    for(i in 1:length(dims)){
      samples <- samples %>% 
        rename_with(~gsub(dims[i],paste0("d",i),.x,fixed = TRUE))
    }
    if(length(dims) == 1){
      sums = samples %>% 
        group_by(d1) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) 
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
    }
    if(length(dims) == 2){
      sums = samples %>% 
        group_by(d1,d2) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))
      
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1,d2) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))
          
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
      
    }
    
    if(length(dims) == 3){
      sums = samples %>% 
        group_by(d1,d2,d3) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d2",dims[2],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))
      
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1,d2,d3) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))
          
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
      
    }
    
    
    if(length(dims) == 4){
      sums = samples %>% 
        group_by(d1,d2,d3,d4) %>% 
        summarise(mean = mean(.value),
                  median = median(.value),
                  sd = sd(.value),
                  lci = quantile(.value,cis[1]),
                  uci = quantile(.value,cis[2]),
                  .groups = "keep") %>% 
        rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d2",dims[2],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d3",dims[3],.x,fixed = TRUE)) %>% 
        rename_with(~gsub("d4",dims[4],.x,fixed = TRUE))
      
      
      if(!is.null(quantiles)){
        qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
        for(i in 1:length(quantiles)){
          qq = quantiles[i]
          qn = qs[i]
          
          sumt = samples %>% 
            group_by(d1,d2,d3,d4) %>% 
            summarise(tt = as.numeric(quantile(.value,qq)),
                      .groups = "keep") %>% 
            rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
            rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))  %>% 
            rename_with(~gsub("d4",dims[4],.x,fixed = TRUE))
          
          
          sums = left_join(sums,sumt,by = dims)
        }
        
      }
      
      
    }
    
    if(length(dims) > 4){stop("ERROR this function cannot handle more than 4 dimensions, but it could be easily modified")}
    
  }else{
    
    
    sums = samples %>% 
      summarise(mean = mean(.value),
                median = median(.value),
                sd = sd(.value),
                lci = quantile(.value,cis[1]),
                uci = quantile(.value,cis[2]),
                .groups = "keep")
    
    if(!is.null(quantiles)){
      qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
      for(i in 1:length(quantiles)){
        qq = quantiles[i]
        qn = qs[i]
        
        sumt = samples %>% 
          summarise(tt = as.numeric(quantile(.value,qq)),
                    .groups = "keep") %>% 
          rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE))
        sums = bind_cols(sums,sumt)
      }
      
    }
    
  }
  return(sums)
  
}
# ------------------------------------------------------------------------------

