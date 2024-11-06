# setup ------------------------------------------------------------------------
library(tidyverse)
library(bbsBayes2)
library(cmdstanr)
library(patchwork)
library(posterior)
library(bayesplot)
library(sf)
library(geofacet)
library(ggrepel)
library(HDInterval)

HDL <- function(x,int,upper = TRUE){
  b <- HDInterval::hdi(x,int)
  return(ifelse(upper,b[["upper"]],b[["lower"]]))
}

# some some helper functions
source("functions/neighbours_define.R")
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("functions/map_trends.R")

species <- "American Dipper"
stratification <- "bbs_usgs"
models <- model <- "first_diff"
model_variants <- model_variant <- "spatial"
data_sets <- data_set <- "cbc"
setwd("C:/Users/tmeehan/Documents/GitHub/CBC_Spatial_Modeling")
# ------------------------------------------------------------------------------



# get fit object and look at summaries -----------------------------------------
# uses data_prep and fit from data prep and model fit script
# bring in model fit object
fit <- readRDS("output/fit_CBC_spatial_first_diff.rds")

# make a stratum df
strat_df <- data_prep %>% 
  select(strata_name,strata_vec,non_zero,area_sq_km) %>% 
  rename(stratum = strata_vec) %>% 
  distinct() %>% 
  mutate(area_weight = area_sq_km/sum(area_sq_km))

# make a vector of years
yrs <- data_prep %>% 
  select(year_vec,count_year) %>% 
  rename(yr = year_vec,
         year = count_year) %>% 
  distinct()

# make a stratum map
map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
  select(-area_sq_km) %>% 
  inner_join(.,strat_df,
             by = "strata_name") %>% 
  arrange(stratum)
# ------------------------------------------------------------------------------



# stratum level indices and trends ---------------------------------------------
# posterior samples of stratum indices
ind_samples <- posterior_samples(fit = fit, parm = "n",
                                 dims = c("stratum","yr")) %>% 
  inner_join(., strat_df, by = "stratum") %>% 
  inner_join(., yrs, by = "yr")

# stratum level index summaries
inds_strat <- 
  ind_samples %>%
  group_by(strata_name,year) %>% 
  summarise(index = median(.value),
            index_q_0.05 = HDL(.value,0.9,upper = FALSE),
            index_q_0.95 = HDL(.value,0.9,upper = TRUE),
            .groups = "keep") %>% 
  mutate(index_type = "full")

# use stratum indices to compute stratum trends of different length
first_years <- c(1966, 1970, 1993, 2009)
trends_strata <- NULL
for(j in 1:length(first_years)){
  ys <- first_years[j]
  ye <- 2019
  nyrs <- ye-ys
  trend_tmp <- ind_samples %>% 
    filter(year %in% c(ys,ye)) %>% 
    ungroup() %>% 
    select(-matches(match = "yr",ignore.case = FALSE)) %>% 
    pivot_wider(names_from = year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",ys),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",ye),.x,
                         fixed = TRUE))%>% 
    group_by(.draw,strata_name) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end/start,ny = nyrs),
              ch = chng(end/start),
              .groups = "keep") %>% 
    group_by(strata_name) %>% 
    summarise(trend = mean(t),
              lci = quantile(t,0.025,names = FALSE),
              uci = quantile(t,0.975,names = FALSE),
              width_CI = uci-lci) %>% 
    mutate(model = model,
           model_variant = model_variant,
           start_year = ys,
           end_year = ye,
           trend_length = nyrs,
           data_set = data_set)
  trends_strata <- bind_rows(trends_strata, trend_tmp)
}
# ------------------------------------------------------------------------------



# continental level indices and trends -----------------------------------------
# make continent scaled index samples by summing across strata
comp_samples <- 
  ind_samples %>%
  mutate(.value = .value * area_weight) %>% 
  group_by(.draw, year) %>% #draw-wise summary of area-weighted indices
  summarise(.value = sum(.value),
            .groups = "drop")

# continent scaled index summaries
# inds and smooth are the same because its first difference
inds_comp <- inds_comp_smooth <-
  comp_samples %>%  
  group_by(year) %>% 
  summarise(index = median(.value),
            index_q_0.05 = HDL(.value,0.9,upper = FALSE),
            index_q_0.95 = HDL(.value,0.9,upper = TRUE),
            .groups = "keep") %>% 
  mutate(strata_name = "continent") %>% 
  mutate(index_type = "full")

# combine stratum and continental index summaries for later
inds_all <- bind_rows(inds_comp,
                      inds_strat) %>% 
  mutate(model = "first_difference",
         model_variant = "spatial",
         data_set = "cbc")

# use indices to make continent scaled trends of different lengths
trends_comp <- NULL
for(j in 1:length(first_years)){
  ys <- first_years[j]
  ye <- 2019
  nyrs <- ye-ys
  trend_tmp <- comp_samples %>% 
    filter(year %in% c(ys,ye)) %>% 
    ungroup() %>% 
    select(-matches(match = "yr",ignore.case = FALSE)) %>% 
    pivot_wider(names_from = year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(., ~gsub(replacement = "start",
                         pattern = paste0("Y",ys),.x,
                         fixed = TRUE))%>% 
    rename_with(., ~gsub(replacement = "end",
                         pattern = paste0("Y",ye),.x,
                         fixed = TRUE))%>% 
    group_by(.draw) %>% 
    summarise(end = sum(end),
              start = sum(start),
              t = texp(end/start,ny = nyrs),
              ch = chng(end/start),
              .groups = "drop") %>% 
    summarise(trend = mean(t),
              lci = quantile(t,0.025,names = FALSE),
              uci = quantile(t,0.975,names = FALSE),
              width_CI = uci-lci) %>% 
    mutate(model = model,
           model_variant = model_variant,
           start_year = ys,
           end_year = ye,
           trend_length = nyrs,
           data_set = data_set,
           strata_name = "continent")
  trends_comp <- bind_rows(trends_comp, trend_tmp)
}
# ------------------------------------------------------------------------------



# index plots ------------------------------------------------------------------
# continent
tmp <- inds_all %>% 
  filter(strata_name == "continent")
tp <- ggplot(data = tmp,
            aes(x = year, y = index)) +
  geom_ribbon(aes(ymin = index_q_0.05,
                  ymax = index_q_0.95),
              alpha = 0.2, fill="darkblue") +
  geom_line(color="darkblue") +
  theme_bw(); tp

# strata
obs_mean <- data_prep %>% 
  select(strata_name, year=count_year, how_many) %>% 
  group_by(strata_name, year) %>% 
  summarise(index = mean(how_many))
tmp <- inds_all %>% 
  filter(strata_name != "continent",
         data_set == "cbc",
         index_type == "full")
tp <- ggplot() +
  geom_ribbon(data = tmp, 
              aes(x = year, ymin = index_q_0.05, ymax = index_q_0.95),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = tmp, 
            aes(x = year, y = index), color="darkblue") +
  geom_point(data = obs_mean, 
             aes(x = year, y = index),
             alpha=0.3, size=2, shape=1) +
  facet_wrap(~ strata_name, scales = "free") +
  theme_bw(); tp
# ------------------------------------------------------------------------------



# map trends -------------------------------------------------------------------
trends_out <- bind_rows(trends_strata, trends_comp)
yrpairs <- trends_out %>% 
  filter(data_set == as.character("cbc")) %>% 
  select(start_year, end_year) %>% 
  distinct() %>% 
  arrange(start_year)
for(j in 1:nrow(yrpairs)){
  sy <- as.integer(yrpairs[j, "start_year"])
  ey <- as.integer(yrpairs[j, "end_year"])
  model_variantsel <- model_variant
  trend_tmp2 <- trends_out %>% 
    filter(model == model,
           model_variant == model_variant,
           data_set == data_set,
           strata_name != "continent",
           start_year == sy,
           end_year == ey)
  m2 <- map_trends(trend_tmp2,
                   base_map_blank = map,
                   title = paste(data_set, model, model_variant),
                   region_name = "strata_name")
  print(m2)
}
# ------------------------------------------------------------------------------



# look into effort correction --------------------------------------------------
# effort par summaries
fit$summary(variables = c("b_raw", "p_raw"), "mean", "sd") %>% View()
mcmc_hist(fit$draws("b_raw"))
mcmc_hist(fit$draws("B"))
mcmc_hist(fit$draws("p_raw"))
mcmc_hist(fit$draws("P"))

# eff correction function
bps <- fit$summary(variables = "p_raw", "mean", "sd") %>% 
  mutate(strat=unlist(str_extract_all(variable, "\\d+"))) %>% 
  select(strat, p=mean, p_sd=sd) %>% 
  left_join(fit$summary(variables = "b_raw", "mean", "sd") %>% 
  mutate(strat=unlist(str_extract_all(variable, "\\d+"))) %>% 
  select(strat, B=mean, B_sd=sd)) %>% 
  select(1,4,5,2,3) %>% 
  mutate(strat=as.numeric(strat))
eff_curve <- expand.grid(strat=as.numeric(bps$strat), eff=seq(1, 100, 3)) %>% 
  arrange(strat) %>% left_join(bps) %>% 
  left_join(strata_information %>% rename(strat=strata_vec)) %>% 
  mutate(sc_eff=eff/mean(eff),
         y=exp(B*((sc_eff^p)-1)/p))

# eff correction plot
ggplot(eff_curve, aes(x=eff, y=y)) + 
  geom_line() + facet_wrap(~strata_name)
eff_curve %>% select(strata_name, B, B_sd, p, p_sd) %>% distinct()

# compare to gam per stratum
library(mgcv)
library(mgcViz)
names(data_prep)
td1 <- data_prep %>% filter(strata_name=="CA-AB-10")
td1 <- data_prep %>% filter(strata_name=="CA-AB-11")
td1 <- data_prep %>% filter(strata_name=="US-AK-4")
td1 <- data_prep %>% filter(strata_name=="US-CA-15")
td1 <- data_prep %>% filter(strata_name=="US-SD-17")
summary(tm1 <- gam(how_many ~ s(count_year, k=12) + s(field_hours, k=4) +
                     s(circle_vec, bs="re"),
                         data=td1, family="nb"))
b <- getViz(tm1)
gridPrint(plot(sm(b, 1)) + l_points(shape = 19, size = 3, alpha = 0.3) + 
            l_fitLine() + l_ciLine(),
          plot(sm(b, 2)) + l_points(shape = 19, size = 3, alpha = 0.3) + 
            l_fitLine() + l_ciLine(), 
          ncol = 2)
# ------------------------------------------------------------------------------

