# setup ------------------------------------------------------------------------
# library(patchwork)
library(posterior)
library(bayesplot)
library(patchwork)
#library(geofacet)
#library(ggrepel)
library(ggpattern)
library(HDInterval)
library(naturecounts)
library(cmdstanr)
library(sf)
library(tidyverse)
library(bbsBayes2)

# some some helper functions
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions/neighbours_define.R")
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("functions/map_trends.R")
HDL <- function(x,int,upper = TRUE){
  b <- HDInterval::hdi(x,int)
  return(ifelse(upper,b[["upper"]],b[["lower"]]))
}

# settings
species <- "American Dipper"
species_l <- "Cinclus_mexicanus"
stratification <- "bbs_usgs"
models <- model <- "first_diff"
model_variants <- model_variant <- "spatial"
data_sets <- data_set <- "cbc"
# ------------------------------------------------------------------------------



# get fit object and other input data ------------------------------------------
# uses data_prep and fit from data prep and model fit script
# bring in model fit object
fit <- readRDS(paste0("output/fit_",species_l,"_CBC_spatial_first_diff.rds"))
data_prep <- readRDS(paste0("data/data_prep_",species_l,".rds"))
stan_data <- readRDS(paste0("output/datalist_",species_l,"_CBC_spatial_first_diff.rds"))

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
strata_map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
  select(-area_sq_km) %>% 
  inner_join(.,strat_df,
             by = "strata_name") %>% 
  mutate(strata_vec = as.integer(factor(strata_name))) %>% 
  arrange(stratum)

# get generation length
sp_id <- meta_species_taxonomy() %>% 
  filter(english_name %in% species) %>% pull(species_id)
gen_years <- nc_query_table(table="SpeciesLifeHistory") %>%
  filter(subcategDescr=="Average generation length (years)") %>% 
  filter(speciesID %in% sp_id) %>% pull(value) %>% as.numeric()
gen_3_years <- round(gen_years * 3)
# ------------------------------------------------------------------------------



# # stratum level indices and trends ---------------------------------------------
# # posterior samples of stratum indices
# ind_samples <- posterior_samples(fit = fit, parm = "n",
#                                  dims = c("stratum","yr")) %>% 
#   inner_join(., strat_df, by = "stratum") %>% 
#   inner_join(., yrs, by = "yr")
# 
# # stratum level index summaries
# inds_strat <- 
#   ind_samples %>%
#   group_by(strata_name,year) %>% 
#   summarise(index = median(.value),
#             index_q_0.05 = HDL(.value,0.9,upper = FALSE),
#             index_q_0.95 = HDL(.value,0.9,upper = TRUE),
#             .groups = "keep") %>% 
#   mutate(index_type = "full")
# 
# # use stratum indices to compute stratum trends of different length
# first_years <- c(1966, 1970, 1993, 2009)
# trends_strata <- NULL
# for(j in 1:length(first_years)){
#   ys <- first_years[j]
#   ye <- 2019
#   nyrs <- ye-ys
#   trend_tmp <- ind_samples %>% 
#     filter(year %in% c(ys,ye)) %>% 
#     ungroup() %>% 
#     select(-matches(match = "yr",ignore.case = FALSE)) %>% 
#     pivot_wider(names_from = year,
#                 values_from = .value,
#                 names_prefix = "Y") %>% 
#     rename_with(., ~gsub(replacement = "start",
#                          pattern = paste0("Y",ys),.x,
#                          fixed = TRUE))%>% 
#     rename_with(., ~gsub(replacement = "end",
#                          pattern = paste0("Y",ye),.x,
#                          fixed = TRUE))%>% 
#     group_by(.draw,strata_name) %>% 
#     summarise(end = sum(end),
#               start = sum(start),
#               t = texp(end/start,ny = nyrs),
#               ch = chng(end/start),
#               .groups = "keep") %>% 
#     group_by(strata_name) %>% 
#     summarise(trend = mean(t),
#               lci = quantile(t,0.025,names = FALSE),
#               uci = quantile(t,0.975,names = FALSE),
#               width_CI = uci-lci) %>% 
#     mutate(model = model,
#            model_variant = model_variant,
#            start_year = ys,
#            end_year = ye,
#            trend_length = nyrs,
#            data_set = data_set)
#   trends_strata <- bind_rows(trends_strata, trend_tmp)
# }
# # ------------------------------------------------------------------------------
# 
# 
# 
# # continental level indices and trends -----------------------------------------
# # make continent scaled index samples by summing across strata
# comp_samples <- 
#   ind_samples %>%
#   mutate(.value = .value * area_weight) %>% 
#   group_by(.draw, year) %>% #draw-wise summary of area-weighted indices
#   summarise(.value = sum(.value),
#             .groups = "drop")
# 
# # continent scaled index summaries
# # inds and smooth are the same because its first difference
# inds_comp <- inds_comp_smooth <-
#   comp_samples %>%  
#   group_by(year) %>% 
#   summarise(index = median(.value),
#             index_q_0.05 = HDL(.value,0.9,upper = FALSE),
#             index_q_0.95 = HDL(.value,0.9,upper = TRUE),
#             .groups = "keep") %>% 
#   mutate(strata_name = "continent") %>% 
#   mutate(index_type = "full")
# 
# # combine stratum and continental index summaries for later
# inds_all <- bind_rows(inds_comp,
#                       inds_strat) %>% 
#   mutate(model = "first_difference",
#          model_variant = "spatial",
#          data_set = "cbc")
# 
# # use indices to make continent scaled trends of different lengths
# trends_comp <- NULL
# for(j in 1:length(first_years)){
#   ys <- first_years[j]
#   ye <- 2019
#   nyrs <- ye-ys
#   trend_tmp <- comp_samples %>% 
#     filter(year %in% c(ys,ye)) %>% 
#     ungroup() %>% 
#     select(-matches(match = "yr",ignore.case = FALSE)) %>% 
#     pivot_wider(names_from = year,
#                 values_from = .value,
#                 names_prefix = "Y") %>% 
#     rename_with(., ~gsub(replacement = "start",
#                          pattern = paste0("Y",ys),.x,
#                          fixed = TRUE))%>% 
#     rename_with(., ~gsub(replacement = "end",
#                          pattern = paste0("Y",ye),.x,
#                          fixed = TRUE))%>% 
#     group_by(.draw) %>% 
#     summarise(end = sum(end),
#               start = sum(start),
#               t = texp(end/start,ny = nyrs),
#               ch = chng(end/start),
#               .groups = "drop") %>% 
#     summarise(trend = mean(t),
#               lci = quantile(t,0.025,names = FALSE),
#               uci = quantile(t,0.975,names = FALSE),
#               width_CI = uci-lci) %>% 
#     mutate(model = model,
#            model_variant = model_variant,
#            start_year = ys,
#            end_year = ye,
#            trend_length = nyrs,
#            data_set = data_set,
#            strata_name = "continent")
#   trends_comp <- bind_rows(trends_comp, trend_tmp)
# }
# # ------------------------------------------------------------------------------
# 
# 
# 
# # index plots ------------------------------------------------------------------
# # continent
# tmp <- inds_all %>% 
#   filter(strata_name == "continent")
# tp <- ggplot(data = tmp,
#             aes(x = year, y = index)) +
#   geom_ribbon(aes(ymin = index_q_0.05,
#                   ymax = index_q_0.95),
#               alpha = 0.2, fill="darkblue") +
#   geom_line(color="darkblue") +
#   theme_bw(); tp
# 
# # strata
# obs_mean <- data_prep %>% 
#   select(strata_name, year=count_year, how_many) %>% 
#   group_by(strata_name, year) %>% 
#   summarise(index = mean(how_many))
# tmp <- inds_all %>% 
#   filter(strata_name != "continent",
#          data_set == "cbc",
#          index_type == "full")
# tp <- ggplot() +
#   geom_ribbon(data = tmp, 
#               aes(x = year, ymin = index_q_0.05, ymax = index_q_0.95),
#               alpha = 0.3, fill="darkblue") +
#   geom_line(data = tmp, 
#             aes(x = year, y = index), color="darkblue") +
#   geom_point(data = obs_mean, 
#              aes(x = year, y = index),
#              alpha=0.3, size=2, shape=1) +
#   facet_wrap(~ strata_name, scales = "free") +
#   theme_bw(); tp
# # ------------------------------------------------------------------------------
# 
# 
# 
# # map trends -------------------------------------------------------------------
# trends_out <- bind_rows(trends_strata, trends_comp)
# yrpairs <- trends_out %>% 
#   filter(data_set == as.character("cbc")) %>% 
#   select(start_year, end_year) %>% 
#   distinct() %>% 
#   arrange(start_year)
# for(j in 1:nrow(yrpairs)){
#   sy <- as.integer(yrpairs[j, "start_year"])
#   ey <- as.integer(yrpairs[j, "end_year"])
#   model_variantsel <- model_variant
#   trend_tmp2 <- trends_out %>% 
#     filter(model == model,
#            model_variant == model_variant,
#            data_set == data_set,
#            strata_name != "continent",
#            start_year == sy,
#            end_year == ey)
#   m2 <- map_trends(trend_tmp2,
#                    base_map_blank = strata_map,
#                    title = paste(data_set, model, model_variant),
#                    region_name = "strata_name")
#   print(m2)
# }
# # ------------------------------------------------------------------------------



# effort correction by stratum -------------------------------------------------

# please help understanding scale of effort on x axes **************************

effort_preds <- data.frame(p_of_mean_effort = stan_data$effort_preds,
                           effort = c(1:stan_data$n_effort_preds))

# effort curve by stratum samples
eff <- posterior_samples(fit = fit, parm = "effort_strata",
                            dims = c("stratum","effort")) %>% 
  posterior_sums(dims = c("stratum","effort")) %>% 
  inner_join(., strat_df, by = "stratum") %>% 
  inner_join(effort_preds, by = "effort")

# plot
vis_eff <- ggplot(data = eff,
                  aes(x = p_of_mean_effort, y = mean))+
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975),
              colour = NA, alpha = 0.2)+
  geom_line()+
  geom_rug(data = data_prep,
             aes(x = scaled_effort),
           inherit.aes = FALSE)+ # rug plot to show observed effort
  scale_x_continuous(limits=c(0, quantile(eff$p_of_mean_effort, probs=0.8))) +
  facet_wrap(vars(strata_name),
             scales = "free_y"); vis_eff
# ------------------------------------------------------------------------------



# effort hyperparameters -------------------------------------------------------
EFF <- posterior_samples(fit = fit, parm = "EFFORT",
                             dims = c("effort")) %>% 
  posterior_sums(dims = c("effort")) %>% 
  inner_join(effort_preds,by = "effort")
vis_EFF <- ggplot(data = EFF,
                  aes(x = p_of_mean_effort, y = mean))+
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975),
              colour = NA, alpha = 0.2)+
  geom_line()+
  geom_rug(data = data_prep,
           aes(x = scaled_effort),
           inherit.aes = FALSE); vis_EFF

# effort par summaries
eff_par_sum <- fit$summary(
  variables = c("b_raw", "p_raw"),
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.025, .975))
) %>% mutate(sig=ifelse(q2.5>0 | q97.5<0, 1, 0))
View(eff_par_sum)

mcmc_hist(fit$draws("b_raw"))
mcmc_hist(fit$draws("B"))
mcmc_hist(fit$draws("p_raw"))
mcmc_hist(fit$draws("P"))
# ------------------------------------------------------------------------------




# define trend periods ---------------------------------------------------------
# need to make sure I get the right YEAR ***************************************
# start year for indices
year_1 <- 1966 # first year
year_exp <- 1993 # year when expanded BBS analyses started
year_N <- 2019 # last year
year_10 <- year_N - 10 + 1
year_3g <- year_N - gen_3_years + 1
# ------------------------------------------------------------------------------



# bcr index summaries ----------------------------------------------------------
# weights table for aggregation during index and trend production
bcr_wts <- strata_map %>% 
  select(Stratum_Factored=strata_vec, Area=area_sq_km, bcr=bcr) %>% 
  st_drop_geometry() %>% as.data.frame()

# bcr indices
bcr_idxs <- index_function(fit = fit,
                           parameter = "n",
                           strat = "Stratum_Factored",
                           year = "Year",
                           first_dim = "s",
                           quant = 0.95,
                           weights_df = bcr_wts,
                           area = "Area",
                           summary_regions = "bcr",
                           year_1 = year_1,
                           to_summarise = T)

# observed counts ********************************** any ideas on getting effort 
# corrected 'observed'?
obs_mean <- data_prep %>% 
  select(bcr, year=count_year, how_many) %>% 
  group_by(bcr, year) %>% 
  summarise(index = mean(how_many)) %>% 
  ungroup()

# plot indices for full bcr time series
bcr_ind_plot <- ggplot() +
  geom_ribbon(data = bcr_idxs$indices, 
              aes(x = Year+year_1-1, ymin = lci, ymax = uci),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = bcr_idxs$indices, 
            aes(x = Year+year_1-1, y = mean), color="darkblue") +
  geom_point(data = obs_mean,
             aes(x = year, y = index),
             alpha=0.3, size=2, shape=1) +
  facet_wrap(~ paste0("BCR", bcr), scales = "free") +
  scale_x_continuous(breaks=seq(min(bcr_idxs$indices$Year+year_1-1), 
                                max(bcr_idxs$indices$Year+year_1-1), 
                                7)) +
  theme_bw() +
  labs(x="Year", y="Abundance index (95% CrI) per year and BCR"); bcr_ind_plot

# turn indices into lt trends
bcr_lt_trds <- trends_function(ind_list = bcr_idxs, start_year = year_1, 
                               end_year = year_N, quant = 0.95) %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
bcr_exp_trds <- trends_function(ind_list = bcr_idxs, start_year = year_exp, 
                                end_year = year_N, quant = 0.95) %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
bcr_10yr_trds <- trends_function(ind_list = bcr_idxs, start_year = year_10, 
                                 end_year = year_N, quant = 0.95) %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_year = gen_3_years)

# turn indices into 3 gen trends
bcr_3gen_trds <- trends_function(ind_list = bcr_idxs, start_year = year_3g, 
                                 end_year = year_N, quant = 0.95) %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# plot bcr trends
tp_lt <- map_function(trds = bcr_lt_trds, spunit = "bcr")
tp_exp <- map_function(trds = bcr_exp_trds, spunit = "bcr")
tp_10yr <- map_function(trds = bcr_10yr_trds, spunit = "bcr")
tp_3gen <- map_function(trds = bcr_3gen_trds, spunit = "bcr")

# all together
ptch1 <- ((tp_lt + tp_exp) / (tp_10yr + tp_3gen))
ptch1 <- ptch1 + plot_layout(guides = 'collect')
cairo_pdf(paste0("./output/", gsub(" ", "_", species), "_bcr_trend_map.pdf"),  
          width = 9.25, height = 10)
ptch1
dev.off()
# ------------------------------------------------------------------------------



# province and state index summaries -------------------------------------------
# weights table for aggregation during index and trend production
ps_wts <- strata_map %>% 
  select(Stratum_Factored=strata_vec, Area=area_sq_km, prov_state=prov_state) %>% 
  st_drop_geometry() %>% as.data.frame()

# province and state indices
ps_idxs <- index_function(fit = fit,
                           parameter = "n",
                           strat = "Stratum_Factored",
                           year = "Year",
                           first_dim = "s",
                           quant = 0.95,
                           weights_df = ps_wts,
                           area = "Area",
                           summary_regions = "prov_state",
                           year_1 = year_1,
                           to_summarise = T)

obs_mean <- data_prep %>% 
  select(prov_state, year=count_year, how_many) %>% 
  group_by(prov_state, year) %>% 
  summarise(index = mean(how_many)) %>% 
  ungroup()

# plot indices for full bcr time series
ps_ind_plot <- ggplot() +
  geom_ribbon(data = ps_idxs$indices, 
              aes(x = Year+year_1-1, ymin = lci, ymax = uci),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = ps_idxs$indices, 
            aes(x = Year+year_1-1, y = mean), color="darkblue") +
  geom_point(data = obs_mean,
             aes(x = year, y = index),
             alpha=0.3, size=2, shape=1) +
  facet_wrap(~ prov_state, scales = "free") +
  scale_x_continuous(breaks=seq(min(ps_idxs$indices$Year+year_1-1), 
                                max(ps_idxs$indices$Year+year_1-1), 
                                7)) +
  theme_bw() +
  labs(x="Year", 
       y="Abundance index (95% CrI) per year and province or state"); ps_ind_plot

# turn indices into lt trends
ps_lt_trds <- trends_function(ind_list = ps_idxs, start_year = year_1, 
                              end_year = year_N, quant = 0.95) %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
ps_exp_trds <- trends_function(ind_list = ps_idxs, start_year = year_exp, 
                               end_year = year_N, quant = 0.95) %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
ps_10yr_trds <- trends_function(ind_list = ps_idxs, start_year = year_10, 
                                end_year = year_N, quant = 0.95) %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
ps_3gen_trds <- trends_function(ind_list = ps_idxs, start_year = year_3g, 
                                end_year=year_N, quant = 0.95) %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# plot province state trends
tp_lt <- map_function(trds = ps_lt_trds, spunit = "prov_state")
tp_exp <- map_function(trds = ps_exp_trds, spunit = "prov_state")
tp_10yr <- map_function(trds = ps_10yr_trds, spunit = "prov_state")
tp_3gen <- map_function(trds = ps_3gen_trds, spunit = "prov_state")

# all together
ptch1 <- ((tp_lt + tp_exp) / (tp_10yr + tp_3gen))
ptch1 <- ptch1 + plot_layout(guides = 'collect')
cairo_pdf(paste0("./output/", gsub(" ", "_", species), "_prov_state_trend_map.pdf"),  
          width = 9.25, height = 10)
ptch1
dev.off()
# ------------------------------------------------------------------------------






# country index summaries ------------------------------------------------------
# country weights table
cntry_wts <- strata_map %>% 
  select(Stratum_Factored=strata_vec, Area=area_sq_km, country=country) %>% 
  st_drop_geometry() %>% as.data.frame()

# country indices
cntry_idxs <- index_function(fit = fit,
                          parameter = "n",
                          strat = "Stratum_Factored",
                          year = "Year",
                          first_dim = "s",
                          quant = 0.95,
                          weights_df = cntry_wts,
                          area = "Area",
                          summary_regions = "country",
                          year_1 = 1966,
                          to_summarise = T)

obs_mean <- data_prep %>% 
  select(country, year=count_year, how_many) %>% 
  group_by(country, year) %>% 
  summarise(index = mean(how_many)) %>% 
  ungroup()

# country time series
cntry_idx_plot <- ggplot() +
  geom_ribbon(data = cntry_idxs$indices, 
              aes(x = Year+year_1-1, ymin = lci, ymax = uci),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = cntry_idxs$indices, 
            aes(x = Year+year_1-1, y = mean), color="darkblue") +
  geom_point(data = obs_mean,
             aes(x = year, y = index),
             alpha=0.3, size=2, shape=1) +
  facet_wrap(~ country, scales = "free") +
  labs(x = "Year", y = "Abundance index (95% CrI)") +
  theme_bw(); cntry_idx_plot

# turn indices into lt trends
cntry_lt_trds <- trends_function(ind_list = cntry_idxs, start_year = year_1, 
                                 end_year = year_N, quant = 0.95) %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
cntry_exp_trds <- trends_function(ind_list = cntry_idxs, start_year = year_exp, 
                               end_year = year_N, quant = 0.95) %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
cntry_10yr_trds <- trends_function(ind_list = cntry_idxs, start_year = year_10, 
                                end_year = year_N, quant = 0.95) %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
cntry_3gen_trds <- trends_function(ind_list = cntry_idxs, start_year = year_3g, 
                                end_year = year_N, quant = 0.95) %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)
# ------------------------------------------------------------------------------



# continent summaries ----------------------------------------------------------
# continent weights table
cont_wts <- strata_map %>% 
  select(Stratum_Factored=strata_vec, Area=area_sq_km) %>% 
  mutate(continent="Canada and United States of America") %>% # add a continent field
  st_drop_geometry() %>% as.data.frame()

# continent indices
cont_idxs <- index_function(fit = fit,
                             parameter = "n",
                             strat = "Stratum_Factored",
                             year = "Year",
                             first_dim = "s",
                             quant = 0.95,
                             weights_df = cont_wts,
                             area = "Area",
                             summary_regions = "continent",
                             year_1 = 1966,
                             to_summarise = T)

obs_mean <- data_prep %>% 
  mutate(continent="Canada and United States of America") %>% # add a continent field
  select(continent, year=count_year, how_many) %>% 
  group_by(continent, year) %>% 
  summarise(index = mean(how_many))

# continent time series
cont_idx_plot <- ggplot() +
  geom_ribbon(data = cont_idxs$indices, 
              aes(x = Year+year_1-1, ymin = lci, ymax = uci),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = cont_idxs$indices, 
            aes(x = Year+year_1-1, y = mean), color="darkblue") +
  geom_point(data = obs_mean,
             aes(x = year, y = index),
             alpha=0.3, size=2, shape=1) +
  facet_wrap(~ continent, scales = "free") +
  labs(x = "Year", y = "Abundance index (95% CrI)") +
  theme_bw(); cont_idx_plot


# turn indices into lt trends
cont_lt_trds <- trends_function(ind_list = cont_idxs, start_year = year_1, 
                                 end_year = year_N, quant = 0.95) %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
cont_exp_trds <- trends_function(ind_list = cont_idxs, start_year = year_exp, 
                                  end_year = year_N, quant = 0.95) %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
cont_10yr_trds <- trends_function(ind_list = cont_idxs, start_year = year_10, 
                                   end_year = year_N, quant = 0.95) %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
cont_3gen_trds <- trends_function(ind_list = cont_idxs, start_year = year_3g, 
                                   end_year = year_N, quant = 0.95) %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)
# ------------------------------------------------------------------------------



