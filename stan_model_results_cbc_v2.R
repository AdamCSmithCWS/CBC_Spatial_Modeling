


# setup ------------------------------------------------------------------------
# library(patchwork)
# library(posterior) # this is called in script but not loaded due to masking
#library(bayesplot)
library(patchwork)
#library(geofacet)
#library(ggrepel)
library(ggforce)
library(ggpattern)
library(naturecounts)
library(cmdstanr)
library(sf)
library(bbsBayes2)
library(tidyverse)

# some some helper functions
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("functions/draft_functions_v2.R")

# settings
species <- "American Dipper"
species_l <- "Cinclus_mexicanus"
stratification <- "bbs_usgs"
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

# start year for indices
year_1 <- 1966 # first year
year_exp <- 1993 # year when expanded BBS analyses started
year_N <- 2019 # last year
year_10 <- year_N - 10 # includes eleven years data, i.e. inclusive
year_3g <- year_N - gen_3_years
pif_quant <- c(0.025, 0.05, 0.165, 0.5, 0.835, 0.95, 0.975)
# ------------------------------------------------------------------------------




# diagnostics ------------------------------------------------------------------
# abundance indices
nit_par_diags <- fit$draws() %>% 
  posterior::subset_draws(variable = c("n")) %>% 
  posterior::as_draws_df() %>% 
  posterior::summarise_draws(mean, sd, 
                  ~quantile(.x, probs = pif_quant),
                  posterior::default_convergence_measures()) %>% 
  rename_at(vars(c(4:8)), .funs = ~paste0("q", .x)) %>% 
  mutate(stratum = as.numeric(str_split(string = gsub("n|\\[|\\]", "", variable),
                           pattern = ",", simplify = T)[,1]),
         year = as.numeric(str_split(string = gsub("n|\\[|\\]", "", variable),
                           pattern = ",", simplify = T)[,2]))

# nit diagnostic warnings
nit_rhat_test <- any(nit_par_diags$rhat > 1.05)
if(nit_rhat_test == T){
  fileConn <- file("./output/nit_rhat_WARNING.txt")
  writeLines("nit parameters did not converge well, rhat > 1.05", fileConn)
  close(fileConn)
}
nit_ess_test <- any(nit_par_diags$ess_bulk<200)
if(nit_ess_test == T){
  fileConn <- file("./output/nit_ess_WARNING.txt")
  writeLines("nit parameters have low sample size, ess < 200", fileConn)
  close(fileConn)
}

# effort parameters
eff_par_diags <- fit$draws() %>% 
  posterior::subset_draws(variable = c("B", "b_raw", "P", "p_raw")) %>% 
  posterior::as_draws_df() %>% 
  mutate_at(vars(starts_with("b_raw")), .funs = ~ .x + B) %>% 
  mutate_at(vars(starts_with("p_raw")), .funs = ~ .x + P) %>% 
  rename_with(~str_replace(., "b_raw", "B_strat"), starts_with("b_raw")) %>% 
  rename_with(~str_replace(., "p_raw", "P_strat"), starts_with("p_raw")) %>% 
  posterior::as_draws_df() %>% 
  posterior::summarise_draws(mean, sd, 
                  ~quantile(.x, probs = pif_quant)) %>% 
  rename_at(vars(c(4:8)), .funs = ~paste0("q", .x)) %>% 
  bind_cols(fit$draws() %>% 
              posterior::subset_draws(variable = c("B", "b_raw", "P", "p_raw")) %>% 
              posterior::as_draws_df() %>% 
              posterior::summarise_draws(posterior::default_convergence_measures()) %>% 
              select(-1)) %>% 
  mutate(stratum = as.numeric(str_extract(string = variable, pattern = "\\d+")),
         year = NA)

# eff diagnostic warnings
eff_rhat_test <- any(eff_par_diags$rhat > 1.5)
if(eff_rhat_test == T){
  fileConn<-file("./output/eff_rhat_WARNING.txt")
  writeLines("effort parameters did not converge well, rhat > 1.5", fileConn)
  close(fileConn)
}
eff_ess_test <- any(eff_par_diags$ess_bulk < 50)
if(eff_ess_test == T){
  fileConn <- file("./output/eff_ess_WARNING.txt")
  writeLines("effort parameters have low sample size, ess < 50", fileConn)
  close(fileConn)
}

# combine diagnostic summaries
par_est_sums <- nit_par_diags %>% bind_rows(eff_par_diags) %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(strata_name, stratum), by="stratum") %>% 
  mutate(real_year = year + year_1 - 1) %>% 
  select(strata_name, real_year, stratum, year, variable, everything())
 
# save diagnostic summaries
write_csv(par_est_sums, 
          paste0("./output/", 
                 gsub(" ", "_", species), "_parameter_estimate_summaries.csv"))
# ------------------------------------------------------------------------------




# effort correction by stratum -------------------------------------------------
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
                  aes(x = p_of_mean_effort, y = mean)) +
  geom_ribbon(aes(ymin = Q_025,ymax = Q_975),
              colour = NA, alpha = 0.2, fill = "darkblue") +
  geom_line(col = "darkblue") +
  geom_rug(data = data_prep, aes(x = scaled_effort), inherit.aes = FALSE) + 
  scale_x_continuous(limits = c(0, quantile(eff$p_of_mean_effort, probs = 0.8))) +
  facet_wrap(vars(strata_name), scales = "free_y") +
  labs(x="Party hours / mean party hours", 
       y="Scaled partial effect per stratum") +
  theme_bw()

# save
cairo_pdf(paste0("./output/", gsub(" ", "_", species), 
                 "_stratum_effort_correction_plots.pdf"),  
          width = 11, height = 8.5)
vis_eff
dev.off()
# ------------------------------------------------------------------------------




# bcr index summaries ----------------------------------------------------------
# weights table for aggregation
bcr_wts <- strata_map %>% 
  select(stratum = strata_vec, Area = area_sq_km, bcr = bcr) %>% 
  st_drop_geometry() %>% as.data.frame()

# bcr indices
bcr_idx_lst <- index_function(fit = fit,
                             parameter = "n",
                             strat = "stratum",
                             year = "Year",
                             first_dim = "s",
                             quant = 0.95,
                             weights_df = bcr_wts,
                             area = "Area",
                             summary_regions = "bcr",
                             year_1 = year_1,
                             to_summarise = T)

# add smooth indices
bcr_smooth_idxs <- get_smooth_ns(idx_list = bcr_idx_lst)
bcr_idx_lst$indices_smooth <- bcr_smooth_idxs$indices_smooth
bcr_idx_lst$samples_smooth <- bcr_smooth_idxs$samples_smooth

# get data to plot indices for full time series
pd1 <- bcr_idx_lst$indices_smooth %>% 
  mutate(bcr_sort = paste0("BCR", str_pad(bcr, pad = "0", side = "left", 
                                          width = 2))) %>% 
  arrange(bcr_sort, Year)
pd2 <- bcr_idx_lst$indices %>% 
  mutate(bcr_sort = paste0("BCR", str_pad(bcr, pad = "0", side = "left", 
                                          width = 2))) %>% 
  arrange(bcr_sort, Year)
regs <- sort(unique(pd2$bcr))
pgs <- 1:ceiling(length(regs) / 9)

# loop through plot pages
dev.off()
for(pg in 1:length(pgs)){
  stt <- pg
  if(stt == 1) stt <- 1
  if(stt > 1) stt <- ((pg - 1) * 9) + 1
  edd <- stt + 9 - 1
  fst <- regs[stt]
  ltt <- regs[edd]
  if(pg == max(pgs)) ltt <- last(regs)
  fn <- paste0("./output/", gsub(" ", "_", species),
         "_annual_indices_bcr_", fst, "-", ltt, ".pdf")
  plt_i <- ggplot() +
    geom_ribbon(data = pd1, 
                aes(x = true_year, ymin = lci, ymax = uci),
                alpha = 0.3, fill="darkblue") +
    geom_line(data = pd1, 
              aes(x = true_year, y = mean), color="darkblue") +
    geom_linerange(data = pd2, 
                   aes(x = true_year, ymin = lci, ymax = uci),
                   alpha = 0.6, linewidth=0.5) +
    geom_point(data = pd2, 
               aes(x = true_year, y = mean), alpha = 0.6, shape=16) +
    ggforce::facet_wrap_paginate(~ bcr_sort, scales="free",
                                 ncol = 3, nrow = 3, page = pg) +
    scale_x_continuous(breaks=seq(min(pd1$true_year), 
                                  max(pd1$true_year), 
                                  7)) +
    labs(x="Year", y="Abundance index (95% CrI) per year and BCR") +
    theme_bw()
  cairo_pdf(fn, width = 11, height = 8.5)
  print(plt_i)
  dev.off()
}
  
# turn indices into lt trends
bcr_lt_trds <- trends_function(ind_list = bcr_idx_lst, start_year = year_1, 
                              end_year = year_N, quant = 0.95)[[
                                "smooth_endpoint_trends"]] %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
bcr_exp_trds <- trends_function(ind_list = bcr_idx_lst, start_year = year_exp, 
                               end_year = year_N, quant = 0.95)[[
                                 "smooth_endpoint_trends"]] %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
bcr_10yr_trds <- trends_function(ind_list = bcr_idx_lst, start_year = year_10, 
                                end_year = year_N, quant = 0.95)[[
                                  "smooth_endpoint_trends"]] %>% 
  rename(stratum = bcr) %>% 
  mutate(strata_name = paste0("BCR", stratum),
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
bcr_3gen_trds <- trends_function(ind_list = bcr_idx_lst, start_year = year_3g, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
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
# weights table for aggregation
ps_wts <- strata_map %>% 
  select(stratum = strata_vec, Area = area_sq_km, prov_state = prov_state) %>% 
  st_drop_geometry() %>% as.data.frame()

# province and state indices
ps_idx_lst <- index_function(fit = fit,
                           parameter = "n",
                           strat = "stratum",
                           year = "Year",
                           first_dim = "s",
                           quant = 0.95,
                           weights_df = ps_wts,
                           area = "Area",
                           summary_regions = "prov_state",
                           year_1 = year_1,
                           to_summarise = T)

# add smooth indices
ps_smooth_idxs <- get_smooth_ns(idx_list = ps_idx_lst)
ps_idx_lst$indices_smooth <- ps_smooth_idxs$indices_smooth
ps_idx_lst$samples_smooth <- ps_smooth_idxs$samples_smooth

# plot indices for full time series
pd1 <- ps_idx_lst$indices_smooth %>% arrange(prov_state, Year)
pd2 <- ps_idx_lst$indices %>% arrange(prov_state, Year)
regs <- sort(unique(pd2$prov_state))
pgs <- 1:ceiling(length(regs) / 9)

# loop through plot pages
dev.off()
for(pg in 1:length(pgs)){
  stt <- pg
  if(stt == 1) stt <- 1
  if(stt > 1) stt <- ((pg - 1) * 9) + 1
  edd <- stt + 9 - 1
  fst <- regs[stt]
  ltt <- regs[edd]
  if(pg == max(pgs)) ltt <- last(regs)
  fn <- paste0("./output/", gsub(" ", "_", species),
               "_annual_indices_prov_state_", fst, "-", ltt, ".pdf")
  plt_i <- ggplot() +
    geom_ribbon(data = pd1, 
                aes(x = true_year, ymin = lci, ymax = uci),
                alpha = 0.3, fill="darkblue") +
    geom_line(data = pd1, 
              aes(x = true_year, y = mean), color="darkblue") +
    geom_linerange(data = pd2, 
                   aes(x = true_year, ymin = lci, ymax = uci),
                   alpha = 0.6, linewidth=0.5) +
    geom_point(data = pd2, 
               aes(x = true_year, y = mean), alpha = 0.6, shape=16) +
    ggforce::facet_wrap_paginate(~ prov_state, scales="free",
                                 ncol = 3, nrow = 3, page = pg) +
    scale_x_continuous(breaks=seq(min(pd1$true_year), 
                                  max(pd1$true_year), 
                                  7)) +
    labs(x="Year", y="Abundance index (95% CrI) per year and province or state") +
    theme_bw()
  cairo_pdf(fn, width = 11, height = 8.5)
  print(plt_i)
  dev.off()
}

# turn indices into lt trends
ps_lt_trds <- trends_function(ind_list = ps_idx_lst, start_year = year_1, 
                              end_year = year_N, quant = 0.95)[[
                                "smooth_endpoint_trends"]] %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
ps_exp_trds <- trends_function(ind_list = ps_idx_lst, start_year = year_exp, 
                               end_year = year_N, quant = 0.95)[[
                                 "smooth_endpoint_trends"]] %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
ps_10yr_trds <- trends_function(ind_list = ps_idx_lst, start_year = year_10, 
                                end_year = year_N, quant = 0.95)[[
                                  "smooth_endpoint_trends"]] %>% 
  rename(stratum = prov_state) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
ps_3gen_trds <- trends_function(ind_list = ps_idx_lst, start_year = year_3g, 
                                end_year=year_N, quant = 0.95)[[
                                  "smooth_endpoint_trends"]] %>% 
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
# weights table for aggregation
cntry_wts <- strata_map %>% 
  select(stratum = strata_vec, Area = area_sq_km, country = country) %>% 
  st_drop_geometry() %>% as.data.frame()

# country indices
cntry_idx_lst <- index_function(fit = fit,
                          parameter = "n",
                          strat = "stratum",
                          year = "Year",
                          first_dim = "s",
                          quant = 0.95,
                          weights_df = cntry_wts,
                          area = "Area",
                          summary_regions = "country",
                          year_1 = 1966,
                          to_summarise = T)

# add smooth indices
cntry_smooth_idxs <- get_smooth_ns(idx_list = cntry_idx_lst)
cntry_idx_lst$indices_smooth <- cntry_smooth_idxs$indices_smooth
cntry_idx_lst$samples_smooth <- cntry_smooth_idxs$samples_smooth

# plot indices for full time series
pd1 <- cntry_idx_lst$indices_smooth %>% arrange(country, Year)
pd2 <- cntry_idx_lst$indices %>% arrange(country, Year)
cntry_ind_plot <- ggplot() +
  geom_ribbon(data = pd1, 
              aes(x = true_year, ymin = lci, ymax = uci),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = pd1, 
            aes(x = true_year, y = mean), color="darkblue") +
  geom_linerange(data = pd2, 
                 aes(x = true_year, ymin = lci, ymax = uci),
                 alpha = 0.6, linewidth=0.5) +
  geom_point(data = pd2, 
             aes(x = true_year, y = mean), alpha = 0.6, shape=16) +
  ggforce::facet_wrap_paginate(~ country, scales="free",
                               ncol = 2, nrow = 1, page = 1) +
  scale_x_continuous(breaks=seq(min(pd1$true_year), 
                                max(pd1$true_year), 
                                7)) +
  labs(x="Year", y="Abundance index (95% CrI) per year and country") +
  theme_bw()

# save index plot
cairo_pdf(paste0("./output/", gsub(" ", "_", species), 
                 "_annual_indices_country.pdf"),  
          width = 11, height = 8.5)
cntry_ind_plot
dev.off()

# turn indices into lt trends
cntry_lt_trds <- trends_function(ind_list = cntry_idx_lst, start_year = year_1, 
                              end_year = year_N, quant = 0.95)[[
                                "smooth_endpoint_trends"]] %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
cntry_exp_trds <- trends_function(ind_list = cntry_idx_lst, start_year = year_exp, 
                                  end_year = year_N, quant = 0.95)[[
                                    "smooth_endpoint_trends"]] %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
cntry_10yr_trds <- trends_function(ind_list = cntry_idx_lst, start_year = year_10, 
                                   end_year = year_N, quant = 0.95)[[
                                     "smooth_endpoint_trends"]] %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
cntry_3gen_trds <- trends_function(ind_list = cntry_idx_lst, start_year = year_3g, 
                                   end_year = year_N, quant = 0.95)[[
                                     "smooth_endpoint_trends"]] %>% 
  rename(stratum = country) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)
# ------------------------------------------------------------------------------




# continent summaries ----------------------------------------------------------
# continent weights table
cont_wts <- strata_map %>% 
  select(stratum=strata_vec, Area=area_sq_km) %>% 
  mutate(continent="Canada and United States of America") %>% # add a continent field
  st_drop_geometry() %>% as.data.frame()

# continent indices
cont_idx_lst <- index_function(fit = fit,
                             parameter = "n",
                             strat = "stratum",
                             year = "Year",
                             first_dim = "s",
                             quant = 0.95,
                             weights_df = cont_wts,
                             area = "Area",
                             summary_regions = "continent",
                             year_1 = year_1,
                             to_summarise = T)

# add smooth indices
cont_smooth_idxs <- get_smooth_ns(idx_list = cont_idx_lst)
cont_idx_lst$indices_smooth <- cont_smooth_idxs$indices_smooth
cont_idx_lst$samples_smooth <- cont_smooth_idxs$samples_smooth

# plot indices for full time series
cont_ind_plot <- ggplot() +
  geom_ribbon(data = cont_idx_lst$indices_smooth, 
              aes(x = true_year, ymin = lci, ymax = uci),
              alpha = 0.3, fill="darkblue") +
  geom_line(data = cont_idx_lst$indices_smooth, 
            aes(x = true_year, y = mean), color="darkblue") +
  geom_linerange(data = cont_idx_lst$indices, 
                 aes(x = true_year, ymin = lci, ymax = uci),
                 alpha = 0.6, linewidth=0.5) +
  geom_point(data = cont_idx_lst$indices, 
             aes(x = true_year, y = mean), alpha = 0.6, shape=16) +
  facet_wrap(~ continent, scales = "free") +
  scale_x_continuous(breaks=seq(min(cont_idx_lst$indices$true_year), 
                                max(cont_idx_lst$indices$true_year), 
                                7)) +
  labs(x="Year", y="Abundance index (95% CrI) per year") +
  theme_bw()

# save index plot
cairo_pdf(paste0("./output/", gsub(" ", "_", species), 
                 "_annual_indices_survey.pdf"),  
          width = 11, height = 8.5)
cont_ind_plot
dev.off()

# turn indices into lt trends
cont_lt_trds <- trends_function(ind_list = cont_idx_lst, start_year = year_1, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
cont_exp_trds <- trends_function(ind_list = cont_idx_lst, start_year = year_exp, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
cont_10yr_trds <- trends_function(ind_list = cont_idx_lst, start_year = year_10, 
                                  end_year = year_N, quant = 0.95)[[
                                    "smooth_endpoint_trends"]] %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
cont_3gen_trds <- trends_function(ind_list = cont_idx_lst, start_year = year_3g, 
                                  end_year = year_N, quant = 0.95)[[
                                    "smooth_endpoint_trends"]] %>% 
  rename(stratum = continent) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)
# ------------------------------------------------------------------------------




# stratum index summaries ------------------------------------------------------
# stratum indices
strat_idx_lst <- index_function(fit = fit,
                                parameter = "n",
                                strat = "stratum",
                                year = "Year",
                                first_dim = "s",
                                quant = 0.95,
                                weights_df = NULL,
                                area = "Area",
                                summary_regions = "stratum",
                                year_1 = 1966,
                                to_summarise = F)

# add smooth indices
strat_smooth_idxs <- get_smooth_ns(idx_list = strat_idx_lst)
strat_idx_lst$indices_smooth <- strat_smooth_idxs$indices_smooth
strat_idx_lst$samples_smooth <- strat_smooth_idxs$samples_smooth

# plot indices for full time series
pd1 <- strat_idx_lst$indices_smooth %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(stratum, strata_name), by="stratum") %>% 
  arrange(stratum, Year)
pd2 <- strat_idx_lst$indices %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(stratum, strata_name), by="stratum") %>% 
  arrange(stratum, Year)
regs <- sort(unique(pd2$stratum))
pgs <- 1:ceiling(length(regs) / 9)

# loop through plot pages
for(pg in 1:length(pgs)){
  stt <- pg
  if(stt == 1) stt <- 1
  if(stt > 1) stt <- ((pg - 1) * 9) + 1
  edd <- stt + 9 - 1
  fst <- regs[stt]
  ltt <- regs[edd]
  if(pg == max(pgs)) ltt <- last(regs)
  fn <- paste0("./output/", gsub(" ", "_", species),
               "_annual_indices_stratum_", fst, "-", ltt, ".pdf")
  plt_i <- ggplot() +
    geom_ribbon(data = pd1, 
                aes(x = true_year, ymin = lci, ymax = uci),
                alpha = 0.3, fill="darkblue") +
    geom_line(data = pd1, 
              aes(x = true_year, y = mean), color="darkblue") +
    geom_linerange(data = pd2, 
                   aes(x = true_year, ymin = lci, ymax = uci),
                   alpha = 0.6, linewidth=0.5) +
    geom_point(data = pd2, 
               aes(x = true_year, y = mean), alpha = 0.6, shape=16) +
    ggforce::facet_wrap_paginate(~ strata_name, scales="free",
                                 ncol = 3, nrow = 3, page = pg) +
    scale_x_continuous(breaks=seq(min(pd1$true_year), 
                                  max(pd1$true_year), 
                                  7)) +
    labs(x="Year", y="Abundance index (95% CrI) per year and stratum") +
    theme_bw()
  cairo_pdf(fn, width = 11, height = 8.5)
  print(plt_i)
  dev.off() 
}

# turn indices into lt trends
strat_lt_trds <- trends_function(ind_list = strat_idx_lst, start_year = year_1, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
  mutate(stratum = as.numeric(stratum)) %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(stratum, strata_name), by="stratum") %>% 
  select(stratum=strata_name, trend, lci, uci, q0.5,
         q0.05, q0.95, q0.165, q0.835,
         percent_change, p_ch_lci, 
         p_ch_uci, prob_decline, prob_decline_gt30, prob_decline_gt50,
         prob_decline_gt70, region_type, start_year, end_year) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into medium term trends
strat_exp_trds <- trends_function(ind_list = strat_idx_lst, start_year = year_exp, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
  mutate(stratum = as.numeric(stratum)) %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(stratum, strata_name), by="stratum") %>% 
  select(stratum=strata_name, trend, lci, uci, q0.5,
         q0.05, q0.95, q0.165, q0.835,
         percent_change, p_ch_lci, 
         p_ch_uci, prob_decline, prob_decline_gt30, prob_decline_gt50,
         prob_decline_gt70, region_type, start_year, end_year) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 10 yr trends
strat_10yr_trds <- trends_function(ind_list = strat_idx_lst, start_year = year_10, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
  mutate(stratum = as.numeric(stratum)) %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(stratum, strata_name), by="stratum") %>% 
  select(stratum=strata_name, trend, lci, uci, q0.5, 
         q0.05, q0.95, q0.165, q0.835,
         percent_change, p_ch_lci, 
         p_ch_uci, prob_decline, prob_decline_gt30, prob_decline_gt50,
         prob_decline_gt70, region_type, start_year, end_year) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# turn indices into 3 gen trends
strat_3gen_trds <- trends_function(ind_list = strat_idx_lst, start_year = year_3g, 
                                 end_year = year_N, quant = 0.95)[[
                                   "smooth_endpoint_trends"]] %>% 
  mutate(stratum = as.numeric(stratum)) %>% 
  left_join(strata_map %>% st_drop_geometry() %>% 
              select(stratum, strata_name), by="stratum") %>% 
  select(stratum=strata_name, trend, lci, uci, q0.5, 
         q0.05, q0.95, q0.165, q0.835,
         percent_change, p_ch_lci, 
         p_ch_uci, prob_decline, prob_decline_gt30, prob_decline_gt50,
         prob_decline_gt70, region_type, start_year, end_year) %>% 
  mutate(strata_name = stratum,
         sig = ifelse(lci>0 | uci<0, 1, 0),
         gen_3_years = gen_3_years)

# plot bcr trends
tp_lt <- map_function(trds = strat_lt_trds, spunit = "bbs_usgs")
tp_exp <- map_function(trds = strat_exp_trds, spunit = "bbs_usgs")
tp_10yr <- map_function(trds = strat_10yr_trds, spunit = "bbs_usgs")
tp_3gen <- map_function(trds = strat_3gen_trds, spunit = "bbs_usgs")

# all together
ptch1 <- ((tp_lt + tp_exp) / (tp_10yr + tp_3gen))
ptch1 <- ptch1 + plot_layout(guides = 'collect')
cairo_pdf(paste0("./output/", gsub(" ", "_", species), "_stratum_trend_map.pdf"),  
          width = 9.25, height = 10)
ptch1
dev.off()
# ------------------------------------------------------------------------------





# aggregate indices and trends -------------------------------------------------
# get all trends
agg_trends <- rbind(cont_lt_trds, cont_exp_trds, cont_10yr_trds, cont_3gen_trds,
      cntry_lt_trds, cntry_exp_trds, cntry_10yr_trds, cntry_3gen_trds,
  bcr_lt_trds, bcr_exp_trds, bcr_10yr_trds, bcr_3gen_trds,
  ps_lt_trds, ps_exp_trds, ps_10yr_trds, ps_3gen_trds,
  strat_lt_trds, strat_exp_trds, strat_10yr_trds, strat_3gen_trds) %>%
  mutate(species = species, 
         trend_type = ifelse(start_year == year_1, "long-term",
                             ifelse(start_year == year_exp, "medium-term",
                                    ifelse(start_year == year_10, "ten-year",
                                           "three-generation")))) %>% 
  select(species, region = strata_name, region_type,
         trend_type, length_3_gens=gen_3_years, 
         year_start = start_year, year_end = end_year,
         trend_mean = trend, trend_median = q0.05, 
         trend_q0.025 = lci, trend_q0.975 = uci, 
         trend_diff_zero_95 = sig, 
         trend_q0.05 = q0.05, trend_q0.95 = q0.95, 
         trend_q0.165 = q0.165, trend_q0.835 = q0.835,
         percent_change_mean = percent_change, 
         percent_change_q0.025 = p_ch_lci, percent_change_q0.975 = p_ch_uci,
         prob_decline, prob_decline_gt30, prob_decline_gt50, prob_decline_gt70)

# get all indices and smooth indices
agg_indices <- rbind(cont_idx_lst$indices %>% rename(region=continent) %>% 
                       mutate(index_type = "first_diff"), 
                     cont_idx_lst$indices_smooth %>% rename(region=continent) %>% 
                       mutate(index_type = "smooth_first_diff"),
                     cntry_idx_lst$indices %>% rename(region=country) %>% 
                       mutate(index_type = "first_diff"), 
                     cntry_idx_lst$indices_smooth %>% rename(region=country) %>% 
                       mutate(index_type = "smooth_first_diff"),
                     bcr_idx_lst$indices %>% rename(region=bcr) %>% 
                       mutate(index_type = "first_diff",
                              region = paste0("BCR", 
                                              str_pad(region, width = 2, 
                                                      side = "left", pad = "0"))), 
                     bcr_idx_lst$indices_smooth %>% rename(region=bcr) %>% 
                       mutate(index_type = "smooth_first_diff",
                              region = paste0("BCR", 
                                              str_pad(region, width = 2, 
                                                      side = "left", pad = "0"))),
                     ps_idx_lst$indices %>% rename(region=prov_state) %>% 
                       mutate(index_type = "first_diff"), 
                     ps_idx_lst$indices_smooth %>% rename(region=prov_state) %>% 
                       mutate(index_type = "smooth_first_diff"),
                     strat_idx_lst$indices %>% 
                       left_join(strata_map %>% st_drop_geometry() %>% 
                                   select(stratum, strata_name), by="stratum") %>% 
                       rename(region=stratum) %>% 
                       mutate(region = strata_name, index_type = "first_diff") %>% 
                       select(-strata_name), 
                     strat_idx_lst$indices_smooth %>% 
                       left_join(strata_map %>% st_drop_geometry() %>% 
                                   select(stratum, strata_name), by="stratum") %>% 
                       rename(region=stratum) %>% 
                       mutate(region = strata_name, index_type = "smooth_first_diff") %>% 
                       select(-strata_name)) %>% 
  rename(count_year = Year) %>% 
  mutate(species = species, 
         count_number = count_year + year_1 - 1900,
         year_end = true_year + 1) %>% 
  select(species, region, region_type, index_type, 
         count_number, year_start = true_year, year_end = year_end, 
         index_mean = mean, index_median = median, 
         index_q0.025 = lci, index_q0.975 = uci, 
         index_q0.05 = q0.05, index_q0.95 = q0.95, index_q0.165 = q0.165, 
         index_q0.835 = q0.835)

# save
write.csv(agg_indices, paste0("./output/", gsub(" ", "_", species), 
                              "_relative_abundance_indices.csv"), 
          na="", row.names=F)
write.csv(agg_trends, paste0("./output/", gsub(" ", "_", species), 
                             "_relative_abundance_trends.csv"),
          na="", row.names=F)
# ------------------------------------------------------------------------------

