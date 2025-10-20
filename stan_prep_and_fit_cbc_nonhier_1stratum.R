# setup ------------------------------------------------------------------------
library(tidyverse)
library(bbsBayes2)
library(cmdstanr)
library(patchwork)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# ------------------------------------------------------------------------------



# select species and load data -------------------------------------------------
species <- "Cinclus_mexicanus"
dir("./data")
species_f <- species
data_1 <- read.csv(paste0("data/",species,"_modeled_records.csv"))
data_1 <- data_1 %>% 
  #filter(scaled_effort < 10) %>% # this is to drop super high effort counts
  mutate(strata_name = paste(country,state,bcr,sep = "-")) # match BBS strata

# optional removal of low abundance strata
strat_means <- data_1 %>% 
  group_by(strata_name) %>% 
  summarise(mean_obs = mean(how_many),
            lmean = log(mean_obs)) %>% 
  arrange(-(lmean))
# ------------------------------------------------------------------------------

data_1 <- data_1 %>% 
  filter(strata_name %in% strat_means$strata_name[1:2])
 
# get, tweak and view strata df ------------------------------------------------
strata_information <- data_1 %>% 
  select(strata_name,
         circles_per_stratum,nonzero_circles) %>% 
  distinct() %>% 
    mutate(non_zero = nonzero_circles/circles_per_stratum,
           strata_vec = as.integer(factor(strata_name)))


nstrata <- 2
nonzeroweight <- as.numeric(strata_information$non_zero)
# ------------------------------------------------------------------------------


# join strata back to data table and drop unnecessary columns
data_prep <- data_1 %>% 
  select(circle,count_year,
         lon,lat,how_many, field_hours,
         scaled_effort,strata_name) %>% 
  left_join(., strata_information,
            by = c("strata_name")) %>% 
  mutate(circle_vec = as.integer(factor(circle)),
         year_vec = count_year - (min(count_year)-1))

saveRDS(data_prep,paste0("data/data_prep_1_stratum_", species, ".rds"))
# ------------------------------------------------------------------------------



# circle per strata work -------------------------------------------------------
nsites <- max(data_prep$circle_vec)

# list of site and strat combos
sites_df <- data_prep %>% 
  select(strata_vec,circle_vec) %>% 
  distinct() %>% 
  arrange(strata_vec,
          circle_vec) 

# number of sites in each stratum
nsites_strata <- data_prep %>% 
  select(strata_vec,circle_vec) %>% 
  distinct() %>% 
  arrange(strata_vec,
          circle_vec) %>% 
  group_by(strata_vec) %>% 
  summarise(nsites = n())

nsites_strata <- as.integer(nsites_strata$nsites)
maxnsites_strata <- max(nsites_strata)

# matrix of which sites are in which strata
ste_mat <- matrix(data = 0,
                  nrow = nstrata,
                  ncol = maxnsites_strata)
for(i in 1:nstrata){
  ste_mat[i,1:nsites_strata[i]] <- sites_df[which(sites_df$strata_vec == i),
                                            "circle_vec"]
}
# ------------------------------------------------------------------------------



# data list for stan -----------------------------------------------------------
ncounts <- nrow(data_prep)
nyears <- max(data_prep$year_vec)
count <- data_prep$how_many
strat <- data_prep$strata_vec
year <- data_prep$year_vec
site <- data_prep$circle_vec
hours <- data_prep$scaled_effort
field_hours <- data_prep$field_hours

stan_data <- list(# scalar indicators
  nsites = nsites,
  nstrata = nstrata,
  ncounts = ncounts,
  nyears = nyears,
  
  # basic data
  count = count,
  strat = strat,
  year = year,
  site = site,
  
  # effort information
  hours = hours,
  
  # ragged array information to link sites to strata
  nsites_strata = nsites_strata,
  maxnsites_strata = maxnsites_strata,
  ste_mat = ste_mat,

  # weights
  nonzeroweight = nonzeroweight
)
# ------------------------------------------------------------------------------




stan_data[["fixed_year"]] <- floor(stan_data$nyears/2)
stan_data[["zero_betas"]] <- rep(0,stan_data$nstrata)
stan_data[["Iy1"]] <- c((stan_data$fixed_year-1):1)
stan_data[["nIy1"]] <- length(stan_data[["Iy1"]])
stan_data[["Iy2"]] <- c((stan_data$fixed_year+1):stan_data$nyears)
stan_data[["nIy2"]] <- length(stan_data[["Iy2"]])


# Add additional effort visualisation variables
stan_data[["neffort_preds"]] <- 100
stan_data[["effort_preds"]] <- c(seq(from = min(stan_data$hours), 
                                     to = max(stan_data$hours),
                                     length.out = 100))

# get model file
#mod.file <- "models/first_diff_nonhier_CBC_nonnonhier_effort.stan"

#mod.file <- "models/first_diff_nonhier_CBC.stan"
mod.file <- "models/first_diff_nonhier_CBC.stan"

# compile model
#set_cmdstan_path(path = "D:/Users/tmeehan/Documents/.cmdstan/cmdstan-2.35.0")
model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))
# ------------------------------------------------------------------------------



# run stan model ---------------------------------------------------------------

# start sampling
stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=4, 
  iter_sampling=200,
  iter_warmup=500,
  parallel_chains = 4,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 11,
  #seed = 123,
  init = 1,
  show_exceptions = TRUE)

# save cmdstan objects
summ <- stanfit$summary()

stanfit$save_object(paste0("output/fit_",species,"_CBC_nonhier_first_diff.rds"))
saveRDS(stan_data, paste0("output/datalist_",species,"_CBC_nonhier_first_diff.rds"))
saveRDS(summ, paste0("output/parameter_summary_",species,"_CBC_nonhier_first_diff.rds"))


tmp <- summ %>% filter(grepl("n[",variable,fixed = TRUE))
tmp$year <- rep(1:stan_data$nyears+(min(data_1$count_year)-1),
                each = stan_data$nstrata)
tmp$strata_vec <- rep(1:stan_data$nstrata,
                times = stan_data$nyears)


obs <- data_prep %>% 
  group_by(count_year,
           strata_vec) %>% 
  summarise(obs_mean = mean(how_many))

tmp <- tmp %>% 
  left_join(obs,by = c("year" = "count_year",
                       "strata_vec"))

tst <- ggplot(data = tmp,
              aes(x = year,y = median))+
  geom_ribbon(aes(ymin = q5,ymax = q95),
              alpha = 0.3)+
  geom_line()+
  geom_point(aes(x = year, y = obs_mean))+
  facet_wrap(vars(strata_vec))
tst



effort <- summ %>% filter(grepl("effort_vis[",variable,fixed = TRUE))
effort$hours <- rep(c(seq(from = min(stan_data$hours), 
                          to = max(stan_data$hours),
                          length.out = 100)),
                each = stan_data$nstrata)
effort$strata_vec <- rep(1:stan_data$nstrata,
                      times = stan_data$neffort_preds)




tst <- ggplot(data = effort,
              aes(x = hours,
                  y = mean))+
  geom_ribbon(aes(ymin = q5,ymax = q95),
              alpha = 0.3)+
  geom_line()+
  facet_wrap(vars(strata_vec))

tst


