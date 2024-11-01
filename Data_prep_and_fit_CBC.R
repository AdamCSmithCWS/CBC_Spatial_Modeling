# setup ------------------------------------------------------------------------
library(tidyverse)
library(bbsBayes2)
library(cmdstanr)
library(patchwork)
setwd("C:/Users/tmeehan/Documents/Github/CBC_Spatial_Modeling")
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
            lmean = log(mean_obs))
# ------------------------------------------------------------------------------


 
# get, tweak and view strata df ------------------------------------------------
strat_df <- data_1 %>% 
  select(strata_name,
         circles_per_stratum,nonzero_circles) %>% 
  distinct() %>% 
    mutate(non_zero = nonzero_circles/circles_per_stratum)

strata_map <- bbsBayes2::load_map(stratify_by = "bbs_usgs") %>% 
  inner_join(.,strat_df,
             by = "strata_name") %>% 
  mutate(strata_vec = as.integer(factor(strata_name))) %>% 
  arrange(strata_vec)

map_view <- ggplot(strata_map)+
  geom_sf(aes(fill = strata_vec))+
  scale_colour_viridis_c(); map_view

nstrata <- max(strata_map$strata_vec)
nonzeroweight <- as.numeric(strata_map$non_zero)
# ------------------------------------------------------------------------------



# build neighbour relationships ------------------------------------------------
source("functions/neighbours_define.R")
neighbours <- neighbours_define(strata_map,
                  species = species,
                  strat_indicator = "strata_vec",
                  plot_dir = "",
                  plot_file = "_strata_map")

N_edges <- neighbours$N_edges
node1 <- neighbours$node1
node2 <- neighbours$node2

strata_information <- strata_map %>% 
  sf::st_set_geometry(.,NULL)

# join strata back to data table and drop unnecessary columns
data_prep <- data_1 %>% 
  select(circle,count_year,
         lon,lat,how_many, field_hours,
         scaled_effort,strata_name) %>% 
  left_join(., strata_information,
            by = c("strata_name")) %>% 
  mutate(circle_vec = as.integer(factor(circle)),
         year_vec = count_year - (min(count_year)-1))
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
  
  # spatial structure
  N_edges = N_edges,
  node1 = node1,
  node2 = node2,
  
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



# set up spatial first difference model ----------------------------------------
stan_data[["N_edges"]] <- N_edges
stan_data[["node1"]] <- node1
stan_data[["node2"]] <- node2
stan_data[["fixed_year"]] <- floor(stan_data$nyears/2)
stan_data[["zero_betas"]] <- rep(0,stan_data$nstrata)
stan_data[["Iy1"]] <- c((stan_data$fixed_year-1):1)
stan_data[["nIy1"]] <- length(stan_data[["Iy1"]])
stan_data[["Iy2"]] <- c((stan_data$fixed_year+1):stan_data$nyears)
stan_data[["nIy2"]] <- length(stan_data[["Iy2"]])

# get model file
mod.file <- "models/first_diff_spatial_CBC.stan"

# compile model
set_cmdstan_path(path = "D:/Users/tmeehan/Documents/.cmdstan/cmdstan-2.35.0")
model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))
# ------------------------------------------------------------------------------



# run stan model ---------------------------------------------------------------
## initial Values (not used?) 
# init_def <- function(){ list(strata_raw = rnorm(nstrata,0,0.1),
#                              STRATA = 0,
#                              sdstrata = runif(1,0.01,0.1),
#                              ste_raw = rnorm(nsites,0,0.1),
#                              sdnoise = runif(1,0.01,0.2),
#                              sdb = runif(1,0.01,0.1),
#                              sdp = runif(1,0.01,0.1),
#                              b_raw = rnorm(nstrata,0,0.01),
#                              p_raw = rnorm(nstrata,0,0.01),
#                              B = 0,
#                              P = 0,
#                              sdste = runif(1,0.01,0.2),
#                              sdbeta = runif(1,0.01,0.1),
#                              sdBETA = runif(1,0.01,0.1),
#                              BETA_raw = rnorm(nyears-1,0,0.1),
#                              beta_raw = matrix(rnorm((nyears-1)*nstrata,0,0.01),nrow = nstrata,ncol = nyears-1))}

# start sampling
stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=4, 
  iter_sampling=500,
  iter_warmup=500,
  parallel_chains = 4,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 11,
  #seed = 123,
  init = 1,
  show_exceptions = FALSE)

# save cmdstan objects
stanfit$save_object("output/fit_CBC_spatial_first_diff.rds")
saveRDS(stan_data, "output/datalist_CBC_spatial_first_diff.rds")
# ------------------------------------------------------------------------------
