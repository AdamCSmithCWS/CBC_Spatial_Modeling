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

# some some helper functions
source("functions/neighbours_define.R")
source("functions/indices_cmdstan.R")
source("functions/posterior_summary_functions.R")
source("functions/map_trends.R")
# ------------------------------------------------------------------------------




# get fit object and look at summaries -----------------------------------------
fit <- readRDS("output/fit_CBC_spatial_first_diff.rds")

# effort par summaries
fit$summary(variables = c("b_raw", "p_raw"), "mean", "sd") %>% View()
mcmc_hist(fit$draws("b_raw"))
mcmc_hist(stanfit$draws("p_raw"))
mcmc_hist(stanfit$draws("sdnoise"))

# plot continent indices
inds <- index_function(fit)
ggplot(data = inds$indices,
                aes(x = Year, y = median))+
  geom_ribbon(aes(ymin = lci, ymax = uci),
              alpha = 0.3)+
  geom_line() +
  labs(y="Estimated annual relative abundance", x="Year")+
  theme_bw()

# plot stratum indices
inds <- index_function(fit, strat = "stratum")
ggplot(data = inds$indices,
                aes(x = Year, y = median))+
  geom_ribbon(aes(ymin = lci, ymax = uci),
              alpha = 0.3)+
  geom_line() +
  facet_wrap(~stratum, scales="free") +
  labs(y="Estimated annual relative abundance", x="Year")+
  theme_bw()

# generate trends DOESN'T WORK
trends <- generate_trends(inds)


