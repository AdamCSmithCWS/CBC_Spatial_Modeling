# some trend models in inla for comparison.
# this script uses data from the Data_prep_and_fit_CBC.R script, so run that
# first to get the stan data list.




# prepare data for inla for comparison -----------------------------------------
# new libraries
library(INLA)
library(spdep)
library(brinla)

# get data
real_dat <- data.frame(count=as.integer(stan_data$count),
                       strat_idx1=as.integer(stan_data$strat),
                       year_idx1=as.integer(stan_data$year),
                       site_idx1=as.integer(stan_data$site),
                       log_hrs=log(field_hours))

# make prediction data with na for count and unimportant variables
pred_dat <- expand_grid(strat_idx1=1:stan_data$nstrata,
                        year_idx1=1:stan_data$nyears) %>% 
  mutate(count=NA, site_idx1=NA, log_hrs=log(mean(field_hours))) %>% 
  select(count, strat_idx1, year_idx1, site_idx1, log_hrs)

# merge data
inla_dat <- bind_rows(real_dat, pred_dat) %>% 
  mutate(strat_idx2=strat_idx1,
         strat_idx3=strat_idx1,
         year_idx2=year_idx1,
         year_idx3=year_idx1,
  )

# get row index for prediction data
pred_idxs <- (nrow(real_dat)+1):(nrow(real_dat)+nrow(pred_dat))

# make stratum key
strat_key <- strata_map %>% select(strat_idx1=strata_vec, everything())

# make neighbors
lw1 <- mat2listw(neighbours$adj_matrix, style="B")
nb1 <- lw1$neighbours; str(nb1)
nb2INLA("map.adj", nb1)
g1 <- inla.read.graph(filename = "map.adj")
plot(g1)
# ------------------------------------------------------------------------------



# hier first diff (rw1) model (log hours) in inla for comparison ---------------
# priors
prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
bym_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), 
                  phi = list(prior = "pc", param = c(0.5, 0.5)))

# formula for non-spatial rw1
fnsrw1 <- count ~ 
  # global intercept
  1 + 
  # random site deviation from global intercept
  f(site_idx1, model="iid", constr=T, hyper=prec_prior) +
  # random stratum deviation from global intercept
  f(strat_idx1, model="iid", constr=T, hyper=prec_prior) +
  # global effort effect
  log_hrs + 
  # random stratum deviation from global effort effect
  f(strat_idx2, log_hrs, model="iid", constr=T, hyper=prec_prior) + 
  # rw1 smooth per stratum
  f(year_idx1, model="rw1", constr=T, scale.model=T, hyper=prec_prior,
    replicate=strat_idx3)

# run model
inla_fit1 <- inla(fnsrw1, family = "nbinomial", data=inla_dat,
                  control.predictor = list(compute=T, link=1),
                  control.compute=list(waic=T),
                  verbose=TRUE)
summary(inla_fit1)
bri.hyperpar.summary(inla_fit1)

# plot trajectory per stratum
inla_dat %>% slice(pred_idxs) %>% 
  mutate(fit=inla_fit1$summary.fitted.values$`0.5quant`[pred_idxs],
         lcl=inla_fit1$summary.fitted.values$`0.025quant`[pred_idxs],
         ucl=inla_fit1$summary.fitted.values$`0.975quant`[pred_idxs]) %>% 
  left_join(strat_key) %>% 
  ggplot(aes(x=year_idx1, y=fit, ymin=lcl, ymax=ucl)) + 
  geom_ribbon(alpha=0.3) +
  geom_line() +
  facet_wrap(~strata_name, scales="free")

# plot effort correction function per stratum
eff_pred <- expand_grid(exponent=inla_fit1$summary.fixed$`0.5quant`[2] + 
                          inla_fit1$summary.random$strat_idx2$`0.5quant`,
                        hours=seq(min(data_1$field_hours), max(data_1$field_hours), 
                                  length.out=50)) %>% 
  mutate(intercept=rep(inla_fit1$summary.fixed$`0.5quant`[1] + 
                         inla_fit1$summary.random$strat_idx1$`0.5quant`, 
                       each=50),
         log_hrs=log(hours),
         effort_correction=exp(intercept+(log_hrs*exponent)),
         stratum=rep(unique(strat_key$strata_name), each=50)) 
eff_pred %>% 
  ggplot(aes(x=hours, y=effort_correction)) +
  geom_line() +
  facet_wrap(~stratum)
# ------------------------------------------------------------------------------



# spatial (bym) first diff (rw1) model (log hours) in inla for comparison ------
# formula for spatial rw1
fsrw1 <- count ~ 
  # global intercept
  1 + 
  # random site deviation from global intercept
  f(site_idx1, model="iid", constr=T, hyper=prec_prior) +
  # global effort effect
  log_hrs + 
  # spatial stratum deviation from global effort effect
  f(strat_idx1, log_hrs, model="bym2", graph=g1, constr=T, hyper=bym_prior) + 
  # spatial rw1 smooth per stratum
  f(strat_idx2, model="bym2", graph=g1, constr=T, scale.model=T, hyper=bym_prior,
    group=year_idx1, control.group=list(model="rw1", hyper=prec_prior))

# run model
inla_fit2 <- inla(fsrw1, family = "nbinomial", data=inla_dat,
                  control.predictor = list(compute=T, link=1),
                  control.compute=list(waic=T),
                  verbose=TRUE)
summary(inla_fit2)
bri.hyperpar.summary(inla_fit2)

# plot trajectory per stratum
inla_dat %>% slice(pred_idxs) %>% 
  mutate(fit1=inla_fit1$summary.fitted.values$`0.5quant`[pred_idxs],
         lcl1=inla_fit1$summary.fitted.values$`0.025quant`[pred_idxs],
         ucl1=inla_fit1$summary.fitted.values$`0.975quant`[pred_idxs]) %>% 
  mutate(fit2=inla_fit2$summary.fitted.values$`0.5quant`[pred_idxs],
         lcl2=inla_fit2$summary.fitted.values$`0.025quant`[pred_idxs],
         ucl2=inla_fit2$summary.fitted.values$`0.975quant`[pred_idxs]) %>% 
  left_join(strat_key) %>% 
  ggplot() + 
  geom_ribbon(aes(x=year_idx1, y=fit1, ymin=lcl1, ymax=ucl1), alpha=0.2, fill="darkred") +
  geom_line(aes(x=year_idx1, y=fit1), col="darkred") +
  geom_ribbon(aes(x=year_idx1, y=fit2, ymin=lcl2, ymax=ucl2), alpha=0.2, fill="darkblue") +
  geom_line(aes(x=year_idx1, y=fit2), col="darkblue") +
  facet_wrap(~strata_name, scales="free") +
  labs(x="Years since 1965", y="Posterior median fit and 95% CrI (red = nonspatial, blue = spatial)")

# waics
inla_fit1$waic$waic
inla_fit2$waic$waic
# ------------------------------------------------------------------------------