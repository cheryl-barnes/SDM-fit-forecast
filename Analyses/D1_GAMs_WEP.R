# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Climate-informed models benefit hindcasting but may present challenges when forecasting species-habitat associations. Ecogr. 

# References: 
# Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.

# Overarching Objective: Identify best-fit generalized additive models (GAMs) for hindcasting (based on conventional statistics) and evaluate forecast skill (based on retrospective skill testing; Thorson 2019) for presence-absence, numerical abundance, and biomass of groundfishes in the Bering Sea (1982-2018). Compare species distribution models (SDMs) with varying degrees of complexity to assess whether the addition of time-varying processes to status quo static SDMs improves hindcast performance and/or forecast skill.

# Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

# Simple Dynamic GAMs, Walleye Pollock (WEP):
setwd("~/Documents/JISAO/SDM-fit-forecast/")

require(tidyverse)
require(tidyr)
require(dplyr)
require(mgcv)
require(MuMIn)
  options(na.action = "na.fail")
require(lmtest)
require(mapdata)
require(rgdal)
require(rgeos)
require(sf)
require(ggplot2)
require(visreg)
require(rnaturalearth)
require(rnaturalearthdata)
require(rnaturalearthhires)
require(raster)
require(beepr)
require(Metrics)
    
# Load survey data linked to model covariates:
load("Data/dynamic_model_data.rda")

#########################################################################
# Presence-Absence:
dynamic.model.data$Year = as.numeric(dynamic.model.data$Year)
WEP = dynamic.model.data %>%
  filter(Species == "Walleye.Pollock")

# Full model, year = fixed effect:
d1GAM.WEPpa.full = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4) + s(phi, bs='tp', k=4) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp")
  summary(d1GAM.WEPpa.full)
  save(d1GAM.WEPpa.full, file = "Analyses/d1GAM_WEPpa_full.rda")
# load("Analyses/d1GAM_WEPpa_full.rda")
  
# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d1GAM.WEPpa.check.full = round(as.data.frame(concurvity(d1GAM.WEPpa.full, full=F)), digits=3); d1GAM.WEPpa.check.full
# bathy-BPI > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
d1GAM.WEPpa.ind.cov.full = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp")  
  summary(d1GAM.WEPpa.ind.cov.full)
# R-sq = 0.428; Dev. = 41.2%; UBRE = -0.52077
  save(d1GAM.WEPpa.ind.cov.full, file = "Analyses/d1GAM_WEPpa_ind_cov.rda")
# load("Analyses/d1GAM_WEPpa_ind_cov.rda")
  
# Generate alternative models using backwards step-wise selection:
# Remove sponges (p > 0.1):
d1GAM.WEPpa.alt.1 = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d1GAM.WEPpa.alt.1)
# R-sq = 0.428; Dev. = 41.2%; UBRE = -0.52089
  save(d1GAM.WEPpa.alt.1, file = "Analyses/d1GAM_WEPpa_alt1.rda")
# load("Analyses/d1GAM_WEPpa_alt1.rda")

# Remove whips (p > 0.1):
d1GAM.WEPpa.alt.2 = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1)+ s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d1GAM.WEPpa.alt.2)
# R-sq = 0.428; Dev. = 41.1%; UBRE = -0.52086
  save(d1GAM.WEPpa.alt.2, file = "Analyses/d1GAM_WEPpa_alt2.rda")
# load("Analyses/d1GAM_WEPpa_alt2.rda")

# Save best-fit model based on GCV/UBRE:  
d1GAM.WEPpa.best = d1GAM.WEPpa.alt.1 
  summary(d1GAM.WEPpa.best)
  save(d1GAM.WEPpa.best, file = "Analyses/d1GAM_WEPpa_best.rda")
# load("Analyses/d1GAM_WEPpa_best.rda")

#########################################################################
# Positive Catches:
WEP_catch = subset(WEP, adjN > 0 & adjWT > 0)

# Abundance #
# Full model:
d1GAM.WEP.N.full = gam(adjN ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4) + s(phi, bs='tp', k=4) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = poisson(link="log"), method = "GCV.Cp")  
  summary(d1GAM.WEP.N.full)
  save(d1GAM.WEP.N.full, file = "Analyses/d1GAM_WEP_N_full.rda")
# load("Analyses/d1GAM_WEP_N_full.rda")

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d1GAM.WEP.N.check = round(as.data.frame(concurvity(d1GAM.WEP.N.full, full=F)), digits=3); d1GAM.WEP.N.check
# bathy-lon/lat > 0.5, exclude bathy
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
d1GAM.WEP.N.ind.cov = gam(adjN ~ te(lon, lat, bs='tp', m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = poisson(link="log"), method = "GCV.Cp") 
  summary(d1GAM.WEP.N.ind.cov)
# adj-R2 = 0.104; Deviance = 27.6%; UBRE = 907.81
  save(d1GAM.WEP.N.ind.cov, file = "Analyses/d1GAM_WEP_N_ind_cov.rda")
# load("Analyses/d1GAM_WEP_N_ind_cov.rda")
  
# Save best-fit model: 
d1GAM.WEP.N.best = d1GAM.WEP.N.ind.cov
  summary(d1GAM.WEP.N.best)
  save(d1GAM.WEP.N.best, file = "Analyses/d1GAM_WEP_N_best.rda")
# load("Analyses/d1GAM_WEP_N_best.rda")
    
#########################################################################
# Biomass #
# Full model:
d1GAM.WEP.WT.full = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4) + s(phi, bs='tp', k=4) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = Gamma(link="log"), method = "GCV.Cp")  
  summary(d1GAM.WEP.WT.full)
  save(d1GAM.WEP.WT.full, file = "Analyses/d1GAM_WEP_WT_full.rda")
# load("Analyses/d1GAM_WEP_WT_full.rda")

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d1GAM.WEP.WT.check = round(as.data.frame(concurvity(d1GAM.WEP.WT.full, full=F)), digits=3); d1GAM.WEP.WT.check
# bathy-lon/lat, exclude BPI
# phi-lon/lat, exclude phi
# phi-CPE, exclude phi
d1GAM.WEP.WT.ind.cov = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = Gamma(link="log"), method = "GCV.Cp") 
  summary(d1GAM.WEP.WT.ind.cov)
# adj-R2 = 0.0838; Deviance = 28.3%; GCV = 2.0123
  save(d1GAM.WEP.WT.ind.cov, file = "Analyses/d1GAM_WEP_WT_ind_cov.rda")
# load("Analyses/d1GAM_WEP_WT_ind_cov.rda")
  
# Save best-fit model: 
d1GAM.WEP.WT.best = d1GAM.WEP.WT.ind.cov
  summary(d1GAM.WEP.WT.best)
  save(d1GAM.WEP.WT.best, file = "Analyses/d1GAM_WEP_WT_best.rda")
# load("Analyses/d1GAM_WEP_WT_best.rda")

#########################################################################
# Predict probability of occurrence, abundance, and biomass:
load("Data/dynamic_model_grid.rda")

# Probability of Occurrence:
WEPpa.pred = predict.gam(d1GAM.WEPpa.best, newdata = dynamic.model.grid,
                         type = "response", se.fit = T) 
WEPpa.predict = cbind(dynamic.model.grid, WEPpa.pred)
WEPpa.predict = WEPpa.predict %>%
  rename(pa.fit = fit, pa.se.fit = se.fit)
  save(WEPpa.predict, file = "Analyses/d1WEP_pa_predict.rda")
# load("Analyses/d1WEP_pa_predict.rda")

# Abundance:
WEP.N.pred = predict.gam(d1GAM.WEP.N.best, newdata = dynamic.model.grid,
                         type = "response", se.fit = T)
WEP.N.predict = cbind(dynamic.model.grid, WEP.N.pred)
WEP.N.predict = WEP.N.predict %>%
  rename(N.fit = fit, N.se.fit = se.fit)
  save(WEP.N.predict, file = "Analyses/d1WEP_N_predict.rda")
# load("Analyses/d1WEP_N_predict.rda")

# Biomass:
WEP.WT.pred = predict.gam(d1GAM.WEP.WT.best, newdata = dynamic.model.grid,
                          type = "response", se.fit = T)
WEP.WT.predict = cbind(dynamic.model.grid,WEP.WT.pred)
WEP.WT.predict = WEP.WT.predict %>%
  rename(WT.fit = fit, WT.se.fit = se.fit)
  save(WEP.WT.predict, file = "Analyses/d1WEP_WT_predict.rda")
# load("Analyses/d1WEP_WT_predict.rda")

# Merge predictions:
WEP.predict = WEPpa.predict %>% 
  full_join(WEP.N.predict) %>%
  full_join(., WEP.WT.predict)

# Estimate abundance and biomass, accounting for probability of occurrence:
WEP.predict$WEP_Abun = with(WEP.predict, pa.fit * N.fit)
WEP.predict$WEP_Bio = with(WEP.predict, pa.fit * WT.fit)
  save(WEP.predict, file = "Analyses/d1WEP_predict.rda")
# load("Analyses/d1WEP_predict.rda")  
  
##################################################################  
# Estimate Spearman's correlation coefficients for hindcasts (tendency to co-vary; range: -1 to 1):

# Merge data from model fitting and predicting:
WEP_sf = st_as_sf(WEP, coords = c("lon", "lat"))
  st_crs(WEP_sf) = 4326
WEP.predict_sf = st_as_sf(WEP.predict, coords = c("lon", "lat"))
  st_crs(WEP.predict_sf) = 4326
  
d1WEP.obs.pred = list()
for(i in unique(WEP_sf$Year)) {
  OBS = subset(WEP_sf, Year == i)
  PRED = subset(WEP.predict_sf, Year == i)
  d1WEP.obs.pred[[i]] = st_join(OBS, PRED, join = st_nearest_feature, left = T)}
d1WEP_full = as.data.frame(dplyr::bind_rows(d1WEP.obs.pred))
d1WEP_full = d1WEP_full %>%
  rename(Year = Year.x)
  save(d1WEP_full, file = "Analyses/d1WEP_obs_pred.rda")
# load("Analyses/d1WEP_obs_pred.rda")

# correlation coefficients (tendency to co-vary, range: -1 to 1):
# presence-absence:
pa.cor.S = cor.test(d1WEP_full$PA, d1WEP_full$pa.fit, method = c("spearman")); pa.cor.S
  # S = 3.0011e+11; p < 0.001; rho = 0.5053377 

pa.cor.yr = NULL
for(i in unique(d1WEP_full$Year)) {
  cor.data = subset(d1WEP_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  pa.cor.yr = rbind(pa.cor.yr, coef) }

pa.cor.yr = as.data.frame(pa.cor.yr)
  colnames(pa.cor.yr) = c("Year", "S", "p-value", "rho")
pa.cor.yr[, c("Year", "S", "rho")] = 
  lapply(pa.cor.yr[, c("Year", "S", "rho")], as.numeric)  

# abundance:
N.cor.S = cor.test(d1WEP_full$adjN, d1WEP_full$WEP_Abun, method = c("spearman")); N.cor.S
# S = 2.0133e+11; p < 0.001; rho = 0.6681532 

N.cor.yr = NULL
for(i in unique(d1WEP_full$Year)) {
  cor.data = subset(d1WEP_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  N.cor.yr = rbind(N.cor.yr, coef) }

N.cor.yr = as.data.frame(N.cor.yr)
  colnames(N.cor.yr) = c("Year", "S", "p-value", "rho")
N.cor.yr[, c("Year", "S", "rho")] = 
  lapply(N.cor.yr[, c("Year", "S", "rho")], as.numeric)  

# biomass:
WT.cor.S = cor.test(d1WEP_full$adjWT, d1WEP_full$WEP_Bio, method = c("spearman")); WT.cor.S
# S = 2.2149e+11; p < 0.001; rho = 0.63492 

WT.cor.yr = NULL
for(i in unique(d1WEP_full$Year)) {
  cor.data = subset(d1WEP_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  WT.cor.yr = rbind(WT.cor.yr, coef) }

WT.cor.yr = as.data.frame(WT.cor.yr)
  colnames(WT.cor.yr) = c("Year", "S", "p-value", "rho")
WT.cor.yr[, c("Year", "S", "rho")] = 
  lapply(WT.cor.yr[, c("Year", "S", "rho")], as.numeric)  

#################################################################
# Retrospective skill testing:
load("Data/dyn_forecasts_stack_10km.rda")

# presence-absence:
WEP = WEP %>%
  arrange(Year)

pa.predictions_all = NULL
for(i in 1992:2018) {
  pa.data = subset(WEP, Year < i) # for model fitting (i = first year forecasted)
  
pa.model = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = pa.data, family = binomial(link="cloglog"), method = "GCV.Cp")   

pa.dyn.cov = subset(d.forecasts.10km, Year >= i & Year <= i + 21)
pa.pred = predict.gam(pa.model, newdata = pa.dyn.cov, 
                      type="response", se.fit=T)
pa.dyn.cov_pers = subset(d.forecasts.10km, Year == i - 1)
    pa.pers = predict.gam(pa.model, newdata = pa.dyn.cov_pers, 
                          type="response", se.fit=T)
pa.pers = as.data.frame(pa.pers)
pa.pers = pa.pers %>% rename(fit_pa.pers = fit, se.fit_pa.pers = se.fit)
  pa.obs.pred = cbind(pa.dyn.cov, pa.pred, pa.pers)
pa.obs.pred$forecast.yr = pa.obs.pred$Year
  pa.obs.pred$fitted.through = i - 1
pa.predictions_all = rbind(pa.predictions_all, pa.obs.pred) }

pa.predictions_all = pa.predictions_all %>%
  rename(pa.fit = fit, pa.se.fit = se.fit)
  save(pa.predictions_all, file = "Analyses/d1WEP_pa_retro.rda") 
# load("Analyses/d1WEP_pa_retro.rda")
  
# positive catches:
WEP_catch = WEP_catch %>%
  arrange(Year)

# abundance:
N.predictions_all = NULL
for(i in 1992:2018) {
  N.data = subset(WEP_catch, Year < i) # for model fitting (i = first year forecasted)
  
N.model = gam(adjN ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = N.data, family = poisson(link="log"), method = "GCV.Cp") 

N.dyn.cov = subset(d.forecasts.10km, Year >= i & Year <= i + 21)
N.pred = predict.gam(N.model, newdata = N.dyn.cov, 
                     type="response", se.fit=T)
  N.dyn.cov_pers = subset(d.forecasts.10km, Year == i - 1)
N.pers = predict.gam(N.model, newdata = N.dyn.cov_pers, 
                     type="response", se.fit=T)
N.pers = as.data.frame(N.pers)
N.pers = N.pers %>% rename(fit_N.pers = fit, se.fit_N.pers = se.fit)
  N.obs.pred = cbind(N.dyn.cov, N.pred, N.pers)
N.obs.pred$forecast.yr = N.obs.pred$Year
  N.obs.pred$fitted.through = i - 1
N.predictions_all = rbind(N.predictions_all, N.obs.pred)}

N.predictions_all = N.predictions_all %>%
  rename(N.fit = fit, N.se.fit = se.fit)
  save(N.predictions_all, file = "Analyses/d1WEP_N_retro.rda") 
# load("Analyses/d1WEP_N_retro.rda")
  
# biomass:
WT.predictions_all = NULL
for(i in 1992:2018) {
  WT.data = subset(WEP_catch, Year < i) # for model fitting (i = first year forecasted)
  
WT.model = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WT.data, family = Gamma(link="log"), method = "GCV.Cp") 

WT.dyn.cov = subset(d.forecasts.10km, Year >= i & Year <= i + 21)
WT.pred = predict.gam(WT.model, newdata = WT.dyn.cov, 
                      type="response", se.fit=T)
  WT.dyn.cov_pers = subset(d.forecasts.10km, Year == i - 1)
WT.pers = predict.gam(WT.model, newdata = WT.dyn.cov_pers, 
                          type="response", se.fit=T)
WT.pers = as.data.frame(WT.pers)
WT.pers = WT.pers %>% rename(fit_WT.pers = fit, se.fit_WT.pers = se.fit)
  WT.obs.pred = cbind(WT.dyn.cov, WT.pred, WT.pers)
WT.obs.pred$forecast.yr = WT.obs.pred$Year
  WT.obs.pred$fitted.through = i - 1
WT.predictions_all = rbind(WT.predictions_all, WT.obs.pred)}

WT.predictions_all = WT.predictions_all %>%
  rename(WT.fit = fit, WT.se.fit = se.fit)
  save(WT.predictions_all, file = "Analyses/d1WEP_WT_retro.rda") 
# load("Analyses/d1WEP_WT_retro.rda")

# all models:
pred.pos_all = N.predictions_all %>% full_join(WT.predictions_all)
  save(pred.pos_all, file = "Analyses/d1WEP_retro_pos.rda") 
# load("Analyses/d1WEP_retro_pos.rda")

pred_all = pa.predictions_all %>% full_join(pred.pos_all)

# Estimate abundance and biomass, accounting for probability of occurrence (estimate / AreaSwept => N or kg per sq. km):
pred_all$WEP_Abun = with(pred_all, (pa.fit * N.fit))
pred_all$WEP_Bio = with(pred_all, (pa.fit * WT.fit))
  pred_all$Abun.thousands = pred_all$WEP_Abun / 1000
  pred_all$Bio.MT = pred_all$WEP_Bio / 1000

# persistence
pred_all$WEP_Abun.pers = with(pred_all, (pa.fit * fit_N.pers))
  pred_all$Abun_pers.thousands = pred_all$WEP_Abun.pers / 1000
pred_all$WEP_Bio.pers = with(pred_all, (pa.fit * fit_WT.pers))
  pred_all$Bio_pers.MT = pred_all$WEP_Bio.pers / 1000 
  
  save(pred_all, file = "Analyses/d1WEP_retro.rda") 
# load("Analyses/d1WEP_retro.rda")
  
#################################################################
# Estimate Spearman's correlation coefficients for forecasts (tendency to co-vary; range: -1 to 1):
WEP_sf = st_as_sf(WEP, coords = c("lon", "lat"))
  st_crs(WEP_sf) = 4326 # WGS84
pred_all_sf = st_as_sf(pred_all, coords = c("lon", "lat"))
  st_crs(pred_all_sf) = 4326 # WGS84
  
# Merge data from model fitting and predicting:
d1WEP.obs.pred_all = NULL
for(model in unique(pred_all_sf$fitted.through)) {
  d1WEP.pred = subset(pred_all_sf, fitted.through == model)
   for(yr in unique(d1WEP.pred$Year)) {
     d1WEP.pred.yr = subset(d1WEP.pred, Year == yr)
     d1WEP.obs.yr = subset(WEP_sf, Year == yr)
     d1WEP.obs.pred_yr = st_join(d1WEP.obs.yr, d1WEP.pred.yr, join = st_nearest_feature, left=T)
d1WEP.obs.pred_all = rbind(d1WEP.obs.pred_all, d1WEP.obs.pred_yr)} }

d1WEP_retro_all = as.data.frame(d1WEP.obs.pred_all)
  save(d1WEP_retro_all, file = "Analyses/d1WEP_obs_pred_retro.rda")
# load("Analyses/d1WEP_obs_pred_retro.rda")

# Estimate Spearman's correlation coefficient and RMSE:
d1WEP_retro_all$fitted.through = as.factor(d1WEP_retro_all$fitted.through)
  format.pval(c(0.1, 0.0001, 1e-27))

d1WEP_retro_all = d1WEP_retro_all %>%
  rename(Year = Year.x)

# presence-absence:
pa.obs.pred.cor_all = NULL
for(model in unique(d1WEP_retro_all$fitted.through)) {
  cor.pa.data = subset(d1WEP_retro_all, fitted.through == model)
for(yr in unique(cor.pa.data$Year)) {
  cor.pa.yr.data = subset(cor.pa.data, Year == yr) 
    corr.pa.test.s = cor.test(cor.pa.yr.data$PA, cor.pa.yr.data$pa.fit, 
                              method = "spearman")
    corr.pa.test.s_pers = cor.test(cor.pa.yr.data$PA, cor.pa.yr.data$fit_pa.pers, 
                              method = "spearman")
  coef.pa = cbind(model, yr, 
               corr.pa.test.s$statistic, 
               format.pval(corr.pa.test.s$p.value, digits=5, eps=0.001), 
                  corr.pa.test.s$estimate, 
               corr.pa.test.s_pers$statistic, 
               format.pval(corr.pa.test.s_pers$p.value, digits=5, eps=0.001), 
                  corr.pa.test.s_pers$estimate)
pa.obs.pred.cor_all = rbind(pa.obs.pred.cor_all, coef.pa) }}

pa.obs.pred.cor_all = as.data.frame(pa.obs.pred.cor_all) 
colnames(pa.obs.pred.cor_all) = c("fitted.through", "Year", "S", "p-value.s", "rho", "S_pers", "p-value.s_pers", "rho_pers")

pa.obs.pred.cor_all = pa.obs.pred.cor_all %>%
  mutate(fitted.through = as.numeric(fitted.through)) %>%
  arrange(desc(fitted.through))
pa.obs.pred.cor_all$Metric = "pa"

# positive catches (accounting for probability of occurrence):
# abundance:
N.obs.pred.cor_all = NULL
for(model in unique(d1WEP_retro_all$fitted.through)) {
  cor.N.data = subset(d1WEP_retro_all, fitted.through == model)
for(yr in unique(cor.N.data$Year)) {
  cor.N.yr.data = subset(cor.N.data, Year == yr) # various correlation coefficients
    corr.N.test.s = cor.test(cor.N.yr.data$adjN, cor.N.yr.data$WEP_Abun, 
                              method = "spearman")
    corr.N.test.s_pers = cor.test(cor.N.yr.data$adjN, cor.N.yr.data$WEP_Abun.pers, 
                              method = "spearman")
  coef.N = cbind(model, yr, 
               corr.N.test.s$statistic,
               format.pval(corr.N.test.s$p.value, digits=5, eps=0.001), 
                  corr.N.test.s$estimate, 
               corr.N.test.s_pers$statistic,
               format.pval(corr.N.test.s_pers$p.value, digits=5, eps=0.001), 
                  corr.N.test.s_pers$estimate)
N.obs.pred.cor_all = rbind(N.obs.pred.cor_all, coef.N) }}

N.obs.pred.cor_all = as.data.frame(N.obs.pred.cor_all) 
colnames(N.obs.pred.cor_all) = c("fitted.through", "Year", "S", "p-value.s", "rho", "S_pers", "p-value.s_pers", "rho_pers")

N.obs.pred.cor_all = N.obs.pred.cor_all %>%
  mutate(fitted.through = as.numeric(fitted.through)) %>%
  arrange(desc(fitted.through))
N.obs.pred.cor_all$Metric = "N"

# biomass:
WT.obs.pred.cor_all = NULL
for(model in unique(d1WEP_retro_all$fitted.through)) {
  cor.WT.data = subset(d1WEP_retro_all, fitted.through == model)
for(yr in unique(cor.WT.data$Year)) {
    cor.WT.yr.data = subset(cor.WT.data, Year == yr) # various correlation coefficients
    corr.WT.test.s = cor.test(cor.WT.yr.data$adjWT, cor.WT.yr.data$WEP_Bio, 
                              method = "spearman")
    corr.WT.test.s_pers = cor.test(cor.WT.yr.data$adjWT, cor.WT.yr.data$WEP_Bio.pers, 
                              method = "spearman")
  coef.WT = cbind(model, yr, 
               corr.WT.test.s$statistic, 
               format.pval(corr.WT.test.s$p.value, digits=5, eps=0.001), 
                  corr.WT.test.s$estimate, 
               corr.WT.test.s_pers$statistic, 
               format.pval(corr.WT.test.s_pers$p.value, digits=5, eps=0.001), 
                  corr.WT.test.s_pers$estimate)
WT.obs.pred.cor_all = rbind(WT.obs.pred.cor_all, coef.WT) }}

WT.obs.pred.cor_all = as.data.frame(WT.obs.pred.cor_all) 
colnames(WT.obs.pred.cor_all) = c("fitted.through", "Year", "S", "p-value.s", "rho", "S_pers", "p-value.s_pers", "rho_pers")

WT.obs.pred.cor_all = WT.obs.pred.cor_all %>%
  mutate(fitted.through = as.numeric(fitted.through)) %>%
  arrange(desc(fitted.through))
WT.obs.pred.cor_all$Metric = "WT"

obs.pred.pa.cor_all = rbind(pa.obs.pred.cor_all, N.obs.pred.cor_all, WT.obs.pred.cor_all)
obs.pred.pa.cor_all$Metric = as.factor(obs.pred.pa.cor_all$Metric)
  obs.pred.pa.cor_all$Metric = ordered(obs.pred.pa.cor_all$Metric, levels = c("pa", "N", "WT"))

obs.pred.pa.cor_all = obs.pred.pa.cor_all %>%
  mutate_at(c("Year", "rho", "rho_pers"), as.numeric)   

obs.pred.pa.cor_all$yr.forecast = obs.pred.pa.cor_all$Year - obs.pred.pa.cor_all$fitted.through
levels(obs.pred.pa.cor_all$Metric) = c("Presence-Absence", "Abundance", "Biomass")
  save(obs.pred.pa.cor_all, file = "Analyses/d1WEP_obs_pred_pa_cor.rda")
# load("Analyses/d1WEP_obs_pred_pa_cor.rda")