# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Model complexity has contrasting benefits for hindcasting and forecasting species responses to climate change. Ecogr. 

# References: 
# Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.

# Overarching Objective: Identify best-fit generalized additive models (GAMs) for hindcasting (based on conventional statistics) and evaluate forecast skill (based on retrospective skill testing; Thorson 2019) for presence-absence, numerical abundance, and biomass of groundfishes in the Bering Sea (1982-2018). Compare species distribution models (SDMs) with varying degrees of complexity to assess whether the addition of time-varying processes to status quo static SDMs improves hindcast performance and/or forecast skill.

# Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

# Static GAMs, Arrowtooth Flounder (ATF):
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
require(sp)
require(ggplot2)
require(visreg)
require(rnaturalearth)
require(rnaturalearthdata)
require(rnaturalearthhires)
require(raster)
require(Metrics)

# Load survey data linked to model covariates:
load("Data/static_model_data.rda")

##################################################################
# Presence-Absence:
ATF = static.model.data %>%
  filter(Species == "Arrowtooth.Flounder") 

# Full model:
sGAM.ATFpa.full = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(sGAM.ATFpa.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
sGAM.ATFpa.check = round(as.data.frame(concurvity(sGAM.ATFpa.full, full=F)), digits=3); sGAM.ATFpa.check
# bathy-BPI > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# ROMS_BT-lon/lat > 0.5, keep both - interested in both for predictions...
sGAM.ATFpa.ind.cov = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")  
  summary(sGAM.ATFpa.ind.cov)
# adj-R2 = 0.612; Deviance = 55.3%; UBRE = -0.39584

# Generate alternative models using backwards step-wise selection:
  # Remove whips (p > 0.1):
sGAM.ATFpa.alt.1 = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")  
  summary(sGAM.ATFpa.alt.1)
# adj-R2 = 0.612; Deviance = 55.3%; UBRE = -0.3959

# Save best-fit model based on GCV/GCV/UBRE: 
sGAM.ATFpa.best = sGAM.ATFpa.ind.cov
  summary(sGAM.ATFpa.best)

##################################################################
# Positive Catches:
ATF_catch = subset(ATF, adjN > 0 & adjWT > 0)

# Abundance #
# Full model:
sGAM.ATF.N.full = gam(adjN ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = poisson(link="log"), method = "GCV.Cp")  
  summary(sGAM.ATF.N.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
sGAM.ATF.N.check = round(as.data.frame(concurvity(sGAM.ATF.N.full, full=F)), digits=3); sGAM.ATF.N.check
# BPI-lon/lat > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# ROMS_BT-lon/lat > 0.5, keep both - interested in both for predictions...
sGAM.ATF.N.ind.cov = gam(adjN ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, k=4) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = poisson(link="log"), method = "GCV.Cp") 
  summary(sGAM.ATF.N.ind.cov)
# adj-R2 = 0.154; Deviance = 34.0%; UBRE = 22.058

# Save best-fit (full) model:
sGAM.ATF.N.best = sGAM.ATF.N.ind.cov 
  summary(sGAM.ATF.N.best); sum(sGAM.ATF.N.best$edf)
# adj-R2 = 0.154; Deviance = 34.0%; UBRE = 22.058

##################################################################
# Biomass #
# Full model:
sGAM.ATF.WT.full = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")  
  summary(sGAM.ATF.WT.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values > 0.5 * significant concurvity
sGAM.ATF.WT.check = round(as.data.frame(concurvity(sGAM.ATF.WT.full, full=F)), digits=3); sGAM.ATF.WT.check
# BPI-bathy, exclude BPI
# phi-lon/lat, exclude phi
# ROMS_BT-lon/lat, keep both - interested in both for predictions...
sGAM.ATF.WT.ind.cov = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")
  summary(sGAM.ATF.WT.ind.cov)
# adj-R2 = 0.148; Deviance = 37.8%; GCV = 0.93465

# Generate alternative models using backwards step-wise selection:
    # Remove sponges (p > 0.1):
sGAM.ATF.WT.alt_1 = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")
  summary(sGAM.ATF.WT.alt_1)
# adj-R2 = 0.148; Deviance = 37.8%; GCV = 0.93435

# Remove whips (p > 0.1):
sGAM.ATF.WT.alt_2 = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + corals + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")  
  summary(sGAM.ATF.WT.alt_2)
# adj-R2 = 0.148; Deviance = 37.8%; GCV = 0.93415

# Save best-fit model:
sGAM.ATF.WT.best = sGAM.ATF.WT.alt_2  
  summary(sGAM.ATF.WT.best)

##################################################################
# Predict probability of occurrence, abundance, and biomass:
load("Data/static_model_grid.rda")

# Probability of Occurrence:
ATFpa.pred = predict.gam(sGAM.ATFpa.best, newdata = static.model.grid, type = "response", se.fit = T)
ATFpa.predict = cbind(static.model.grid, ATFpa.pred)
ATFpa.predict = ATFpa.predict %>%
  rename(pa.fit = fit, pa.se.fit = se.fit)

# Abundance:
ATF.N.pred = predict.gam(sGAM.ATF.N.best, newdata = static.model.grid, type = "response", se.fit = T)
ATF.N.predict = cbind(static.model.grid, ATF.N.pred)
ATF.N.predict = ATF.N.predict %>%
  rename(N.fit = fit, N.se.fit = se.fit)

# Biomass:
ATF.WT.pred = predict.gam(sGAM.ATF.WT.best, newdata = static.model.grid, type = "response", se.fit = T)
ATF.WT.predict = cbind(static.model.grid, ATF.WT.pred)
ATF.WT.predict = ATF.WT.predict %>%
  rename(WT.fit = fit, WT.se.fit = se.fit)

# Merge predictions:
ATF.predict = ATFpa.predict %>% 
  full_join(ATF.N.predict) %>%
  full_join(., ATF.WT.predict)

# Estimate abundance and biomass, accounting for probability of occurrence (estimate / AreaSwept => N or kg per sq. km):
ATF.predict$ATF_Abun = (with(ATF.predict, pa.fit * N.fit) / ATF.predict$AreaSwept)
ATF.predict$ATF_Bio = (with(ATF.predict, pa.fit * WT.fit) / ATF.predict$AreaSwept)

##################################################################  
# Estimate Spearman's correlation coefficients for hindcasts (tendency to co-vary; range: -1 to 1):

# Merge data from model fitting and predicting:
ATF_sf = st_as_sf(ATF, coords = c("lon", "lat"))
  st_crs(ATF_sf) = 4326
ATF.predict_sf = st_as_sf(ATF.predict, coords = c("lon", "lat"))
  st_crs(ATF.predict_sf) = 4326
sATF_full = st_join(ATF_sf, ATF.predict_sf, join = st_nearest_feature, left = T)

# presence-absence:
pa.cor.S = cor.test(as.numeric(sATF_full$PA), 
                    sATF_full$pa.fit, method = c("spearman")); pa.cor.S
  # S = 1.4054e+11; p < 0.001; rho = 0.7509014 

# abundance:
N.cor.S = cor.test(sATF_full$adjN, 
                   sATF_full$ATF_Abun, method = c("spearman")); N.cor.S
  # S = 1.3019e+11; p < 0.001; rho = 0.7692519 

# biomass:
WT.cor.S = cor.test(sATF_full$adjWT, 
                    sATF_full$ATF_Bio, method = c("spearman")); WT.cor.S
  # S = 1.285e+11; p < 0.001; rho = 0.7722422 

#################################################################
# Retrospective skill testing:
# Load dynamic survey data linked to model covariates (for averaging each sequential time series):
load("Data/dynamic_model_data.rda")
  dynamic.model.data$Year = as.numeric(dynamic.model.data$Year)
load("Data/s_forecasts_stack_10km.rda")

# presence-absence:
ATF.retro = dynamic.model.data %>%
  filter(Species == "Arrowtooth.Flounder") %>%
  arrange(Year)

pa.predictions_all = NULL
for(i in 1992:2018) {
  
  pa.data = ATF.retro %>%
    group_by(lon, lat) %>%
    mutate_at(c("sponges", "corals", "whips"), as.character) %>%
    mutate_at(c("sponges", "corals", "whips"), as.numeric) %>%
    mutate_at(c("bathy", "slope", "ROMS_BT", "sponges", "corals", "whips"), mean) %>%
    as.data.frame() %>%
    mutate_at(c("sponges", "corals", "whips"), as.factor)
  
pa.model = gam(PA ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = pa.data, family = binomial(link="cloglog"), method = "GCV.Cp") 

pa.s.cov_pers = s.forecasts.10km
    pa.pers = predict.gam(pa.model, newdata = pa.s.cov_pers, 
                          type="response", se.fit=T)
    pa.pers = as.data.frame(pa.pers)
pa.pers = pa.pers %>% 
  rename(fit_pa.pers = fit, 
         se.fit_pa.pers = se.fit)
  pa.obs.pred = cbind(pa.s.cov_pers, pa.pers)
pa.obs.pred$forecast.yr = pa.obs.pred$Year
  pa.obs.pred$fitted.through = i - 1
pa.predictions_all = rbind(pa.predictions_all, pa.obs.pred)}

# positive catches:
# abundance:
catch.data = subset(pa.data, PA == 1)
N.predictions_all = NULL
for(i in 1992:2018) {
  
  N.data = catch.data %>%
    group_by(lon, lat) %>%
    mutate_at(c("sponges", "corals", "whips"), as.character) %>%
    mutate_at(c("sponges", "corals", "whips"), as.numeric) %>%
    mutate_at(c("bathy", "slope", "ROMS_BT", "sponges", "corals", "whips"), mean) %>%
    as.data.frame() %>%
    mutate_at(c("sponges", "corals", "whips"), as.factor)
  
N.model = gam(adjN ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, k=4) + s(ROMS_BT, bs='tp', k=4, m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = catch.data, family = poisson(link="log"), method = "GCV.Cp")  

N.s.cov_pers = s.forecasts.10km
    N.pers = predict.gam(N.model, newdata = N.s.cov_pers, 
                         type="response", se.fit=T)
    N.pers = as.data.frame(N.pers)
N.pers = N.pers %>% 
  rename(fit_N.pers = fit, 
         se.fit_N.pers = se.fit)
  N.obs.pred = cbind(N.s.cov_pers, N.pers)
N.obs.pred$forecast.yr = N.obs.pred$Year
  N.obs.pred$fitted.through = i - 1
N.predictions_all = rbind(N.predictions_all, N.obs.pred)}

# biomass:
WT.predictions_all = NULL
for(i in 1992:2018) {
  
  WT.data = catch.data %>%
    group_by(lon, lat) %>%
    mutate_at(c("sponges", "corals", "whips"), as.character) %>%
    mutate_at(c("sponges", "corals", "whips"), as.numeric) %>%
    mutate_at(c("bathy", "slope", "ROMS_BT", "sponges", "corals", "whips"), mean) %>%
    as.data.frame() %>%
    mutate_at(c("sponges", "corals", "whips"), as.factor)
  
WT.model = gam(adjWT ~ te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + corals + offset(log(AreaSwept)), data = catch.data, family = Gamma(link="log"), method = "GCV.Cp")  

WT.s.cov_pers = s.forecasts.10km
    WT.pers = predict.gam(WT.model, newdata = WT.s.cov_pers, 
                          type="response", se.fit=T)
    WT.pers = as.data.frame(WT.pers)
WT.pers = WT.pers %>% 
  rename(fit_WT.pers = fit, 
         se.fit_WT.pers = se.fit)
  WT.obs.pred = cbind(WT.s.cov_pers, WT.pers)
WT.obs.pred$forecast.yr = WT.obs.pred$Year
  WT.obs.pred$fitted.through = i - 1
WT.predictions_all = rbind(WT.predictions_all, WT.obs.pred)}

# all models:
pred.pos = N.predictions_all %>% full_join(WT.predictions_all)
pred.all = pa.predictions_all %>% full_join(pred.pos)

# Estimate abundance and biomass, accounting for probability of occurrence (estimate / AreaSwept => N or kg per sq. km):
pred.all$ATF_Abun = with(pred.all, (fit_pa.pers * fit_N.pers))
pred.all$ATF_Bio = with(pred.all, (fit_pa.pers * fit_WT.pers))
  pred.all$Abun.thousands = pred.all$ATF_Abun / 1000
  pred.all$Bio.MT = pred.all$ATF_Bio / 1000
  
  save(pred.all, file = "Analyses/sATF_retro.rda") 
# load("Analyses/sATF_retro.rda")
  
#################################################################
# Estimate Spearman's correlation coefficients for forecasts (tendency to co-vary; range: -1 to 1):
ATF_sf = st_as_sf(ATF.retro, coords = c("lon", "lat"))
  st_crs(ATF_sf) = 4326 # WGS84
pred.all_sf = st_as_sf(pred.all, coords = c("lon", "lat"))
  st_crs(pred.all_sf) = 4326 # WGS84
  
# Merge data from model fitting and predicting:
sATF.obs.pred = NULL
for(model in unique(pred.all_sf$fitted.through)) {
  sATF.pred = subset(pred.all_sf, fitted.through == model)
sATF = st_join(ATF_sf, sATF.pred, 
                   join = st_nearest_feature, left = T) 
sATF.obs.pred = rbind(sATF.obs.pred, sATF)} 

# Remove "hindcasts":
sATF_retro = as.data.frame(sATF.obs.pred)
sATF_retro = sATF_retro %>%
  filter(Year > fitted.through)

# presence-absence:
pa.coef = NULL
for(model in unique(sATF_retro$fitted.through)) {
  cor.pa.data = subset(sATF_retro, fitted.through == model)
    corr.pa.test.s = cor.test(cor.pa.data$PA, 
                              cor.pa.data$fit_pa.pers, 
                              method = "spearman")
  coef.pa = cbind(model, corr.pa.test.s$statistic, 
               format.pval(corr.pa.test.s$p.value, digits=5, eps=0.001), 
                  corr.pa.test.s$estimate) 
pa.coef = rbind(pa.coef, coef.pa) }

pa.obs.pred.cor = as.data.frame(pa.coef) 
colnames(pa.obs.pred.cor) = c("fitted.through", "S", "p-value.s", "rho")

pa.obs.pred.cor = pa.obs.pred.cor %>%
  mutate(fitted.through = as.numeric(fitted.through)) %>%
  arrange(desc(fitted.through))
pa.obs.pred.cor$Metric = "pa"
  
# positive catches (accounting for probability of occurrence):
# abundance:
N.coef = NULL
for(model in unique(sATF_retro$fitted.through)) {
  cor.N.data = subset(sATF_retro, fitted.through == model)
    corr.N.test.s = cor.test(cor.N.data$adjN, 
                             cor.N.data$ATF_Abun, 
                             method = "spearman")
  coef.N = cbind(model, corr.N.test.s$statistic, 
               format.pval(corr.N.test.s$p.value, digits=5, eps=0.001), 
                  corr.N.test.s$estimate)
N.coef = rbind(N.coef, coef.N) }

N.obs.pred.cor = as.data.frame(N.coef) 
colnames(N.obs.pred.cor) = c("fitted.through", "S", "p-value.s", "rho")

N.obs.pred.cor = N.obs.pred.cor %>%
  mutate(fitted.through = as.numeric(fitted.through)) %>%
  arrange(desc(fitted.through))
N.obs.pred.cor$Metric = "N"

# biomass:
WT.coef = NULL
for(model in unique(sATF_retro$fitted.through)) {
  cor.WT.data = subset(sATF_retro, fitted.through == model)
    corr.WT.test.s = cor.test(cor.WT.data$adjWT, 
                              cor.WT.data$ATF_Bio, 
                              method = "spearman")
  coef.WT = cbind(model, corr.WT.test.s$statistic, 
               format.pval(corr.WT.test.s$p.value, digits=5, eps=0.001), 
                  corr.WT.test.s$estimate)
WT.coef = rbind(WT.coef, coef.WT) }

WT.obs.pred.cor = as.data.frame(WT.coef) 
colnames(WT.obs.pred.cor) = c("fitted.through", "S", "p-value.s", "rho")

WT.obs.pred.cor = WT.obs.pred.cor %>%
  mutate(fitted.through = as.numeric(fitted.through)) %>%
  arrange(desc(fitted.through))
WT.obs.pred.cor$Metric = "WT"

obs.pred.pa.cor = rbind(pa.obs.pred.cor, N.obs.pred.cor, WT.obs.pred.cor)
obs.pred.pa.cor$Metric = as.factor(obs.pred.pa.cor$Metric)
  obs.pred.pa.cor$Metric = ordered(obs.pred.pa.cor$Metric, levels = c("pa", "N", "WT"))

obs.pred.pa.cor = obs.pred.pa.cor %>%
  mutate_at(c("rho"), as.numeric)  

levels(obs.pred.pa.cor$Metric) = c("Presence-Absence", "Abundance", "Biomass")
  save(obs.pred.pa.cor, file = "Analyses/sATF_obs_pred_pa_cor.rda")
# load("Analyses/sATF_obs_pred_pa_cor.rda")