# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Climate-informed models benefit hindcasting but may present challenges when forecasting species-habitat associations. Ecogr.  

# References: 
# Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.

# Overarching Objective: Identify best-fit generalized additive models (GAMs) for hindcasting (based on conventional statistics) and evaluate forecast skill (based on retrospective skill testing; Thorson 2019) for presence-absence, numerical abundance, and biomass of groundfishes in the Bering Sea (1982-2018). Compare species distribution models (SDMs) with varying degrees of complexity to assess whether the addition of time-varying processes to status quo static SDMs improves hindcast performance and/or forecast skill.

# Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

# Intermediate Dynamic GAMs, Arrowtooth Flounder (ATF):
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
require(Metrics)
    
# Load survey data linked to model covariates:
load("Data/dynamic_model_data.rda")

##################################################################
# Presence-Absence:
dynamic.model.data$Year = as.numeric(dynamic.model.data$Year)
ATF = dynamic.model.data %>%
  filter(Species == "Arrowtooth.Flounder")

# Full model, year = fixed effect:
d2GAM.ATFpa.full = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")
  summary(d2GAM.ATFpa.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d2GAM.ATFpa.check.full = round(as.data.frame(concurvity(d2GAM.ATFpa.full, full=F)), digits=3); d2GAM.ATFpa.check.full
# bathy-BPI > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
d2GAM.ATFpa.ind.cov.full = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")  
  summary(d2GAM.ATFpa.ind.cov.full)
# R-sq = 0.721; Dev. = 67.3%; UBRE = -0.55095
  
# Generate alternative models using backwards step-wise selection:
# Remove whips (p > 0.1):
d2GAM.ATFpa.alt.1 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d2GAM.ATFpa.alt.1)
# R-sq = 0.721; Dev. = 67.3%; UBRE = -0.55096

# Remove slope (p > 0.1):
d2GAM.ATFpa.alt.2 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d2GAM.ATFpa.alt.2)
# R-sq = 0.721; Dev. = 67.2%; UBRE = -0.5509

# Remove sponges (p > 0.1):
d2GAM.ATFpa.alt.3 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d2GAM.ATFpa.alt.3)
# R-sq = 0.721; Dev. = 67.2%; UBRE = -0.551
  
# Save best-fit model based on GCV/UBRE:  
d2GAM.ATFpa.best = d2GAM.ATFpa.alt.3
  summary(d2GAM.ATFpa.best)

#########################################################################
# Positive Catches:
ATF_catch = subset(ATF, adjN > 0 & adjWT > 0)

# Abundance #
# Full model:
d2GAM.ATF.N.full = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = poisson(link="log"), method = "GCV.Cp")  
  summary(d2GAM.ATF.N.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d2GAM.ATF.N.check = round(as.data.frame(concurvity(d2GAM.ATF.N.full, full=F)), digits=3); d2GAM.ATF.N.check
# BPI-bathy > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
d2GAM.ATF.N.ind.cov = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = poisson(link="log"), method = "GCV.Cp")
  summary(d2GAM.ATF.N.ind.cov)
# adj-R2 = 0.193; Deviance = 41.3%; UBRE = 19.523
  
# Save best-fit model: 
d2GAM.ATF.N.best = d2GAM.ATF.N.ind.cov
  summary(d2GAM.ATF.N.best); sum(d2GAM.ATF.N.best$edf)

#########################################################################
# Biomass #
# Full model:
d2GAM.ATF.WT.full = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")  
  summary(d2GAM.ATF.WT.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d2GAM.ATF.WT.check = round(as.data.frame(concurvity(d2GAM.ATF.WT.full, full=F)), digits=3); d2GAM.ATF.WT.check
# BPI-bathy, exclude BPI
# phi-lon/lat, exclude phi
# phi-CPE, exclude phi
d2GAM.ATF.WT.ind.cov = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")
  summary(d2GAM.ATF.WT.ind.cov)
# adj-R2 = 0.160; Deviance = 47.0%; GCV = 0.80703

# Generate alternative models using backwards step-wise selection:
# Remove sponges (p > 0.1):
d2GAM.ATF.WT.alt.1 = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")
  summary(d2GAM.ATF.WT.alt.1)
# adj-R2 = 0.161; Deviance = 46.9%; GCV = 0.80721
  
# Save best-fit model: 
d2GAM.ATF.WT.best = d2GAM.ATF.WT.ind.cov
  summary(d2GAM.ATF.WT.best)
    
#########################################################################
# Predict probability of occurrence, abundance, and biomass:
load("Data/dynamic_model_grid.rda")

# Probability of Occurrence:
ATFpa.pred = predict.gam(d2GAM.ATFpa.best, newdata = dynamic.model.grid, 
                         type = "response", se.fit = T) 
ATFpa.predict = cbind(dynamic.model.grid, ATFpa.pred)
ATFpa.predict = ATFpa.predict %>%
  rename(pa.fit = fit, pa.se.fit = se.fit)

# Abundance:
ATF.N.pred = predict.gam(d2GAM.ATF.N.best, newdata = dynamic.model.grid, 
                         type = "response", se.fit = T)
ATF.N.predict = cbind(dynamic.model.grid, ATF.N.pred)
ATF.N.predict = ATF.N.predict %>%
  rename(N.fit = fit, N.se.fit = se.fit)

# Biomass:
ATF.WT.pred = predict.gam(d2GAM.ATF.WT.best, newdata = dynamic.model.grid, 
                          type = "response", se.fit = T)
ATF.WT.predict = cbind(dynamic.model.grid, ATF.WT.pred)
ATF.WT.predict = ATF.WT.predict %>%
  rename(WT.fit = fit, WT.se.fit = se.fit)

# Merge predictions:
ATF.predict = ATFpa.predict %>% 
  full_join(ATF.N.predict) %>%
  full_join(., ATF.WT.predict)

# Estimate abundance and biomass, accounting for probability of occurrence:
ATF.predict$ATF_Abun = with(ATF.predict, pa.fit * N.fit)
ATF.predict$ATF_Bio = with(ATF.predict, pa.fit * WT.fit)
  
##################################################################  
# Estimate Spearman's correlation coefficients for hindcasts (tendency to co-vary; range: -1 to 1):

# Merge data from model fitting and predicting:
ATF_sf = st_as_sf(ATF, coords = c("lon", "lat"))
  st_crs(ATF_sf) = 4326
ATF.predict_sf = st_as_sf(ATF.predict, coords = c("lon", "lat"))
  st_crs(ATF.predict_sf) = 4326
  
d2ATF.obs.pred = list()
for(i in unique(ATF_sf$Year)) {
  OBS = subset(ATF_sf, Year == i)
  PRED = subset(ATF.predict_sf, Year == i)
  d2ATF.obs.pred[[i]] = st_join(OBS, PRED, join = st_nearest_feature, left = T)}
d2ATF_full = as.data.frame(dplyr::bind_rows(d2ATF.obs.pred))
d2ATF_full = d2ATF_full %>%
  rename(Year = Year.x)

# correlation coefficients (tendency to co-vary, range: -1 to 1):
# presence-absence:
pa.cor.S = cor.test(d2ATF_full$PA, d2ATF_full$pa.fit, method = c("spearman")); pa.cor.S
  # S = 1.1244e+11; p < 0.001; rho = 0.800705 

# abundance:
N.cor.S = cor.test(d2ATF_full$adjN, d2ATF_full$ATF_Abun, method = c("spearman")); N.cor.S
# S = 1.0432e+11; p < 0.001; rho = 0.815092

# biomass:
WT.cor.S = cor.test(d2ATF_full$adjWT, d2ATF_full$ATF_Bio, method = c("spearman")); WT.cor.S
# S = 1.0107e+11; p < 0.001; rho = 0.8208627

#################################################################
# Retrospective skill testing:
load("Data/dyn_forecasts_stack_10km.rda")

# presence-absence:
ATF = ATF %>%
  arrange(Year)

pa.predictions_all = NULL
for(i in 1992:2018) {
  pa.data = subset(ATF, Year < i) # for model fitting (i = first year forecasted)
  
pa.model = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = pa.data, family = binomial(link="cloglog"), method = "GCV.Cp") 

pa.dyn.cov = subset(d.forecasts.10km, Year >= i)
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
  
# positive catches:
ATF_catch = ATF_catch %>%
  arrange(Year)

# abundance:
N.predictions_all = NULL
for(i in 1992:2018) {
  N.data = subset(ATF_catch, Year < i) # for model fitting (i = first year forecasted)
  
N.model = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = N.data, family = poisson(link="log"), method = "GCV.Cp")

N.dyn.cov = subset(d.forecasts.10km, Year >= i)
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
  
# biomass:
WT.predictions_all = NULL
for(i in 1992:2018) {
  WT.data = subset(ATF_catch, Year < i) # for model fitting (i = first year forecasted)
  
WT.model = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WT.data , family = Gamma(link="log"), method = "GCV.Cp")

WT.dyn.cov = subset(d.forecasts.10km, Year >= i)
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

# all models:
pred.pos_all = N.predictions_all %>% full_join(WT.predictions_all)
pred_all = pa.predictions_all %>% full_join(pred.pos_all)

# Estimate abundance and biomass, accounting for probability of occurrence (estimate / AreaSwept => N or kg per sq. km):
pred_all$ATF_Abun = with(pred_all, (pa.fit * N.fit))
pred_all$ATF_Bio = with(pred_all, (pa.fit * WT.fit))
  pred_all$Abun.thousands = pred_all$ATF_Abun / 1000
  pred_all$Bio.MT = pred_all$ATF_Bio / 1000

# persistence
pred_all$ATF_Abun.pers = with(pred_all, (pa.fit * fit_N.pers))
  pred_all$Abun_pers.thousands = pred_all$ATF_Abun.pers / 1000
pred_all$ATF_Bio.pers = with(pred_all, (pa.fit * fit_WT.pers))
  pred_all$Bio_pers.MT = pred_all$ATF_Bio.pers / 1000 
  
  save(pred_all, file = "Analyses/d2ATF_retro.rda") 
# load("Analyses/d2ATF_retro.rda")
  
#################################################################
# Estimate Spearman's correlation coefficients for forecasts (tendency to co-vary; range: -1 to 1):
ATF_sf = st_as_sf(ATF, coords = c("lon", "lat"))
  st_crs(ATF_sf) = 4326 # WGS84
pred_all_sf = st_as_sf(pred_all, coords = c("lon", "lat"))
  st_crs(pred_all_sf) = 4326 # WGS84
  
# Merge data from model fitting and predicting:
d2ATF.obs.pred_all = NULL
for(model in unique(pred_all_sf$fitted.through)) {
  d2ATF.pred = subset(pred_all_sf, fitted.through == model)
   for(yr in unique(d2ATF.pred$Year)) {
     d2ATF.pred.yr = subset(d2ATF.pred, Year == yr)
     d2ATF.obs.yr = subset(ATF_sf, Year == yr)
     d2ATF.obs.pred_yr = st_join(d2ATF.obs.yr, d2ATF.pred.yr, join = st_nearest_feature, left=T)
d2ATF.obs.pred_all = rbind(d2ATF.obs.pred_all, d2ATF.obs.pred_yr)} }

d2ATF_retro_all = as.data.frame(d2ATF.obs.pred_all)

# Estimate Spearman's correlation coefficient and RMSE:
d2ATF_retro_all$fitted.through = as.factor(d2ATF_retro_all$fitted.through)
  format.pval(c(0.1, 0.0001, 1e-27))

d2ATF_retro_all = d2ATF_retro_all %>%
  rename(Year = Year.x)

# presence-absence:
pa.obs.pred.cor_all = NULL
for(model in unique(d2ATF_retro_all$fitted.through)) {
  cor.pa.data = subset(d2ATF_retro_all, fitted.through == model)
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
for(model in unique(d2ATF_retro_all$fitted.through)) {
  cor.N.data = subset(d2ATF_retro_all, fitted.through == model)
for(yr in unique(cor.N.data$Year)) {
  cor.N.yr.data = subset(cor.N.data, Year == yr) # various correlation coefficients
    corr.N.test.s = cor.test(cor.N.yr.data$adjN, cor.N.yr.data$ATF_Abun, 
                              method = "spearman")
    corr.N.test.s_pers = cor.test(cor.N.yr.data$adjN, cor.N.yr.data$ATF_Abun.pers, 
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
for(model in unique(d2ATF_retro_all$fitted.through)) {
  cor.WT.data = subset(d2ATF_retro_all, fitted.through == model)
for(yr in unique(cor.WT.data$Year)) {
    cor.WT.yr.data = subset(cor.WT.data, Year == yr) # various correlation coefficients
    corr.WT.test.s = cor.test(cor.WT.yr.data$adjWT, cor.WT.yr.data$ATF_Bio, 
                              method = "spearman")
    corr.WT.test.s_pers = cor.test(cor.WT.yr.data$adjWT, cor.WT.yr.data$ATF_Bio.pers, 
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
  save(obs.pred.pa.cor_all, file = "Analyses/d2ATF_obs_pred_pa_cor.rda")
# load("Analyses/d2ATF_obs_pred_pa_cor.rda")