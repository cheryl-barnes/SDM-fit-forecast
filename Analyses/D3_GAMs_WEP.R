# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Climate-informed models benefit hindcasting but may present challenges when forecasting species-habitat associations. Ecogr. 

# References: 
# Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.

# Overarching Objective: Identify best-fit generalized additive models (GAMs) for hindcasting (based on conventional statistics) and evaluate forecast skill (based on retrospective skill testing; Thorson 2019) for presence-absence, numerical abundance, and biomass of groundfishes in the Bering Sea (1982-2018). Compare species distribution models (SDMs) with varying degrees of complexity to assess whether the addition of time-varying processes to status quo static SDMs improves hindcast performance and/or forecast skill.

# Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

# Complex Dynamic GAMs, Walleye Pollock (WEP):
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

#########################################################################
# Presence-Absence:
dynamic.model.data$Year = as.numeric(dynamic.model.data$Year)
WEP = dynamic.model.data %>%
  filter(Species == "Walleye.Pollock")

# Full model, year = fixed effect:
d3GAM.WEPpa.full = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp")
  summary(d3GAM.WEPpa.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d3GAM.WEPpa.check.full = round(as.data.frame(concurvity(d3GAM.WEPpa.full, full=F)), digits=3); d3GAM.WEPpa.check.full
# bathy-BPI > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
d3GAM.WEPpa.ind.cov.full = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp")  
  summary(d3GAM.WEPpa.ind.cov.full)
# R-sq = 0.496; Dev. = 48.5%; UBRE = -0.57047
  
# Generate alternative models using backwards step-wise selection:
# Remove sponges (p > 0.1):
d3GAM.WEPpa.alt.1 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d3GAM.WEPpa.alt.1)
# R-sq = 0.496; Dev. = 48.5%; UBRE = -0.5706

# Remove whips (p > 0.1):
d3GAM.WEPpa.alt.2 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = WEP, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d3GAM.WEPpa.alt.2)
# R-sq = 0.496; Dev. = 48.5%; UBRE = -0.5706

# Save best-fit model based on GCV/UBRE:  
d3GAM.WEPpa.best = d3GAM.WEPpa.alt.2 
  # same GCV, more parsimonious
  summary(d3GAM.WEPpa.best)
      
######################################################################### 
# Plot partial effects of model covariates (Fig. S4):
plot.visreg = function() {
  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, family="Arial", size=12),
        panel.background = element_rect(fill=NA, colour=NA, size=0.01, line="solid"),
        panel.grid = element_blank(), 
        axis.text = element_text(family="Arial", size=12),
        axis.title.x = element_text(vjust=-0.13, size=12),
        axis.title.y = element_text(vjust=2.0, size=12),
        legend.background = element_rect(fill="transparent")) }

lon_min = -178.65
lon_max = -156.5
lat_min = 50
lat_max = 65.9

data(worldHiresMapEnv)
  d3GAM.WEPpa.best$data = WEP

# position:
  pdf("Plots/FigS4_d3GAM_WEP_PA_lon_lat.pdf", width=6, height=8)
vis.gam(d3GAM.WEPpa.best, view=c("lon", "lat"), plot.type="contour", type="response", contour.col="black", color="heat", xlab="", ylab="", main="", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")
  dev.off()

# year:
d3GAM.WEPpa.best_year = visreg(d3GAM.WEPpa.best, xvar="Year", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(limits=c(1982,2019), breaks = seq(1983,2019,5), expand = c(0,0.03), labels=c("1983","","1993","","2003","","2013","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0"))
ggsave(filename="Plots/FigS4_d3GAM_WEP_PA_year.png", plot=d3GAM.WEPpa.best_year, dpi=500, width=2.5, height=3, units="in")

# bathy:
range(na.omit(WEP$GEAR_DEPTH))
d3GAM.WEPpa.best_bathy = visreg(d3GAM.WEPpa.best, xvar="bathy", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,1750), breaks=seq(0,1000,200), labels=c("0","","400","","","1000")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  coord_cartesian(xlim=c(0,1200)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_WEP_PA_bathy.png", plot=d3GAM.WEPpa.best_bathy, dpi=500, width=2.5, height=3, units="in")

# slope:
range(na.omit(WEP$slope))
d3GAM.WEPpa.best_slope = visreg(d3GAM.WEPpa.best, xvar="slope", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,26), breaks=seq(0,25,5), labels=c("0","","10","","20","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_WEP_PA_slope.png", plot=d3GAM.WEPpa.best_slope, dpi=500, width=2.5, height=3, units="in")

# bottom temperature (deg C):
range(WEP$ROMS_BT)
d3GAM.WEPpa.best_BT = visreg(d3GAM.WEPpa.best, xvar="ROMS_BT", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-2.1,13), breaks=seq(-2,13,2), labels=c("","0","","","6","","","12")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  coord_cartesian(xlim=c(-2.1,13)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_WEP_PA_BT.png", plot=d3GAM.WEPpa.best_BT, dpi=500, width=2.5, height=3, units="in")

# CPE:
range(WEP$CPE)
d3GAM.WEPpa.best_CPE = visreg(d3GAM.WEPpa.best, xvar="CPE", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,345), breaks=seq(0,300,50), labels=c("0","","100","","200","", "300")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_WEP_PA_CPE.png", plot=d3GAM.WEPpa.best_CPE, dpi=500, width=2.5, height=3, units="in")

# corals (presence-absence):
d3GAM.WEPpa.best_corals = visreg(d3GAM.WEPpa.best, xvar="corals", fill=list(col="red", fill="red", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0"))
ggsave(filename="Plots/FigS4_d3GAM_WEP_PA_corals.png", plot=d3GAM.WEPpa.best_corals, dpi=500, width=2.5, height=3, units="in")

#########################################################################
# Positive Catches:
WEP_catch = subset(WEP, adjN > 0 & adjWT > 0)

# Abundance #
# Full model:
d3GAM.WEP.N.full = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = poisson(link="log"), method = "GCV.Cp")  
  summary(d3GAM.WEP.N.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d3GAM.WEP.N.check = round(as.data.frame(concurvity(d3GAM.WEP.N.full, full=F)), digits=3); d3GAM.WEP.N.check
# bathy-lon/lat > 0.5, exclude bathy
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
d3GAM.WEP.N.ind.cov = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = poisson(link="log"), method = "GCV.Cp") 
  summary(d3GAM.WEP.N.ind.cov)
# adj-R2 = 0.147; Deviance = 35.0%; UBRE = 814.93
  
# Save best-fit model: 
d3GAM.WEP.N.best = d3GAM.WEP.N.ind.cov
  summary(d3GAM.WEP.N.best)
    
#########################################################################
# Plot partial effects of model covariates (Fig. S5):  
d3GAM.WEP.N.best$data = WEP_catch
    
# position:
  pdf("Plots/FigS5_d3GAM_WEP_N_lon_lat.pdf", width=6, height=8)
vis.gam(d3GAM.WEP.N.best, view=c("lon", "lat"), plot.type="contour", type="link", contour.col="black", color="heat", xlab="", ylab="", main="", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")
  dev.off()

# year:
d3GAM.WEP.N.best_year = visreg(d3GAM.WEP.N.best, xvar="Year", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(limits=c(1982,2019), breaks = seq(1983,2019,5), expand = c(0,0.03), labels=c("1983","","1993","","2003","","2013","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200")) 
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_year.png", plot=d3GAM.WEP.N.best_year, dpi=500, width=2.5, height=3, units="in")

# slope:
range(na.omit(WEP$slope))
d3GAM.WEP.N.best_slope = visreg(d3GAM.WEP.N.best, xvar="slope", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,26), breaks=seq(0,25,5), labels=c("0","","10","","20","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_slope.png", plot=d3GAM.WEP.N.best_slope, dpi=500, width=2.5, height=3, units="in")

# BPI:
range(na.omit(WEP$BPI))
d3GAM.WEP.N.best_BPI = visreg(d3GAM.WEP.N.best, xvar="BPI", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-200,1000), breaks=seq(-200,1000,200), labels=c("","0","","400","","800","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_BPI.png", plot=d3GAM.WEP.N.best_BPI, dpi=500, width=2.5, height=3, units="in")

# bottom temperature (deg C):
range(WEP$ROMS_BT)
d3GAM.WEP.N.best_BT = visreg(d3GAM.WEP.N.best, xvar="ROMS_BT", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-2.1,13), breaks=seq(-2,13,2), labels=c("","0","","","6","","","12")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_BT.png", plot=d3GAM.WEP.N.best_BT, dpi=500, width=2.5, height=3, units="in")

# CPE:
range(WEP$CPE)
d3GAM.WEP.N.best_CPE = visreg(d3GAM.WEP.N.best, xvar="CPE", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,345), breaks=seq(0,300,50), labels=c("0","","100","","200","", "300")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_CPE.png", plot=d3GAM.WEP.N.best_CPE, dpi=500, width=2.5, height=3, units="in")

# sponges (presence-absence):
d3GAM.WEP.N.best_sponges = visreg(d3GAM.WEP.N.best, xvar="sponges", fill=list(col="orange", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="orange"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) + 
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_sponges.png", plot=d3GAM.WEP.N.best_sponges, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# corals (presence-absence):
d3GAM.WEP.N.best_corals = visreg(d3GAM.WEP.N.best, xvar="corals", fill=list(col="red", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) + 
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_corals.png", plot=d3GAM.WEP.N.best_corals, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# whips (presence-absence):
d3GAM.WEP.N.best_whips = visreg(d3GAM.WEP.N.best, xvar="whips", fill=list(col="blue", alpha=0.1), points=list(col="cadet1", alpha=0.5, cex=0.25), line=list(col="blue"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) + 
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1300,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS5_d3GAM_WEP_N_whips.png", plot=d3GAM.WEP.N.best_whips, dpi=500, width=2.5, height=3, bg="transparent", units="in")

#########################################################################
# Biomass #
# Full model:
d3GAM.WEP.WT.full = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = Gamma(link="log"), method = "GCV.Cp")  
  summary(d3GAM.WEP.WT.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d3GAM.WEP.WT.check = round(as.data.frame(concurvity(d3GAM.WEP.WT.full, full=F)), digits=3); d3GAM.WEP.WT.check
# bathy-lon/lat, exclude bathy
# phi-lon/lat, exclude phi
# phi-CPE, exclude phi
d3GAM.WEP.WT.ind.cov = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WEP_catch, family = Gamma(link="log"), method = "GCV.Cp") 
  summary(d3GAM.WEP.WT.ind.cov)
# adj-R2 = 0.0926; Deviance = 35.9%; GCV = 1.8164

# Save best-fit model: 
d3GAM.WEP.WT.best = d3GAM.WEP.WT.ind.cov
  summary(d3GAM.WEP.WT.best)
    
#########################################################################
# Plot partial effects of model covariates (Fig. S6):  
d3GAM.WEP.WT.best$data = WEP_catch

# position:
  pdf("Plots/FigS6_d3GAM_WEP_WT_lon_lat.pdf", width=6, height=8)
vis.gam(d3GAM.WEP.WT.best, view=c("lon", "lat"), plot.type="contour", type="link", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")
  dev.off()

# year:
d3GAM.WEP.WT.best_year = visreg(d3GAM.WEP.WT.best, xvar="Year", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(limits=c(1982,2019), breaks = seq(1983,2019,5), expand = c(0,0.03), labels=c("1983","","1993","","2003","","2013","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_year.png", plot=d3GAM.WEP.WT.best_year, dpi=500, width=2.5, height=3, units="in")

# slope:
range(na.omit(WEP$slope))
d3GAM.WEP.WT.best_slope = visreg(d3GAM.WEP.WT.best, xvar="slope", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,26), breaks=seq(0,25,5), labels=c("0","","10","","20","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_slope.png", plot=d3GAM.WEP.WT.best_slope, dpi=500, width=2.5, height=3, units="in")

# BPI:
range(na.omit(WEP$BPI))
d3GAM.WEP.WT.best_BPI = visreg(d3GAM.WEP.WT.best, xvar="BPI", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-200,1000), breaks=seq(-200,1000,200), labels=c("","0","","400","","800","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_BPI.png", plot=d3GAM.WEP.WT.best_BPI, dpi=500, width=2.5, height=3, units="in")

# bottom temperature (deg C):
range(WEP$ROMS_BT)
d3GAM.WEP.WT.best_BT = visreg(d3GAM.WEP.WT.best, xvar="ROMS_BT", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-2.1,13), breaks=seq(-2,13,2), labels=c("","0","","","6","","","12")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_BT.png", plot=d3GAM.WEP.WT.best_BT, dpi=500, width=2.5, height=3, units="in")

# CPE:
range(WEP$CPE)
d3GAM.WEP.WT.best_CPE = visreg(d3GAM.WEP.WT.best, xvar="CPE", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,345), breaks=seq(0,300,50), labels=c("0","","100","","200","", "300")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_CPE.png", plot=d3GAM.WEP.WT.best_CPE, dpi=500, width=2.5, height=3, units="in")

# sponges (presence-absence):
d3GAM.WEP.WT.best_sponges = visreg(d3GAM.WEP.WT.best, xvar="sponges", fill=list(col="orange", fill="orange", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="orange"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_sponges.png", plot=d3GAM.WEP.WT.best_sponges, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# corals (presence-absence):
d3GAM.WEP.WT.best_corals = visreg(d3GAM.WEP.WT.best, xvar="corals", fill=list(col="red", fill="red", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_corals.png", plot=d3GAM.WEP.WT.best_corals, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# whips (presence-absence):
d3GAM.WEP.WT.best_whips = visreg(d3GAM.WEP.WT.best, xvar="whips", fill=list(col="blue", fill="blue", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="blue"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1350), breaks=seq(0,1200,300), labels=c("0","","600","","1200"))
ggsave(filename="Plots/FigS6_d3GAM_WEP_WT_whips.png", plot=d3GAM.WEP.WT.best_whips, dpi=500, width=2.5, height=3, bg="transparent", units="in")

#########################################################################
# Predict probability of occurrence, abundance, and biomass:
load("Data/dynamic_model_grid.rda")

# Probability of Occurrence:
WEPpa.pred = predict.gam(d3GAM.WEPpa.best, newdata = dynamic.model.grid, 
                         type = "response", se.fit = T) 
WEPpa.predict = cbind(dynamic.model.grid, WEPpa.pred)
WEPpa.predict = WEPpa.predict %>%
  rename(pa.fit = fit, pa.se.fit = se.fit)

# Abundance:
WEP.N.pred = predict.gam(d3GAM.WEP.N.best, newdata = dynamic.model.grid, 
                         type = "response", se.fit = T)
WEP.N.predict = cbind(dynamic.model.grid, WEP.N.pred)
WEP.N.predict = WEP.N.predict %>%
  rename(N.fit = fit, N.se.fit = se.fit)

# Biomass:
WEP.WT.pred = predict.gam(d3GAM.WEP.WT.best, newdata = dynamic.model.grid, 
                          type = "response", se.fit = T)
WEP.WT.predict = cbind(dynamic.model.grid, WEP.WT.pred)
WEP.WT.predict = WEP.WT.predict %>%
  rename(WT.fit = fit, WT.se.fit = se.fit)

# Merge predictions:
WEP.predict = WEPpa.predict %>% 
  full_join(WEP.N.predict) %>%
  full_join(., WEP.WT.predict)

# Estimate abundance and biomass, accounting for probability of occurrence:
WEP.predict$WEP_Abun = with(WEP.predict, pa.fit * N.fit)
WEP.predict$WEP_Bio = with(WEP.predict, pa.fit * WT.fit)
    
##################################################################  
# Estimate Spearman's correlation coefficients for hindcasts (tendency to co-vary; range: -1 to 1):

# Merge data from model fitting and predicting:
WEP_sf = st_as_sf(WEP, coords = c("lon", "lat"))
  st_crs(WEP_sf) = 4326
WEP.predict_sf = st_as_sf(WEP.predict, coords = c("lon", "lat"))
  st_crs(WEP.predict_sf) = 4326
  
d3WEP.obs.pred = list()
for(i in unique(WEP_sf$Year)) {
  OBS = subset(WEP_sf, Year == i)
  PRED = subset(WEP.predict_sf, Year == i)
  d3WEP.obs.pred[[i]] = st_join(OBS, PRED, join = st_nearest_feature, left = T)}
d3WEP_full = as.data.frame(dplyr::bind_rows(d3WEP.obs.pred))
d3WEP_full = d3WEP_full %>%
  rename(Year = Year.x)

# presence-absence:
pa.cor.S = cor.test(d3WEP_full$PA, d3WEP_full$pa.fit, method = c("spearman")); pa.cor.S
  # S = 2.8368e+11; p < 0.001; rho = 0.5324179

pa.cor.yr = NULL
for(i in unique(d3WEP_full$Year)) {
  cor.data = subset(d3WEP_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  pa.cor.yr = rbind(pa.cor.yr, coef) }

pa.cor.yr = as.data.frame(pa.cor.yr)
  colnames(pa.cor.yr) = c("Year", "S", "p-value", "rho")
pa.cor.yr[, c("Year", "S", "rho")] = 
  lapply(pa.cor.yr[, c("Year", "S", "rho")], as.numeric)  

# abundance:
N.cor.S = cor.test(d3WEP_full$adjN, d3WEP_full$WEP_Abun, method = c("spearman")); N.cor.S
# S = 1.7558e+11; p < 0.001; rho = 0.7106028 

N.cor.yr = NULL
for(i in unique(d3WEP_full$Year)) {
  cor.data = subset(d3WEP_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  N.cor.yr = rbind(N.cor.yr, coef) }

N.cor.yr = as.data.frame(N.cor.yr)
  colnames(N.cor.yr) = c("Year", "S", "p-value", "rho")
N.cor.yr[, c("Year", "S", "rho")] = 
  lapply(N.cor.yr[, c("Year", "S", "rho")], as.numeric)  

# biomass:
WT.cor.S = cor.test(d3WEP_full$adjWT, d3WEP_full$WEP_Bio, method = c("spearman")); WT.cor.S
# S = 1.884e+11; p < 0.001; rho = 0.6894704

WT.cor.yr = NULL
for(i in unique(d3WEP_full$Year)) {
  cor.data = subset(d3WEP_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  WT.cor.yr = rbind(WT.cor.yr, coef) }

WT.cor.yr = as.data.frame(WT.cor.yr)
  colnames(WT.cor.yr) = c("Year", "S", "p-value", "rho")
WT.cor.yr[, c("Year", "S", "rho")] = 
  lapply(WT.cor.yr[, c("Year", "S", "rho")], as.numeric)  

#################################################################
# Plot model predictions as population percentiles (Figs. 2, S10-S12):
world = ne_countries(scale = "large", returnclass = "sf")

plot.theme = function() {
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(), 
        legend.title = element_text(family="Arial", size=11, hjust=0),
        legend.text = element_text(family="Arial", size=10),
        legend.text.align = 0,
        legend.key.width = unit(4.0, "mm"),
        legend.key.height = unit(5.0, "mm"),
        legend.background = element_rect(fill="transparent"), 
        axis.title.x = element_text(vjust=-0.1, hjust=0.51, family="Arial", size=11),
        axis.title.y = element_text(vjust=2.0, hjust=0.485, family="Arial", size=11),
        axis.text = element_text(family="Arial", size=10),
        strip.background = element_blank(),
        strip.text = element_text(family="Arial", size=10.5, hjust=0.5)) }

# Northern Bering Sea boundary layers:
setwd("~/Documents/JISAO/SDM-fit-forecast/Data/shapefiles/")
  NBS = readOGR(".", "northern_bs")
NBS.proj = spTransform(NBS, CRS("+init=epsg:4326"))

# Dissolve boundaries within the NBS:
NBS.bound = data.frame()
  NBS.data = rbind(NBS.bound, NBS.proj@data)
NBS.data$Region = "NBS"

NBS.proj@data = full_join(NBS.proj@data, NBS.data)
  row.names(NBS.proj) = row.names(NBS.proj@data)
NBS.proj = spChFIDs(NBS.proj, row.names(NBS.proj))
NBS.proj = gUnaryUnion(NBS.proj, id = NBS.proj@data$Region)
  row.names(NBS.proj) = as.character(1:length(NBS.proj))
NBS.area = subset(fortify(NBS.proj), group != "1.2")
setwd("~/Documents/JISAO/SDM-fit-forecast/")

# Presence-absence - predicted percentiles:
WEP.pa.predict.yr = WEP.predict %>%
  group_by(Year) %>% 
  arrange(desc(pa.fit)) %>%
  mutate(pa.percent.yr = round((100 * cumsum(pa.fit) / sum(pa.fit)), 
                            digits=3)) %>%
  mutate(pa.EFH.perc.yr = ifelse(pa.percent.yr < 25, 25, # hot spot
                          ifelse(pa.percent.yr < 50, 50, # core habitat
                          ifelse(pa.percent.yr < 75, 75, # suitable habitat
                          ifelse(pa.percent.yr < 95, 95, # EFH
                                                      0))))) 
WEP.pa.predict.yr$pa.EFH.perc.yr = as.factor(WEP.pa.predict.yr$pa.EFH.perc.yr)
  WEP.EFH.pa.predict.yr = subset(WEP.pa.predict.yr, pa.EFH.perc.yr != "0")
WEP.EFH.pa.predict.yr$pa.EFH.perc.yr = ordered(WEP.EFH.pa.predict.yr$pa.EFH.perc.yr, levels=c("95","75","50","25"))


# Fig. 2 (select years):
WEP.EFH.pa.select = subset(WEP.EFH.pa.predict.yr, 
                           Year == 1986 | 
                           Year == 1994 |
                           Year == 2002 |
                           Year == 2010 |
                           Year == 2018)

WEPpa.cont_select = ggplot(data=world) +
  geom_tile(data=WEP.EFH.pa.select, aes(x=lon, y=lat, fill=pa.EFH.perc.yr)) +
  scale_fill_manual(values = c("midnightblue", "deepskyblue4", "seagreen4", "gold1"), na.value = NA, name = "") +
  geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.width = unit(1.5, "mm")) +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60","")) +
  facet_wrap(~ Year, nrow=1)
ggsave(WEPpa.cont_select, filename="Plots/Fig2_WEP_pa.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")


# Fig. S10 (all years):
WEP.EFH.pa.predict_yr = WEP.EFH.pa.predict.yr[,c("Year", "lon", "lat", "pa.EFH.perc.yr")] %>%
  as.data.frame()

WEP.EFH.pa.predict_yr = add_row(WEP.EFH.pa.predict_yr)
WEP.EFH.pa.predict_yr[nrow(WEP.EFH.pa.predict_yr),] = list(1980, NA, NA, NA)
WEP.EFH.pa.predict_yr = add_row(WEP.EFH.pa.predict_yr)
WEP.EFH.pa.predict_yr[nrow(WEP.EFH.pa.predict_yr),] = list(1981, NA, NA, NA)
WEP.EFH.pa.predict_yr = add_row(WEP.EFH.pa.predict_yr)
WEP.EFH.pa.predict_yr[nrow(WEP.EFH.pa.predict_yr),] = list(2019, NA, NA, NA)

WEPpa.cont_yr = ggplot(data=world) +
  geom_tile(data=WEP.EFH.pa.predict_yr, aes(x=lon, y=lat, fill=pa.EFH.perc.yr)) +
  scale_fill_manual(values = c("midnightblue", "deepskyblue4", "seagreen4", "gold1"), na.value = NA, name = "") +
  geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.width = unit(1.5, "mm")) +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60","")) +
  facet_wrap(~ Year, ncol=10, drop=F)
ggsave(WEPpa.cont_yr, filename="Plots/FigS10.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")
  
# Abundance - predicted percentiles:
WEP.Abun.predict.yr = WEP.predict %>%
  group_by(Year) %>% 
  arrange(desc(WEP_Abun)) %>%
  mutate(Abun.percent.yr = round((100 * cumsum(WEP_Abun) / sum(WEP_Abun)), digits=3)) %>%
  mutate(Abun.EFH.perc.yr = ifelse(Abun.percent.yr < 25, 25, # hot spot
                            ifelse(Abun.percent.yr < 50, 50, # core habitat
                            ifelse(Abun.percent.yr < 75, 75, # suitable habitat
                            ifelse(Abun.percent.yr < 95, 95, # EFH
                                                          0))))) 
WEP.Abun.predict.yr$Abun.EFH.perc.yr = as.factor(WEP.Abun.predict.yr$Abun.EFH.perc.yr)
  WEP.EFH.Abun.predict.yr = subset(WEP.Abun.predict.yr, Abun.EFH.perc.yr != "0")
WEP.EFH.Abun.predict.yr$Abun.EFH.perc.yr = ordered(WEP.EFH.Abun.predict.yr$Abun.EFH.perc.yr, levels=c("95","75","50","25"))


# Fig. S11 (all years):
WEP.EFH.Abun.predict_yr = WEP.EFH.Abun.predict.yr[,c("Year", "lon", "lat", "Abun.EFH.perc.yr")] %>%
  as.data.frame()

WEP.EFH.Abun.predict_yr = add_row(WEP.EFH.Abun.predict_yr)
WEP.EFH.Abun.predict_yr[nrow(WEP.EFH.Abun.predict_yr),] = list(1980, NA, NA, NA)
WEP.EFH.Abun.predict_yr = add_row(WEP.EFH.Abun.predict_yr)
WEP.EFH.Abun.predict_yr[nrow(WEP.EFH.Abun.predict_yr),] = list(1981, NA, NA, NA)
WEP.EFH.Abun.predict_yr = add_row(WEP.EFH.Abun.predict_yr)
WEP.EFH.Abun.predict_yr[nrow(WEP.EFH.Abun.predict_yr),] = list(2019, NA, NA, NA)

WEP.Abun.cont_yr = ggplot(data=world) +
  geom_tile(data=WEP.EFH.Abun.predict_yr, aes(x=lon, y=lat, fill=Abun.EFH.perc.yr)) +
  scale_fill_manual(values = c("midnightblue", "deepskyblue4", "seagreen4", "gold1"), na.value = NA, name = "") +
  geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.width = unit(1.5, "mm")) +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60","")) +
  facet_wrap(~ Year, ncol=10, drop=F)
ggsave(WEP.Abun.cont_yr, filename="Plots/FigS11.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")

# Biomass - predicted percentiles:
WEP.Bio.predict.yr = WEP.predict %>%
  group_by(Year) %>% 
  arrange(desc(WEP_Bio)) %>%
  mutate(Bio.percent.yr = round((100 * cumsum(WEP_Bio) / sum(WEP_Bio)), 
                            digits=3)) %>%
  mutate(Bio.EFH.perc.yr = ifelse(Bio.percent.yr < 25, 25, # hot spot
                           ifelse(Bio.percent.yr < 50, 50, # core habitat
                           ifelse(Bio.percent.yr < 75, 75, # suitable habitat
                           ifelse(Bio.percent.yr < 95, 95, # EFH
                                                      0))))) 
WEP.Bio.predict.yr$Bio.EFH.perc.yr = as.factor(WEP.Bio.predict.yr$Bio.EFH.perc.yr)
  WEP.EFH.Bio.predict.yr = subset(WEP.Bio.predict.yr, Bio.EFH.perc.yr != "0")
WEP.EFH.Bio.predict.yr$Bio.EFH.perc.yr = ordered(WEP.EFH.Bio.predict.yr$Bio.EFH.perc.yr, levels=c("95","75","50","25"))

  
# Fig. 2 (select years):
WEP.EFH.WT.select = subset(WEP.EFH.Bio.predict.yr, 
                           Year == 1986 | 
                           Year == 1994 |
                           Year == 2002 |
                           Year == 2010 |
                           Year == 2018)

WEP.WT.cont_select = ggplot(data=world) +
  geom_tile(data=WEP.EFH.WT.select, aes(x=lon, y=lat, fill=Bio.EFH.perc.yr)) +
  scale_fill_manual(values = c("midnightblue", "deepskyblue4", "seagreen4", "gold1"), na.value = NA, name = "") +
  geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.width = unit(1.5, "mm")) +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60","")) +
  facet_wrap(~ Year, nrow=1)
ggsave(WEP.WT.cont_select, filename="Plots/Fig2_WEP_WT.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")

# Fig. S12 (all years):
WEP.EFH.Bio.predict_yr = WEP.EFH.Bio.predict.yr[,c("Year", "lon", "lat", "Bio.EFH.perc.yr")] %>%
  as.data.frame()

WEP.EFH.Bio.predict_yr = add_row(WEP.EFH.Bio.predict_yr)
WEP.EFH.Bio.predict_yr[nrow(WEP.EFH.Bio.predict_yr),] = list(1980, NA, NA, NA)
WEP.EFH.Bio.predict_yr = add_row(WEP.EFH.Bio.predict_yr)
WEP.EFH.Bio.predict_yr[nrow(WEP.EFH.Bio.predict_yr),] = list(1981, NA, NA, NA)
WEP.EFH.Bio.predict_yr = add_row(WEP.EFH.Bio.predict_yr)
WEP.EFH.Bio.predict_yr[nrow(WEP.EFH.Bio.predict_yr),] = list(2019, NA, NA, NA)

WEP.Bio.cont_yr = ggplot(data=world) +
  geom_tile(data=WEP.EFH.Bio.predict_yr, aes(x=lon, y=lat, fill=Bio.EFH.perc.yr)) +
  scale_fill_manual(values = c("midnightblue", "deepskyblue4", "seagreen4", "gold1"), na.value = NA, name = "") +
  geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.key.width = unit(1.5, "mm")) +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60","")) +
  facet_wrap(~ Year, ncol=10, drop=F)
ggsave(WEP.Bio.cont_yr, filename="Plots/FigS12.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")
        
#################################################################
# Retrospective skill testing:
load("Data/dyn_forecasts_stack_10km.rda")

# presence-absence:
WEP = WEP %>%
  arrange(Year)

pa.predictions_all = NULL
for(i in 1992:2018) {
  pa.data = subset(WEP, Year < i) # for model fitting (i = first year forecasted)
  
pa.model = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = pa.data, family = binomial(link="cloglog"), method = "GCV.Cp") 

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
  
# positive catches:
WEP_catch = WEP_catch %>%
  arrange(Year)

# abundance:
N.predictions_all = NULL
for(i in 1992:2018) {
  N.data = subset(WEP_catch, Year < i) # for model fitting (i = first year forecasted)
  
N.model = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = N.data, family = poisson(link="log"), method = "GCV.Cp") 

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
  
# biomass:
WT.predictions_all = NULL
for(i in 1992:2018) {
  WT.data = subset(WEP_catch, Year < i) # for model fitting (i = first year forecasted)
  
WT.model = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = WT.data, family = Gamma(link="log"), method = "GCV.Cp") 

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

# all models:
pred.pos_all = N.predictions_all %>% full_join(WT.predictions_all)
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
  
  save(pred_all, file = "Analyses/d3WEP_retro.rda") 
# load("Analyses/d3WEP_retro.rda")
  
#################################################################
# Estimate Spearman's correlation coefficients for forecasts (tendency to co-vary; range: -1 to 1):
WEP_sf = st_as_sf(WEP, coords = c("lon", "lat"))
  st_crs(WEP_sf) = 4326 # WGS84
pred_all_sf = st_as_sf(pred_all, coords = c("lon", "lat"))
  st_crs(pred_all_sf) = 4326 # WGS84
  
# Merge data from model fitting and predicting:
d3WEP.obs.pred_all = NULL
for(model in unique(pred_all_sf$fitted.through)) {
  d3WEP.pred = subset(pred_all_sf, fitted.through == model)
   for(yr in unique(d3WEP.pred$Year)) {
     d3WEP.pred.yr = subset(d3WEP.pred, Year == yr)
     d3WEP.obs.yr = subset(WEP_sf, Year == yr)
     d3WEP.obs.pred_yr = st_join(d3WEP.obs.yr, d3WEP.pred.yr, join = st_nearest_feature, left=T)
d3WEP.obs.pred_all = rbind(d3WEP.obs.pred_all, d3WEP.obs.pred_yr)} }

d3WEP_retro_all = as.data.frame(d3WEP.obs.pred_all)

# Estimate Spearman's correlation coefficient and RMSE:
d3WEP_retro_all$fitted.through = as.factor(d3WEP_retro_all$fitted.through)
  format.pval(c(0.1, 0.0001, 1e-27))

d3WEP_retro_all = d3WEP_retro_all %>%
  rename(Year = Year.x)

  save(d3WEP_retro_all, file = "Analyses/d3WEP_obs_pred_retro.rda")
# load("Analyses/d3WEP_obs_pred_retro.rda")  
  
# presence-absence:
pa.obs.pred.cor_all = NULL
for(model in unique(d3WEP_retro_all$fitted.through)) {
  cor.pa.data = subset(d3WEP_retro_all, fitted.through == model)
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
for(model in unique(d3WEP_retro_all$fitted.through)) {
  cor.N.data = subset(d3WEP_retro_all, fitted.through == model)
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
for(model in unique(d3WEP_retro_all$fitted.through)) {
  cor.WT.data = subset(d3WEP_retro_all, fitted.through == model)
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
  save(obs.pred.pa.cor_all, file = "Analyses/d3WEP_obs_pred_pa_cor.rda")
# load("Analyses/d3WEP_obs_pred_pa_cor.rda")