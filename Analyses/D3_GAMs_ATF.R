# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Climate-informed models benefit hindcasting but may present challenges when forecasting species-habitat associations. Ecogr. 

# References: 
# Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.

# Overarching Objective: Identify best-fit generalized additive models (GAMs) for hindcasting (based on conventional statistics) and evaluate forecast skill (based on retrospective skill testing; Thorson 2019) for presence-absence, numerical abundance, and biomass of groundfishes in the Bering Sea (1982-2018). Compare species distribution models (SDMs) with varying degrees of complexity to assess whether the addition of time-varying processes to status quo static SDMs improves hindcast performance and/or forecast skill.

# Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

# Complex Dynamic GAMs, Arrowtooth Flounder (ATF):
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
d3GAM.ATFpa.full = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")
  summary(d3GAM.ATFpa.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d3GAM.ATFpa.check.full = round(as.data.frame(concurvity(d3GAM.ATFpa.full, full=F)), digits=3); d3GAM.ATFpa.check.full
# bathy-BPI > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
d3GAM.ATFpa.ind.cov.full = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")  
  summary(d3GAM.ATFpa.ind.cov.full)
# R-sq = 0.737; Dev. = 69.5%; UBRE = -0.57304
  
# Generate alternative models using backwards step-wise selection:
# Remove sponges (p > 0.1):
d3GAM.ATFpa.alt.1 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d3GAM.ATFpa.alt.1)
# R-sq = 0.736; Dev. = 69.5%; UBRE = -0.57313

# Remove whips (p > 0.1):
d3GAM.ATFpa.alt.2 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp") 
  summary(d3GAM.ATFpa.alt.2)
# R-sq = 0.736; Dev. = 69.5%; UBRE = -0.57319

# Remove slope (p > 0.1):
d3GAM.ATFpa.alt.3 = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = ATF, family = binomial(link="cloglog"), method = "GCV.Cp")
  summary(d3GAM.ATFpa.alt.3)
# R-sq = 0.736; Dev. = 69.4%; UBRE = -0.57312

# Save best-fit model based on GCV/UBRE:  
d3GAM.ATFpa.best = d3GAM.ATFpa.alt.2 
  summary(d3GAM.ATFpa.best)
      
  
# Plot partial effects of model covariates (Fig. S4):
plot.visreg = function() {
  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, family="Arial", size=12),
        panel.background = element_rect(fill=NA, colour=NA, size=0.01, line="solid"),         panel.grid = element_blank(), 
        axis.text = element_text(family="Arial", size=12),
        axis.title.x = element_text(vjust=-0.13, size=12),
        axis.title.y = element_text(vjust=2.0, size=12),
        legend.background = element_rect(fill="transparent")) }

lon_min = -178.65
lon_max = -156.5
lat_min = 50
lat_max = 65.9

data(worldHiresMapEnv)
  d3GAM.ATFpa.best$data = ATF

# position:
  pdf("Plots/FigS4_d3GAM_ATF_PA_lon_lat.pdf", width=6, height=8)
vis.gam(d3GAM.ATFpa.best, view=c("lon", "lat"), plot.type="contour", type="response", contour.col="black", color="heat", xlab="", ylab="", main="", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")
  dev.off()

# year:
d3GAM.ATFpa.best_year = visreg(d3GAM.ATFpa.best, xvar="Year", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(limits=c(1982,2019), breaks = seq(1983,2019,5), expand = c(0,0.03), labels=c("1983","","1993","","2003","","2013","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0"))
ggsave(filename="Plots/FigS4_d3GAM_ATF_PA_year.png", plot=d3GAM.ATFpa.best_year, dpi=500, width=2.5, height=3, units="in")

# bathy:
range(na.omit(ATF$GEAR_DEPTH))
d3GAM.ATFpa.best_bathy = visreg(d3GAM.ATFpa.best, xvar="bathy", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,1750), breaks=seq(0,1000,200), labels=c("0","","400","","800","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  coord_cartesian(xlim=c(0,1200)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_ATF_PA_bathy.png", plot=d3GAM.ATFpa.best_bathy, dpi=500, width=2.5, height=3, units="in")

# slope:
range(na.omit(ATF$slope))
d3GAM.ATFpa.best_slope = visreg(d3GAM.ATFpa.best, xvar="slope", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,26), breaks=seq(0,25,5), labels=c("0","","10","","20","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_ATF_PA_slope.png", plot=d3GAM.ATFpa.best_slope, dpi=500, width=2.5, height=3, units="in")

# bottom temperature (deg C):
range(ATF$ROMS_BT)
d3GAM.ATFpa.best_BT = visreg(d3GAM.ATFpa.best, xvar="ROMS_BT", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-2.1,13), breaks=seq(-2,13,2), labels=c("","0","","","6","","","12")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  coord_cartesian(xlim=c(-2.1,13)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_ATF_PA_BT.png", plot=d3GAM.ATFpa.best_BT, dpi=500, width=2.5, height=3, units="in")

# CPE:
range(ATF$CPE)
d3GAM.ATFpa.best_CPE = visreg(d3GAM.ATFpa.best, xvar="CPE", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,345), breaks=seq(0,300,50), labels=c("0","","100","","200","", "300")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS4_d3GAM_ATF_PA_CPE.png", plot=d3GAM.ATFpa.best_CPE, dpi=500, width=2.5, height=3, units="in")

# corals (presence-absence):
d3GAM.ATFpa.best_corals = visreg(d3GAM.ATFpa.best, xvar="corals", fill=list(col="red", fill="red", alpha=0.1), points=list(col="red", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.25), labels=c("0.0","","0.5","","1.0"))
ggsave(filename="Plots/FigS4_d3GAM_ATF_PA_corals.png", plot=d3GAM.ATFpa.best_corals, dpi=500, width=2.5, height=3, units="in")

#########################################################################
# Positive Catches:
ATF_catch = subset(ATF, adjN > 0 & adjWT > 0)

# Abundance #
# Full model: 
d3GAM.ATF.N.full = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = poisson(link="log"), method = "GCV.Cp")  
  summary(d3GAM.ATF.N.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d3GAM.ATF.N.check = round(as.data.frame(concurvity(d3GAM.ATF.N.full, full=F)), digits=3); d3GAM.ATF.N.check
# BPI-bathy > 0.5, exclude BPI
# phi-lon/lat > 0.5, exclude phi
# phi-CPE > 0.5, exclude phi
# ROMS_BT-lon/lat > 0.5, keep both - interested in both for predictions...
d3GAM.ATF.N.ind.cov = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = poisson(link="log"), method = "GCV.Cp") 
  summary(d3GAM.ATF.N.ind.cov)
# adj-R2 = 0.224; Deviance = 46.5%; UBRE = 17.737
  
# Save best-fit model: 
d3GAM.ATF.N.best = d3GAM.ATF.N.ind.cov
  summary(d3GAM.ATF.N.best)
    
#########################################################################
# Plot partial effects of model covariates (Fig. S5):  
d3GAM.ATF.N.best$data = ATF_catch
    
# position:
  pdf("Plots/FigS5_d3GAM_ATF_N_lon_lat.pdf", width=6, height=8)
vis.gam(d3GAM.ATF.N.best, view=c("lon", "lat"), plot.type="contour", type="link", contour.col="black", color="heat", xlab="", ylab="", main="", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")
  dev.off()

# year:
d3GAM.ATF.N.best_year = visreg(d3GAM.ATF.N.best, xvar="Year", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(limits=c(1982,2019), breaks = seq(1983,2019,5), expand = c(0,0.03), labels=c("1983","","1993","","2003","","2013","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80"))
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_year.png", plot=d3GAM.ATF.N.best_year, dpi=500, width=2.5, height=3, units="in")

# bathy:
range(na.omit(ATF$GEAR_DEPTH))
d3GAM.ATF.N.best_bathy = visreg(d3GAM.ATF.N.best, xvar="bathy", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,1750), breaks=seq(0,1000,200), labels=c("0","","400","","","1000")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80")) +
  coord_cartesian(xlim=c(0,1200)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_bathy.png", plot=d3GAM.ATF.N.best_bathy, dpi=500, width=2.5, height=3, units="in")

# slope:
range(na.omit(ATF$slope))
d3GAM.ATF.N.best_slope = visreg(d3GAM.ATF.N.best, xvar="slope", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,26), breaks=seq(0,25,5), labels=c("0","","10","","20","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_slope.png", plot=d3GAM.ATF.N.best_slope, dpi=500, width=2.5, height=3, units="in")

# bottom temperature (deg C):
range(ATF$ROMS_BT)
d3GAM.ATF.N.best_BT = visreg(d3GAM.ATF.N.best, xvar="ROMS_BT", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-2.1,13), breaks=seq(-2,13,2), labels=c("","0","","","6","","","12")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80")) +
  coord_cartesian(xlim=c(-2.1,13)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_BT.png", plot=d3GAM.ATF.N.best_BT, dpi=500, width=2.5, height=3, units="in")

# CPE:
range(ATF$CPE)
d3GAM.ATF.N.best_CPE = visreg(d3GAM.ATF.N.best, xvar="CPE", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,345), breaks=seq(0,300,50), labels=c("0","","100","","200","", "300")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_CPE.png", plot=d3GAM.ATF.N.best_CPE, dpi=500, width=2.5, height=3, units="in")

# sponges (presence-absence):
d3GAM.ATF.N.best_sponges = visreg(d3GAM.ATF.N.best, xvar="sponges", fill=list(col="red", fill="red", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80"))
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_sponges.png", plot=d3GAM.ATF.N.best_sponges, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# corals (presence-absence):
d3GAM.ATF.N.best_corals = visreg(d3GAM.ATF.N.best, xvar="corals", fill=list(col="red", fill="red", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80"))
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_corals.png", plot=d3GAM.ATF.N.best_corals, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# whips (presence-absence):
d3GAM.ATF.N.best_whips = visreg(d3GAM.ATF.N.best, xvar="whips", fill=list(col="blue", fill="blue", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="blue"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,90), breaks=seq(0,90,20), labels=c("0","","40","", "80"))
ggsave(filename="Plots/FigS5_d3GAM_ATF_N_whips.png", plot=d3GAM.ATF.N.best_whips, dpi=500, width=2.5, height=3, bg="transparent", units="in")

#########################################################################
# Biomass #
# Full model:
d3GAM.ATF.WT.full = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(BPI, bs='tp', k=4, m=1) + s(phi, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp")  
  summary(d3GAM.ATF.WT.full)

# Re-run the full model without covariates that exhibit high levels of concurvity:
  # 'estimate' values >= 0.5 * significant concurvity
d3GAM.ATF.WT.check = round(as.data.frame(concurvity(d3GAM.ATF.WT.full, full=F)), digits=3); d3GAM.ATF.WT.check
# BPI-bathy, exclude BPI
# phi-lon/lat, exclude phi
# phi-CPE, exclude phi
# ROMS_BT-lon/lat, keep both - interested in both for predictions...
d3GAM.ATF.WT.ind.cov = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp") 
  summary(d3GAM.ATF.WT.ind.cov)
# adj-R2 = 0.185; Deviance = 51.5%; GCV = 0.75375

# Generate alternative models using backwards step-wise selection:
# Remove sponges (p > 0.1):
d3GAM.ATF.WT.alt.1 = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = ATF_catch, family = Gamma(link="log"), method = "GCV.Cp") 
  summary(d3GAM.ATF.WT.alt.1)
# adj-R2 = 0.185; Deviance = 51.5%; GCV = 0.75352
  
# Save best-fit model: 
d3GAM.ATF.WT.best = d3GAM.ATF.WT.alt.1
  summary(d3GAM.ATF.WT.best)
    
#########################################################################
# Plot partial effects of model covariates (Fig. S6):  
  d3GAM.ATF.WT.best$data = ATF_catch

# position:
  pdf("Plots/FigS6_d3GAM_ATF_WT_lon_lat.pdf", width=6, height=8)
vis.gam(d3GAM.ATF.WT.best, view=c("lon", "lat"), plot.type="contour", type="link", contour.col="black", color="heat", xlab="", ylab="", main="", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")
  dev.off()

# year:
d3GAM.ATF.WT.best_year = visreg(d3GAM.ATF.WT.best, xvar="Year", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(limits=c(1982,2019), breaks = seq(1983,2019,5), expand = c(0,0.03), labels=c("1983","","19
                                                                                                  93","","2003","","2013","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200"))
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_year.png", plot=d3GAM.ATF.WT.best_year, dpi=500, width=2.5, height=3, units="in")

# bathy:
range(na.omit(ATF$GEAR_DEPTH))
d3GAM.ATF.WT.best_bathy = visreg(d3GAM.ATF.WT.best, xvar="bathy", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,1750), breaks=seq(0,1000,200), labels=c("0","","400","","","1000")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200")) +
  coord_cartesian(xlim=c(0,1200)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_bathy.png", plot=d3GAM.ATF.WT.best_bathy, dpi=500, width=2.5, height=3, units="in")

# slope:
range(na.omit(ATF$slope))
d3GAM.ATF.WT.best_slope = visreg(d3GAM.ATF.WT.best, xvar="slope", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,26), breaks=seq(0,25,5), labels=c("0","","10","","20","")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_slope.png", plot=d3GAM.ATF.WT.best_slope, dpi=500, width=2.5, height=3, units="in")

# bottom temperature (deg C):
range(ATF$ROMS_BT)
d3GAM.ATF.WT.best_BT = visreg(d3GAM.ATF.WT.best, xvar="ROMS_BT", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(-2.1,13), breaks=seq(-2,13,2), labels=c("","0","","","6","","","12")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200")) +
  coord_cartesian(xlim=c(-2.1,13)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_BT.png", plot=d3GAM.ATF.WT.best_BT, dpi=500, width=2.5, height=3, units="in")

# CPE:
range(ATF$CPE)
d3GAM.ATF.WT.best_CPE = visreg(d3GAM.ATF.WT.best, xvar="CPE", fill=list(col="gray", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  labs(x="", y="") +
  scale_x_continuous(expand=c(0,0.1), limits=c(0,345), breaks=seq(0,300,50), labels=c("0","","100","","200","", "300")) +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_CPE.png", plot=d3GAM.ATF.WT.best_CPE, dpi=500, width=2.5, height=3, units="in")

# corals (presence-absence):
d3GAM.ATF.WT.best_corals = visreg(d3GAM.ATF.WT.best, xvar="corals", fill=list(col="red", fill="red", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200"))
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_corals.png", plot=d3GAM.ATF.WT.best_corals, dpi=500, width=2.5, height=3, bg="transparent", units="in")

# whips (presence-absence):
d3GAM.ATF.WT.best_whips = visreg(d3GAM.ATF.WT.best, xvar="whips", fill=list(col="blue", fill="blue", alpha=0.1), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="blue"), band=T, partial=F, rug=F, scale="response", type="conditional", nn=500, gg=T) +
  plot.visreg() +
  theme(plot.background = element_blank()) +
  labs(x="", y="") +
  scale_y_continuous(expand=c(0,0.1), limits=c(0,220), breaks=seq(0,200,50), labels=c("0","","100","","200"))
ggsave(filename="Plots/FigS6_d3GAM_ATF_WT_whips.png", plot=d3GAM.ATF.WT.best_whips, dpi=500, width=2.5, height=3, bg="transparent", units="in")

#########################################################################
# Predict probability of occurrence, abundance, and biomass:
load("Data/dynamic_model_grid.rda")

# Probability of Occurrence:
ATFpa.pred = predict.gam(d3GAM.ATFpa.best, newdata = dynamic.model.grid, 
                         type = "response", se.fit = T) 
ATFpa.predict = cbind(dynamic.model.grid, ATFpa.pred)
ATFpa.predict = ATFpa.predict %>%
  rename(pa.fit = fit, pa.se.fit = se.fit)

# Abundance:
ATF.N.pred = predict.gam(d3GAM.ATF.N.best, newdata = dynamic.model.grid, 
                         type = "response", se.fit = T)
ATF.N.predict = cbind(dynamic.model.grid, ATF.N.pred)
ATF.N.predict = ATF.N.predict %>%
  rename(N.fit = fit, N.se.fit = se.fit)

# Biomass:
ATF.WT.pred = predict.gam(d3GAM.ATF.WT.best, newdata = dynamic.model.grid, 
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
  
d3ATF.obs.pred = list()
for(i in unique(ATF_sf$Year)) {
  OBS = subset(ATF_sf, Year == i)
  PRED = subset(ATF.predict_sf, Year == i)
  d3ATF.obs.pred[[i]] = st_join(OBS, PRED, join = st_nearest_feature, left = T)}
d3ATF_full = as.data.frame(dplyr::bind_rows(d3ATF.obs.pred))
d3ATF_full = d3ATF_full %>%
  rename(Year = Year.x)

# presence-absence:
pa.cor.S = cor.test(d3ATF_full$PA, d3ATF_full$pa.fit, method = c("spearman")); pa.cor.S
  # S = 1.0941e+11; p < 0.001; rho = 0.806078

pa.cor.yr = NULL
for(i in unique(d3ATF_full$Year)) {
  cor.data = subset(d3ATF_full, Year == i)
  corr.test = cor.test(cor.data$PA, cor.data$pa.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  pa.cor.yr = rbind(pa.cor.yr, coef) }

pa.cor.yr = as.data.frame(pa.cor.yr)
  colnames(pa.cor.yr) = c("Year", "S", "p-value", "rho")
pa.cor.yr[, c("Year", "S", "rho")] = 
  lapply(pa.cor.yr[, c("Year", "S", "rho")], as.numeric)  

# abundance:
N.cor.S = cor.test(d3ATF_full$adjN, d3ATF_full$ATF_Abun, method = c("spearman")); N.cor.S
# S = 1.0041e+11; p < 0.001; rho = 0.82202467

N.cor.yr = NULL
for(i in unique(d3ATF_full$Year)) {
  cor.data = subset(d3ATF_full, Year == i)
  corr.test = cor.test(cor.data$adjN, cor.data$N.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  N.cor.yr = rbind(N.cor.yr, coef) }

N.cor.yr = as.data.frame(N.cor.yr)
  colnames(N.cor.yr) = c("Year", "S", "p-value", "rho")
N.cor.yr[, c("Year", "S", "rho")] = 
  lapply(N.cor.yr[, c("Year", "S", "rho")], as.numeric)  

# biomass:
WT.cor.S = cor.test(d3ATF_full$adjWT, d3ATF_full$ATF_Bio, method = c("spearman")); WT.cor.S
# S = 9.7687e+10; p < 0.001; rho = 0.8268537

WT.cor.yr = NULL
for(i in unique(d3ATF_full$Year)) {
  cor.data = subset(d3ATF_full, Year == i)
  corr.test = cor.test(cor.data$adjWT, cor.data$WT.fit, method = c("spearman"))
  coef = cbind(i, corr.test$statistic, format.pval(corr.test$p.value, digits=5, eps=0.001), corr.test$estimate)
  WT.cor.yr = rbind(WT.cor.yr, coef) }

WT.cor.yr = as.data.frame(WT.cor.yr)
  colnames(WT.cor.yr) = c("Year", "S", "p-value", "rho")
WT.cor.yr[, c("Year", "S", "rho")] = 
  lapply(WT.cor.yr[, c("Year", "S", "rho")], as.numeric)  

#################################################################
# Plot model predictions as population percentiles (Figs. 2, S7-S9):
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
ATF.pa.predict.yr = ATF.predict %>%
  group_by(Year) %>% 
  arrange(desc(pa.fit)) %>%
  mutate(pa.percent.yr = round((100 * cumsum(pa.fit) / sum(pa.fit)), 
                            digits=3)) %>%
  mutate(pa.EFH.perc.yr = ifelse(pa.percent.yr < 25, 25, # hot spot
                          ifelse(pa.percent.yr < 50, 50, # core habitat
                          ifelse(pa.percent.yr < 75, 75, # suitable habitat
                          ifelse(pa.percent.yr < 95, 95, # EFH
                                                      0))))) 
ATF.pa.predict.yr$pa.EFH.perc.yr = as.factor(ATF.pa.predict.yr$pa.EFH.perc.yr)
  ATF.EFH.pa.predict.yr = subset(ATF.pa.predict.yr, pa.EFH.perc.yr != "0")
ATF.EFH.pa.predict.yr$pa.EFH.perc.yr = ordered(ATF.EFH.pa.predict.yr$pa.EFH.perc.yr, levels=c("95","75","50","25"))
  
  
# Fig. 2 (select years):
ATF.EFH.pa.select = subset(ATF.EFH.pa.predict.yr, 
                           Year == 1986 | 
                           Year == 1994 |
                           Year == 2002 |
                           Year == 2010 |
                           Year == 2018)

ATFpa.cont_select = ggplot(data=world) +
  geom_tile(data=ATF.EFH.pa.select, aes(x=lon, y=lat, fill=pa.EFH.perc.yr)) +
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
ggsave(ATFpa.cont_select, filename="Plots/Fig2_ATF_pa.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")

# Fig. S7 (all years):
ATF.EFH.pa.predict_yr = ATF.EFH.pa.predict.yr[,c("Year", "lon", "lat", "pa.EFH.perc.yr")] %>%
  as.data.frame()

ATF.EFH.pa.predict_yr = add_row(ATF.EFH.pa.predict_yr)
ATF.EFH.pa.predict_yr[nrow(ATF.EFH.pa.predict_yr),] = list(1980, NA, NA, NA)
ATF.EFH.pa.predict_yr = add_row(ATF.EFH.pa.predict_yr)
ATF.EFH.pa.predict_yr[nrow(ATF.EFH.pa.predict_yr),] = list(1981, NA, NA, NA)
ATF.EFH.pa.predict_yr = add_row(ATF.EFH.pa.predict_yr)
ATF.EFH.pa.predict_yr[nrow(ATF.EFH.pa.predict_yr),] = list(2019, NA, NA, NA)

ATFpa.cont_yr = ggplot(data=world) +
  geom_tile(data=ATF.EFH.pa.predict_yr, aes(x=lon, y=lat, fill=pa.EFH.perc.yr)) +
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
ggsave(ATFpa.cont_yr, filename="Plots/FigS7.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")

# Abundance - predicted percentiles:
ATF.Abun.predict.yr = ATF.predict %>%
  group_by(Year) %>% 
  arrange(desc(ATF_Abun)) %>%
  mutate(Abun.percent.yr = round((100 * cumsum(ATF_Abun) / sum(ATF_Abun)), digits=3)) %>%
  mutate(Abun.EFH.perc.yr = ifelse(Abun.percent.yr < 25, 25, # hot spot
                            ifelse(Abun.percent.yr < 50, 50, # core habitat
                            ifelse(Abun.percent.yr < 75, 75, # suitable habitat
                            ifelse(Abun.percent.yr < 95, 95, # EFH
                                                          0))))) 
ATF.Abun.predict.yr$Abun.EFH.perc.yr = as.factor(ATF.Abun.predict.yr$Abun.EFH.perc.yr)
  ATF.EFH.Abun.predict.yr = subset(ATF.Abun.predict.yr, Abun.EFH.perc.yr != "0")
ATF.EFH.Abun.predict.yr$Abun.EFH.perc.yr = ordered(ATF.EFH.Abun.predict.yr$Abun.EFH.perc.yr, levels=c("95","75","50","25"))


# Fig. S8 (all years):
ATF.EFH.Abun.predict_yr = ATF.EFH.Abun.predict.yr[,c("Year", "lon", "lat", "Abun.EFH.perc.yr")] %>%
  as.data.frame()

ATF.EFH.Abun.predict_yr = add_row(ATF.EFH.Abun.predict_yr)
ATF.EFH.Abun.predict_yr[nrow(ATF.EFH.Abun.predict_yr),] = list(1980, NA, NA, NA)
ATF.EFH.Abun.predict_yr = add_row(ATF.EFH.Abun.predict_yr)
ATF.EFH.Abun.predict_yr[nrow(ATF.EFH.Abun.predict_yr),] = list(1981, NA, NA, NA)
ATF.EFH.Abun.predict_yr = add_row(ATF.EFH.Abun.predict_yr)
ATF.EFH.Abun.predict_yr[nrow(ATF.EFH.Abun.predict_yr),] = list(2019, NA, NA, NA)

ATF.Abun.cont_yr = ggplot(data=world) +
  geom_tile(data=ATF.EFH.Abun.predict_yr, aes(x=lon, y=lat, fill=Abun.EFH.perc.yr)) +
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
ggsave(ATF.Abun.cont_yr, filename="Plots/FigS8.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")

# Biomass - predicted percentiles:
ATF.Bio.predict.yr = ATF.predict %>%
  group_by(Year) %>% 
  arrange(desc(ATF_Bio)) %>%
  mutate(Bio.percent.yr = round((100 * cumsum(ATF_Bio) / sum(ATF_Bio)), 
                            digits=3)) %>%
  mutate(Bio.EFH.perc.yr = ifelse(Bio.percent.yr < 25, 25, # hot spot
                           ifelse(Bio.percent.yr < 50, 50, # core habitat
                           ifelse(Bio.percent.yr < 75, 75, # suitable habitat
                           ifelse(Bio.percent.yr < 95, 95, # EFH
                                                      0))))) 
ATF.Bio.predict.yr$Bio.EFH.perc.yr = as.factor(ATF.Bio.predict.yr$Bio.EFH.perc.yr)
  ATF.EFH.Bio.predict.yr = subset(ATF.Bio.predict.yr, Bio.EFH.perc.yr != "0")
ATF.EFH.Bio.predict.yr$Bio.EFH.perc.yr = ordered(ATF.EFH.Bio.predict.yr$Bio.EFH.perc.yr, levels=c("95","75","50","25"))

  
# Fig. 2 (select years):
ATF.EFH.WT.select = subset(ATF.EFH.Bio.predict.yr, 
                           Year == 1986 | 
                           Year == 1994 |
                           Year == 2002 |
                           Year == 2010 |
                           Year == 2018)

ATF.WT.cont_select = ggplot(data=world) +
  geom_tile(data=ATF.EFH.WT.select, aes(x=lon, y=lat, fill=Bio.EFH.perc.yr)) +
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
ggsave(ATF.WT.cont_select, filename="Plots/Fig2_ATF_WT.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")


# Fig. S9 (all years):
ATF.EFH.Bio.predict_yr = ATF.EFH.Bio.predict.yr[,c("Year", "lon", "lat", "Bio.EFH.perc.yr")] %>%
  as.data.frame()

ATF.EFH.Bio.predict_yr = add_row(ATF.EFH.Bio.predict_yr)
ATF.EFH.Bio.predict_yr[nrow(ATF.EFH.Bio.predict_yr),] = list(1980, NA, NA, NA)
ATF.EFH.Bio.predict_yr = add_row(ATF.EFH.Bio.predict_yr)
ATF.EFH.Bio.predict_yr[nrow(ATF.EFH.Bio.predict_yr),] = list(1981, NA, NA, NA)
ATF.EFH.Bio.predict_yr = add_row(ATF.EFH.Bio.predict_yr)
ATF.EFH.Bio.predict_yr[nrow(ATF.EFH.Bio.predict_yr),] = list(2019, NA, NA, NA)

ATF.Bio.cont_yr = ggplot(data=world) +
  geom_tile(data=ATF.EFH.Bio.predict_yr, aes(x=lon, y=lat, fill=Bio.EFH.perc.yr)) +
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
ggsave(ATF.Bio.cont_yr, filename="Plots/FigS9.png", dpi=500, width=9, height=6.5, units="in", bg="transparent")

#################################################################
# Retrospective skill testing:
load("Data/dyn_forecasts_stack_10km.rda")

# presence-absence:
ATF = ATF %>%
  arrange(Year)

pa.predictions_all = NULL
for(i in 1992:2018) {
  pa.data = subset(ATF, Year < i) # for model fitting (i = first year forecasted)
  
pa.model = gam(PA ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + offset(log(AreaSwept)), data = pa.data, family = binomial(link="cloglog"), method = "GCV.Cp") 

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
  
N.model = gam(adjN ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + sponges + corals + whips + offset(log(AreaSwept)), data = N.data, family = poisson(link="log"), method = "GCV.Cp") 

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
  
WT.model = gam(adjWT ~ s(Year, bs="cr", m=1) + te(lon, lat, bs='tp', m=1) + ti(lon, lat, Year, bs='tp', m=1) + s(bathy, bs='tp', k=6, m=1) + s(slope, bs='tp', k=4, m=1) + s(ROMS_BT, bs='tp', k=4, m=1) + s(lon, lat, by=CPE, bs='tp', m=1) + corals + whips + offset(log(AreaSwept)), data = WT.data, family = Gamma(link="log"), method = "GCV.Cp")  

  WT.dyn.cov = subset(d.forecasts.10km, Year >= i)
WT.pred = predict.gam(WT.model, newdata = WT.dyn.cov, 
                      type="response", se.fit=T)
  WT.dyn.cov_pers = subset(d.forecasts.10km, Year == i - 1)
    WT.pers = predict.gam(WT.model, newdata = WT.dyn.cov_pers, type="response", se.fit=T)
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
  
  save(pred_all, file = "Analyses/d3ATF_retro.rda") 
# load("Analyses/d3ATF_retro.rda")
  
#################################################################
# Estimate Spearman's correlation coefficients for forecasts (tendency to co-vary; range: -1 to 1):
ATF_sf = st_as_sf(ATF, coords = c("lon", "lat"))
  st_crs(ATF_sf) = 4326 # WGS84
pred_all_sf = st_as_sf(pred_all, coords = c("lon", "lat"))
  st_crs(pred_all_sf) = 4326 # WGS84
  
# Merge data from model fitting and predicting:
d3ATF.obs.pred_all = NULL
for(model in unique(pred_all_sf$fitted.through)) {
  d3ATF.pred = subset(pred_all_sf, fitted.through == model)
   for(yr in unique(d3ATF.pred$Year)) {
     d3ATF.pred.yr = subset(d3ATF.pred, Year == yr)
     d3ATF.obs.yr = subset(ATF_sf, Year == yr)
     d3ATF.obs.pred_yr = st_join(d3ATF.obs.yr, d3ATF.pred.yr, join = st_nearest_feature, left=T)
d3ATF.obs.pred_all = rbind(d3ATF.obs.pred_all, d3ATF.obs.pred_yr)} }

d3ATF_retro_all = as.data.frame(d3ATF.obs.pred_all)

# Estimate Spearman's correlation coefficient and RMSE:
d3ATF_retro_all$fitted.through = as.factor(d3ATF_retro_all$fitted.through)
  format.pval(c(0.1, 0.0001, 1e-27))

d3ATF_retro_all = d3ATF_retro_all %>%
  rename(Year = Year.x)

# presence-absence:
pa.obs.pred.cor_all = NULL
for(model in unique(d3ATF_retro_all$fitted.through)) {
  cor.pa.data = subset(d3ATF_retro_all, fitted.through == model)
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
for(model in unique(d3ATF_retro_all$fitted.through)) {
  cor.N.data = subset(d3ATF_retro_all, fitted.through == model)
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
for(model in unique(d3ATF_retro_all$fitted.through)) {
  cor.WT.data = subset(d3ATF_retro_all, fitted.through == model)
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
  save(obs.pred.pa.cor_all, file = "Analyses/d3ATF_obs_pred_pa_cor.rda")
# load("Analyses/d3ATF_obs_pred_pa_cor.rda")