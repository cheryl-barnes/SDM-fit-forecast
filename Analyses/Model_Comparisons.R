# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Climate-informed models benefit hindcasting but may present challenges when forecasting species-habitat associations. Ecogr. 

# Script Objectives: 1) Compare Spearman's rho (rank-based correlation coefficient; measure of forecast skill) among species, model types, and population metrics. 2) Plot standardized residuals for static and complex dynamic models to assess spatiotemporal variation in forecast bias. 3) Evaluate differences in bottom temperatures used to fit and forecast from each model as a potential explanation for differences in forecast skill.

# Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

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
require(sp)
require(sf)
require(ggplot2)
require(visreg)
require(rnaturalearth)
require(rnaturalearthdata)
require(rnaturalearthhires)
require(raster)
require(zoo)
require(RColorBrewer)

myPalette = rev(brewer.pal(11, "RdBu"))
myPalette.rev = brewer.pal(11, "RdBu")

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
        axis.title.x = element_text(vjust=-0.1, hjust=0.51, family="Arial", size=12),
        axis.title.y = element_text(vjust=2.0, hjust=0.485, family="Arial", size=12),
        axis.text = element_text(family="Arial", size=12),
        strip.background = element_blank(),
        strip.text = element_text(family="Arial", size=12, hjust=0.5)) }

world = ne_countries(scale = "large", returnclass = "sf")

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
##################################################################
# Input Spearman's rho from retrospective skill testing (exploratory models without bottom temperature):   
load("Analyses/sATF_obs_pred_pa_cor.rda")
  sRetro.ATF_cor = obs.pred.pa.cor
  sRetro.ATF_cor$model = "S"
  sRetro.ATF_cor$species = "ATF"
load("Analyses/d1ATF_obs_pred_pa_cor.rda")
  dRetro.ATF_cor_1 = obs.pred.pa.cor_all
  dRetro.ATF_cor_1$model = "D1"
  dRetro.ATF_cor_1$species = "ATF"
load("Analyses/d2ATF_obs_pred_pa_cor.rda")
  dRetro.ATF_cor_2 = obs.pred.pa.cor_all
  dRetro.ATF_cor_2$model = "D2"
  dRetro.ATF_cor_2$species = "ATF"
load("Analyses/d3ATF_obs_pred_pa_cor.rda")
  dRetro.ATF_cor_3 = obs.pred.pa.cor_all
  dRetro.ATF_cor_3$model = "D3"
  dRetro.ATF_cor_3$species = "ATF"
  
load("Analyses/sWEP_obs_pred_pa_cor.rda")
  sRetro.WEP_cor = obs.pred.pa.cor
  sRetro.WEP_cor$model = "S"
  sRetro.WEP_cor$species = "WEP"
load("Analyses/d1WEP_obs_pred_pa_cor.rda")
  dRetro.WEP_cor_1 = obs.pred.pa.cor_all
  dRetro.WEP_cor_1$model = "D1"
  dRetro.WEP_cor_1$species = "WEP"
load("Analyses/d2WEP_obs_pred_pa_cor.rda")
  dRetro.WEP_cor_2 = obs.pred.pa.cor_all
  dRetro.WEP_cor_2$model = "D2"
  dRetro.WEP_cor_2$species = "WEP"
load("Analyses/d3WEP_obs_pred_pa_cor.rda")
  dRetro.WEP_cor_3 = obs.pred.pa.cor_all
  dRetro.WEP_cor_3$model = "D3"
  dRetro.WEP_cor_3$species = "WEP"
  
obs.forecast_cor = sRetro.WEP_cor %>% 
  full_join(., dRetro.WEP_cor_1) %>%
  full_join(., dRetro.WEP_cor_2) %>%
  full_join(., dRetro.WEP_cor_3) %>%
  full_join(., sRetro.ATF_cor) %>%
  full_join(., dRetro.ATF_cor_1) %>%
  full_join(., dRetro.ATF_cor_2) %>%
  full_join(., dRetro.ATF_cor_3)

obs.forecast_cor$model = as.factor(obs.forecast_cor$model)
obs.forecast_cor$model = ordered(obs.forecast_cor$model, 
            levels = c("S", "D1", "D2", "D3")) 

obs.forecast_cor$species = as.factor(obs.forecast_cor$species)
  levels(obs.forecast_cor$species) = 
    c("Arrowtooth Flounder", "Walleye Pollock")
  
# Estimate mean and sd Spearman's rho (observed vs. forecasted) using 10-Year moving windows:
obs.forecast_s.cor = subset(obs.forecast_cor, model == "S")
obs.forecast_s.block = bind_rows(replicate(15, obs.forecast_s.cor, simplify = F))
  obs.forecast_s.block$yr.forecast = rep(1:15, 
                        each=nrow(obs.forecast_s.cor))

obs.forecast_dyn.cor = subset(obs.forecast_cor, model != "S")

obs.forecast_dyn.cor$window = with(obs.forecast_dyn.cor, 
                      ifelse(yr.forecast < 11, 5,
                      ifelse(yr.forecast >= 6 & yr.forecast < 16, 10,
                      ifelse(yr.forecast >= 11 & yr.forecast < 21, 15, NA))))

obs.forecast_block = obs.forecast_dyn.cor %>%
  group_by(species, Metric, model, window) %>%
  mutate(mean.rho = mean(rho, na.rm = T)) %>%
  mutate(sd.rho = sd(rho, na.rm = T))
obs.forecast_block = as.data.frame(obs.forecast_block)
obs.forecast_block = subset(obs.forecast_block, !is.na(window))

obs.forecast_block_plot = ggplot() +
  geom_smooth(data=obs.forecast_block, aes(x=yr.forecast, y=rho, group=model, col=model, fill=model, linetype=model), lwd=0.75, method="loess", span=2, se=T) +
  geom_smooth(data=obs.forecast_s.block, aes(x=yr.forecast, y=rho, group=model, col=model, fill=model, linetype=model), lwd=0.75, method="loess", span=2, se=T) +
  scale_color_manual(values = c("gold1", "darkorange1", "firebrick1", "firebrick4"), breaks = c("S", "D1", "D2", "D3"), name = "Model") +
  scale_fill_manual(values = c("gold1", "darkorange1", "firebrick1", "firebrick4"), breaks = c("S", "D1", "D2", "D3"), name = "Model") +
  scale_linetype_manual(values = c("dotted", "dotdash", "longdash", "solid"), breaks = c("S", "D1", "D2", "D3"),name = "Model") +
  plot.theme() +
  theme(legend.position = c(0.075,0.405), 
        legend.key.height = unit(3.3, "mm"),
        legend.key.width = unit(8.5, "mm"),
        strip.text.y = element_blank()) +
  facet_grid(species ~ Metric) +
  labs(x = "Forecast Year", y = expression("Spearman's Correlation Coefficient ("~rho~")")) +
  scale_x_continuous(limits=c(0,15))
ggsave(filename="Plots/Fig3.png", plot=obs.forecast_block_plot, dpi=500, width=6.5, height=6.5, units="in") 

##################################################################
# Plot standardized residuals for static and complex dynamic models:
load("Analyses/sATF_obs_pred_retro.rda")
  sATF.retro = sATF_retro
  sATF.retro$model = as.factor("S")
  sATF.retro$species = "ATF"
load("Analyses/d3ATF_obs_pred_retro.rda")
  d3ATF.retro = d3ATF_retro_all
  d3ATF.retro$model = as.factor("D3")
  d3ATF.retro$species = "ATF"
load("Analyses/sWEP_obs_pred_retro.rda")
  sWEP.retro = sWEP_retro
  sWEP.retro$model = as.factor("S")
  sWEP.retro$species = "WEP"
load("Analyses/d3WEP_obs_pred_retro.rda")
  d3WEP.retro = d3WEP_retro_all
  d3WEP.retro$model = as.factor("D3")
  d3WEP.retro$species = "WEP"
  
obs.pred_retro = sATF.retro %>% 
  full_join(., d3ATF.retro) %>%
  full_join(., sWEP.retro) %>% 
  full_join(., d3WEP.retro)

# Biomass, select years (Fig. 4):
sATF.retro_sf = st_as_sf(sATF.retro)  
  sATF_retro = sATF.retro_sf %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    as.data.frame()
sATF_retro$resid.WT = with(sATF_retro, 
                          log(adjWT+0.0001) - log(ATF_Bio+0.0001))
sATF_retro = sATF_retro[,c("model", "Year", "ID", "lon", "lat", "fitted.through", "resid.WT")]

d3ATF.retro_sf = st_as_sf(d3ATF.retro)  
  d3ATF_retro = d3ATF.retro_sf %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    as.data.frame()
d3ATF_retro$resid.WT = with(d3ATF_retro, 
                          log(adjWT+0.0001) - log(ATF_Bio+0.0001))
d3ATF_retro = d3ATF_retro %>%
  rename(Year = Year.x)
d3ATF_retro = d3ATF_retro[,c("model", "Year", "ID", "lon", "lat", "fitted.through", "resid.WT")]

ATF.resid = sATF_retro %>% full_join(d3ATF_retro)
ATF.resid_multi = subset(ATF.resid, 
                         fitted.through == 2002 & Year == 2003 | 
                         fitted.through == 2009 & Year == 2010 | 
                         fitted.through == 2016 & Year == 2017)  

sWEP.retro_sf = st_as_sf(sWEP.retro)  
  sWEP_retro = sWEP.retro_sf %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    as.data.frame() 
sWEP_retro$resid.WT = with(sWEP_retro, 
                          log(adjWT+0.0001) - log(WEP_Bio+0.0001))
sWEP_retro = sWEP_retro[,c("model", "Year", "ID", "lon", "lat", "fitted.through", "resid.WT")]

d3WEP.retro_sf = st_as_sf(d3WEP.retro)  
  d3WEP_retro = d3WEP.retro_sf %>%
    dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                  lat = sf::st_coordinates(.)[,2]) %>%
    as.data.frame()
d3WEP_retro$resid.WT = with(d3WEP_retro, 
                          log(adjWT+0.0001) - log(WEP_Bio+0.0001))
d3WEP_retro = d3WEP_retro %>%
  rename(Year = Year.x)
d3WEP_retro = d3WEP_retro[,c("model", "Year", "ID", "lon", "lat", "fitted.through", "resid.WT")]

WEP.resid = sWEP_retro %>% full_join(d3WEP_retro)
WEP.resid_multi = subset(WEP.resid, 
                         fitted.through == 2002 & Year == 2003 | 
                         fitted.through == 2009 & Year == 2010 | 
                         fitted.through == 2016 & Year == 2017) 

# Plot standardized residuals, S & D3 biomass in select years (Fig. 4):
ATF.resid_plot = ggplot(data=world) +
    geom_point(data=ATF.resid_multi, aes(x=lon, y=lat, group=model, col=resid.WT), position = "jitter") +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(legend.position = c(0.916, 0.068), 
          legend.key.height = unit(2.5, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.text.align = 1) +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) +  
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks=c(-175,-165))
ggsave(filename="Plots/Fig4_ATF.png", plot=ATF.resid_plot, dpi=500, width=5.4, height=7.4, units="in", bg="transparent")  
  
WEP.resid_plot = ggplot(data=world) +
    geom_point(data=WEP.resid_multi, aes(x=lon, y=lat, group=model, col=resid.WT), position = "jitter") +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(legend.position = c(0.916, 0.068), 
          legend.key.height = unit(2.5, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.text.align = 1) +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) +  
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks=c(-175,-165))
ggsave(filename="Plots/Fig4_WEP.png", plot=WEP.resid_plot, dpi=500, width=5.4, height=7.4, units="in", bg="transparent") 


# Plot standardized residuals, biomass for all years (Fig. S15):
ATF.resid_early = subset(ATF.resid, Year < 2001)
WT.resid.plot_ATF.early = ggplot(data=world) +
    geom_point(data=ATF.resid_early, aes(x=lon, y=lat, col=resid.WT), position = "jitter", size=0.4) +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=9),
          strip.text = element_text(size=10.5),
          legend.position = "none") +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) + 
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60",""))
ggsave(filename="Plots/FigS15_ATF_early.png", plot=WT.resid.plot_ATF.early, dpi=500, width=4.5, height=9, units="in", bg="transparent") 
  
WEP.resid_early = subset(WEP.resid, Year < 2001)
WT.resid.plot_WEP.early = ggplot(data=world) +
    geom_point(data=WEP.resid_early, aes(x=lon, y=lat, col=resid.WT), position = "jitter", size=0.4) +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=9),
          strip.text = element_text(size=10.5),
          legend.position = "none") +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) + 
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60",""))
ggsave(filename="Plots/FigS15_WEP_early.png", plot=WT.resid.plot_WEP.early, dpi=500, width=4.5, height=9, units="in", bg="transparent")  

ATF.resid_mid = subset(ATF.resid, Year >= 2001 & Year < 2010)
WT.resid.plot_ATF.mid = ggplot(data=world) +
    geom_point(data=ATF.resid_mid, aes(x=lon, y=lat, col=resid.WT), position = "jitter", size=0.4) +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=9),
          strip.text = element_text(size=10.5),
          legend.position = "none") +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) + 
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60",""))
ggsave(filename="Plots/FigS15_ATF_mid.png", plot=WT.resid.plot_ATF.mid, dpi=500, width=4.5, height=9, units="in", bg="transparent")  

WEP.resid_mid = subset(WEP.resid, Year >= 2001 & Year < 2010)
WT.resid.plot_WEP.mid = ggplot(data=world) +
    geom_point(data=WEP.resid_mid, aes(x=lon, y=lat, col=resid.WT), position = "jitter", size=0.4) +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=9),
          strip.text = element_text(size=10.5),
          legend.position = "none") +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) + 
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60",""))
ggsave(filename="Plots/FigS15_WEP_mid.png", plot=WT.resid.plot_WEP.mid, dpi=500, width=4.5, height=9, units="in", bg="transparent") 

ATF.resid_late = subset(ATF.resid, Year >= 2010)
WT.resid.plot_ATF.late = ggplot(data=world) +
    geom_point(data=ATF.resid_late, aes(x=lon, y=lat, col=resid.WT), position = "jitter", size=0.4) +
    scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=9),
          strip.text = element_text(size=10.5),
          legend.position = "none") +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) + 
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60",""))
ggsave(filename="Plots/FigS15_ATF_late.png", plot=WT.resid.plot_ATF.late, dpi=500, width=4.5, height=9, units="in", bg="transparent")  

WEP.resid_late = subset(WEP.resid, Year >= 2010)
WT.resid.plot_WEP.late = ggplot(data=world) +
    geom_point(data=WEP.resid_late, aes(x=lon, y=lat, col=resid.WT), position = "jitter", size=0.4) +
   scale_color_gradientn(colors=myPalette.rev, limits=c(-14,14), breaks=c(-12,-6,0,6,12), labels=c("-12", "", "0", "", "12"), na.value="#67001F", name = "") +
    plot.theme() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size=9),
          strip.text = element_text(size=10.5),
          legend.position = "none") +
    geom_path(data=NBS.proj, aes(x=long, y=lat, group=group), col="black", lwd=0.25) +
  facet_grid(Year ~ model) + 
  geom_sf(size=0.2) +    
    coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
    labs(x="Longitude", y="Latitude")  +
  scale_x_continuous(breaks = c(-175,-170,-165,-160), labels=c("","-170","","-160")) +
  scale_y_continuous(breaks = c(50,55,60,65), labels=c("50","","60",""))
ggsave(filename="Plots/FigS15_WEP_late.png", plot=WT.resid.plot_WEP.late, dpi=500, width=4.5, height=9, units="in", bg="transparent") 
##################################################################
# Bottom temperatures used to fit and forecast from each model type (Fig. S4):

# Data used to fit:
load("Data/dynamic_model_data.rda")
  dynamic.model.data$Year = as.numeric(dynamic.model.data$Year)

dynamic.model.data %>% 
  summarise(mean(ROMS_BT))
dynamic.model.data %>% 
  summarise(sd(ROMS_BT))
range(dynamic.model.data$ROMS_BT)

dynamic.model.data = dynamic.model.data[,c("Year", "lon", "lat", "Species", "adjN", "adjWT", "PA", "ROMS_BT")]

# ATF
ATF = dynamic.model.data %>%
  filter(Species == "Arrowtooth.Flounder")
ATF_catch = subset(ATF, adjN > 0 & adjWT > 0)

# WEP
WEP = dynamic.model.data %>%
  filter(Species == "Walleye.Pollock")
WEP_catch = subset(WEP, adjN > 0 & adjWT > 0)


# Data used to predict/forecast:
load("Data/s_forecasts_stack_10km.rda")
s.forecasts.10km %>% 
  summarise(mean(ROMS_BT))
s.forecasts.10km %>% 
  summarise(sd(ROMS_BT))
range(s.forecasts.10km$ROMS_BT)

load("Data/dyn_forecasts_stack_10km.rda")
d.forecasts.10km %>% 
  summarise(mean(ROMS_BT))
d.forecasts.10km %>% 
  summarise(sd(ROMS_BT))
range(d.forecasts.10km$ROMS_BT)

# presence-absence:
ATF.pa = NULL
for(i in 1992:2018) {
  ATF.pa.data = subset(ATF, Year < i)
    ATF.pa.data$fitted.through = i - 1
  ATF.pa = rbind(ATF.pa, ATF.pa.data) }
ATF.pa$modeled = "presence-absence"

# positive catches:
ATF.pos = NULL
for(i in 1992:2018) {
  ATF.pos.data = subset(ATF_catch, Year < i)
  ATF.pos.data$fitted.through = i - 1
ATF.pos = rbind(ATF.pos, ATF.pos.data) }
ATF.pos$modeled = "positive.catches"

ATF.fit = ATF.pa %>% full_join(ATF.pos)

ATF.fit %>% 
  group_by(modeled) %>%
  summarise(mean(ROMS_BT))
ATF.fit %>% 
  group_by(modeled) %>%
  summarise(sd(ROMS_BT))
ATF.fit %>% 
  group_by(modeled) %>%
  summarise(range(ROMS_BT))

ATF.fit$modeled = as.factor(ATF.fit$modeled)
ATF.fit$modeled = ordered(ATF.fit$modeled, levels = c("presence-absence", "positive.catches"))
levels(ATF.fit$modeled) = c("Presence / Absence", "Positive Catches")

# presence-absence:
WEP.pa = NULL
for(i in 1992:2018) {
  WEP.pa.data = subset(WEP, Year < i)
  WEP.pa.data$fitted.through = i - 1
WEP.pa = rbind(WEP.pa, WEP.pa.data) }
WEP.pa$modeled = "presence-absence"

# positive catches:
WEP.pos = NULL
for(i in 1992:2018) {
  WEP.pos.data = subset(WEP_catch, Year < i)
  WEP.pos.data$fitted.through = i - 1
WEP.pos = rbind(WEP.pos, WEP.pos.data) }
WEP.pos$modeled = "positive.catches"

WEP.fit = WEP.pa %>% full_join(WEP.pos)

WEP.fit %>% 
  group_by(modeled) %>%
  summarise(mean(ROMS_BT))
WEP.fit %>% 
  group_by(modeled) %>%
  summarise(sd(ROMS_BT))
WEP.fit %>% 
  group_by(modeled) %>%
  summarise(range(ROMS_BT))

WEP.fit$modeled = as.factor(WEP.fit$modeled)
WEP.fit$modeled = ordered(WEP.fit$modeled, levels = c("presence-absence", "positive.catches"))
levels(WEP.fit$modeled) = c("Presence / Absence", "Positive Catches")

# predicted/forecasted:
static.cov = bind_rows(replicate(2, s.forecasts.10km, simplify = F))
  static.cov$modeled = rep(c("presence-absence", "positive.catches"), 
                           each=nrow(s.forecasts.10km))
  static.cov$modeled = as.factor(static.cov$modeled)
static.cov$modeled = ordered(static.cov$modeled, 
                             levels = c("presence-absence", "positive.catches"))
levels(static.cov$modeled) = c("Presence / Absence", "Positive Catches")

dyn.cov = bind_rows(replicate(2, d.forecasts.10km, simplify = F))
  dyn.cov$modeled = rep(c("presence-absence", "positive.catches"), 
                        each=nrow(d.forecasts.10km))
  dyn.cov$modeled = as.factor(dyn.cov$modeled)
dyn.cov$modeled = ordered(dyn.cov$modeled, 
                          levels = c("presence-absence", "positive.catches"))

dyn.cov$fitted.through = dyn.cov$Year
  dyn.cov_red = subset(dyn.cov, Year > 1991)
  dyn.cov_red$fitted.through = dyn.cov_red$Year - 1
levels(dyn.cov_red$modeled) = c("Presence / Absence", "Positive Catches")

ATF.fit$spp = "ATF"
WEP.fit$spp = "WEP"
static.cov$spp = "Static"
dyn.cov_red$spp = "Dynamic"

temp.df = ATF.fit %>%
  full_join(WEP.fit) %>%
  full_join(., static.cov) %>%
  full_join(., dyn.cov_red)

temp.df$spp = as.factor(temp.df$spp)
  temp.df$spp = ordered(temp.df$spp, levels = c("ATF", "WEP", "Static", "Dynamic"))
  
temp.df$facet = as.factor(paste(temp.df$modeled, temp.df$spp, sep="."))
temp.data = temp.df[temp.df$facet %in% c("Presence / Absence.WEP",
                                         "Positive Catches.ATF", 
                                         "Positive Catches.WEP", 
                                         "Presence / Absence.Static",
                                         "Presence / Absence.Dynamic"), ]
temp.data$facet = droplevels(temp.data$facet) 
levels(temp.data$facet) = c("ATF", "WEP", "Dynamic", "Static", "All Tows")
  temp.data$facet = ordered(temp.data$facet, levels = c("All Tows", "WEP", "ATF", "Static", "Dynamic"))
  
# Plot fitted and forecasted bottom temperatures (Fig. 5):
obs_pred_temp_vio = ggplot(data = temp.data) +   
  geom_jitter(aes(x=facet, y=ROMS_BT, color=facet, alpha=facet), size=0.05, stroke=0.2, show.legend = F) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,0.85,0.5)) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values = c("darkgreen", "cadetblue3", "blue4", "gold", "red4"), name="") +
  geom_violin(aes(x=facet, y=ROMS_BT, color=facet, fill=facet), stat = "ydensity", alpha=0.3, adjust=2.2, trim=F) +
  scale_fill_manual(values = c("chartreuse3", "darkslategray1", "dodgerblue2", "yellow", "firebrick1"), name="") +
  plot.theme() +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(size=10),
        legend.position = "none") +
  scale_y_continuous(breaks = c(-2.5,0,2.5,5,7.5,10,12.5,15), labels=c("","0","","5","","10","","15")) +
  labs(x = "", y = expression(paste("Bottom Temperature (",degree,"C)")))
ggsave(filename="Plots/Fig5.png", plot=obs_pred_temp_vio, dpi=500, width=6.5, height=4.5, units="in")