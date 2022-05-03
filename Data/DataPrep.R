# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)

# Citation: 
# Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Model complexity has contrasting benefits for hindcasting and forecasting species responses to climate change. Ecogr. 

# References: 
# Hermann AJ, GA Gibson, NA Bond, EN Curchitser, K Hedstrom, W Cheng, M Wang, PJ Stabeno, L Eisner, and KD Cieciel. 2013. A multivariate analysis of observed and modeled biophysical variability on the Bering Sea shelf: multidecadal hindcasts (1970-2009) and forecasts (2010-2040). Deep-Sea Res II: Top Stud Oceanogr. 94:121–139.
# Hoff GR. 2016. Results of the 2016 eastern Bering Sea upper continental slope survey of groundfish and invertebrate resources. US Dep Comm NOAA Tech Memo. NMFS AFSC-339. 272 pp.
# Kearney K, A Hermann, W Cheng, I Ortiz, and K Aydin. 2020. A coupled pelagic-benthic-sympagic biogeochemical model for the Bering Sea: documentation and validation of the BESTNPZ model (v2019.08.23) within a high-resolution regional ocean model. Geosci Model Dev. 13:597-650.
# Laman EA, CN Rooper, K Turner, S Rooney, DW Cooper, and M Zimmerman. 2018. Using species distribution models to describe essential fish habitat in Alaska. Can J Fish Aquat Sci. 75:1230–1255. 
# Lauth R and E Acuna. 2007. 2005 bottom trawl survey of the eastern Bering Sea continental shelf. AFSC Processed Rep. US Dept Comm. 2007-1. 164 pp. 
# Lauth RR, EJ Dawson, and J Connor. 2019. Results of the 2017 eastern and northern Bering Sea continental shelf bottom trawl survey of groundfish and invertebrate fauna. NOAA Tech Mem 396. 270 pp.
# Pirtle JL, SK Shotwell, M Zimmermann, JA Reid, and N Golden. 2019. Habitat suitability models for groundfish in the Gulf of Alaska. Deep-Sea Res. Pt. II. 165: 303–321
# Rooper CN, MF Sigler, P Goddard, P Malecha, R Towler, K Williams, R Wilborn, and M Zimmermann. 2016. Validation and improvement of species distribution models for structure-forming invertebrates in the eastern Bering Sea with an independent survey. Mar Ecol Prog Ser. 551:117-130.
# Sigler MF, CN Rooper, GR Hoff, RP Stone, RA McConnaughey, and TK Wilderbuer. 2015. Faunal features of submarine canyons on the eastern Bering Sea slope. Mar Ecol Prog Ser. 526:21-40.
# Stahl JP and GH Kruse. 2008. Spatial and temporal variability in size at maturity of Walleye Pollock in the Eastern Bering Sea. Trans Amer Fish Soc. 137:1543–1557.
# Stevenson DE and RR Lauth. 2019. Bottom trawl surveys in the northern Bering Sea indicate recent shifts in the distribution of marine species. Polar Biol. 42:407–421.
# Zimmermann M, MM Prescott, and CN Rooper. 2013. Smooth Sheet Bathymetry of the Aleutian Islands. US Dep Comm NOAA Tech Memo. NMFS-AFSC-250. 43pp.
# Zimmermann M, and M Prescott. 2018. Bathymetry and Canyons of the Eastern Bering Sea Slope. Geosciences. 8:184.


# Script Objectives: 1) Adjust survey data to reflect region (Bering Sea), species (Arrowtooth Flounder and Walleye Pollock), and life stage (adults) of interest. 2) Prepare all model covariates (static and/or dynamic) to join with bottom trawl survey data (for model fitting) and uniform grid (for model predictions). Spatially-explicit but time-invariant habitat variables were provided by the Alaska Regional Office, NOAA Fisheries. Regional Ocean Modeling System (ROMS) hindcasts of bottom temperature (BT; degrees C) cold pool extent (CPE; sq. km) were provided by the Alaska Climate Integrated Modeling (ACLIM) Project (see https://github.com/kholsman/ACLIM2 for more information). 

setwd("~/Documents/JISAO/GitHub/")

require(tidyverse)
require(tidyr)
require(dplyr)
require(sf)
require(sp)
require(rgeos)
require(raster)
require(rgdal)
require(gstat)
require(rasterVis)
require(mapdata)
require(rnaturalearth)
require(rnaturalearthdata)
require(rnaturalearthhires)
require(ggplot2)
require(Hmisc)
require(corrplot)
require(mgcv)

# Plotting info:
lon_min = -178.65
lon_max = -156.5
lat_min = 50
lat_max = 65.9

data(worldHiresMapEnv)
world = ne_countries(scale = "large", returnclass = "sf")

plot.theme = function() {
  theme_bw() +
  theme(panel.border = element_rect(color="black"),
        panel.grid = element_blank(),
        legend.key = element_blank(), 
        legend.title = element_text(family="Arial", size=12),
        legend.text = element_text(family="Arial", size=10),
        legend.text.align = 1,
        legend.key.width = unit(1.5, "mm"),
        legend.key.height = unit(2.5, "mm"),
        legend.background = element_rect(fill="transparent"), 
        axis.title.x = element_text(vjust=-0.1, hjust=0.51, family="Arial", size=12),
        axis.title.y = element_text(vjust=2.0, hjust=0.485, family="Arial", size=12),
        axis.text = element_text(family="Arial", size=11, color="black"),
        strip.background = element_blank(),
        strip.text = element_text(family="Arial", size=12)) }  

##################################################################
# Import and prepare survey data (for model fitting):
haul.SBS = read.csv("Data/RACEBase_HAUL.csv") # SE Bering Sea
haul.SBS = haul.SBS %>%
  mutate(across(c(DURATION, DISTANCE_FISHED, END_LATITUDE, END_LONGITUDE, GEAR_DEPTH, BOTTOM_DEPTH, WIRE_LENGTH, ACCESSORIES), as.numeric)) %>%
    filter(ABUNDANCE_HAUL == "Y")
haul.SBS = subset(haul.SBS, select = -c(ABUNDANCE_HAUL))

haul.NBS = read.csv("Data/NBS_HAULS_not_incl.csv") # N Bering Sea
  haul.NBS = subset(haul.NBS, select = -c(YEAR))

haul = rbind(haul.SBS, haul.NBS) 

catch.SBS = read.csv(unz("Data/RACEBase_CATCH.csv.zip", 
                  "RACEBase_CATCH.csv"), header=T) 
catch.NBS = read.csv("Data/NBS_CATCH_not_incl.csv") # N Bering Sea
  catch = rbind(catch.SBS, catch.NBS)

species = read.csv("Data/RACEBase_SPECIES.csv")

# Create a unique identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to length data below):
haul$ID = paste(haul$VESSEL, 
                haul$CRUISE, 
                haul$HAUL, sep="")
catch$ID = paste(catch$VESSEL, 
                catch$CRUISE, 
                catch$HAUL, sep="")

# Create year column:
haul$Year = haul$CRUISE
haul = haul %>%
  tidyr::separate(Year, into = c("Year", "Cruise"), sep = 4)

catch$Year = catch$CRUISE
catch = catch %>%
  tidyr::separate(Year, into = c("Year", "Cruise"), sep = 4)

# Join and subset for region, tow type, and years of interest:
trawl = haul %>% 
  left_join(catch, by = c("ID", "REGION", "Year", "CRUISEJOIN", "HAULJOIN", "VESSEL", "CRUISE", "HAUL", "Cruise")) %>%
  filter(REGION == "BS" & 
         Year >= 1982 & Year < 2019) %>%
  rename(lon = START_LONGITUDE, lat = START_LATITUDE)

# Set correct classes:
trawl$REGION = as.factor(trawl$REGION)
trawl$Year = as.integer(trawl$Year)
trawl$NET_WIDTH = as.numeric(trawl$NET_WIDTH)
trawl$NUMBER_FISH = as.numeric(trawl$NUMBER_FISH)

# Estimate area swept (sq. km) for model offset:
  # Predict net width from wire length where NA:
net.sub = subset(trawl, NET_WIDTH > 0 & WIRE_LENGTH > 0)
net.fit = nls(NET_WIDTH ~ a - (b/WIRE_LENGTH), data = net.sub, start=list(a=19.56, b=605.23)); summary(net.fit); net.fit$m$getAllPars() 
  # starting values from 2005 cruise report (see https://www.arlis.org/docs/vol1/81343364.pdf for details)
net.sub$net.pred = fitted.values(net.fit)

trawl$net.pred = with(trawl, 
        net.fit$m$getAllPars()[1] - (net.fit$m$getAllPars()[2]/WIRE_LENGTH))
trawl = trawl %>%
  mutate(Net.Width = ifelse(is.na(NET_WIDTH), net.pred, NET_WIDTH))
trawl$AreaSwept = with(trawl, (Net.Width / 1000) * DISTANCE_FISHED)
  # net width (m); distance fished (km)

# Join haul and catch data with species codes:
survey = trawl %>% left_join(species)

# Relabel species of interest and group all others (for presence-absence):
survey$Species = with(survey,
  ifelse(SPECIES_CODE == 10112, "Arrowtooth.Flounder", # Kamchatka Flounder
  ifelse(SPECIES_CODE == 10111, "Arrowtooth.Flounder", # Atheresthes spp.
  ifelse(SPECIES_CODE == 10110, "Arrowtooth.Flounder", # Arrowtooth Flounder
  ifelse(SPECIES_CODE == 21740, "Walleye.Pollock", "Other")))))

# Expand data frame to make implicit missing values explicit:
all.hauls = survey %>% 
  tidyr::expand(nesting(Year, ID, lon, lat, AreaSwept), Species)
survey.hauls = survey %>% right_join(all.hauls, all.y = T)

# Combine WEIGHTS and NO. FISH by haul, year, and species:
survey.all = survey.hauls %>%
  filter(Species != "Other") %>%
  group_by(ID, Year, Species) %>%
    mutate(WEIGHT.all = sum(WEIGHT)) %>%
    mutate(NUMBER.FISH.all = sum(NUMBER_FISH)) 

survey.data = survey.all %>%
  dplyr::select(-c(WEIGHT, NUMBER_FISH, SPECIES_CODE, COMMON_NAME, SPECIES_NAME, AUDITJOIN.x, AUDITJOIN.y, CATCHJOIN, VOUCHER))
haul.catch.data = distinct(survey.data)

# Fill NA values with 0:
haul.catch.data$WEIGHT.all[is.na(haul.catch.data$WEIGHT.all)] = 0
haul.catch.data$NUMBER.FISH.all[is.na(haul.catch.data$NUMBER.FISH.all)] = 0

##################################################################
# Adjust catches to reflect adult fish only (Bering Sea):
  # Arrowtooth Flounder > 480 mm and Walleye Pollock > 381 mm
       
# Estimate length-weight relationships from subsampled fishes:
  # Read in specimen data (i.e., length-weight-age tables):
LWA.SBS = read.csv(unz("Data/racebase_specimen.csv.zip", 
                  "racebase_specimen.csv"), header=T) # standard data
LWA.NBS = read.csv("Data/NBS_SPECIMEN_1982_2019.csv") # NBS data
  LWA.NBS = LWA.NBS %>%
    rename(SPECIMEN_ID = SPECIMENID) #match SE Bering Sea
  LWA.all = rbind(LWA.SBS, LWA.NBS) # join SE and N Bering Sea data

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to survey data above):
LWA.all$ID = paste(LWA.all$VESSEL, 
                   LWA.all$CRUISE, 
                   LWA.all$HAUL, sep="")

# Create year column:
LWA.all$Year = LWA.all$CRUISE
LWA.all = LWA.all %>%
  tidyr::separate(Year, into = c("Year", "Cruise"), sep = 4)
LWA.all$Year = as.integer(LWA.all$Year)

# Subset for region and years of interest:
LWA.data = LWA.all %>% 
  filter(REGION == "BS" & Year >= 1982 & Year < 2019)

# Join length-weight-age data and species codes:
LWA = LWA.data %>% 
  left_join(species) 

# Relabel species of interest and group all others:
LWA$Species = with(LWA,
  ifelse(SPECIES_CODE == 10112, "Arrowtooth.Flounder", # Kamchatka Flounder
  ifelse(SPECIES_CODE == 10111, "Arrowtooth.Flounder", # Atheresthes spp.
  ifelse(SPECIES_CODE == 10110, "Arrowtooth.Flounder", # Arrowtooth Flounder
  ifelse(SPECIES_CODE == 21740, "Walleye.Pollock", "Other")))))

# Remove species not of interest:
LWA.spp = subset(LWA, Species != "Other")

# Estimate weight (g)-length (cm) parameters using allometric growth model with bias-correction (Brodziak 2012):
LWA.spp$WEIGHT = as.numeric(LWA.spp$WEIGHT)
  LWA.spp$logW = log(as.numeric(LWA.spp$WEIGHT))
LWA.spp$logL = log(LWA.spp$LENGTH)
  LWA.spp = subset(LWA.spp, logW != "NA")

# Remove weights and lengths with NAs or zeros:
LWA.spp = subset(LWA.spp, LENGTH > 0 & WEIGHT > 0)

for(i in unique(LWA.spp$Species)) {
  sp = subset(LWA.spp, Species == i)
L_W = lm(logW ~ logL, data = sp)
    a = L_W$coefficients[1]
    b = L_W$coefficients[2]
    syx = summary(L_W)$sigma
    cf = exp((syx^2)/2) # correction factor for predicting on orig. scale

sp$pred.WT = predict(L_W, data.frame(logL=sp$logL), interval="c") 
sp$pred.WT.corr = cf *(exp(sp$pred.WT))
  assign(paste(i, "_", "pars", sep=""), c(a, b, cf))
  assign(paste(i), sp) }

##################################################################
# Read in and format haul-specific length compositions:
lengths.SBS = read.csv(unz("Data/RACEBase_length.csv.zip", 
                  "RACEBase_length.csv"), header=T)  # standard survey data
  colnames(lengths.SBS) = toupper(colnames(lengths.SBS))
lengths.NBS = read.csv("Data/NBS_LENGTH_1982_2019.csv") # NBS data
  lengths.all = rbind(lengths.SBS, lengths.NBS) # join standard and NBS data

# Create a unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to survey data above):
lengths.all$ID = paste(lengths.all$VESSEL, 
                       lengths.all$CRUISE, 
                       lengths.all$HAUL, sep="")

# Create year column:
lengths.all$Year = lengths.all$CRUISE
lengths.all = lengths.all %>%
  tidyr::separate(Year, into = c("Year", "Cruise"), sep = 4)

# Subset for region and years of interest:
lengths.all = lengths.all %>% 
  filter(REGION == "BS" & Year >= 1982 & Year < 2019)

# Join length-weight-age data and species codes:
length.data = lengths.all %>% 
  left_join(species) 

# Relabel species of interest and group all others:
length.data$Species = with(length.data,
  ifelse(SPECIES_CODE == 10112, "Arrowtooth.Flounder", # Kamchatka Flounder
  ifelse(SPECIES_CODE == 10111, "Arrowtooth.Flounder", # Atheresthes spp.
  ifelse(SPECIES_CODE == 10110, "Arrowtooth.Flounder", # Arrowtooth Flounder
  ifelse(SPECIES_CODE == 21740, "Walleye.Pollock", "Other")))))
length.data$Species = as.factor(length.data$Species)

# Remove species not of interest:
length.data = subset(length.data, Species != "Other")

# Create fork length bins (mm) separating adults from other life stages:
length.data$LifeStage = with(length.data, 
ifelse(Species == "Walleye.Pollock" & LENGTH > 381, "adult",
ifelse(Species == "Arrowtooth.Flounder" & LENGTH > 480, "adult", "other"))) 

# Predict weight (g) from length (mm):
length.data$WT.g = with(length.data, 
ifelse(Species == "Walleye.Pollock", 
       Walleye.Pollock_pars[3] * (exp(Walleye.Pollock_pars[1] + (Walleye.Pollock_pars[2]*(log(LENGTH))))),
ifelse(Species == "Arrowtooth.Flounder", 
       Arrowtooth.Flounder_pars[3] * (exp(Arrowtooth.Flounder_pars[1] + (Arrowtooth.Flounder_pars[2]*(log(LENGTH))))), NA)))

##################################################################
# Calculate proportions of fish (by weight) sampled in each life stage and haul:

# Convert WT to kg and calculate total weight per haul:
length.data$WT.kg = length.data$WT / 1000
length.data$WT.tot = with(length.data, (FREQUENCY * WT.kg))


# Total weight (all stages), by haul:
catch.all.WT = length.data %>%
  group_by(Year, ID, Species) %>%
  summarise(Total.WT = sum(WT.tot)) %>%
  as.data.frame()

# Total weight (adults only), by haul:
catch.adults.WT = length.data %>%
  filter(LifeStage == "adult") %>%
  group_by(Year, ID, Species) %>%
  summarise(Adult.WT = sum(WT.tot)) %>%
  as.data.frame()

# Calculate proportion of haul made up of adults:
catch.prop.WT = catch.all.WT %>% 
  full_join(catch.adults.WT) %>%
  as.data.frame()

catch.prop.WT$Adult.WT[is.na(catch.prop.WT$Adult.WT)] = 0
catch.prop.WT$Adult.WT.prop = with(catch.prop.WT, (Adult.WT / Total.WT))
catch.prop.WT$Year = as.integer(catch.prop.WT$Year)
  

# Total N (all stages), by haul:
catch.all.N = length.data %>%
  group_by(Year, ID, Species) %>%
  summarise(Total.N = sum(FREQUENCY)) %>%
  as.data.frame()

# Total N (adults only), by haul:
catch.adults.N = length.data %>%
  filter(LifeStage == "adult") %>%
  group_by(Year, ID, Species) %>%
  summarise(Adult.N = sum(FREQUENCY)) %>%
  as.data.frame()

# Calculate proportion of haul made up of adults:
catch.prop.N = catch.all.N %>% 
  full_join(catch.adults.N) %>%
  as.data.frame()

catch.prop.N$Adult.N[is.na(catch.prop.N$Adult.N)] = 0
catch.prop.N$Adult.N.prop = with(catch.prop.N, (Adult.N / Total.N))
catch.prop.N$Year = as.integer(catch.prop.N$Year)

##################################################################
# Join bottom trawl survey and proportional length and weight data:
size.haul.catch = haul.catch.data %>% 
    left_join(catch.prop.WT, all=T) %>% 
    left_join(.,catch.prop.N, all=T) %>% 
    as.data.frame()
table(size.haul.catch$Species) # check

# Assign zeros to adjusted values where no fish were caught:
size.haul.catch$Adult.WT.prop = with(size.haul.catch,
  ifelse(WEIGHT.all == 0 & NUMBER.FISH.all == 0, 0, Adult.WT.prop))
size.haul.catch$Adult.N.prop = with(size.haul.catch,
  ifelse(NUMBER.FISH.all == 0 & WEIGHT.all == 0, 0, Adult.N.prop))

# Remove tows without fish measurements (can't estimate proportion of adults):
size.catch.adults = size.haul.catch %>%
  filter(!is.na(Adult.WT.prop)) %>%
  filter(!is.na(Adult.N.prop))
size.catch.adults$Species = as.factor(size.catch.adults$Species)
table(size.catch.adults$Species) # check

# Adjust catch based on proportions of adults subsampled for lengths:
size.catch.adults$adjWT = with(size.catch.adults, 
            round(WEIGHT.all * Adult.WT.prop, digits=3))
size.catch.adults$adjN = with(size.catch.adults, 
            round(NUMBER.FISH.all * Adult.N.prop, digits=0))

# Label presence-absence for binomial models:
size.catch.adults$PA = as.numeric(size.catch.adults$adjN > 0)
    save(size.catch.adults, file = "Data/size_catch_data.rda")
  # load("Data/size_catch_data.rda")
    
##################################################################
# Import and prepare covariate data:
setwd("~/Documents/JISAO/GitHub/Data/")
bathy = raster::raster("bathy", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") 
slope = raster::raster("slope", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
BPI = raster::raster("BPI", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") 
phi = raster::raster("phi", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
sponges = raster::raster("sponge_factor", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
corals = raster::raster("coral_factor", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
whips = raster::raster("whips_factor", 
  crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

lat = raster::init(bathy, v ='y')
lat = raster::mask(lat, bathy, filename = "Lat",
     overwrite = T)
lon = raster::init(bathy, v ='x')
lon = raster::mask(lon, bathy, filename = "Lon",
     overwrite = T)

bathy_proj = as.data.frame(rasterToPoints(projectRaster(bathy, crs="+init=epsg:4326"))); colnames(bathy_proj) = c("lon", "lat", "bathy")
slope_proj = as.data.frame(rasterToPoints(projectRaster(slope, crs="+init=epsg:4326"))); colnames(slope_proj) = c("lon", "lat", "slope")
BPI_proj = as.data.frame(rasterToPoints(projectRaster(BPI, crs="+init=epsg:4326"))); colnames(BPI_proj) = c("lon", "lat", "BPI")
phi_proj = as.data.frame(rasterToPoints(projectRaster(phi, crs="+init=epsg:4326"))); colnames(phi_proj) = c("lon", "lat", "phi")
sponges_proj = as.data.frame(rasterToPoints(projectRaster(sponges, crs="+init=epsg:4326"))); colnames(sponges_proj) = c("lon", "lat", "sponges")
corals_proj = as.data.frame(rasterToPoints(projectRaster(corals, crs="+init=epsg:4326"))); colnames(corals_proj) = c("lon", "lat", "corals")
whips_proj = as.data.frame(rasterToPoints(projectRaster(whips, crs="+init=epsg:4326"))); colnames(whips_proj) = c("lon", "lat", "whips")

# Combine rasters with same exact locations:
bathy.phi = bathy_proj %>% full_join(phi_proj) 
SFI = cbind(sponges_proj, corals_proj, whips_proj) 
  SFI = SFI[,c("lon", "lat", "sponges", "corals", "whips")]

# Incrementally combine rasters with slightly different locations:
  bathy.phi_sf = st_as_sf(bathy.phi, coords = c("lon", "lat"))
  slope_sf = st_as_sf(slope_proj, coords = c("lon", "lat"))
terrain.cov = st_join(bathy.phi_sf, slope_sf, 
                      join = st_nearest_feature, left = T)

  BPI_sf = st_as_sf(BPI_proj, coords = c("lon", "lat"))
s.covariates = st_join(terrain.cov, BPI_sf, 
                         join = st_nearest_feature, left = T)

  SFI_sf = st_as_sf(SFI, coords = c("lon", "lat"))
covariate.rasters = st_join(s.covariates, SFI_sf, 
                            join = st_nearest_feature, left = T)

# Reformat SFI covariates as factors:
covariate.rasters$sponges = as.factor(as.character(as.numeric(covariate.rasters$sponges > 0)))
covariate.rasters$corals = as.factor(as.character(as.numeric(covariate.rasters$corals > 0)))
covariate.rasters$whips = as.factor(as.character(as.numeric(covariate.rasters$whips > 0)))

  save(covariate.rasters, file = "Data/covariate_rasters.rda")
# load("Data/covariate_rasters.rda")

# Plot (Fig. S1):
setwd("~/Documents/JISAO/GitHub/")  

bathy_proj$depth = ifelse(bathy_proj$bathy > 1000, 1000, bathy_proj$bathy)
Bathy = ggplot(data=world) +
  geom_tile(data=bathy_proj, aes(x=lon, y=lat, fill=depth)) +
  scale_fill_gradientn(colors = c("black", "darkorchid4", "darkslateblue", "dodgerblue1", "cyan3", "darkslategray1", "aliceblue"), breaks = c(250,500,750,1000), labels = c("250","","750", ""), trans="reverse") + 
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim = c(-178.65, -156.5), ylim = c(50, 65.9)) +
  plot.theme() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank()) +
  labs(x = "", y = "")
ggsave(filename="Plots/covar_depth.png", plot=Bathy, dpi=500, width=6.5, height=9, units="in")

slope_proj$slope_mod = ifelse(slope_proj$slope > 3, 3, slope_proj$slope)
Slope = ggplot(data=world) +
  geom_tile(data=slope_proj, aes(x=lon, y=lat, fill=slope_mod)) +
  scale_fill_gradientn(colors = c("lavender", "mediumpurple1", "mediumpurple4"), limits=c(0,3), breaks = c(0,1,2,3), labels=c("","1", "", "3+")) + 
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim = c(-178.65, -156.5), ylim = c(50, 65.9)) +
  plot.theme() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank()) +
  labs(x = "", y = "")
ggsave(filename="Plots/covar_slope.png", plot=Slope, dpi=500, width=6.5, height=9, units="in")

BPI_proj$bpi = ifelse(BPI_proj$BPI > 500, 500, BPI_proj$BPI)
bpi = ggplot(data=world) +
  geom_tile(data=BPI_proj, aes(x=lon, y=lat, fill=bpi)) +
  scale_fill_gradient2(low="mediumblue", mid="seashell", high="firebrick3", midpoint=0, limits = c(-345,500), breaks = c(-250,0,250,500), labels = c("","0","","500+")) + 
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim = c(-178.65, -156.5), ylim = c(50, 65.9)) +
  plot.theme() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank()) +
  labs(x = "", y = "")
ggsave(filename="Plots/covar_BPI.png", plot=bpi, dpi=500, width=6.5, height=9, units="in")

sponges_proj$Sponges = ifelse(sponges_proj$sponges > 0, 1, 0)
sponges_proj$Sponges = as.factor(as.character(sponges_proj$Sponges))

Sponges = ggplot(data=world) +
  geom_tile(data=sponges_proj, aes(x=lon, y=lat, fill=Sponges)) +
  scale_fill_manual(values = c("cornsilk", "chocolate3"), labels = c("absent", "present")) +  
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim = c(-178.65, -156.5), ylim = c(50, 65.9)) +
  plot.theme() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank()) +
  labs(x = "", y = "")
ggsave(filename="Plots/covar_sponges.png", plot=Sponges, dpi=500, width=6.5, height=9, units="in")

corals_proj$Corals = ifelse(corals_proj$corals > 0, 1, 0)
corals_proj$Corals = as.factor(as.character(corals_proj$Corals))

Corals = ggplot(data=world) +
  geom_tile(data=corals_proj, aes(x=lon, y=lat, fill=Corals)) +
  scale_fill_manual(values = c("cornsilk", "chocolate3"), labels = c("absent", "present")) +  
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim = c(-178.65, -156.5), ylim = c(50, 65.9)) +
  plot.theme() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank()) +
  labs(x = "", y = "")
ggsave(filename="Plots/covar_corals.png", plot=Corals, dpi=500, width=6.5, height=9, units="in")

whips_proj$Whips = ifelse(whips_proj$whips > 0, 1, 0)
whips_proj$Whips = as.factor(as.character(whips_proj$Whips))

Whips = ggplot(data=world) +
  geom_tile(data=whips_proj, aes(x=lon, y=lat, fill=Whips)) +
  scale_fill_manual(values = c("cornsilk", "chocolate3"), labels = c("absent", "present")) +  
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim = c(-178.65, -156.5), ylim = c(50, 65.9)) +
  plot.theme() +
    theme(legend.position = "none",
          axis.title = element_blank(), 
          axis.text = element_blank()) +
  labs(x = "", y = "")
ggsave(filename="Plots/covar_whips.png", plot=Whips, dpi=500, width=6.5, height=9, units="in")

##################################################################
# Import ROMS hindcasts:
load("Data/ROMS_NPZ_output.Rdata")
load("Data/ROMS_station_metadata.Rdata")

# Extract and format bottom trawl survey-replicated station metadata:
meta = ROMSNPZ_output$aclim_hindcast_201904$station_metadata
  station_meta = meta$surveyRep_StationsOnly
all_stations_meta$SRVY = "NEBS"
all_stations_meta$SRVY[all_stations_meta$station %in% 
      as.character(station_meta$station)] = "SEBS"
all_stations_meta$SRVY[all_stations_meta$SRVY == "NEBS" & all_stations_meta$lat < 60.5] = "SEBS_Resamp"
all_stations_meta$SRVY = factor(all_stations_meta$SRVY, 
                                levels = c("SEBS","NEBS","SEBS_Resamp"))
    
# Select the most updated hindcast run at the time of analyses (i.e., 2019):
hind = ROMSNPZ_output$aclim_hindcast_201904


# Extract and format BT and CPE hindcasts:
BT = hind$station$BottomTemp
  BT = as_tibble(cbind(all_stations_meta, BT))
CPE = hind$station$ColdPool
  CPE = as_tibble(cbind(all_stations_meta, CPE))
  
# Remove erroneous southern point:
BT = BT[- which(BT$SRVY == "SEBS_Resamp" & BT$lat < 56), ]
CPE = CPE[- which(CPE$SRVY == "SEBS_Resamp" & CPE$lat < 56), ]

# Change data from wide to long format and exclude "survey resampling" data:
BT.all = BT %>% 
  pivot_longer(names_to = "YEAR", values_to = "Temp", "1970":"2018")
BT.surveyRep = BT.all %>%
  filter(SRVY %in% c("NEBS","SEBS"))

CPE.all = CPE %>% 
  pivot_longer(names_to = "YEAR", values_to = "ColdPool", "1970":"2018")
CPE.surveyRep = CPE.all %>%
  filter(SRVY %in% c("NEBS","SEBS"))

# Select and order years of interest (i.e., those with standardized survey methods):
BT.surveyRep$Year = as.integer(BT.surveyRep$YEAR)
  BT.surveyRep = subset(BT.surveyRep, Year > 1981)
BT.surveyRep = BT.surveyRep[order(as.numeric(BT.surveyRep$Year)), ]

CPE.surveyRep$Year = as.integer(CPE.surveyRep$YEAR)
  CPE_surveyRep = subset(CPE.surveyRep, Year > 1981)
CPE.surveyRep = CPE.surveyRep[order(as.numeric(CPE.surveyRep$Year)), ]

# Reformat longitude values:
BT.surveyRep$lon = BT.surveyRep$lon - 360
CPE.surveyRep$lon = CPE.surveyRep$lon - 360

################################################################## 
# BT, Static Models:

# Estimate long-term mean BT by location:
ROMS.BT.mean = BT.surveyRep %>%
  group_by(lon, lat) %>%
  summarise(BT = mean(Temp)) 
ROMS.BT.mean = as.data.frame(ROMS.BT.mean)

ROMS.BT.sd = BT.surveyRep %>%
  group_by(lon, lat) %>%
  summarise(BT = sd(Temp)) 
ROMS.BT.sd = as.data.frame(ROMS.BT.sd)

# Plot (Fig. S2):
BT.mean = ggplot(data=world) +
  geom_point(data=ROMS.BT.mean, aes(x=lon, y=lat, color=BT), size=2) +
  scale_color_gradientn(colors=c("darkslateblue", "dodgerblue1", "green3", "yellow1", "chocolate1", "firebrick3"), name = expression(bar(BT)~degree*C), limits=c(-1.4,11.5), breaks=c(0,2.5,5,7.5,10), labels=c("0","","5","","10")) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  labs(x="Longitude", y="Latitude") +
  theme(legend.position = c(0.875, 0.14))
ggsave(BT.mean, filename="Plots/ROMS_BT_mean.png", dpi=500, width=3.2, height=4.5, units="in", bg="transparent")

BT.sd = ggplot(data=world) +
  geom_point(data=ROMS.BT.sd, aes(x=lon, y=lat, color=BT), size=2) +
  scale_color_gradientn(colors=c("darkslateblue", "dodgerblue1", "green3", "yellow1", "chocolate1", "firebrick3"), name = expression(BT^~sd~degree*C), limits=c(-0.05,3.05), breaks=c(0,1,2,3), labels=c("0","","2","")) +
  geom_sf(size=0.2) +    
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  plot.theme() +
  labs(x="Longitude", y="Latitude", color="Prob.") +
  theme(legend.position = c(0.846, 0.143),
        legend.title.align = 0)
ggsave(BT.sd, filename="Plots/ROMS_BT_sd.png", dpi=500, width=3.2, height=4.5, units="in", bg="transparent")


# Interpolate BT using ordinary kriging and uniform prediction grid:
grid = read.csv("Data/prediction_grid.csv")

# Format as a spatial data frame:
grid_sf = st_as_sf(grid, coords = c("lon", "lat"))
  coordinates(grid) = ~ lon + lat
crs(grid) = "+init=epsg:4326" # WGS84
  coordinates(ROMS.BT.mean) = ~ lon + lat # format as spatial data frame
crs(ROMS.BT.mean) = "+init=epsg:4326" # WGS84

BT.vgm = variogram(BT ~ 1, data = ROMS.BT.mean)
BT.fit = fit.variogram(BT.vgm, model=vgm(psill=8, model="Gau", 
                          range=400, nugget=0.25)); plot(BT.vgm, BT.fit)
ROMS.BT.krig = as.data.frame(krige(BT ~ 1, ROMS.BT.mean, grid, model=BT.fit))

# Rename interpolated BT:
ROMS.BT.krig = ROMS.BT.krig %>%
  rename(ROMS_BT = var1.pred, 
         ROMS_BTvar = var1.var)

  save(ROMS.BT.krig, file = "Data/ROMS_BT_krig_mean.rda")
# load("Data/ROMS_BT_krig_mean.rda") 
        
# Plot:
BT_ROMS_krig = ggplot(data=world) +
  geom_tile(data=ROMS.BT.krig, aes(x=lon, y=lat, fill=ROMS_BT)) +
  scale_fill_gradientn(colors=c("darkslateblue", "dodgerblue1", "green3", "yellow1", "chocolate1", "firebrick3"), limits=c(-2.5,15.4), breaks=c(0,2,4,6,8,10,12), labels=c("0","","4","","8","","12")) +
    geom_sf(size=0.2, color="black") +
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  labs(x="", y="") +
  plot.theme() +
  theme(legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()) 
ggsave(BT_ROMS_krig, filename="Plots/ROMS_BT_krig_mean.png", dpi=500, width=6.5, height=9, units="in")


# Link interpolated long-term mean BT to static covariate rasters:
  # Note: No index of cold pool extent for static models. 
cov.rasters_sf = st_as_sf(covariate.rasters, coords = c("lon", "lat"))
  st_crs(cov.rasters_sf) = 4326 # WGS84
ROMS.BT.krig_sf = st_as_sf(ROMS.BT.krig, coords = c("lon", "lat"))
  st_crs(ROMS.BT.krig_sf) = 4326 # WGS84

static.cov_sf = st_join(cov.rasters_sf, ROMS.BT.krig_sf, 
                               join = st_nearest_feature, left = T)

# Join covariate and survey data (for model fitting):
haul_sf = st_as_sf(size.catch.adults, coords = c("lon", "lat"))
  st_crs(haul_sf) = 4326 # WGS84
  
haul.static.cov_sf = st_join(haul_sf, static.cov_sf, 
                               join = st_nearest_feature, left = T)
static.model.data = as.data.frame(as_Spatial(haul.static.cov_sf))

# Rename coordinates:
static.model.data = static.model.data %>% 
  rename(lon = coords.x1, lat = coords.x2)

  save(static.model.data, file = "Data/static_model_data.rda")
# load("Data/static_model_data.rda")
  
#################################################################
# BT and CPE, Dynamic Models:
  
# Format and convert ROMS hindcasts to a spatial data frame:
ROMS.BT_yr = BT.surveyRep
ROMS.BT_yr$BT = ROMS.BT_yr$Temp
coordinates(ROMS.BT_yr) = ~ lon + lat
  crs(ROMS.BT_yr) = "+init=epsg:4326" # WGS84

# Interpolate year-specific BT:
tmp = list()
for(i in unique(ROMS.BT_yr$Year)) {
  BTi_df = subset(ROMS.BT_yr, Year == i)
BT.yr.vgm = variogram(BT ~ 1, data = BTi_df)
BT.yr.fit = fit.variogram(BT.yr.vgm, model=vgm(psill=8, model="Gau", range=400, nugget=0.25))
BT.yr_krig = krige(BT ~ 1, BTi_df, grid, model=BT.yr.fit)
  BT.yr_krig = as.data.frame(BT.yr_krig)
BT.yr_krig$Year = i
tmp[[i]] = BT.yr_krig }

ROMS.BT.krig_yr = dplyr::bind_rows(tmp)
  
# Rename interpolated BT:
ROMS.BT.krig_yr = ROMS.BT.krig_yr %>% 
  rename(ROMS_BT = var1.pred, 
         ROMS_BTvar = var1.var)

  save(ROMS.BT.krig_yr, file = "Data/ROMS_BT_krig_yr.rda")
# load("Data/ROMS_BT_krig_yr.rda")

# Plot:
for(i in unique(ROMS.BT.krig_yr$Year)) {
  df = subset(ROMS.BT.krig_yr, Year == i)  
  
BT_ROMS_krig.yr = ggplot(data=world) +
  geom_tile(data=df, aes(x=lon, y=lat, fill=ROMS_BT)) +
  scale_fill_gradientn(colors=c("darkslateblue", "dodgerblue1", "green3", "yellow1", "chocolate1", "firebrick3"), 
                       limits=c(-2.5,15.4), breaks=c(0,2,4,6,8,10,12), labels=c("0","","4","","8","","12")) +
  geom_sf(size=0.2, color="black") +
  coord_sf(xlim=c(-178.65, -156.5), ylim=c(50, 65.9)) +
  labs(x="", y="") +
  plot.theme() +
  theme(legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()) 
ggsave(filename=paste("Plots/ROMS_BT_krig_", i, ".png", sep=""), plot=BT_ROMS_krig.yr, dpi=500, width=6.5, height=9, units="in") }

# Interannual index of CPE (sq. km):  
ROMS.CPE = CPE_surveyRep %>%
  group_by(Year) %>%
  summarise(CPE = sum(ColdPool)) 

# Add estimated year-specific CPE to BT data (for model fitting):
ROMS.BT.krig_yr$Year = as.integer(ROMS.BT.krig_yr$Year)
ROMS.BT.CPE = ROMS.BT.krig_yr %>% 
  left_join(ROMS.CPE)
  
# Quantify the smoothed relationship between year and CPE (for model predictions):
ROMS.CPE.gam = gam(CPE ~ s(Year), data=ROMS.CPE, 
                   family=gaussian, method = "GCV.Cp"); summary(ROMS.CPE.gam)

CPE.gam_data = as.data.frame(ROMS.CPE$Year)
  colnames(CPE.gam_data) = c("Year")
ROMS.CPE_pred = predict.gam(ROMS.CPE.gam, newdata = CPE.gam_data)
  ROMS.CPE_pred = as.data.frame(ROMS.CPE_pred)
ROMS.CPE_pred$Year = CPE.gam_data$Year
  colnames(ROMS.CPE_pred) = c("CPE", "Year")
ROMS.CPE_pred$CPE = as.numeric(ROMS.CPE_pred$CPE)

# Plot CPE index and smoothed predictions (Fig. S2):
CPE = ggplot() +
  geom_line(data=ROMS.CPE, aes(x=Year, y=CPE), lwd=1, color="gray15") +
  geom_line(data=ROMS.CPE_pred, aes(x=Year, y=CPE), lwd=1, color="blue1") +
  plot.theme() +
  theme(plot.margin = unit(c(0.25,1.0,0.25,0.25),"cm")) +
  labs(x="Survey Year", y=expression(Cold~Pool~Extent~(km^2)), size=12) +
  scale_x_continuous(expand=c(0,0.1), breaks=seq(1982,2019, by=4))
ggsave(CPE, filename="Plots/ROMS_CPE.png", dpi=500, width=6.3, height=2.8, units="in", bg="transparent")

# Add smoothed CPE predictions to interpolated year-specific BT:
ROMS.BT.CPE = ROMS.BT.CPE %>% 
  left_join(ROMS.CPE_pred)

  save(ROMS.BT.CPE, file = "Data/ROMS_BT_CPE.rda")
# load("Data/ROMS_BT_CPE.rda")

# Link interpolated long-term mean BT to static covariate rasters:
  # Note: No index of cold pool extent for static models. 
ROMS.BT.CPE_sf = st_as_sf(ROMS.BT.CPE, coords = c("lon", "lat"))
  st_crs(ROMS.BT.CPE_sf) = 4326 # WGS84

dynamic.cov_sf = st_join(ROMS.BT.CPE_sf, cov.rasters_sf,
                               join = st_nearest_feature, left = T)

# Join covariate and survey data:
  # Separate into 5-yr chunks to stay within computer programming limits.

# 1982 to 1986:
haul_82_86 = subset(haul_sf, Year > 1981 & Year < 1987)
dyn.cov_82_86 = subset(dynamic.cov_sf, Year > 1981 & Year < 1987)

link_82_86 = list()
for(i in unique(haul_82_86$Year)) {
  HAUL_82_86 = subset(haul_82_86, Year == i)
  COV_82_86 = subset(dyn.cov_82_86, Year == i)
link_82_86[[i]] = st_join(HAUL_82_86, COV_82_86, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_82_86 = as.data.frame(dplyr::bind_rows(link_82_86))

  save(haul.dyn.cov_82_86, file = "Data/haul_dyn_cov_82_86.rda")
# load("Data/haul_dyn_cov_82_86.rda")

# 1987 to 1991:
haul_87_91 = subset(haul_sf, Year > 1986 & Year < 1992)
dyn.cov_87_91 = subset(dynamic.cov_sf, Year > 1986 & Year < 1992)

link_87_91 = list()
for(i in unique(haul_87_91$Year)) {
  HAUL_87_91 = subset(haul_87_91, Year == i)
  COV_87_91 = subset(dyn.cov_87_91, Year == i)
link_87_91[[i]] = st_join(HAUL_87_91, COV_87_91, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_87_91 = as.data.frame(dplyr::bind_rows(link_87_91))

  save(haul.dyn.cov_87_91, file = "Data/haul_dyn_cov_87_91.rda")
# load("Data/haul_dyn_cov_87_91.rda")

# 1992 to 1996:
haul_92_96 = subset(haul_sf, Year > 1991 & Year < 1997)
dyn.cov_92_96 = subset(dynamic.cov_sf, Year > 1991 & Year < 1997)

link_92_96 = list()
for(i in unique(haul_92_96$Year)) {
  HAUL_92_96 = subset(haul_92_96, Year == i)
  COV_92_96 = subset(dyn.cov_92_96, Year == i)
link_92_96[[i]] = st_join(HAUL_92_96, COV_92_96, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_92_96 = as.data.frame(dplyr::bind_rows(link_92_96))

  save(haul.dyn.cov_92_96, file = "Data/haul_dyn_cov_92_96.rda")
# load("Data/haul_dyn_cov_92_96.rda")

# 1997 to 2001:
haul_97_01 = subset(haul_sf, Year > 1996 & Year < 2002)
dyn.cov_97_01 = subset(dynamic.cov_sf, Year > 1996 & Year < 2002)

link_97_01 = list()
for(i in unique(haul_97_01$Year)) {
  HAUL_97_01 = subset(haul_97_01, Year == i)
  COV_97_01 = subset(dyn.cov_97_01, Year == i)
link_97_01[[i]] = st_join(HAUL_97_01, COV_97_01, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_97_01 = as.data.frame(dplyr::bind_rows(link_97_01))

  save(haul.dyn.cov_97_01, file = "Data/haul_dyn_cov_97_01.rda")
# load("Data/haul_dyn_cov_97_01.rda")

# 2002 to 2006:
haul_02_06 = subset(haul_sf, Year > 2001 & Year < 2007)
dyn.cov_02_06 = subset(dynamic.cov_sf, Year > 2001 & Year < 2007)

link_02_06 = list()
for(i in unique(haul_02_06$Year)) {
  HAUL_02_06 = subset(haul_02_06, Year == i)
  COV_02_06 = subset(dyn.cov_02_06, Year == i)
link_02_06[[i]] = st_join(HAUL_02_06, COV_02_06, 
                          join = st_nearest_feature, left = T)} 
haul.dyn.cov_02_06 = as.data.frame(dplyr::bind_rows(link_02_06))

  save(haul.dyn.cov_02_06, file = "Data/haul_dyn_cov_02_06.rda")
# load("Data/haul_dyn_cov_02_06.rda")

# 2007 to 2011:
haul_07_11 = subset(haul_sf, Year > 2006 & Year < 2012)
dyn.cov_07_11 = subset(dynamic.cov_sf, Year > 2006 & Year < 2012)

link_07_11 = list()
for(i in unique(haul_07_11$Year)) {
  HAUL_07_11 = subset(haul_07_11, Year == i)
  COV_07_11 = subset(dyn.cov_07_11, Year == i)
link_07_11[[i]] = st_join(HAUL_07_11,COV_07_11, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_07_11 = as.data.frame(dplyr::bind_rows(link_07_11))

  save(haul.dyn.cov_07_11, file = "Data/haul_dyn_cov_07_11.rda")
# load("Data/haul_dyn_cov_07_11.rda")

# 2012 to 2016:
haul_12_16 = subset(haul_sf, Year > 2011 & Year < 2017)
dyn.cov_12_16 = subset(dynamic.cov_sf, Year > 2011 & Year < 2017)

link_12_16 = list()
for(i in unique(haul_12_16$Year)) {
  HAUL_12_16 = subset(haul_12_16, Year == i)
  COV_12_16 = subset(dyn.cov_12_16, Year == i)
link_12_16[[i]] = st_join(HAUL_12_16, COV_12_16, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_12_16 = dplyr::bind_rows(link_12_16)

  save(haul.dyn.cov_12_16, file = "Data/haul_dyn_cov_12_16.rda")
# load("Data/haul_dyn_cov_12_16.rda")

# 2017 and 2018:
haul_17_18 = subset(haul_sf, Year > 2016 & Year < 2019)
dyn.cov_17_18 = subset(dynamic.cov_sf, Year > 2016 & Year < 2019)

link_17_18 = list() 
for(i in unique(haul_17_18$Year)) {
  HAUL_17_18 = subset(haul_17_18, Year == i)
  COV_17_18 = subset(dyn.cov_17_18, Year == i)
link_17_18[[i]] = st_join(HAUL_17_18, COV_17_18, 
                          join = st_nearest_feature, left = T)}
haul.dyn.cov_17_18 = as.data.frame(dplyr::bind_rows(link_17_18))

  save(haul.dyn.cov_17_18, file = "Data/haul_dyn_cov_17_18.rda")
# load("Data/haul_dyn_cov_17_18.rda")

# Combine 5-yr chunks into a single data frame:
haul.dyn.cov_82_18 = dplyr::bind_rows(
    haul.dyn.cov_82_86, 
    haul.dyn.cov_87_91, 
    haul.dyn.cov_92_96, 
    haul.dyn.cov_97_01, 
    haul.dyn.cov_02_06, 
    haul.dyn.cov_07_11, 
    haul.dyn.cov_12_16, 
    haul.dyn.cov_17_18)

coords = st_coordinates(haul.dyn.cov_82_18$geometry)
dynamic.model.data = cbind(haul.dyn.cov_82_18, coords)
dynamic.model.data = dynamic.model.data %>% 
  rename(Year = Year.x, lon = X, lat = Y)

  save(dynamic.model.data, file = "Data/dynamic_model_data.rda")
# load("Data/dynamic_model_data.rda")
  
#################################################################
# Estimate mean station-specific area swept (for static model predictions):
s.effort_mean = trawl %>%
  summarise(AreaSwept = mean(AreaSwept, na.rm = T)) 
s.effort_sd = trawl %>%
  summarise(AreaSwept = sd(AreaSwept, na.rm = T)) 

# Add to covariate data:
static.cov_sf$AreaSwept = as.numeric(s.effort_mean)

static.model.grid = as.data.frame(as_Spatial(static.cov_sf))
static.model.grid = static.model.grid %>% 
  rename(lon = coords.x1, lat = coords.x2)
  save(static.model.grid, file = "Data/static_model_grid.rda")
# load("Data/static_model_grid.rda")    
  
# Estimate mean area swept by station and year (for dynamic model predictions; Table S1):
d.effort_mean = trawl %>%
  group_by(Year) %>%
  summarise(AreaSwept = mean(AreaSwept, na.rm = T)) 
d.effort_sd = trawl %>%
  group_by(Year) %>%
  summarise(AreaSwept = mean(AreaSwept, na.rm = T))

# Add to covariate data:
dynamic.model.grid = dynamic.cov_sf %>% left_join(d.effort_mean)
dynamic.model.grid = as.data.frame(as_Spatial(dynamic.model.grid))
dynamic.model.grid = dynamic.model.grid %>% 
  rename(lon = coords.x1, lat = coords.x2)
  save(dynamic.model.grid, file = "Data/dynamic_model_grid.rda")
# load("Data/dynamic_model_grid.rda")  

##################################################################
# Illustrate covariate relationships (Fig. S3):
    # Assumes linear relationships, for illustration purposes only.
dyn.covariates = dynamic.model.data %>%
  rename(Lon = lon, Lat =lat, BT = ROMS_BT, Depth = bathy, Phi = phi, Slope = slope, SFI.sponges = sponges, SFI.corals = corals, SFI.whips = whips) %>%
  mutate(Location = Lon * Lat) %>%
  dplyr::select(c(Location, Year, Depth, Phi, Slope, BPI, SFI.sponges, SFI.corals, SFI.whips, BT, CPE))

CorrMat = rcorr(as.matrix(dyn.covariates), type="spearman")
corrplot(CorrMat$r, type="upper", order="original", diag=F, 
                    tl.pos="td", tl.cex=1, cl.cex=1, tl.col="black", 
                    p.mat = CorrMat$P, sig.level = 0.1, insig = "blank", 
                    method="number") # Plots/pairs_plot.png

##################################################################
# Prepare covariates for retrospective skill testing:

# Static covariates:
# Format prediction grid as a series of rasters (to decrease resolution from 1 sq.km to 10 sq.km [to accoutn for processing limitations], retrospective skill testing only):
s.forecast.data = static.model.grid
coordinates(s.forecast.data) = ~ lon + lat
  gridded(s.forecast.data) = T

# Convert each variable to a raster layer:
crs(s.forecast.data) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    s.forecasts.bathy_1km = raster(s.forecast.data, layer="bathy")
    s.forecasts.phi_1km = raster(s.forecast.data, layer="phi")
    s.forecasts.slope_1km = raster(s.forecast.data, layer="slope")
    s.forecasts.BPI_1km = raster(s.forecast.data, layer="BPI")
    s.forecasts.BT_1km = raster(s.forecast.data, layer="ROMS_BT")
    s.forecasts.sponges_1km = raster(s.forecast.data, layer="sponges")
    s.forecasts.corals_1km = raster(s.forecast.data, layer="corals")
    s.forecasts.whips_1km = raster(s.forecast.data, layer="whips")
    
# join static covariate rasters:
s.forecasts.stack_1km = stack(s.forecasts.bathy_1km,
                            s.forecasts.phi_1km,
                            s.forecasts.slope_1km, 
                            s.forecasts.BPI_1km, 
                            s.forecasts.BT_1km,
                            s.forecasts.sponges_1km, 
                            s.forecasts.corals_1km, 
                            s.forecasts.whips_1km)
  save(s.forecasts.stack_1km, file="Data/s_forecasts_stack_1km.rda")
# load("Data/s_forecasts_stack_1km.rda")    

res(s.forecasts.stack_1km)
  nrow(s.forecasts.stack_1km) # original resolution
s.forecasts.stack_10km = aggregate(s.forecasts.stack_1km, fact=10, FUN=mean)
  res(s.forecasts.stack_10km); nrow(s.forecasts.stack_10km) # new resolution
    
# Reformat as data frame and add area swept back in:
s.forecasts.10km = as.data.frame(rasterToPoints(s.forecasts.stack_10km))
s.forecasts.10km = s.forecasts.10km %>%
  rename(lon = x, lat = y) %>%
  mutate(sponges = as.factor(ifelse(sponges > 1, 1, 0))) %>%
  mutate(corals = as.factor(ifelse(corals > 1, 1, 0))) %>%
  mutate(whips = as.factor(ifelse(whips > 1, 1, 0)))
s.forecasts.10km$AreaSwept = s.effort_mean$AreaSwept

  save(s.forecasts.10km, file="Data/s_forecasts_stack_10km.rda")
# load("Data/s_forecasts_stack_10km.rda")

  
# Dynamic covariates:
# Format prediction grid as a series of rasters (to decrease resolution from 1 sq.km to 10 sq.km [to accoutn for processing limitations], retrospective skill testing only):
d.forecast.data = dynamic.model.grid
coordinates(d.forecast.data) = ~ lon + lat
  gridded(d.forecast.data) = T

# Convert each variable to a raster layer:
  # year-specific bottom temperature:
for(i in unique(d.forecast.data$Year)) {
dyn.BT = subset(d.forecast.data, Year == i)
  dyn.BT_1km = raster(dyn.BT, layer="ROMS_BT") 
assign(paste("dyn.BT_1km_", i, sep=""), dyn.BT_1km) }
forecasts.BT_1km = stack(mget(ls(pattern="dyn.BT_1km_")))
  crs(forecasts.BT_1km) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  save(forecasts.BT_1km, file = "Data/dyn_forecasts_BT_1km.rda")
# load("Data/dyn_forecasts_BT_1km.rda")

# static covariates:
crs(d.forecast.data) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    d.forecasts.bathy_1km = raster(d.forecast.data, layer="bathy")
    d.forecasts.phi_1km = raster(d.forecast.data, layer="phi")
    d.forecasts.slope_1km = raster(d.forecast.data, layer="slope")
    d.forecasts.BPI_1km = raster(d.forecast.data, layer="BPI")
    d.forecasts.sponges_1km = raster(d.forecast.data, layer="sponges")
    d.forecasts.corals_1km = raster(d.forecast.data, layer="corals")
    d.forecasts.whips_1km = raster(d.forecast.data, layer="whips")
    
# join dynamic and static covariate rasters:
d.forecasts.stack_1km = stack(d.forecasts.BT_1km,
                            d.forecasts.bathy_1km,
                            d.forecasts.phi_1km,
                            d.forecasts.slope_1km, 
                            d.forecasts.BPI_1km, 
                            d.forecasts.sponges_1km, 
                            d.forecasts.corals_1km, 
                            d.forecasts.whips_1km)
  save(d.forecasts.stack_1km, file="Data/dyn_forecasts_stack_1km.rda")
# load("Data/dyn_forecasts_stack_1km.rda")    

res(d.forecasts.stack_1km)
  nrow(d.forecasts.stack_1km) # original resolution
d.forecasts.stack_10km = aggregate(d.forecasts.stack_1km, fact=10, FUN=mean)
  res(d.forecasts.stack_10km); nrow(d.forecasts.stack_10km) # new resolution
    
# Reformat as data frame and add year back in:
d.forecasts.df_10km = as.data.frame(rasterToPoints(d.forecasts.stack_10km))
forecasts_10km = gather(d.forecasts.df_10km, key = "Year", 
            value = "ROMS_BT", starts_with("dyn.BT_1km_"))
    
d.forecasts_10km$Year = d.forecasts_10km$Year %>% 
  str_replace("dyn.BT_1km_", "")
d.forecasts_10km$Year = as.numeric(d.forecasts_10km$Year)
    
# Add year-specific CPE and area swept back in:
d.forecasts.10km = d.forecasts_10km %>% 
   left_join(ROMS.CPE) %>%
   left_join(., d.effort) %>%
   rename(lon = x, lat = y) %>%
   mutate(sponges = as.factor(ifelse(sponges > 1, 1, 0))) %>%
   mutate(corals = as.factor(ifelse(corals > 1, 1, 0))) %>%
   mutate(whips = as.factor(ifelse(whips > 1, 1, 0)))
  save(d.forecasts.10km, file="Data/dyn_forecasts_stack_10km.rda")
# load("Data/dyn_forecasts_stack_10km.rda")

#################################################################  
# Calculate min, mean, max, and sd for BT (Table S4):

# Presence-absence:
dynamic.model.data %>%
    summarise(min(ROMS_BT))
dynamic.model.data %>%
    summarise(mean(ROMS_BT))
dynamic.model.data %>%
    summarise(sd(ROMS_BT))
dynamic.model.data %>%
    summarise(max(ROMS_BT))

# Positive Catch, ATF:
dynamic.model.data %>%
    filter(Species == "Arrowtooth.Flounder" & adjN > 0) %>%
    summarise(min(ROMS_BT))
dynamic.model.data %>%
    filter(Species == "Arrowtooth.Flounder" & adjN > 0) %>%
    summarise(mean(ROMS_BT))
dynamic.model.data %>%
    filter(Species == "Arrowtooth.Flounder" & adjN > 0) %>%
    summarise(sd(ROMS_BT))
dynamic.model.data %>%
    filter(Species == "Arrowtooth.Flounder" & adjN > 0) %>%
    summarise(max(ROMS_BT))
    
# Positive Catch, ATF:
dynamic.model.data %>%
    filter(Species == "Walleye.Pollock" & adjN > 0) %>%
    summarise(min(ROMS_BT))
dynamic.model.data %>%
    filter(Species == "Walleye.Pollock" & adjN > 0) %>%
    summarise(mean(ROMS_BT))
dynamic.model.data %>%
    filter(Species == "Walleye.Pollock" & adjN > 0) %>%
    summarise(sd(ROMS_BT))
dynamic.model.data %>%
    filter(Species == "Walleye.Pollock" & adjN > 0) %>%
    summarise(max(ROMS_BT))

# Static model predictions:
static.model.grid %>%
    summarise(min(ROMS_BT))
static.model.grid %>%
    summarise(mean(ROMS_BT))
static.model.grid %>%
    summarise(sd(ROMS_BT))
static.model.grid %>%
    summarise(max(ROMS_BT))

# Dynamic model predictions:
dynamic.model.grid %>%
    summarise(min(ROMS_BT))
dynamic.model.grid %>%
    summarise(mean(ROMS_BT))
dynamic.model.grid %>%
    summarise(sd(ROMS_BT))
dynamic.model.grid %>%
    summarise(max(ROMS_BT))
