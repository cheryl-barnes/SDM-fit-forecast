## Climate-informed models benefit hindcasting but present challenges when forecasting species-habitat associations

#### Citation: Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. Accepted. Climate-informed models benefit hindcasting but present challenges when forecasting species-habitat associations. Ecogr. Forthcoming. <br><br>

<b> Cheryl L. Barnes </b><br>
School of Aquatic and Fishery Sciences, University of Washington <br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br>
Alaska Regional Office, National Marine Fisheries Service (NOAA) <br>

<b> Timothy E. Essington </b><br>
School of Aquatic and Fishery Sciences, University of Washington <br>

<b> Jodi L. Pirtle </b><br>
Alaska Regional Office, National Marine Fisheries Service (NOAA) <br>

<b> Christopher N. Rooper </b><br>
Fisheries and Oceans Canada <br>

<b> Edward A. Laman </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br>

<b> Kirstin K. Holsman </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br>

<b> Kerim Y. Aydin </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br>

<b> James T. Thorson </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br><br>

## Overview
This repository details methods used to compare <i> status quo </i> "static" species distribution models (SDMs), which rely on spatially-explicit but time-invariant environmental covariates, with more climate-informed "dynamic" SDMs that account for both spatially-explicit and time-varying covariates. We assesed the distributions and densities of Walleye Pollock (<i>Gadus chalcogrammus</i>) and Arrowtooth Flounder (<i>Atheresthes stomias</i>) throughout the Bering Sea, from 1982 to 2018. We hypothesized that dynamic SDMs would outperform static SDMs for both hindcasting and forecasting species-habitat associations. <br>

## Data
Predictor variables included presence-absence and positive catches (represented as numerical abundance and biomass) from the Alaska Fisheries Science Center's bottom trawl survey. Environmental covariates were treated differently among four model specifications, as follows: 1) S (static) models accounted for spatial variation, long-term mean bottom temperature (degrees C), and static habitat covariates; 2) D1 (simple dynamic) models accounted for spatial variation, location- and year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; 3) D2 (intermediate dynamic) resembled D1 covariates with the addition of a year effect (<i>i.e.</i>, latent temporal variation; and 4) D3 (complex dynamic) models resemebled D2 covariates with the addition of an interaction of year and location (<i>i.e.</i>, latent spatiotemporal variation). Static habitat covariates included depth (m), slope (degrees), bathymetric position index (BPI; unitless), sediment grain size (phi), and the occurrence of structure-forming invertebrates (sponges, corals, and sea whips/pens). Bottom temperature and cold pool extent were derived from the BERING10K Regional Ocean Modeling System (ROMS). Each of the environmental covariates listed are known to impact the distributions and densities of groundfishes in the Bering Sea and elsewhere. We limited concurvity by removing one of two covariates with associations ≥ 0.5. Please refer to the publication for additional details and references. <br>

## Model Fitting
We specified presence-absence generalized additive models (GAMs) with a binomial distribution and complementary log-log link, numerical abundance GAMs with a Poisson distribution and log link, and biomass GAM with a Gamma distribution and log link. We used a tensor product smooth for location (longitude and latitude) and low rank isotropic smoothers for depth, slope, bathymetric position index, sediment grain size, and bottom temperature. We quantified effects of structure-forming invertebrates as factors and included log-transformed area swept as an offset. All smoothed covariates were estimated using thin plate regression splines and generalized cross-validation (Wood 2003, 2004). We limited substantial changes in predictions across similar values of a particular covariate by restricting the effective degrees of freedom to six (k=6) for depth and four (k=4) for all other variables. Effects of location, however, were unconstrained to allow for reasonably complex patterns in habitat use. We also limited the extrapolation of large probabilities outside the range of covariate data by specifying that smoothing penalties should affect the first derivative (setting m=1 in ‘mgcv’). In doing so, estimated responses would remain constant beyond observed values. <br>

## Assessment and Prediction
We used conventional statistics (R squared, % Deviance Explained, and UBRE/GCV) to identify best-fit generalized additive models (GAMs) for hindcasting distributions and densities of our focal species (1982 to 2018). We used retrospective skill testing to evaluate forecast skill (<i>sensu</i> Thorson 2019). Retrospective skill testing involved using the best-fit SDMs that we identified during hindcasting to fit a sequence of nested submodels. We required a minimum of 10 survey years to fit each nested submodel and iteratively increased the length of the time series until a single year remained for forecasting. We then estimates Spearman's correlation coefficient (rho) for each nested submodel to assess correlations between forecasts and observations. Finally, we used a 10-yr moving window to estimate mean (and variance) in Spearman's rho, a relative measure of forecast skill. <br>

## Financial and Logistical Support
This project was funded through the Joint Institute for the Study of the Atmosphere and Ocean (JISAO) under NOAA Cooperative Agreement NA15OAR4320063 and through the Cooperative Institute for Climate, Ocean, & Ecosystem Studies (CIOCES) under NOAA Cooperative Agreement NA20OAR4320271: Contribution No. 2021-1170. Financial support originated from NOAA’s internal request for Magnuson-Stevens Act Implementation proposals. Facilities, equipment, and in-kind support were provided by the University of Washington, the Alaska Fisheries Science Center (NOAA), and Alaska Regional Office (NOAA). <br>

## Acknowledgments
The AFSC’s Resource Assessment and Conservation Engineering Division collected all bottom trawl survey data used in this study. Static model covariates were provided by the Alaska Regional Office’s Habitat Conservation Division. The Alaska Climate Integrated Modeling (ACLIM) project provided survey-replicated hindcasts of bottom temperature and cold pool extent. We thank Lyle Britt and Stan Kotwicki for assisting with data acquisition and interpretation. This manuscript was improved upon by comments from Margaret Siple, Andrew Allyn, and three anonymous reviewers. The findings and conclusions in this paper are solely those of the authors and do not necessarily represent the views of any affiliation previously listed. <br>

## References
Hermann AJ, GA Gibson, NA Bond, EN Curchitser, K Hedstrom, W Cheng, M Wang, PJ Stabeno, L Eisner, and KD Cieciel. 2013. A multivariate analysis of observed and modeled biophysical variability on the Bering Sea shelf: multidecadal hindcasts (1970-2009) and forecasts (2010-2040). Deep-Sea Res II: Top Stud Oceanogr. 94:121–139.<br><br>
Hoff GR. 2016. Results of the 2016 eastern Bering Sea upper continental slope survey of groundfish and invertebrate resources. US Dep Comm NOAA Tech Memo. NMFS AFSC-339. 272 pp.<br><br>
Kearney K, A Hermann, W Cheng, I Ortiz, and K Aydin. 2020. A coupled pelagic-benthic-sympagic biogeochemical model for the Bering Sea: documentation and validation of the BESTNPZ model (v2019.08.23) within a high-resolution regional ocean model. Geosci Model Dev. 13:597-650.<br><br>
Laman EA, CN Rooper, K Turner, S Rooney, DW Cooper, and M Zimmerman. 2018. Using species distribution models to describe essential fish habitat in Alaska. Can J Fish Aquat Sci. 75:1230–1255. <br><br>
Lauth R and E Acuna. 2007. 2005 bottom trawl survey of the eastern Bering Sea continental shelf. AFSC Processed Rep. US Dept Comm. 2007-1. 164 pp. 
Lauth RR, EJ Dawson, and J Connor. 2019. Results of the 2017 eastern and northern Bering Sea continental shelf bottom trawl survey of groundfish and invertebrate fauna. NOAA Tech Mem 396. 270 pp. <br><br>
Pirtle JL, SK Shotwell, M Zimmermann, JA Reid, and N Golden. 2019. Habitat suitability models for groundfish in the Gulf of Alaska. Deep-Sea Res. Pt. II. 165: 303–321 <br><br>
Rooper CN, MF Sigler, P Goddard, P Malecha, R Towler, K Williams, R Wilborn, and M Zimmermann. 2016. Validation and improvement of species distribution models for structure-forming invertebrates in the eastern Bering Sea with an independent survey. Mar Ecol Prog Ser. 551:117-130. <br><br>
Sigler MF, CN Rooper, GR Hoff, RP Stone, RA McConnaughey, and TK Wilderbuer. 2015. Faunal features of submarine canyons on the eastern Bering Sea slope. Mar Ecol Prog Ser. 526:21-40. <br><br>
Stahl JP and GH Kruse. 2008. Spatial and temporal variability in size at maturity of Walleye Pollock in the Eastern Bering Sea. Trans Amer Fish Soc. 137:1543–1557. <br><br>
Stevenson DE and RR Lauth. 2019. Bottom trawl surveys in the northern Bering Sea indicate recent shifts in the distribution of marine species. Polar Biol. 42:407–421. <br><br>
Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.
Zimmermann M, MM Prescott, and CN Rooper. 2013. Smooth Sheet Bathymetry of the Aleutian Islands. US Dep Comm NOAA Tech Memo. NMFS-AFSC-250. 43pp. <br><br>
Zimmermann M, and M Prescott. 2018. Bathymetry and Canyons of the Eastern Bering Sea Slope. Geosciences. 8:184.
