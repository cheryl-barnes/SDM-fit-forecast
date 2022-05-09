## Climate-informed models benefit hindcasting but present challenges when forecasting species-habitat associations

# Citation
Barnes CL, Essington TE, Pirtle JP, Rooper CN, Laman EA, Holsman KK, Aydin KY, and Thorson JT. In revision. Climate-informed models benefit hindcasting but present challenges when forecasting species-habitat associations. Ecogr. Forthcoming. <br><br>

<b> Cheryl L. Barnes </b><br>
School of Aquatic and Fishery Sciences, University of Washington <br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br>
Alaska Regional Office, National Marine Fisheries Service (NOAA) <br><br>

<b> Timothy E. Essington </b><br>
School of Aquatic and Fishery Sciences, University of Washington <br><br>

<b> Jodi L. Pirtle </b><br>
Alaska Regional Office, National Marine Fisheries Service (NOAA) <br><br>

<b> Christopher N. Rooper </b><br>
Fisheries and Oceans Canada <br><br>

<b> Edward A. Laman </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br><br>

<b> Kirstin K. Holsman </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br><br>

<b> Kerim Y. Aydin </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br><br>

<b> James T. Thorson </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br><br>

## Overview
This repository details the methods used to:
a) Identify best-fit generalized additive models (GAMs) for hindcasting (based on conventional statistics) and evaluate forecast skill (based on retrospective skill testing; Thorson 2019) for presence-absence, numerical abundance, and biomass of groundfishes in the Bering Sea (1982-2018); and
b) Compare species distribution models (SDMs) with varying degrees of complexity to assess whether the addition of time-varying processes to status quo static SDMs improves hindcast performance and/or forecast skill.

## Case Study:
Arrowtooth Flounder (<i>Atheresthes stomias</i>) and Walleye Pollock (<i>Gadus chalcogrammus</i>) in the Bering Sea (1990 to 2018). 

# Model Structure:
Static models rely on spatially-explicit but time-invariant environmental conditions whereas dynamic refer to those that account for spatial and temporal variation in select model covariates. Covariates included in each model type: S (static) – spatial variation, long-term mean bottom temperature, and static habitat covariates; D1 (simple dynamic) – spatial variation, location- and Year-specific bottom temperature, interannual index of cold pool extent, and static habitat covariates; D2 (intermediate dynamic) – D1 covariates plus temporal variation; D3 (complex dynamic) – D2 covariates plus spatiotemporal variation. 

## Financial and Logistical Support
This project was funded through the Joint Institute for the Study of the Atmosphere and Ocean (JISAO) under NOAA Cooperative Agreement NA15OAR4320063 and through the Cooperative Institute for Climate, Ocean, & Ecosystem Studies (CIOCES) under NOAA Cooperative Agreement NA20OAR4320271: Contribution No. 2021-1170. Financial support originated from NOAA’s internal request for Magnuson-Stevens Act Implementation proposals. Facilities, equipment, and in-kind support were provided by the University of Washington, the Alaska Fisheries Science Center (NOAA), and Alaska Regional Office (NOAA).

## Acknowledgments
The AFSC’s Resource Assessment and Conservation Engineering Division collected all bottom trawl survey data used in this study. Static model covariates were provided by the Alaska Regional Office’s Habitat Conservation Division. The Alaska Climate Integrated Modeling (ACLIM) project provided survey-replicated hindcasts of bottom temperature and cold pool extent. We thank Lyle Britt and Stan Kotwicki for assisting with data acquisition and interpretation. This manuscript was improved upon by comments from Margaret Siple, Andrew Allyn, and three anonymous reviewers. The findings and conclusions in this paper are solely those of the authors and do not necessarily represent the views of any affiliation previously listed.  <br>

## References
Hermann AJ, GA Gibson, NA Bond, EN Curchitser, K Hedstrom, W Cheng, M Wang, PJ Stabeno, L Eisner, and KD Cieciel. 2013. A multivariate analysis of observed and modeled biophysical variability on the Bering Sea shelf: multidecadal hindcasts (1970-2009) and forecasts (2010-2040). Deep-Sea Res II: Top Stud Oceanogr. 94:121–139.
Hoff GR. 2016. Results of the 2016 eastern Bering Sea upper continental slope survey of groundfish and invertebrate resources. US Dep Comm NOAA Tech Memo. NMFS AFSC-339. 272 pp.
Kearney K, A Hermann, W Cheng, I Ortiz, and K Aydin. 2020. A coupled pelagic-benthic-sympagic biogeochemical model for the Bering Sea: documentation and validation of the BESTNPZ model (v2019.08.23) within a high-resolution regional ocean model. Geosci Model Dev. 13:597-650.
Laman EA, CN Rooper, K Turner, S Rooney, DW Cooper, and M Zimmerman. 2018. Using species distribution models to describe essential fish habitat in Alaska. Can J Fish Aquat Sci. 75:1230–1255. 
Lauth R and E Acuna. 2007. 2005 bottom trawl survey of the eastern Bering Sea continental shelf. AFSC Processed Rep. US Dept Comm. 2007-1. 164 pp. 
Lauth RR, EJ Dawson, and J Connor. 2019. Results of the 2017 eastern and northern Bering Sea continental shelf bottom trawl survey of groundfish and invertebrate fauna. NOAA Tech Mem 396. 270 pp.
Pirtle JL, SK Shotwell, M Zimmermann, JA Reid, and N Golden. 2019. Habitat suitability models for groundfish in the Gulf of Alaska. Deep-Sea Res. Pt. II. 165: 303–321
Rooper CN, MF Sigler, P Goddard, P Malecha, R Towler, K Williams, R Wilborn, and M Zimmermann. 2016. Validation and improvement of species distribution models for structure-forming invertebrates in the eastern Bering Sea with an independent survey. Mar Ecol Prog Ser. 551:117-130.
Sigler MF, CN Rooper, GR Hoff, RP Stone, RA McConnaughey, and TK Wilderbuer. 2015. Faunal features of submarine canyons on the eastern Bering Sea slope. Mar Ecol Prog Ser. 526:21-40.
Stahl JP and GH Kruse. 2008. Spatial and temporal variability in size at maturity of Walleye Pollock in the Eastern Bering Sea. Trans Amer Fish Soc. 137:1543–1557.
Stevenson DE and RR Lauth. 2019. Bottom trawl surveys in the northern Bering Sea indicate recent shifts in the distribution of marine species. Polar Biol. 42:407–421.
Thorson JT. 2019. Forecast skill for predicting distribution shifts: a retrospective experiment for marine fishes in the Eastern Bering Sea. Fish Fish. 20:159–173.
Zimmermann M, MM Prescott, and CN Rooper. 2013. Smooth Sheet Bathymetry of the Aleutian Islands. US Dep Comm NOAA Tech Memo. NMFS-AFSC-250. 43pp.
Zimmermann M, and M Prescott. 2018. Bathymetry and Canyons of the Eastern Bering Sea Slope. Geosciences. 8:184.
