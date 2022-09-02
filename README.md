## overview:
This repo contains all the code used to produce the analysis in
"Global phenology maps reveal patterns of regional divergence,
intercontinental convergence, and climate-independent tropical asynchrony"
(Terasaki Hart et al. 2022), chapter 3 of Drew E. Terasaki Hart's PhD dissertation.

# < PUT CITATION HERE >

Code was written by Drew Ellison Terasaki Hart,
Thao-Nguyen Bui, and Lauren Di Maggio.
It is made freely available under the MIT License,
to be distributed and/or modified with proper attribution.


Any questions, concerns, or other interest can be directed
to drew <dot> hart <at> berkeley <dot> edu. 


### contents

Analyses were run in two major stages, so the two major directories of content are organized to reflect this:
  1. `/phen/`: calculation of global maps of characteristic seasonal phenology, and associated analysis
  2. `/asynch/`: calculation of global maps of phenological asynchrony, and associated analysis
Other directories include:
  - `/data_prep`: code used to download the SIF dataset used, convert to Geotiff, and upload to Google Earth Enginge (GEE)
  - `/data`: ancillary data (i.e., not our main datasets) that lives in this repo and is used in analysis
  - `/etc`: other odds and ends

***NOTE:*** **All GEE Javascript code must be run in GEE. Other code was designed either to be run on UC Berkeley's Savio compute cluster or on a local machine.**
