-------------------
GENERAL INFORMATION
-------------------


1. Title of Dataset:  FIRED Western Hemisphere s1 t5 January 2017-June 2020

2. Authors: Jennifer K. Balch, Lise A. St. Denis, Adam L. Mahood, Nathan P.  Mietkiewicz, Travis Williams, Joe McGlinchy, Maxwell C. Cook.


3. Contact information: jennifer.balch@colorado.edu; adam.mahood@colorado.edu


4. Date of data collection: January 2017 - June 2020


--------------------------
SHARING/ACCESS INFORMATION
-------------------------- 


1. Licenses/restrictions placed on the data: MIT


2. Links to publications that cite or use the data: https://github.com/earthlab/warming-weakens-the-nighttime-barrier-to-global-fire


3. Links to other publicly accessible locations of the data: None


4. Recommended citation for the data: TBD


---------------------
DATA & FILE OVERVIEW
---------------------


1. File List:


   A. Filename: fired_wh_s1_t5_2017-.gpkg

--------------------------
METHODOLOGICAL INFORMATION
--------------------------

See Balch et al 2020. DOI: https://doi.org/10.3390/rs12213498

Code: http://www.github.com/firedpy

-----------------------------------------
DATA-SPECIFIC INFORMATION: 
-----------------------------------------

1. Number of variables: 24

2. Number of cases/rows: 281,615

3. Projection information (proj4 string): "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

3.1. The projection is the native projetion from the MODIS MCD64A1 burned area product from which this dataset is derived. The MCD64A1 product is a raster grid with a resultion of 463 meters. More info at https://lpdaac.usgs.gov/products/mcd64a1v006/

4. Variable List:
	A. Name: id
	   Description: Unique identifier of the fire event 

        B. Name: date
       	   Description: The date that the polygon burned

	C. Name: pixels
	   Description: Total number of pixels burned that day

	D. Name: l1_eco
	   Description: numeric code for the level 1 ecoregion