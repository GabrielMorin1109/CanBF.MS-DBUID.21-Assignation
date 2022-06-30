# Assign Dissemination Block Unique Identifier (DBUID) to each Building Footprints (BFs) produced by Microsoft either by direct intersections or by nearest features
In order to speed up spatial operation on the BFs, I have assigned to each BF polygon either the DBUID that has an intersection with the polygon or the nearest DBUID feature. DBUID is a geographical information code and is described as follows :
> "Each dissemination block is assigned a three-digit code. In order to uniquely identify each dissemination block in Canada, the two-digit province/territory (PR) code, the two-digit census division (CD) code and the four-digit dissemination area (DA) code must precede the dissemination block (DB) code. " (see : https://www12.statcan.gc.ca/census-recensement/2021/ref/dict/az/Definition-eng.cfm?ID=geo014)

In order to run the procedure on your machine, you must run the scripts in the following order : 
1. [download_database_to_run_Microsoft_building_footprint.R](https://github.com/GabrielMorin1109/CanBF.MS-DBUID.21-Assignation/blob/main/R/download_database_to_run_Microsoft_building_footprint.R)
2. [Microsoft_building_footprint.R](https://github.com/GabrielMorin1109/CanBF.MS-DBUID.21-Assignation/blob/main/R/Microsoft_building_footprint.R)

