# Author : Gabriel Morin
# Date : 30 Juin 2022
# Script that download all the database needed to run the Microsoft Building Footprints script.

# LIST OF THE NEEDED DATABASES:
# -> Microsoft Building Footprints
#   ~ Link: (https://github.com/microsoft/CanadianBuildingFootprints)
# -> 2021 Census – Digital Boundary Files (DBF) with Statistical boundaries selected at the Dissemination Blocks level shapefile
#   ~ Link : (https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?Year=21)
#   ~ Direct Download Link : (https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/ldb_000a21a_e.zip)

# Other : 
# -> PRNAME_dic  [Dictionnary of the Province Name]
# ==============================================================================
# Directory of where the data will be save
HDD <- "/media/gabriel/HDD/Data_memoire/"
setwd(HDD)

# File structure
{
  dir.create("Microsoft")
  dir.create("Microsoft/Data_Microsoft")
  dir.create("Microsoft/Microsoft_DBUID/")
  dir.create("Microsoft/Microsoft_DBUID_nearest")
  dir.create("Microsoft/shp_Microsoft_DBUID_all")
  
  dir.create("dissemination-blocks")
  dir.create("dissemination-blocks/Digital_Boundary_Files")
}
# ==============================================================================
# increase download timeout option in R (base of 60 seconds which is insufficient)
options(timeout = max(300, getOption("timeout")))


# MICROSOFT :
MS_URL <-
  list(
    url   = "https://usbuildingdata.blob.core.windows.net/canadian-buildings-v2/",
    files = 
      c(
        "Alberta.zip",
        "BritishColumbia.zip",
        "Manitoba.zip",
        "NewBrunswick.zip",
        "NewfoundlandAndLabrador.zip",
        "NorthwestTerritories.zip",
        "NovaScotia.zip",
        "Nunavut.zip",
        "Ontario.zip",
        "PrinceEdwardIsland.zip",
        "Quebec.zip",
        "Saskatchewan.zip",
        "YukonTerritory.zip"
      )
  )

# Go to Data_Microsoft directory
setwd("Microsoft_test/Data_Microsoft/")
# Loop that download all Microsoft Building Footprints files
for(file in MS_URL$files){
  download.file(
    url      = paste0(MS_URL$url, file), 
    destfile = file
  )
}

# unzip MS file 
sapply(MS_URL$files, utils::unzip, overwrite = TRUE)

# Remove zip files
unlink(MS_URL$files)
# Return to the original working directory
setwd("../../")


# DISSEMINATION BLOCKS :
# Go to the dissemination block folder from the Data_Microsoft folder
setwd("dissemination-blocks/Digital_Boundary_Files")

DB_URL = list(
  url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/",
  file = "ldb_000a21a_e.zip"
)
# Download the shapefile 
download.file(
  url      = paste0(DB_URL$url, DB_URL$file), 
  destfile = DB_URL$file
)

# Unzip file 
utils::unzip(DB_URL$file, overwrite = TRUE)

# Remove zip files
unlink(DB_URL$file)

# Return to dissemination-blocks folder
setwd("../")

# PROVINCE NAME DICTIONARY
download.file(
  url  = "https://raw.githubusercontent.com/GabrielMorin1109/CanBF.MS-DBUID.21-Assignation/main/PRNAME_dic.csv",
  destfile = "PRNAME_dic.csv"
)

# ==============================================================================
# Microsoft Building Footprints should be download and unzip for each province 
# 2021 Census – Digital Boundary Files (DBF) should be download and unzip
# PRNAME_dic.csv should be download
# ==============================================================================
# The folder structure should be as follow :
# /media/gabriel/HDD/Data_memoire/
# ├── [4.0K]  dissemination-blocks
# │   ├── [4.0K]  Digital_Boundary_Files
# │   │   ├── [   9]  ldb_000a21a_e.cpg
# │   │   ├── [ 48M]  ldb_000a21a_e.dbf
# │   │   ├── [ 524]  ldb_000a21a_e.prj
# │   │   ├── [4.1M]  ldb_000a21a_e.sbn
# │   │   ├── [ 55K]  ldb_000a21a_e.sbx
# │   │   ├── [448M]  ldb_000a21a_e.shp
# │   │   ├── [3.8M]  ldb_000a21a_e.shx
# │   │   └── [ 47K]  ldb_000a21a_e.xml
# │   └── [ 615]  PRNAME_dic.csv
# └── [4.0K]  Microsoft
#     ├── [4.0K]  Data_Microsoft
#     │   ├── [390M]  Alberta.geojson
#     │   ├── [302M]  BritishColumbia.geojson
#     │   ├── [135M]  Manitoba.geojson
#     │   ├── [ 71M]  NewBrunswick.geojson
#     │   ├── [ 51M]  NewfoundlandAndLabrador.geojson
#     │   ├── [3.0M]  NorthwestTerritories.geojson
#     │   ├── [ 81M]  NovaScotia.geojson
#     │   ├── [708K]  Nunavut.geojson
#     │   ├── [809M]  Ontario.geojson
#     │   ├── [ 16M]  PrinceEdwardIsland.geojson
#     │   ├── [512M]  Quebec.geojson
#     │   ├── [146M]  Saskatchewan.geojson
#     │   └── [2.6M]  YukonTerritory.geojson
#     ├── [4.0K]  Microsoft_DBUID
#     ├── [4.0K]  Microsoft_DBUID_nearest
#     └── [4.0K]  shp_Microsoft_DBUID_all
# 
# 7 directories, 22 files
# ==============================================================================
# END SCRIPT
# [May not work on your system] (work in Ubuntu 20.04.4 LTS. You need to have the tree package install.)
if(.Platform$OS.type != "unix"){
  # Print directory tree structure (package info : https://packages.ubuntu.com/jammy/tree)
  system(paste("tree -h", HDD, sep = " "))
}