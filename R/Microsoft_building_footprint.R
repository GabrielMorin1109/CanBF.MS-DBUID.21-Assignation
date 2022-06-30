# Author : Gabriel Morin
# Date : 30 Juin 2022
#'* download_database_to_run_Microsoft_building_footprint.R had to be run before running this script *
#'* The code needs some cleaning, but should work. *

# Assignation of each Building Footprint with the corrected Dissemination Block
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)
library(tools)
library(sf)
sf::sf_use_s2(FALSE) # easy correction for st_centroid problem
library(geojsonsf)
library(magrittr)
library(nngeo) # For computing the nearest neighbour between two sf object
library(progress) # Implanting a progress bar to monitoring the progress of the loop
library(parallel)
library(pbapply)
library(dplyr)
# ==============================================================================
# VARIABLES THAT WILL CHANGE THE BEHAVIOR OF THE SCRIPT
# Do we want to go trough all forloop? (FALSE for debugging)
loop.Q <- FALSE
# set Working Directory
HDD <- "/media/gabriel/HDD/Data_memoire/"; setwd(HDD)
# ==============================================================================

if(FALSE){
  # Read Microsoft Building Footprint dataset ----
  # MS <- 
  #   "Microsoft/Data_Microsoft" %>%                             # select the corrected folder where the data is stored
  #   list.files(full.names = TRUE) %>%                           # list all shapefile into 
  #   lapply(geojsonsf::geojson_sf) %>%                           # read all shapefile into a list of sf object
  #   do.call(what =  dplyr::bind_rows, args = .) %>%             # merged the list of dataset into one sf object
  #   sf::st_make_valid() %>%                                     # Correct the geometry
  #   sf::st_transform(sf::st_crs(3347))                          # used of an CRS that give distance measure into meters (in Canada). CRS is the same as StatsCan Census 2022
  
  # Read DMTI building Footprint (not used) ----
  # DMTI.FP <- 
  #   paste0("DMTI_BuildingFootprints_2021/dmti_buildingfootprints_2020_s_poly") %>% # go to the correct folder
  #   list.files(full.names = TRUE) %>%                           # list all the files 
  #   {.[tools::file_ext(.) == "shp"]} %>%                        # select only the shapefiles
  #   lapply(sf::st_read) %>%                                     # read all shapefiles within a list of sf objects
  #   do.call(what =  dplyr::bind_rows, args = .) %>%             # merged the list of dataset into one sf object
  #   sf::st_make_valid() %>%                                     # Correct the geometry
  #   sf::st_transform(sf::st_crs(3347))                          # used of an CRS that give distance measure into meters (in Canada). CRS is the same as StatsCan Census 2022
  
  # Read DMTI AP ----
  # DMTI.AP <- 
  #   paste0("DMTI_AddressPoints_2021/ExposureWorkflow/", "AddressPointsDBUID.csv") %>% 
  #   data.table::fread(integer64 = "numeric") %>%  
  #   st_as_sf(wkt = 'geometry', crs = st_crs(4326)) %>%          # convert as sf object since the data is a csv with point into the 4326 ESPG
  #   sf::st_transform(sf::st_crs(3347)) %>%                      # used of an CRS that give distance measure into meters (in Canada). CRS is the same as StatsCan Census 2022
  #   {
  #     .[which(!st_is_empty(st_geometry(.))), ]                  # Select non-empty point
  #   }
}



# Read Dissemination Blocks shape file ----
DB.sf <-
  paste0("dissemination-blocks/Digital_Boundary_Files",
         "ldb_000b21a_e.shp") %>% 
  sf::st_read() %>% 
  sf::st_make_valid()
DB.sf$DBUID <- as.numeric(DB.sf$DBUID)                        # to match the format of Address Points
if(sf::st_crs(DB.sf) != sf::st_crs(3347)){
  DB.sf <- st_transform(DB.sf, sf::st_crs(3347))
}

#...............................................................................
# Import PRNAME dictionary
PRNAME.dic <- data.table::fread("dissemination-blocks/PRNAME_dic.csv")
# PRNAME from Microsoft
MS.PRNAME <- "Microsoft/Data_Microsoft" %>%                             # select the corrected folder where the data is stored
  list.files() %>% 
  tools::file_path_sans_ext()
# Convert PRUID column’s class to be in the same form as the dictionary.
DB.sf$PRUID <- as.numeric(DB.sf$PRUID)

if(isTRUE(loop.Q)){
  # Read data in sequence rather than in one big chunk
  for(
    pr in seq_len(nrow(PRNAME.dic))
  ){
    MS.pr <-
      PRNAME.dic[pr, MS_name] %>%
      paste0("Microsoft/Data_Microsoft/", ., ".geojson") %>%
      geojsonsf::geojson_sf() %>%                                 # read geojson file into an sf object
      sf::st_make_valid() %>%                                     # Correct the geometry
      sf::st_transform(sf::st_crs(3347))                          # used of an CRS that give distance measure into meters (in Canada). CRS is the same as StatsCan Census 2022
    
    # FOR LOOP (intersects) 
    DB.sf.pr <- DB.sf[DB.sf$PRUID == PRNAME.dic[pr, PRNAME],] # choose the Dissemination Block that are in the current iteration province 
    DB.dt.pr <- as.data.table(DB.sf.pr)
    # compute the bbox ON CD (to be sure that bbox contain but not limit the amount of the real intersection between the two algorithm)
    DB.dt.pr[
      ,c("xmin","ymin","xmax","ymax") := {
        DB.dt.pr$geometry %>%
          sapply(st_bbox) %>%                                     # apply function to each geometry
          t() %>%                                                 # transpose to have xmin (and other) in column insted of row
          apply(2,as.list)                                        # convert back to a list
      }
    ]
    # select the correct values, store information into the DB database
    {
      MS.pr.poly <- st_cast(MS.pr, "POLYGON") # Can multipolygon be remove?? See : MS.pr[which(st_geometry_type(MS.pr) != "POLYGON"),] %>% mapview::mapview()
      mp_crds <- as.data.table(sf::st_coordinates(MS.pr.poly))
      
      # Test which point is within the bbox of the Dissemination Block
      # Cluster: 
      numCores <- detectCores()
      cl <- makeCluster(numCores-2)
      parallel::clusterExport(cl= cl,varlist = c("DB.dt.pr", "mp_crds")) # export the needed variables into each cluster
      parallel::clusterEvalQ(cl = cl, {                                  # export library into each cluster
        library(data.table)
        library(dplyr)
      })
      
      # NOT OPTIMAL PARALLEL STRUCTURE, NEEDED TO BE REWRITE IF REUSED. 
      # In parallel, find in sequence which point of a polygon is into the bbox. The code may be speed-up because of the rowid structure of data.table
      DB.dt.pr[,MS.pr.select:= NA]
      DB.dt.pr$MS.pr.select <- 
        pbapply::pblapply(
          cl = cl,
          X = seq_len(nrow(DB.dt.pr)), # THE VECTOR SHOULD BE SPLIT INTO THE NUMBER OF OPENED'S CLUSTER
          FUN = function(ii){
            # If one point is within the bbox interval, then keep the polygon ID to test the intersect
            mp_crds[
              X >= unlist(DB.dt.pr[ii]$xmin) & 
                X <= unlist(DB.dt.pr[ii]$xmax) &
                Y >= unlist(DB.dt.pr[ii]$ymin) &
                Y <= unlist(DB.dt.pr[ii]$ymax),
              unique(L2)                                                 # L2 is an indicator of which row of MS is within the bbox of the DB. L1 always equal to 1 since we do not have more than two dimensions polygons.
            ]
          }
        )
      stopCluster(cl)
    } 
    
    # Vector of DBUID is used in writing the file that return the intersects match.
    DBUID.c <- DB.sf.pr$DBUID 
    
    # Creating the progress bar for writing the intersepts results for the current processed province.
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta",
                           total = nrow(DB.sf.pr))
      
    i <- 0 # Increment to track which DBUID is currently being written
    for(DBUID.i in DB.sf.pr$geometry){
      i <- i + 1                                                         # increment +1 for writing the next file
      intersept.i <- sf::st_intersects(                                  # Compute the intersect
        DBUID.i,
        MS.pr[DB.dt.pr[i]$MS.pr.select[[1]],],
        sparse = FALSE
      ) %>% t()
      
      # Write the DBUID and the id on my HDD
      which(intersept.i) %>% 
        {DB.dt.pr[i]$MS.pr.select[[1]][.]} %>%                           # extract the ID in MS from the PR file meaning that each province has it own ID. Thus, there is multiple 1.2.etc across the province
        write.csv(
          paste0("Microsoft/Microsoft_DBUID/",DBUID.c[i],".csv")         # will store DB name to minimize the amount of data to be written
        )
      pb$tick()                                                          # update progressbar
    }                                                                    # next province
  }
}
# Lecture of the intersects between microsoft and DBUID data base and correction of 
# polygon in MS that had not been assigned with a DBUID
{
  # PRNAME dictionary
  PRNAME.dic <- fread("dissemination-blocks/PRNAME_dic.csv")
  
  # Convert PRUID column’s class to be in the same form as the dictionary.
  DB.sf$PRUID <- as.numeric(DB.sf$PRUID)
  
  MS_DBUID.path <- "Microsoft/Microsoft_DBUID" %>% 
    list.files(full.names = TRUE)
  
  MS_DBUID.DBUID <- "Microsoft/Microsoft_DBUID" %>% 
    list.files() %>% 
    tools::file_path_sans_ext()
  
  MS_DBUID.ls <- lapply(MS_DBUID.path, fread)
  for(i in seq_along(MS_DBUID.ls)){
    MS_DBUID.ls[[i]][, DBUID := as.numeric(MS_DBUID.DBUID[i])]
  }
  MS_DBUID <- data.table::rbindlist(MS_DBUID.ls, fill = TRUE)[,V2:=NULL]
  
  # to read the number of row for each geojson file in order to match the intersects id in MS_DBUID
  MS.file.name <- 
    data.table(
      PRNAME_path = list.files("Microsoft/Data_Microsoft/", full.names = TRUE),
      MS_name = list.files("Microsoft/Data_Microsoft/") %>% tools::file_path_sans_ext(),
      max_poly = NA
    )
  
  for(i in seq_len(nrow(MS.file.name))){
    MS.file.name$max_poly[i] <- MS.file.name$PRNAME_path[i] %>%
      geojsonsf::geojson_sf() %>% 
      nrow()
  }
  
  MS.file.name[
    PRNAME.dic,
    PRNAME := PRNAME,
    on = .(MS_name)
  ]
  
  MS_DBUID[,PRNAME := floor(DBUID / (10**(9)))]
  MS_DBUID[
    MS.file.name,
    max_poly := max_poly,
    on = .(PRNAME)
  ]
  setindexv(MS_DBUID,cols = c("PRNAME", "x"))
}
# Correct the non intersects polygon 
if(isTRUE(loop.Q)){
  # Extract where there is no match between footprint and DBUID 
  not.in.seq_poly <- lapply(MS.file.name$PRNAME, function(pr){
    data.table::data.table(
      PRNAME = pr,
      no.match = which(!(
        seq_len(MS.file.name[PRNAME == pr,max_poly]) %in%
        unique(MS_DBUID[(!is.na(x)) & (PRNAME == pr), x])
      ))
    )
  }) %>%
    data.table::rbindlist()
  
  # % of matches (diagnosis)
  MS.file.name[
    na.omit(MS_DBUID)[,unique(x), by = .(PRNAME)][,.(num.match = .N), by = .(PRNAME)], # count the number of matches
    match.p100 := (num.match)/max_poly, 
    on = .(PRNAME)
  ][]
  not.in.seq_poly[,.N, by = .(PRNAME)]
  
  # NEAREST MATCH WITH DBUID AND FOOTPRINT
  PRNAME.to.correct <- not.in.seq_poly[!is.na(no.match), unique(PRNAME)]
  
  # Creating the progress bar for writing the intersepts results for the current processed province.
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta",
                         total = length(PRNAME.to.correct))
  for(pr in PRNAME.to.correct){
    pb$tick()
    MS.pr <-
      PRNAME.dic[PRNAME == pr, MS_name] %>%
      paste0("Microsoft/Data_Microsoft/", ., ".geojson") %>%
      geojsonsf::geojson_sf() %>%                                 # read geojson file into an sf object
      sf::st_make_valid() %>%                                     # Correct the geometry
      sf::st_transform(sf::st_crs(3347))                          # used of an CRS that give distance measure into meters (in Canada). CRS is the same as StatsCan Census 2022
    
    # FOR LOOP (intersects) 
    DB.sf.pr <- DB.sf[DB.sf$PRUID == pr,]                         # Choose the Dissemination Block that are in the current iteration province 
    # Variable that identify which footprint has not been match with a DBUID
    no.match.id <- not.in.seq_poly[PRNAME == pr, no.match]
    
    nearest.match.id <- sf::st_nearest_feature(
      MS.pr[no.match.id,],
      DB.sf.pr
    )
    DBUID.match <- nearest.match.id %>% 
      {DB.sf.pr[.,]$DBUID}
    nearest.match.id.dt <- data.table::data.table(
      no.match.id,
      DBUID.match
    )
    # write correction
    nearest.match.id.dt %>% 
      write.csv(
        paste0("Microsoft/Microsoft_DBUID_nearest/",
               PRNAME.dic[PRNAME == pr, MS_name],".csv") # Will store DB name to minimize the amount of data to be written
      )
  }
  # END FOR LOOP ~~~~
}
# Merge all information ----
{
  MS_DBUID.corrected <- 
    "Microsoft/Microsoft_DBUID_nearest/" %>% 
    list.files(full.names = TRUE) %>% 
    lapply(fread, integer64 = "numeric") %>% 
    data.table::rbindlist() %>% 
    data.table::setnames(c("no.match.id","DBUID.match"), c("x","DBUID"))
  
  
  MS_DBUID.corrected[,PRNAME:=floor(DBUID / (10**(9)))]
  MS_DBUID[,max_poly:=NULL]
  # bind all information
  MS.DBUID <- data.table::rbindlist(
    list(
      MS_DBUID,
      MS_DBUID.corrected
    )
  ) %>% na.omit %>% 
    {.[order(x),.SD,by =.(PRNAME)]}                             # sort by PRNAME
  # COMPUTE THE AREA OF EACH FOOTPRINT POLYGON 
  # Read Microsoft Building Footprint dataset ----
  MS.ls <- 
    "Microsoft/Data_Microsoft" %>%                              # select the corrected folder where the data is stored
    list.files(full.names = TRUE) %>%                           # list all shapefile into 
    lapply(geojsonsf::geojson_sf)                               # read all shapefile into a list of sf object
    
  # order of which provinces in MS had been treated
  MS.PRNAME.order <- PRNAME.dic[sapply(MS.PRNAME, function(x) which(MS_name %in% x)), PRNAME]
  names(MS.ls) <- MS.PRNAME.order
  
  # order the match
  MS.DBUID <- MS.DBUID[
    unlist(lapply(MS.PRNAME.order, function(x) which(PRNAME %in% x))),
  ]
  MS.DBUID[,id:=x]
  
  # create an id to help the merge
  for(i in seq_along(MS.ls)){
    MS.ls[[i]]$id <- seq_len(nrow(MS.ls[[i]]))
  }
  
  # Merge DBUID information with the geometry
  MS.all.ls <- list()
  for(pr in MS.PRNAME.order){
    i <- which(names(MS.ls) %in% as.character(pr))
    MS.all.ls[[i]] <- dplyr::left_join(
      MS.ls[[i]],
      MS.DBUID[PRNAME == pr, .(DBUID, id)],
      by = c('id' = 'id')
    )
  }
  MS.all <- MS.all.ls %>% 
    do.call(what =  dplyr::bind_rows, args = .) %>%             # merged the list of dataset into one sf object
    sf::st_make_valid() %>%                                     # Correct the geometry
    sf::st_transform(sf::st_crs(3347))                          # used of an CRS that give distance measure into meters (in Canada). CRS is the same as StatsCan Census 2022
}
# Write MS into an shp
sf::st_write(
  MS.all,
  "Microsoft/shp_Microsoft_DBUID_all/Microsoft_DBUID.geojson"
)
#...............................................................................
#  THE FILE HAS BEEN SAVED IN THE FOLDER 
# /media/gabriel/HDD/Data_memoire/Microsoft/shp_Microsoft_DBUID_all/
# WITH THE FOLLOWING NAME : Microsoft_DBUID.geojson