library('move2')
library('lubridate')
library('dplyr')
library('magrittr')


source("makeEventClusters.R")
source("makeclustertable.R")
source("clustering.R")

`%!in%` <- Negate(`%in%`)


## The parameter "data" is reserved for the data object passed on from the previous app

# to display messages to the user in the log file of the App in MoveApps
# one can use the function from the logger.R file:
# logger.fatal(), logger.error(), logger.warn(), logger.info(), logger.debug(), logger.trace()

# Showcase injecting app setting (parameter `year`)
rFunction = function(data, clusterstart, clusterend, clusterstep = 1, clusterwindow = 7, clustexpiration = 14) {
  
  # Issues with the 'between' function within this MoveApp
  # Doesn't seem to recognise the POSIX as desired
  # Why?
  
  
  # Check suitability of inputs
  if(!is.instant(as.Date(clusterstart)) | is.null(clusterstart)) {
    logger.error(paste0("Start date of clustering ", clusterstart, " is not a valid date Defaulting to start date 2 weeks before final date."))
    clusterstart <- max(mt_time(data), na.rm = TRUE) - days(14)
  }
  
  if("behav" %!in% colnames(data)) {
    logger.fatal("Classified behaviour column 'behav' is not contained by input data. Unable to perform clustering. Please use classification MoveApp prior to this stage in the workflow")
    stop()
  }
  
  if(!is.instant(as.Date(clusterend)) | is.null(clusterend)) {
    logger.error(paste0("End date of clustering ", clusterend, " is not a valid date Defaulting to end date as most recent timestamp."))
    clusterend <- max(mt_time(data), na.rm = TRUE)
  }
  
  
  logger.info(paste0("Clustering between ", clusterstart, " and ", clusterend))
  
  # Performing clustering
  clusteredData <- clustering(data, as.Date(clusterstart), as.Date(clusterend), clusterstep, clusterwindow, clustexpiration)
  
  # Retrieve tagdata output
  clusteredTagData <- clusteredData$clustereventdata
  
  # Retrieving clustertable and releasing as artefact
  clustertable <- clusteredData$clustereventtable %>%
    mutate(xy.clust = ifelse(!is.na(xy.clust), paste0("A", xy.clust), NA)) %>%
    mt_as_move2(time_column = "firstdatetime", track_id_column = "xy.clust")
  
  # Create path and save
  dir.create(targetDirFiles <- tempdir())
  saveRDS(clustertable, file = paste0(targetDirFiles, "\\clustertable.rds")) 
  
  # Save artifact as .zip
  zip_file <- appArtifactPath(paste0("myfiles.zip"))
  zip::zip(zip_file,
           files = list.files(targetDirFiles, full.names = TRUE),
           mode = "cherry-pick")
  
  
  # Pass tag data onto next MoveApp
  return(clusteredTagData)
}
