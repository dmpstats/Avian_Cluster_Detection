library('move2')
library('lubridate')
library('dplyr')
library('magrittr')
library('tidyr')
library('Gmedian')
library('sf')
library('stringr')
library('units')

# Shortened 'not in':
`%!in%` <- Negate(`%in%`)

# Function to calculate geometric medians:
calcGMedianSF <- function(data) {
  med <- Gmedian::Weiszfeld(st_coordinates(data))$median %>% as.data.frame() %>%
    rename(x = V1, y = V2) %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(data)) %>%
    st_geometry()
  return(med)
}

rFunction <- function(data,
                      clusterstart,
                      clusterend,
                      clusterstep = 1, 
                      clusterwindow = 7, 
                      clustexpiration = 14, 
                      behavsystem = TRUE, 
                      d,
                      clustercode = "") {
  
  # --------------------------------------------------------------
  # Setup & Input Checks -----------------------------------------

  # Check clustercode
  if (clustercode != "") {
    logger.trace(paste0("Provided clustercode is ", clustercode))
    clustercode <- paste0(clustercode, ".")
  } else {
    logger.warn("No clustercode provided. Defaulting no clustercode. ")
    clustercode <- ""
  }
  
  # Check suitability of inputs
  if(!is.instant(as.Date(clusterstart)) | is.null(clusterstart)) {
    logger.error(paste0("Start date of clustering ", clusterstart, " is not a valid date Defaulting to start date 2 weeks before final date."))
    clusterstart <- max(mt_time(data), na.rm = TRUE) - days(14)
  }
  
  if(behavsystem == TRUE & "behav" %!in% colnames(data)) {
    logger.fatal("Classified behaviour column 'behav' is not contained by input data. Unable to perform clustering. Please use classification MoveApp prior to this stage in the workflow")
    stop()
  }
  
  if(!is.instant(as.Date(clusterend)) | is.null(clusterend)) {
    logger.error(paste0("End date of clustering ", clusterend, " is not a valid date Defaulting to end date as most recent timestamp."))
    clusterend <- max(mt_time(data), na.rm = TRUE)
  }
  
  
  logger.info(paste0("Clustering between ", clusterstart, " and ", clusterend))
  
  
  # Filter to relevant time window for overall clustering:
  clusterstart %<>% as_date()
  clusterend %<>% as_date()
  eventdata <- data %>%
    filter(between(mt_time(data), clusterstart - days(clusterwindow), clusterend + days(1)))
  
  
  
  # ------------------------------------------------------------------
  # Clustering Loop
  #'  Operate on a rolling-window basis starting from clusterstart
  #'  We increment by clusterstep (in days) and re-cluster until we reach the clusterend
  #'  
  #'  'clusterdate' parameter will define the final day of data we are currently clustering up to
  # -----------------------------------------------------------------
  
  clusterdate <- floor_date(clusterstart, unit = "days")
  
  # Set up the initial case, in which we have no roll-over of cluster data:
  clusterDataDwnld <- NULL
  tempclustertable <- NULL
  laststep <- FALSE
  rollingstarttime <- Sys.time()
  
  # Begin loop:
  while (
    #clusterdate < floor_date(clusterend, unit = "days") + days(clusterwindow)
    laststep == FALSE # testing alternative 
    ) {
    
    
    # If we're on the final step, set the clusterdate equal to final day:
    if (clusterdate >= floor_date(clusterend, unit = "days")) {
      logger.trace(paste0("Current clusterdate ", as.Date(clusterdate), " is beyond final date. Assigning final date, ", as.Date(clusterend)))
      clusterdate <- floor_date(clusterend, unit = "days")
      laststep <- TRUE
    }

    #' -----------------------------------------------------------------------
    #' 1. Data Setup and Import ----------------------------------------------
    #'  Here, filter to the relevant window and import
    #'  data from the last rolling-window (if available)
    
    # Import cluster data from previous run
    # This will contain only the key information for the clustering rolling-window 
    # i.e. location data and timestamp
    if (!is.null(tempclustertable)) {
      clusterDataDwnld <- tempclustertable
    }
    
    # Filter the location data down to our clustering window
    clusteringData <- filter(eventdata,
                             # Take the number of days we need for clustering:
                             between(mt_time(eventdata), 
                                     clusterdate - days(clusterwindow), 
                                     clusterdate 
                             ))
    
    skipToNext <- FALSE # Determines whether to jump this clusterstep
    
    # Check that there is enough tracking data to cluster, and skip if not:
    if (nrow(clusteringData) < 3) {
      logger.trace(paste0(as.Date(clusterdate), ":      Clustering complete - not enough data within clustering period"))
      skipToNext <- TRUE}
    
    # And if using the behavioural classification system, 
    # check there is enough non-travelling behaviour:
    if (behavsystem == TRUE) {
      if (sum(clusteringData$behav != "STravelling") < 2) { 
        logger.trace(paste0(as.Date(clusterdate), ":      Clustering complete - not enough stationary behaviour within clustering period"))
        skipToNext <- TRUE}}
    
    # If either of these conditions are met, skip ahead to the next possible clusterdate:
    if (skipToNext == TRUE) {
      # Find minimum following date on which we have data
      clusterdate <- filter(eventdata, mt_time(eventdata) > clusterdate) %>%
        mt_time() %>%
        min() + days(clusterwindow)
      logger.trace(paste0("     Skipping to next date with data to cluster: ", as.Date(clusterdate)))
      next}
    
    
    
    # ----------------------------------------------------------------------
    #' 2. Perform Clustering ----------------------------------------------
    #' Now that we know there is enough data, generate the new clusters
    #' 

    
    # If using behavioural classification, filter to stationary behaaviours:
    if (behavsystem == TRUE) {clusterpoints <- filter(clusteringData,
                                                      behav %!in% c("STravelling", "Unknown"))
    } else {clusterpoints <- clusteringData} 
    
    logger.trace(paste0(as.Date(clusterdate), ":     Creating new clusters for ", nrow(clusterpoints), " locations"))
    
    # Generate clusters:
    clusterpoints %<>% mutate(
      ID = 1:nrow(.),
      X = st_coordinates(.)[, 1],
      Y = st_coordinates(.)[, 2]
    )
    
    # Build distance matrix (no SF method):
    mdist <- dist(cbind(clusterpoints$X, clusterpoints$Y))
    hc <- hclust(as.dist(mdist), method = "complete")
    
    # Add cluster information to data:
    clusterpoints$clust <- cutree(hc, h = d)
    genclustertable <- clusterpoints
    
    # Filter out clusters with only one location: (increasing to 3)
    cid <- as.vector(which(table(clusterpoints$clust) < 3))
    genclustertable %<>% filter(clust %!in% cid)
    
    
    # Append cluster data ---------------------------------------------
    # Generate a median location for each centroid:
    clusts <- genclustertable %>%
      group_by(clust) %>%
      summarise(geometry = st_combine(geometry))
    for (j in 1:nrow(clusts)) {
      clustid <- clusts$clust[j]
      clusts$geometry[j] <- calcGMedianSF(clusts[j,])
    }

    # Create output table:
    clusts %<>% 
      mutate(xy.clust = paste0("up", clust)) %>%
      dplyr::select(-clust)
    
    # Output location data:
    xydata <- clusterpoints %>% ungroup() %>%
      mutate(xy.clust = paste0("up", clust)) %>%
      dplyr::select(-clust) %>%
      
      # Remove the non-significant clusters:
      mutate(xy.clust = if_else(xy.clust %in% unique(clusts$xy.clust), xy.clust, "upNA"))
    
    

    # -----------------------------------------------
    # 3. Identify Matching Clusters ----------------------
  
    newclust <- clusts
    logger.trace(paste0(as.Date(clusterdate), 
                        ":       ", 
                        nrow(newclust),
                        " new clusters generated"))
    
    if (is.null(clusterDataDwnld) | nrow(newclust) == 0) {
      # In this case, no matches are possible
      matchingclustermap <- data.frame(existID = NULL, updID = NULL)
      if (is.null(clusterDataDwnld)) {
        existingclust <- NULL
      }
      
    } else {

      existingclust <- clusterDataDwnld %>%
        filter(lastdatetime + days(clustexpiration) > clusterdate - days(clusterwindow)) %>%
        arrange(xy.clust) 
      
      # Generate distance matrix:
      # Rows are existing clusts, columns are new clusts
      dists <- st_distance(existingclust, newclust) %>% units::drop_units()
      rownames(dists) <- existingclust$xy.clust
      colnames(dists) <- newclust$xy.clust
      
      # Identify close clusters (within 175m) and arrange in table:
      closeClusterIndices  <- try(as_tibble(
        which(dists < 175, arr.ind = T, useNames = T)))
      closeClusterIndices$row <- rownames(dists)[closeClusterIndices$row]
      closeClusterIndices$col <-  colnames(dists)[closeClusterIndices$col]
      
      if (nrow(closeClusterIndices) > 0) {
        matchingclustermap <- closeClusterIndices %>%
          as.data.frame() %>% 
          rename(existID = row,
                 updID = col)
      } else {
        matchingclustermap <- data.frame(
          existID = NULL, updID = NULL # nullify if no clusters matched
        )
      }

    }
    
    
    # -----------------------------------------------
    # Merging Clusters: Case 1 ----------------------
    #'
    #' Multiple old clusters -> 1 new cluster
    #' In this case, we want to re-allocate the constituent locations
    #' of the new to the nearest 'old' cluster
    
    # Check to see if this is true:
    if (length(unique(matchingclustermap$updID)) != nrow(matchingclustermap)) {
      
      # Retrieve their IDs
      ids <- which(table(matchingclustermap$updID) > 1)
      for (u in 1:length(ids)) {
        
        # 1: Get merge, location, and old cluster data
        tempdat <- filter(matchingclustermap, updID == names(ids)[u])
        allocpoints <- filter(xydata, xy.clust %in% tempdat$updID)
        existingclustpoints <- filter(clusterDataDwnld, xy.clust %in% tempdat$existID)
        
        # Update cluster names:
        tempdat$updID <- paste0(tempdat$updID, ".", tempdat$existID)
        
        # Generate distances between new cluster locations and old cluster medians:
        pdists <- st_distance(allocpoints, existingclustpoints) %>% 
          units::drop_units() %>%
          as.data.frame() %>%
          mutate(
            # Get index of nearest cluster
            nearind = apply(., 1, which.min),
            allocclust = existingclustpoints$xy.clust[nearind]
          )
        
        # Match to nearest:
        allocpoints$xy.clust <- tempdat$updID[match(pdists$allocclust, tempdat$existID)]
        
        # Remove from clustermap:
        matchingclustermap %<>% filter(updID != names(ids)[u]) %>%
          bind_rows(tempdat)
        }
      }
    
    # -----------------------------------------------
    # Merging Clusters: Case 2 ----------------------
    #'
    #' 1 old cluster -> multiple new clusters
    #' In this case, we just want to keep the old cluster
    #' and we reallocate all new clusters to the same ID
    
    # Check to see if this is true:
    if (length(unique(matchingclustermap$existID)) != nrow(matchingclustermap)) {
      
      # Retrieve their IDs
      ids <- which(table(matchingclustermap$existID) > 1)
      for (u in 1:length(ids)) {
        
        # 1: Deal with cluster-matching table
        tempdat <- filter(matchingclustermap, existID == names(ids)[u])
        newname <- paste0(tempdat$updID[1], ".", tempdat$existID[1]) # new cluster ID
        newmatch <- tempdat[1,]
        newmatch$updID <- newname # this is the match we want
        matchingclustermap %<>% bind_rows(newmatch) # add to the cluster-match table
        
        # 2: Deal with clustertable just by cutting the cluster out
        # They'll be re-introduced with the matching cluster-map
        newclust <- filter(newclust, xy.clust %!in% tempdat$updID)
        
        # 3: Deal with tracking data by assigning all to new cluster
        xydata$xy.clust[xydata$xy.clust %in% tempdat$updID] <- newname
      }
    }
    
    logger.trace(paste0(as.Date(clusterdate), ":       ", nrow(matchingclustermap), "  matched to existing clusters"))
    
    
    # -----------------------------------------------
    # Merging Clusters: Case 3 ----------------------
    #'
    #' Pre-existing clusters with no  merges

    # Retrieve the existing IDs not being merged:
    existIDs <- existingclust$xy.clust[existingclust$xy.clust %!in% matchingclustermap$existID]
    
    # Generate table matching them to nothing
    if (length(existIDs) != 0) {
      existclustermap <- data.frame(
        existID = existIDs,
        updID = NA
      )
    } else {
      existclustermap <- data.frame(existID = NULL, upID = NULL)
    }
    logger.trace(paste0(as.Date(clusterdate), ":       ", nrow(existclustermap), "  clusters not matched"))
    
    
    # -----------------------------------------------
    # Merging Clusters: Case 4 ----------------------
    #'
    #' New clusters with no previous match (didn't previously exist)
    
    # Retrieve the IDs of these clusters:
    newIDs <- newclust$xy.clust[newclust$xy.clust %!in% matchingclustermap$updID]
    
    # We find the lowest 'free' value at which we can start their cluster IDs:
    if (length(newIDs) != 0) {
      if (!is.null(existingclust)) {
        startid <- max(clusterDataDwnld$xy.clust) + 1
      } else {startid <- 1}
      
      # Generate table with these new IDs:
      newclustermap <- data.frame(existID = startid:(startid - 1 + length(newIDs)), 
                                  updID = newIDs)
    } else {
      newclustermap <- data.frame(existID = NULL, upID = NULL)
    }
    logger.trace(paste0(as.Date(clusterdate), ":       ", nrow(newclustermap), "  new clusters have no match"))
    
    
    # ----------------------------------------------
    # 4. Perform Matching ----------------------------
    
    # Combine the three merge tables (and an option for NA clusters):
    updatedclustermap <- rbind(
      matchingclustermap,
      newclustermap,
      existclustermap,
      c(NA, "upNA")
    )
    colnames(updatedclustermap) <- c("existID", "updID")
    
    # First, match location data:    
    xytagdata <- xydata %>%
      rename(updID = xy.clust) %>%
      left_join(updatedclustermap) %>% # Match clusters
      rename(xy.clust = existID) %>%
      mutate(xy.clust = as.numeric(xy.clust)) %>%
      dplyr::select(-updID) %>% # Drop update column
      suppressMessages()

    # Secondly, cluster data:
    newclust <- newclust %>%
      rename(updID = xy.clust) %>%
      left_join(updatedclustermap) %>%
      rename(xy.clust = existID) %>%
      mutate(xy.clust = as.numeric(xy.clust)) %>%
      dplyr::select(-updID) %>%
      suppressMessages()

    # Combine previous and updated clusterdata:
    if (is.null(existingclust) | is.null(newclust)) {
      updatedClusters <- newclust # nothing to add in this case
    } else {
      updatedClusters <- rbind(existingclust[, c("xy.clust", "geometry")], 
                               newclust) %>%
        group_by(xy.clust) %>%
        summarise(geometry = st_union(geometry)) %>% # Join geometries together into multipoint
        st_centroid() # Find mean location of two centroids (we should change this later)
    }

    # Add cluster data into main location data:
    data %<>%
      filter(index %!in% xytagdata$index) %>% # remove rows whose clusters need to be updated
      bind_rows(., xytagdata) %>% # add them back in with the new update
      arrange(mt_track_id(.), mt_time(.)) # reorder

    
    # ----------------------------------------------------
    # 5. Create TEMPORARY clustertable ----------------------
    # This will contain only the essential information for
    # cluster-updating: location, time, clustID
    
    logger.trace(paste0(as.Date(clusterdate), ":     Generating temporary clustertable"))
    
    # Perform this step only if there are updates to perform on the tag data:
    if (nrow(filter(data, xy.clust %in% updatedClusters$xy.clust)) != 0) {
      
      # Bind necessary data on updated clusters:
      tempclusts <- filter(data, xy.clust %in% updatedClusters$xy.clust) %>% mutate(
        datetime = mt_time(.)
      ) 
      tempclusts %<>%
        group_by(xy.clust) %>%
        summarise(
          firstdatetime = min(datetime),
          lastdatetime = max(datetime),
          geometry = st_combine(geometry)
        ) 
      
      # Generate median location:
      for (j in 1:nrow(tempclusts)) {
        tempclusts$geometry[j] <- calcGMedianSF(tempclusts[j,])
      }
      clustertable_update <- tempclusts
    } else {
      
      # If impossible (i.e. no clusters exist), nullify this
      clustertable_update <- NULL
    }
    
    
    # -----------------------------------------------------------------------
    # 6. Return other clusters to the dataset ----------------------------------
    logger.trace(paste0(as.Date(clusterdate), ":     Merging back into clusterDataDwnld"))
    
    if (!is.null(clusterDataDwnld)) {
      
      clusterDataDwnld <- bind_rows(clustertable_update,
                                    filter(clusterDataDwnld, 
                                           xy.clust %!in% clustertable_update$xy.clust)) %>%
        arrange(desc(firstdatetime))
    } else {
      
      if (!is.null(clustertable_update)) {
        # If there is an update to the clustertable 
        # but no original data:
        clusterDataDwnld <- clustertable_update %>%
          arrange(desc(firstdatetime))
      } else {
        # Otherwise, no clusters whatsoever
        clusterDataDwnld <- NULL
      }
    }
    

    
    # ----------------------------------------------------
    # 7. Generate output & move on rolling window -----------
    
    # Log progress:
    
    logger.trace(paste0(as.Date(clusterdate), ":     COMPLETE. Total of ", nrow(clusterDataDwnld), " cumulative clusters generated"))
    if (laststep == FALSE) {
      logger.trace(paste0(as.Date(clusterdate), ": Increasing clusterdate to ", clusterdate + days(clusterstep)))
      clusterdate <- clusterdate + days(clusterstep)
    }
    

  
    # End of clustering loop ---------
  }
  rollingendtime <- Sys.time()
  logger.trace(paste0("CLUSTERING COMPLETE. Time taken: ", 
               difftime(rollingendtime, rollingstarttime, units = "mins"),  " mins. Generating clustertable for all clusters"))
  
  
  # ----------------------------------------------------
  # 8. Generate FULL clustertable ----------------------
  # This will contain output cluster information
  
  tablestarttime <- Sys.time()
  clustertable <- data %>%
    filter(!is.na(xy.clust)) %>%
    
    # Essential data for later computation:
    mutate(
      tag = mt_track_id(.),
      datetime = mt_time(.),
      timediff_hrs = mt_time_lags(.) %>% 
        units::set_units("minutes") %>% 
        units::drop_units()/60,
      hour = mt_time(.) %>% 
        lubridate::hour()
    ) %>%
    group_by(xy.clust) %>%
    summarise(geometry = st_combine(geometry),
              
              # Time data:
              firstdatetime = min(datetime),
              lastdatetime = max(datetime),
              days = length(unique(date(datetime))),
              totalduration = (ceiling_date(lastdatetime, unit = "days") - floor_date(firstdatetime, unit = "days")) %>% as.integer(),
              daysempty = as.numeric(totalduration - days),
              TimeTotal = sum(timediff_hrs, na.rm = T),
              TimeDay = sum(timediff_hrs[hour > 5 & hour < 17], na.rm = T),
              TimeFeed = sum(timediff_hrs[behav == "SFeeding"], na.rm = T),
              MedianHourFeed = median(hour[behav == "SFeeding"], na.rm = T),
              MedianHourDay = median(hour[hour > 5 & hour < 17], na.rm = T),
              
              # Behavioural data:
              nbirds = length(unique(tag)),
              birds = paste(unique(tag), collapse = ", "),
              Total = n(),
              SFeeding = sum(behav == "SFeeding"),
              SRoosting = sum(behav == "SRoosting"),
              SResting = sum(behav == "SResting"),

              DistMedian = NA,
              DistSD = NA,
              dist_to_bird_km = NA,
              within_25k = NA,
              within_50k = NA,
              nightdist_mean = NA, 
              nightdist_med = NA, 
              nightdist_sd = NA,
              arrivaldist_mean = NA,
              arrivaldist_med = NA) 

  rm(list = c("clusterDataDwnld", "clusteringData", "clusterpoints", "eventdata", "genclustertable", "xydata", "xytagdata", "mdist")) # clear some storage for next operations
  clustertable %<>% mt_as_move2(time_column = "firstdatetime", track_id_column = "xy.clust") # switch to move2 for ease
  
  # Generate distance data:
  logger.trace(paste0("Generating distance data for all clusters. This may run slowly"))


  # CLUSTERTABLE LOOPING -------------------------------------------
  
  # Convert data to matrix format for easy searching
  mat.data <- data %>%
    as.data.frame() %>%
    rename(trackID = mt_track_id_column(data),
           timestamp = mt_time_column(data)) %>%
    dplyr::select(all_of(c("trackID", "timestamp", "xy.clust"))) %>%
    cbind(st_coordinates(data)) %>%
    as.matrix()
  
  tempdat <- NULL
  for (k in 1:nrow(clustertable)) {
    
    if (mod(k, 100) == 1) {logger.trace(paste0("    ", round(100 * (k / nrow(clustertable)), 3), "% complete"))} # Log progress
    
    # Convert cluster to its constituent points: 
    clustdat <- st_geometry(clustertable[k,]) %>%
      st_cast("POINT") 
    dists <- clustdat %>%
      st_distance() %>% 
      units::set_units("metres") %>% # to be safe
      units::drop_units()
    
    # Median and SD distance between points:
    dists <- dists[upper.tri(dists)] %>%
      as.vector() 
    tempmed <- median(dists)
    tempsd <- sd(dists)
    clustertable$DistMedian[k] <- tempmed
    clustertable$DistSD[k] <- tempsd
    
    # Generate final median location:
    clustertable$geometry[k] <- calcGMedianSF(clustdat)
    
    # Generate revisit data:
    birdsinclust <- strsplit(clustertable[k,]$birds, split = ", ") %>% unlist()
    tempdat <- data %>%
      ungroup() %>%
      filter(between(
        mt_time(.), 
        clustertable$firstdatetime[k],
        clustertable$lastdatetime[k]),
      mt_track_id(.) %in% birdsinclust) %>%
      mutate(incluster = ifelse(xy.clust != clustertable$xy.clust[k] | is.na(xy.clust), 0, 1),
             indaycluster = ifelse(incluster == 1 & between(hour(mt_time(.)), 10, 15), 1, 0)) %>%
      mutate(temptag = mt_track_id(.),
             tempdate = ifelse(incluster == 0, NA, as_date(mt_time(.)))) %>%
      st_drop_geometry() %>%
      as.data.frame() %>%
      group_by(temptag) %>%
      summarise(
        xy.clust = clustertable$xy.clust[k],
        visitsinevent = sum(rle(incluster)$values),
        dayvisits = sum(rle(indaycluster)$values),
        dayvisits = pmin(dayvisits, visitsinevent), # Fix case where dayvisits > totalvisits 
        ndays = n_distinct(tempdate),
        meanvisits = visitsinevent / ndays,
        meandayvisits = dayvisits / ndays) %>%
      # Take means across birds:
      summarise(xy.clust = clustertable$xy.clust[k],
                #visitsinevent_tot = sum(visitsinevent),
                visitsinevent_mean_pday = mean(meanvisits),
                #dayvisits_tot = sum(dayvisits),
                dayvisits_mean_pday = mean(meandayvisits)) %>%
      bind_rows(tempdat, .)
    
    # Generate data on nearest birds/ids -----------------------
    
    clustdat <- clustertable[k,]
    locavailable <- FALSE
    while (locavailable == FALSE) {
      # Filter to time window of cluster and remove all birds in the cluster (travelling + feeding):
      neartags <- data %>% 
        filter(between(timestamp, clustdat$firstdatetime, clustdat$lastdatetime),
               individual_local_identifier %!in% unlist(stringr::str_split(clustdat$birds, pattern = ", "))
        )
      if (nrow(neartags) != 0) {
        locavailable <- TRUE
      } else {
        # If no locations in window, we extend it 1hr each way and check again
        clustdat$firstdatetime <- clustdat$firstdatetime - lubridate::hours(1)
        clustdat$lastdatetime <- clustdat$lastdatetime + lubridate::hours(1)
      }
    }
    # Get the nearest location
    dist <- neartags[st_nearest_feature(st_geometry(clustdat),
                                             neartags),] %>%
      st_distance(clustdat)
    clustertable$dist_to_bird_km[k] <- units::set_units(dist, "metres") %>%
      units::drop_units() %>%
      divide_by(1000)
    
    # Identify active tags
    activetags <- data %>% 
      filter(between(mt_time(.), clustdat$firstdatetime, clustdat$lastdatetime)) %>%
      mt_track_id(.) %>%
      unique()
    
    # Start with 50km buffer:
    rad50 <- st_buffer(clustdat, dist = 50000)
    neartags2 <- data %>% 
      mutate(X = st_coordinates(.)[, 1],
             Y = st_coordinates(.)[, 2]) %>%
      filter(between(mt_time(.), clustdat$lastdatetime - days(30), clustdat$lastdatetime),
             between(X, st_coordinates(clustdat)[, 1] - 50000, st_coordinates(clustdat)[, 1] + 50000),
             between(Y, st_coordinates(clustdat)[, 2] - 50000, st_coordinates(clustdat)[, 2] + 50000),
             mt_track_id(.) %in% activetags)
    nearpoints50 <- neartags2[st_contains(rad50, neartags2) %>% unlist,]
    nearbirds50 <- mt_track_id(nearpoints50) %>%
      unique() %>% length()
    
    # Use this subset for 25km buffer:
    rad25 <- st_buffer(clustdat, dist = 25000)
    nearpoints25 <- neartags2[st_contains(rad25, neartags2) %>% unlist,]
    nearbirds25 <- mt_track_id(nearpoints25) %>%
      unique() %>% length()
    
    clustertable$within_50k[k] <- nearbirds50
    clustertable$within_25k[k] <- nearbirds25
    
    # Generate night-distance data -------------------------------------
    

    # Isolate locations of cluster-involved birds over its full timespan
    clustpoints <- data %>% 
      filter(
        between(mt_time(.),
          clustdat$firstdatetime - days(1),
          clustdat$lastdatetime),
        mt_track_id(.) %in% (strsplit(clustdat$birds, split = ", ") %>% unlist()))
    atevent <- clustpoints %>% # select only inside-cluster points
      filter(xy.clust == clustertable$xy.clust[k])
    
    # Extract the days on which each bird was at the event:
    daysbybird <- table(mt_track_id(atevent), as_date(mt_time(atevent))) %>%
      as.data.frame() %>%
      filter(Freq > 0)
    nightdat <- atevent[0,]
        # Loop through to create nightdist dataset:
    
    nightdist_temp <- data.frame(id = unique(daysbybird$Var1), nightdist_mean = NA, nightdist_med = NA, nightdist_sd = NA)
    
    # Loop by bird to get individual distances
    for (j in 1:length(unique(daysbybird$Var1))) {
      id <- unique(daysbybird$Var1)[j] # extract id
      daysbybird_temp <- filter(daysbybird, Var1 == id) # extract visits data
      
      for (m in 1:nrow(daysbybird_temp)) {
        newdat <- clustpoints %>%
          filter(mt_track_id(.) == daysbybird_temp[m, 1],
                 date(mt_time(.)) == as_date(daysbybird_temp[m, 2]))
        nightdat <- mt_stack(nightdat, newdat, .track_combine = "merge")
      }
      
      nightdat %<>% mutate(day = case_when(
        hour(mt_time(.)) > 21 | hour(mt_time(.)) < 3 ~ 0,
        TRUE ~ 1
      ))
      if (sum(nightdat$day) != nrow(nightdat)) {
        nightdists <- st_distance(
          nightdat[nightdat$day == 1,],
          nightdat[nightdat$day == 0,]
        ) %>% units::drop_units()
        nightdist_temp$nightdist_mean[j] <- mean(nightdists, na.rm = T)
        nightdist_temp$nightdist_med[j] <- median(nightdists, na.rm = T)
        nightdist_temp$nightdist_sd[j] <- sd(nightdists, na.rm = T)
      }
    }
    

    # Find mean of by-bird operations and add to clustertable
    clustertable$nightdist_mean[k] <- mean(nightdist_temp$nightdist_mean, na.rm = T)
    clustertable$nightdist_med[k] <- mean(nightdist_temp$nightdist_med, na.rm = T)
    clustertable$nightdist_sd[k] <- mean(nightdist_temp$nightdist_sd, na.rm = T)
    
    # Generate night-before distance data --------------------------
    
    # Identify each ID's first date of arrival
    firstarrivals <- daysbybird %>%
      as.data.frame() %>%
      group_by(Var2) %>%
      slice_min(order_by = Var2) %>%
      ungroup()
    
    for (m in nrow(firstarrivals)) {
      nightdat <- clustpoints %>%
        filter(
          # Case 1: Day before arrival, late at night
          (date(mt_time(.)) == as_date(toString(firstarrivals$Var2[m])) - days(1)) &
               (hour(mt_time(.)) > 21) &
               (mt_track_id(.) == firstarrivals$Var1[m]) |
            # Case 2: day of arrival, early morning
          (date(mt_time(.)) == as_date(toString(firstarrivals$Var2[m]))) &
           (hour(mt_time(.)) < 5) &
           (mt_track_id(.) == firstarrivals$Var1[m])
          )
      
      arrivaldists <- st_distance(nightdat, clustdat)
      clustertable$arrivaldist_mean[k] <- mean(arrivaldists, na.rm = T)
      clustertable$arrivaldist_med[k] <- median(arrivaldists, na.rm = T)
      
    }
    
    nightbefore <- clustpoints %>%
      filter(between(mt_time(.),
                     clustdat$firstdatetime - days(1),
                     clustdat$lastdatetime
                     ))
    
    
    

  }
  tableendtime <- Sys.time()
  logger.trace(paste0("Clustertable generation completed. Time taken: ", 
                      difftime(tableendtime, tablestarttime, units = "mins"), " mins."))
  

    clustertable %<>% left_join(tempdat, by = "xy.clust") %>%
      as.data.frame() # temporarily convert to DF to add clustercodes
    logger.trace(paste0("Clustertable is size ", object.size(clustertable) %>% format(units = "Mb")))
  
  # Fix to remove overwritten clusters from clustertable:
  logger.trace(paste0("Removing clusters that have since been overwritten in the tagdata: ", 
                     toString(
                       clustertable$xy.clust[which(clustertable$xy.clust %!in% data$xy.clust)]
                     )))
  
  # Add clustercode:
  clustertable %<>% mutate(xy.clust = ifelse(!is.na(xy.clust), paste0(clustercode, xy.clust), NA)) %>%
    mt_as_move2(time_column = "firstdatetime", track_id_column = "xy.clust")
  data %<>% mutate(xy.clust = ifelse(!is.na(xy.clust), paste0(clustercode, xy.clust), NA))
  
  # Release outputs
  saveRDS(clustertable, file = appArtifactPath("clustertable.rds")) 
  return(data)
}

