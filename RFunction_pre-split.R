library('move2')
library('lubridate')
library('dplyr')
library('magrittr')
library('tidyr')
library('Gmedian')
library('sf')
library('stringr')
library('units')
library('pbapply')

#' TODO
#' - Add dependency on (or compute and bind locally) `nightpoint`column

# Shortened 'not in':
`%!in%` <- Negate(`%in%`)

# Function to calculate geometric medians:
calcGMedianSF <- function(data) {
  
  if (st_geometry_type(data[1,]) == "POINT") {
    med <- data %>% 
      st_coordinates()
    med <- Gmedian::Weiszfeld(st_coordinates(data))$median %>% as.data.frame() %>%
      rename(x = V1, y = V2) %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(data)) %>%
      st_geometry()
  }
  

  if (st_geometry_type(data[1,]) == "MULTIPOINT") {
    med <- data %>% 
      st_coordinates() %>%
      as.data.frame() %>%
      group_by(L1) %>%
      group_map(
        ~ Gmedian::Weiszfeld(.)$median 
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      st_as_sf(coords = colnames(.), crs = st_crs(data)) %>%
      st_geometry()
  }
  
  return(med)
  
}

rFunction <- function(data,
                      clusterstart = NULL,
                      clusterend = NULL,
                      clusterstep = 1, 
                      clusterwindow = 7, 
                      clustexpiration = 14, 
                      behavsystem = TRUE, 
                      d = 500,
                      clustercode = "") {
  
  #' --------------------------------------------------------------
  # 0. Setup & Input Checks -----------------------------------------

  # Check clustercode
  if (clustercode != "") {
    logger.trace(paste0("Provided clustercode is ", clustercode))
    clustercode <- paste0(clustercode, ".")
  } else {
    logger.warn("No clustercode provided. Defaulting no clustercode. ")
    clustercode <- ""
  }
  
  # Check suitability of start/end inputs
  if(!is.instant(as.Date(clusterstart)) | is.null(clusterstart)) {
    logger.error(paste0("Start date of clustering ", clusterstart, " is not a valid date Defaulting to start date 2 weeks before final date."))
    clusterstart <- max(mt_time(data), na.rm = TRUE) - days(14)
  }
  if(!is.instant(as.Date(clusterend)) | is.null(clusterend)) {
    logger.error(paste0("End date of clustering ", clusterend, " is not a valid date Defaulting to end date as most recent timestamp."))
    clusterend <- max(mt_time(data), na.rm = TRUE)
  }
  
  # Check behaviour present, if needed
  if(behavsystem == TRUE & "behav" %!in% colnames(data)) {
    logger.fatal("Classified behaviour column 'behav' is not contained by input data. Unable to perform clustering. Please use classification MoveApp prior to this stage in the workflow")
    stop()
  }
  
  # Check if sunset columns are present
  if ("sunset_timestamp" %in% colnames(data)) {
    suntimes <- TRUE
  } else {
    suntimes <- FALSE
  }
  
  logger.info(paste0("Clustering between ", clusterstart, " and ", clusterend))
  
  # Filter to relevant time window for overall clustering:
  clusterstart %<>% as_date()
  clusterend %<>% as_date()
  eventdata <- data %>%
    filter(between(mt_time(data), clusterstart - days(clusterwindow), clusterend + days(1)))
  
  
  #' ------------------------------------------------------------------
  # Clustering Loop
  #'  Operate on a rolling-window basis starting from clusterstart
  #'  We increment by clusterstep (in days) and re-cluster until we reach the clusterend
  #'  
  #'  'clusterdate' parameter will define the final day of data we are currently clustering up to

  
  # Set up the initial case, in which we have no roll-over of cluster data:
  clusterdate <- floor_date(clusterstart, unit = "days")
  clusterDataDwnld <- NULL
  tempclustertable <- NULL
  laststep <- FALSE
  rollingstarttime <- Sys.time()
  
  
  # Begin loop:
  while (laststep == FALSE) {
    
    
    # If we're on the final step, set the clusterdate equal to final day:
    if (clusterdate >= floor_date(clusterend, unit = "days")) {
      logger.trace(paste0("Current clusterdate ", as.Date(clusterdate), " is beyond final date. Assigning final date, ", as.Date(clusterend)))
      clusterdate <- floor_date(clusterend, unit = "days")
      laststep <- TRUE
    }
    
    
    #' -----------------------------------------------------------------------
    #' 1. Data Setup and Import ----------------------------------------------
    #'  Here, filter to the relevant window and import
    #'  data from the previous rolling-window (if available)
    
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
    skipToNext <- FALSE # Determines whether to jump over this clusterstep
    
    # Check that there is enough tracking data to cluster, and skip if not:
    if (nrow(clusteringData) < 3) {
      logger.trace(paste0(as.Date(clusterdate), ":      Clustering complete - not enough data within clustering period"))
      skipToNext <- TRUE}
    
    # And if using the behavioural classification system, check there is enough non-travelling behaviour:
    #browser()
    if (behavsystem == TRUE) {
      if (sum(clusteringData$behav != "STravelling", na.rm = TRUE) < 2) { 
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
    
    
    
    #' ----------------------------------------------------------------------
    # 2. Perform Clustering ----------------------------------------------
    #' Now that we know there is enough data, generate the new clusters

    
    # If using behavioural classification, filter to stationary behaviours:
    if (behavsystem == TRUE) {
      clusterpoints <- filter(clusteringData, behav %!in% c("STravelling", "Unknown"))
    } else {
      clusterpoints <- clusteringData
      } 
    
    logger.trace(paste0(as.Date(clusterdate), ":     Creating new clusters for ", nrow(clusterpoints), " locations"))
    
    # Generate clusters:
    clusterpoints %<>% mutate(
      ID = 1:nrow(.),
      X = st_coordinates(.)[, 1],
      Y = st_coordinates(.)[, 2]
    )
    
    # Build distance matrix (no SF method available for the hclust step):
    mdist <- dist(cbind(clusterpoints$X, clusterpoints$Y))
    hc <- hclust(as.dist(mdist), method = "complete")
    
    # Add cluster information to data:
    clusterpoints$clust <- cutree(hc, h = d)
    genclustertable <- clusterpoints
    
    # Filter out clusters with only one location: (increasing to 3)
    cid <- as.vector(which(table(clusterpoints$clust) < 3))
    genclustertable %<>% filter(clust %!in% cid)
    
    
    #' -----------------------------------------------------------------
    ## Append cluster data ---------------------------------------------
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
      mutate(xy.clust = if_else(xy.clust %in% unique(clusts$xy.clust), xy.clust, "upNA"))
    
    
    #' ----------------------------------------------------
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
    
    
    #' -----------------------------------------------
    ## Merging Clusters: Case 1 ----------------------
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
    
    #' -----------------------------------------------
    ## Merging Clusters: Case 2 ----------------------
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
    
    
    #' -----------------------------------------------
    ## Merging Clusters: Case 3 ----------------------
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
    
    
    #' -----------------------------------------------
    ## Merging Clusters: Case 4 ----------------------
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
    
    
    #' ----------------------------------------------
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
      bind_rows(., xytagdata) %>%
      #mt_stack(., xytagdata, .track_combine = "merge") %>% # add them back in with the new update
      arrange(mt_track_id(.), mt_time(.)) # reorder

    
    #' ----------------------------------------------------
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
    
    
    #' -----------------------------------------------------------------------
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
    

    
    #' ----------------------------------------------------
    # 7. Generate output & move on rolling window -----------
    
    # Log progress:
    logger.trace(paste0(as.Date(clusterdate), ":     COMPLETE. Total of ", nrow(clusterDataDwnld), " cumulative clusters generated"))
    if (laststep == FALSE) {
      logger.trace(paste0(as.Date(clusterdate), ": Increasing clusterdate to ", clusterdate + days(clusterstep)))
      clusterdate <- clusterdate + days(clusterstep)
    }
    
    # End of clustering loop
  }
  
  rollingendtime <- Sys.time()
  logger.trace(paste0("CLUSTERING COMPLETE. Time taken: ", 
               difftime(rollingendtime, rollingstarttime, units = "mins"),  " mins. Generating clustertable for all clusters"))
  
  
  #' ----------------------------------------------------
  # 8. Generate final clustertable ----------------------
  # This will contain output cluster information
  
  tablestarttime <- Sys.time()
  clustertable <- data %>%
    filter(!is.na(xy.clust)) %>%
    mutate(ID = mt_track_id(.), 
           timestamp = mt_time(.),
           timediff_hrs = mt_time_lags(.) %>% 
             units::set_units("minutes") %>% 
             units::drop_units()/60,) %>%
    #as.data.frame() %>%
    group_by(xy.clust, ID) %>%
    summarise(
      geometry = st_combine(geometry),
      all_cluster_points = geometry,
      
      # Time calculations:
      firstdatetime = min(timestamp),
      lastdatetime = max(timestamp),
      days = length(unique(lubridate::date(timestamp))),
      totalduration = (ceiling_date(lastdatetime, unit = "days") - floor_date(firstdatetime, unit = "days")) %>% as.integer(),
      daysempty = as.numeric(totalduration - days),
      daysemptyprop = daysempty / totalduration,
      no_nightpoints = all(nightpoint == 0),
      
      # Behavioural data:
      Total = n(),
      SFeeding = sum(behav == "SFeeding"),
      SRoosting = sum(behav == "SRoosting"),
      SResting = sum(behav == "SResting"),
      MedianHourFeed = median(hour[behav == "SFeeding"], na.rm = T),
      MedianHourDay = median(hour[hour > 5 & hour < 17], na.rm = T),
      
      # Accelerometer calculations
      # med_var_x = case_when(
      #   "var_acc_x" %in% colnames(data) ~ median(var_acc_x, na.rm = T),
      #   TRUE ~ NA
      # ),
      # med_var_y = case_when(
      #   "var_acc_y" %in% colnames(data) ~ median(var_acc_y, na.rm = T),
      #   TRUE ~ NA
      # ),
      # med_var_z = case_when(
      #   "var_acc_z" %in% colnames(data) ~ median(var_acc_z, na.rm = T),
      #   TRUE ~ NA
      # ),
      # sd_var_x = case_when(
      #   "var_acc_x" %in% colnames(data) ~ sd(var_acc_x, na.rm = T),
      #   TRUE ~ NA
      # ),
      # sd_var_y = case_when(
      #   "var_acc_y" %in% colnames(data) ~ sd(var_acc_y, na.rm = T),
      #   TRUE ~ NA
      # ),
      # sd_var_z = case_when(
      #   "var_acc_z" %in% colnames(data) ~ sd(var_acc_z, na.rm = T),
      #   TRUE ~ NA
      # )
      
      .groups = "keep"
    ) %>%
    st_as_sf(crs = st_crs(data))
  


  rm(list = c("clusterDataDwnld", "clusteringData", "clusterpoints", "eventdata", "genclustertable", "xydata", "xytagdata", "mdist")) # clear some storage for next operations
  clustertable %<>% mt_as_move2(time_column = "firstdatetime", track_id_column = "xy.clust") # switch to move2 for ease
  
  # Generate distance data:
  logger.trace(paste0("Generating distance data for all clusters. This may run slowly"))
  

  ## CLUSTERTABLE DATA PREP ---------------------------------------------------
  
  
  # Reformat data for clustertable calculations
  
  # Update geometry column to be geometric medians 
  # These are per-bird, not per-cluster
  clustertable %<>% st_set_geometry(
    calcGMedianSF(.) %>%
      st_geometry()
  )
  
  # Add second geometry column for cluster-wide centroids
  # These are per-cluster, not per-bird
  wholeclusts <- clustertable %>%
    group_by(xy.clust) %>%
    summarise(
      geometry = st_combine(geometry),
      .groups = "keep"
    ) %>%
    st_set_geometry(
      calcGMedianSF(.) %>%
        st_geometry() 
    ) %>%
    rename(wholeclust_geometry = geometry) %>%
    as.data.frame()
  clustertable %<>% left_join(wholeclusts, by = "xy.clust")
  
  # Create matrix data for filtering speed
  mat.data <- data %>%
    mutate(ID = mt_track_id(.),
           timestamp = mt_time(.)) %>%
    as.data.frame() %>%
    dplyr::select(any_of(c("ID",
                           "geometry",
                           "timestamp",
                           "hour",
                           "dist_m",
                           "behav",
                           "sunrise_timestamp",
                           "sunset_timestamp",
                           "timediff_hrs",
                           "xy.clust",
                           "nightpoint"))) %>%
    st_as_sf(crs = st_crs(data))
  

  ## FURTHER CLUSTERTABLE ATTRIBUTES -------------------------------------------
  
  ### a. Calculate accelerometer data if ACC is available
  if ("var_acc_x" %in% colnames(data)) {
    
    logger.trace("   Accelerometer columns identified. Calculating ACC summaries")
    
    clustertable_acc <- data %>%
      filter(!is.na(xy.clust)) %>%
      mutate(ID = mt_track_id(.), 
             timestamp = mt_time(.),
             timediff_hrs = mt_time_lags(.) %>% 
               units::set_units("minutes") %>% 
               units::drop_units()/60,) %>%
      #as.data.frame() %>%
      group_by(xy.clust, ID) %>%
      summarise(
        # Accelerometer calculations
        med_var_x = median(var_acc_x, na.rm = T),
        med_var_y = median(var_acc_y, na.rm = T),
        med_var_z = median(var_acc_z, na.rm = T),
        
        sd_var_x = sd(var_acc_x, na.rm = T),
        sd_var_y = sd(var_acc_y, na.rm = T),
        sd_var_z = sd(var_acc_z, na.rm = T),
        
        .groups = "keep"
      )
  }
  
  ### b. Time-at-Carcass Calculations -----------------------------------------
  logger.trace("   Generating [a. Time-at-Carcass data]")
  
  timeAtCarcTab <- function(clustdat) {
    
    carctime <- clustdat %>%
      filter(!is.na(xy.clust)) %>%
      as.data.frame() %>%
      group_by(xy.clust, ID, date(timestamp)) %>%
      # Remove all final points to prevent overnight locations messing things up:
      filter(row_number() != n()) %>%
      summarise(count = n(),
                time_spent = sum(timediff_hrs),
                time_spent_day = sum(timediff_hrs * !nightpoint),
                .groups = "keep"
                ) %>%
      
      group_by(xy.clust, ID) %>%
      summarise(
        meanvisit_time = mean(time_spent, na.rm = T),
        meanvisit_time_day = mean(time_spent_day, na.rm = T),
        .groups = "keep") 
    return(carctime)
  }
  carctime <- timeAtCarcTab(mat.data)
  
  
  ### c. Revisitation Calculations --------------------------------------------
  logger.trace("   Generating [b. Revisitation data]")
  
  
  revisitTab <- function(clustdat, clustertable) {
    
    revisit_calc <- function(row) {
      
      clust <- row$xy.clust
      bird <- row$ID
      firstdate <- as_datetime(row$firstdatetime) 
      lastdate <- as_datetime(row$lastdatetime)
      
      nearpoints <- clustdat %>% 
        filter(
          between(timestamp, firstdate, lastdate),
          ID == bird
        ) %>%
        mutate(incluster = ifelse(xy.clust == clust & !is.na(xy.clust), 1, 0),
               indaycluster = ifelse(incluster == 1 & nightpoint == 0, 1, 0)
        ) %>%
        group_by(date(timestamp)) %>%
        summarise(
          visits = sum(rle(incluster)$values),
          dayvisits = sum(rle(indaycluster)$values),
          dayvisits = pmin(dayvisits, visits), # fix case where dayvisits > visits
          .groups = "keep" 
        ) %>%
        ungroup() %>%
        summarise(
          meanvisits = mean(visits, na.rm = T),
          meandayvisits = mean(dayvisits, na.rm = T),
          .groups = "keep" 
        )
      
      return(c(
        nearpoints$meanvisits,
        nearpoints$meandayvisits
      ))
    }
    
    revdat <- pbapply(clustertable, 1, revisit_calc) %>% 
      t() %>%
      as.data.frame() %>%
      rename(meanvisits = V1, meanvisits_daytime = V2)
    outclusts <- cbind(clustertable, revdat) %>%
      as.data.frame() %>%
      dplyr::select(c("xy.clust", "ID", "meanvisits", "meanvisits_daytime"))
    return(outclusts)
  }
  
  # FOLLOWING LINE FAILS
  revisits <- revisitTab(mat.data, clustertable) 
  
  
  
  ### d. Night-Distance Calculations ------------------------------------------
  
  
  
  logger.trace("   Generating [c. Night-distance data]")

  nightTab <- function(clustdat, clustertable) {
    
    # ALTERNATIVE METHOD TEST
    # Firstly, filter to all night locations
    nightpts <- clustdat %>% 
      dplyr::select(-xy.clust) %>%
      as.data.frame() %>%
      filter(nightpoint == 1) %>%
      mutate(date = case_when(
        # Fix night locations being 'split' by midnight:
        # If before midnight, associate it with that same date
        hour(timestamp) > 12 ~ date(timestamp),
        # But if in the morning, associate it with the night before
        hour(timestamp) < 12 ~ date(timestamp) - days(1)
      ))
    
    # Generate second table:
    # clust-bird-date, one entry for each clust visited by a bird on each date
    clustdays <- clustdat %>%
      as.data.frame() %>%
      filter(!is.na(xy.clust)) %>%
      group_by(xy.clust, ID, date(timestamp)) %>%
      summarise() %>%
      rename(date = `date(timestamp)`) %>%
      
      # Bind clust centroid data:
      left_join(
        clustertable %>%
          as.data.frame() %>%
          dplyr::select(c("xy.clust", "wholeclust_geometry")) %>%
          .[!duplicated(.),] %>%
          st_as_sf(crs = st_crs(data)),
        by = "xy.clust", relationship = "many-to-many")  %>%
      .[!duplicated(.),] 
    
    # The following table contains all night locations 
    # and has matched them to a cluster visited by a bird on that same day.
    # Where more than 1 cluster is visited by a bird within a day, the
    # night location has been duplicated (once for each cluster) so that it can be grouped more than once.
    night_table <- left_join(nightpts, clustdays, by = c("ID", "date"), relationship = "many-to-many") %>% 
      filter(!is.na(xy.clust)) 
    
    # Now we introduce a distance column:
    dists <- pbapply::pbmapply(st_distance, night_table$geometry, night_table$wholeclust_geometry)
    nightdists <- cbind(night_table, dists)
    
    # Testing a new variable: proportion of nearby night points 
    # This is the proportion of night points on the same day as this cluster 
    # within 250m
    nearnights <- nightdists %>%
      group_by(xy.clust, ID) %>%
      summarise(near_night_prop = sum(dists < 250) / n())
    
    
    # Group by clust-ID-date and summarise, taking median first
    nightdists_by_day <- nightdists %>%
      group_by(ID, date, xy.clust) %>%
      summarise(nightdist = median(dists, na.rm = T), .groups = "keep")
    
    # Finally, take mean per bird across several days
    nightdists_bird_clust <- nightdists_by_day %>%
      group_by(xy.clust, ID) %>%
      summarise(nightdist_med = mean(nightdist, na.rm = T), .groups = "keep") %>%
      left_join(nearnights, by = c("xy.clust", "ID"))
    
    return(nightdists_bird_clust)
    
  }
  
  nightdists <- nightTab(mat.data, clustertable)
 

  
  ### e. Arrival-Distance Calculations ---------------------------------------
  logger.trace("   Generating [d. Arrival-Distance data]")
  
  arrivalTab <- function(clustdat, clustertable) {

    clustarrivals <- clustdat %>%
      st_drop_geometry() %>%
      as.data.frame() %>%
      group_by(xy.clust, ID) %>%
      summarise(day_of_arrival = min(date(timestamp)), .groups = "keep") %>%
      rename(birdID = ID) %>%
      left_join(
        clustertable %>%
          as.data.frame() %>%
          dplyr::select(c("xy.clust", "wholeclust_geometry")) %>%
          .[!duplicated(.),] %>%
          st_as_sf(crs = st_crs(data)),
        by = "xy.clust", relationship = "many-to-many") %>%
      .[!duplicated(.),] %>%
      st_set_geometry(.$wholeclust_geometry) %>%
      dplyr::select(-wholeclust_geometry)
    
    genClustDists <- function(row) {
      
      clust <- row$xy.clust
      bird <- row$birdID
      arrivaldate <- row$day_of_arrival %>% as_date()
      clustgeometry <- row$wholeclust_geometry %>%
        st_sfc(crs = st_crs(clustdat))
      
      arrivaldist <- clustdat %>%
        filter(between(
          timestamp,
          arrivaldate - hours(12), 
          arrivaldate + hours(12)
        ),
        nightpoint == 1,
        ID == bird) %>%
        mutate(
          dist = st_distance(geometry, clustgeometry)
        ) %>%
        .$dist %>%
        mean(na.rm = T)
      
      return(arrivaldist)
    }
    arrivaldat <- pbapply::pbapply(clustarrivals, 1, genClustDists)
    
    outdat <- clustarrivals %>%
      ungroup() %>%
      mutate(arrivaldists = arrivaldat) %>%
      st_drop_geometry() %>%
      rename(ID = birdID) %>%
      dplyr::select(-c("day_of_arrival"))
    
    return(outdat)
  }
  arrivaldists <- arrivalTab(mat.data, clustertable)
  
  
  ### f. Nearest-Tag Calculations --------------------------------------------
  logger.trace("   Generating [d. Nearest-Tag data]")
  
  nearBirdsTab <- function(clustdat, clustertable) {
    
    # List all clusters
    topclusters <- clustertable %>% 
      ungroup() %>%
      st_set_geometry(.$wholeclust_geometry) %>%
      group_by(xy.clust) %>%
      summarise(
        birds = paste(unique(ID), collapse = ", "),
        firstdatetime = min(firstdatetime, na.rm = T),
        lastdatetime = max(lastdatetime, na.rm = T),
        nbirds = length(unique(ID)),
        .groups = "keep" 
      )
    
    distvals <- function(row) {
      
      clustID <- row$xy.clust
      firstdatetime <- row$firstdatetime %>% as_datetime()
      lastdatetime <- row$lastdatetime %>% as_datetime()
      clustgeometry <- row$geometry %>%
        st_sfc(crs = st_crs(clustdat))
      
      nearpoints <- clustdat %>%
        filter(between(timestamp, firstdatetime - days(14), lastdatetime)) 
      
      activetags <- nearpoints %>%
        filter(between(timestamp, firstdatetime, lastdatetime)) %>%
        .$ID %>%
        unique()
      
      nearpoints2 <- nearpoints %>% 
        mutate(
          dist = st_distance(geometry, clustgeometry)
        ) %>% 
        filter(ID %in% activetags,
               ID %!in% row$birds
        )
      
      mindist_m <- min(nearpoints2$dist, na.rm = T) %>%
        units::set_units("metres") %>%
        units::drop_units()
      within_25k <- length(unique(nearpoints2$ID[nearpoints2$dist < units::set_units(25, "kilometres")]))
      within_50k <- length(unique(nearpoints2$ID[nearpoints2$dist < units::set_units(50, "kilometres")]))
      
      return(c(mindist_m, within_25k, within_50k))
    }
    
    birddat <- pbapply::pbapply(topclusters, 1, distvals) %>%
      t() %>%
      as.data.frame() %>%
      rename(
        mindist_m = V1,
        within_25k = V2, 
        within_50k = V3
      )
    
    topclusters_final <- cbind(topclusters, birddat) %>%
      as.data.frame() %>%
      dplyr::select(c("xy.clust", "birds", "mindist_m", "within_25k", "within_50k"))
    return(topclusters_final)
  }
  nearbirds <- nearBirdsTab(mat.data, clustertable)
  
  ### g. Stack outputs into final clustertable ---------------------------------
  logger.trace("   Merging [a-f] into primary clustertable")
  
  clustertable %<>%
    left_join(carctime, by = c("xy.clust", "ID")) %>%
    left_join(revisits, by = c("xy.clust", "ID")) %>%
    left_join(nightdists, by = c("xy.clust", "ID")) %>%
    left_join(arrivaldists, by = c("xy.clust", "ID")) %>%
    left_join(nearbirds, by = "xy.clust") %>%
    mt_as_move2(time_column = "firstdatetime", track_id_column = "xy.clust")
  
  if ("var_acc_x" %in% colnames(data)) {
    clustertable %<>% 
      left_join(st_drop_geometry(clustertable_acc) %>% 
                  dplyr::select(c("xy.clust", "med_var_x", "med_var_y", "med_var_z", "sd_var_x", "sd_var_y", "sd_var_z")), 
                by = c("xy.clust", "ID"))
  }
  
  
  ## 9. Finalise Outputs ---------------------------------------------------------
  
  # Finally, remove 1-location clusters from tagdata and clustertable
  rem <- clustertable$xy.clust[clustertable$Total == 1]
  data %<>% mutate(
    xy.clust = case_when(
      xy.clust %in% rem ~ NA,
      TRUE ~ xy.clust
    )
  ) %>%
    dplyr::select(-c("ID", "X", "Y"))
  clustertable %<>% filter(xy.clust %!in% rem)
  
  
  # Looping complete - log time and release outputs
  #tableendtime <- Sys.time()
  #logger.trace(paste0("Clustertable generation completed. Time taken: ", 
  #                    difftime(tableendtime, tablestarttime, units = "mins"), " mins."))

    clustertable %<>% as.data.frame() %>% st_as_sf() # temporarily convert to DF to add clustercodes (Move object creates errors)
    logger.trace(paste0("Clustertable is size ", object.size(clustertable) %>% format(units = "Mb")))
  
  # Fix to remove overwritten clusters from clustertable:
  logger.trace(paste0("Removing clusters that have since been overwritten in the tagdata: ", 
                     toString(
                       clustertable$xy.clust[which(clustertable$xy.clust %!in% data$xy.clust)]
                     )))
  
  # Add clustercode:
  clustertable %<>% mutate(xy.clust = ifelse(!is.na(xy.clust), paste0(clustercode, xy.clust), NA)) %>%
    mt_as_move2(time_column = "firstdatetime", track_id_column = "xy.clust")  %>%
    mt_as_track_attribute(c("birds", "mindist_m", "within_25k", "within_50k")) %>%
    st_set_geometry(.$wholeclust_geometry) %>% # change geometry to be whole-cluster centroid (the same location will be shared by several rows)
    dplyr::select(-"wholeclust_geometry")
  data %<>% mutate(xy.clust = ifelse(!is.na(xy.clust), paste0(clustercode, xy.clust), NA))
  
  
  # Add clustertable as an attribute of the main dataset
  attr(data, "cluster_tbl") <- clustertable
  
  # Release outputs
  saveRDS(clustertable, file = appArtifactPath("clustertable.rds")) 
  return(data)
  
}

