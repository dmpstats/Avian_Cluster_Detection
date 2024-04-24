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
                      clustercode = "A") {
  
  
  #' --------------------------------------------------------------
  # 0. Setup & Input Checks -----------------------------------------

  # Check clustercode
  if (not_null(clustercode)) {
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
  
  
  # data %<>% mutate(xy.clust = ifelse(!is.na(xy.clust), paste0(clustercode, xy.clust), NA))
  
  # Pass cluster-appended movement data onto next MoveApp:
  return(data)
  
}


