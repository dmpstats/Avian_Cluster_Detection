library('move2')
library('lubridate')
library('dplyr')
library('Gmedian')
library('sf')
library('units')
library('magrittr')



rFunction <- function(data,
                      clusterstart = NULL,
                      clusterend = NULL,
                      clusterstep = 1L, 
                      clusterwindow = 7L, 
                      d = 500L,
                      match_thresh = 175L,
                      clustexpiration = 14L, 
                      behavsystem = TRUE, 
                      clustercode = "A") {
  
  
  ## 1. Setup & Input Checks -----------------------------------------
  
  max_tm <- max(mt_time(data), na.rm = TRUE)
  min_tm <- min(mt_time(data), na.rm = TRUE) 

  # check on positive integer inputs 
  is_positive_integer(clusterstep)
  is_positive_integer(clusterwindow)
  is_positive_integer(clustexpiration)
  is_positive_integer(match_thresh)
  is_positive_integer(d)
  
  # clusterstart check
  if(is.null(clusterstart)) {
    # if clusterstart is NULL, use the first timepoint as start of clustering period
    logger.warn(paste0("`clusterstart` is set to NULL. Defaulting to clustering start timepoint ",
                       "as the earliest timestamp in input data."))
    clusterstart <- min_tm
    
  } else{
    # parse string input as POSIXct
    clusterstart <- parse_mvapps_datetime(clusterstart)
    # check start-point is within time-span covered by data
    is_clusterbound_within_dt_period(clusterstart, min_tm, max_tm)
  }
  
  # clusterend check
  if(is.null(clusterend)) {
    logger.warn(paste0("`clusterend` is set to NULL. Defaulting to end date as ",
                       "the latest timestamp in input data."))
    clusterend <- max_tm
  } else {
    # parse string input as POSIXct
    clusterend <- parse_mvapps_datetime(clusterend)
    # check end-point is within time-span covered by data
    is_clusterbound_within_dt_period(clusterend, min_tm, max_tm)
  }
  
  
  
  # check clusterstart vs clusterend
  if(clusterstart >= clusterend){
    msg <- paste0(
      "Clustering Start Date-time (`clusterstart`) '", clusterstart,"' is equal, or exceeds, the ",
      "Clustering End Date-time (`clusterend`) '", clusterend, "'. Please ensure ",
      "that Start Date-time < End Date-time.")

    logger.fatal(msg)
    stop(msg, call. = FALSE)
  }

  
  # check clusterwindow Vs clustering timespan
  clustering_timespan <- round(difftime(clusterend, clusterstart, units = "days"), 2)

  if(clusterwindow > clustering_timespan){
    logger.warn(paste0(
      "Specified `clusterwindow` of ", clusterwindow, " days exceeds the available ", 
      clustering_timespan, "-day timespan of data being clustered. Therefore, clustering will ",
      "occur only once across the whole selected clustering period, without a rolling window."
    ))
  }
  
  
  
  
  # check clusterwindow Vs clusterstep
  if(clusterstep > clusterwindow){
    msg <- paste0(
      "Specified Clustering Time-step (`clusterstep`) is larger than ",
      "Clustering Window (`clusterwindow`), which would lead to loss of data. ",
      "Please make sure time-step <= time-window.")
    
    logger.fatal(msg)
    stop(msg, call. = FALSE) 
  }
  
  
  # Check behaviour present, if needed
  if(behavsystem == TRUE & "behav" %!in% colnames(data)) {
    msg <- paste0(
      "Classified behaviour column 'behav' is not contained in the input data. ",
      "Unable to perform clustering. Please deploy the App 'Behaviour Classification",
      "for Vultures' earlier in workflow, before the current App.")
    
    logger.fatal(msg)
    stop(msg, call. = FALSE) 
  }
  
  
  # check if index column is present (from standardization app), generating it otherwise
  if("index" %!in% names(data)){
    logger.info("'index' column not present in data - adding 'index' based on 'track_id' and 'timestamp'.")
    data <- mutate(data, index = paste0(mt_track_id(data), " ", mt_time(data)))
  }
  


  # Check clustercode
  if (not_null(clustercode)) {
    logger.trace(paste0("Provided clustercode is ", clustercode))
    clustercode <- paste0(clustercode, ".")
  } else {
    logger.warn("No clustercode provided. Defaulting no clustercode. ")
    clustercode <- ""
  }
  
  
  #' Check if projection is UTM, modifying it to UTM (ESPG 32733) if not.
  #' Requirement for correct application of the dendrogram's height cut-off
  #' value (`d`) in the hierarchical clustering step below, which is assumed to
  #' be in meters
  transf_to_utm <- FALSE
  if(sf::st_crs(data)$proj != "utm"){

    logger.info(
      paste0("Input data is not projected in a UTM system. Temporarily reprojecting ",
             "locations to 'WGS 84 / UTM zone 33S' for clustering purposes."))

    # store original projection for later back-transformation
    orig_crs <- sf::st_crs(data)
    # transform to 'WGS 84 / UTM zone 33S'
    data <- data |> sf::st_transform(crs = 32733)
    # Update flag
    transf_to_utm <- TRUE
  }
  
  
  logger.info(paste0(
    "Proceeding to clustering between ", clusterstart, " and ", clusterend, 
    ", for a total of ", nrow(data), " locations from ", mt_n_tracks(data),  
    " tracks."))
  
  # Subset main data to period chosen for clustering
  eventdata <- data %>% filter(between(mt_time(data), clusterstart, clusterend))
  
  
  #' ------------------------------------------------------------------
  ##  2. Clustering Loop -----
  #' -----------------------------------------------------------------
  #'  Operate on a iterative rolling-window basis starting from clusterstart
  #'  We increment by clusterstep (in days) and re-cluster until we reach the clusterend
  #'  
  #' 'clusterdate' denotes the last timestamp in the data subset captured by the
  #' rolling time-window at each step that is being submitted to clustering

  # Set up the initial case, in which we have no roll-over of cluster data, that
  # is the timestamp from start date-time that will yield the time-interval
  # satisfying the length of the clustering time-window
  clusterdate <- clusterstart + days(clusterwindow)
  
  clusterDataDwnld <- NULL
  laststep <- FALSE
  rollingstarttime <- Sys.time()
  
  
  # Begin loop:
  while (laststep == FALSE) {
    
    # If we're on the final step, set the clusterdate equal to final day:
    if (clusterdate >= clusterend) {
      logger.trace(paste0("Current clusterdate ", as.Date(clusterdate), " is beyond final date. Assigning final date, ", as.Date(clusterend)))
      clusterdate <- clusterend
      laststep <- TRUE
    }
    
    #' -----------------------------------------------------------------------
    ### 2.1. Data Setup and Import                        ---------------------
  
    #'  Here, filter to the relevant window and import
    #'  data from the previous rolling-window (if available)
    
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
    if (behavsystem == TRUE) {
      if (sum(clusteringData$behav != "STravelling", na.rm = TRUE) < 2) { 
        logger.trace(paste0(as.Date(clusterdate), ":      Clustering complete - not enough stationary behaviour within clustering period"))
        skipToNext <- TRUE}}
    
    # If either of these conditions are met skip ahead to the next iteration
    if (skipToNext == TRUE) {
      # if not the last step, set next possible clusterdate 
      if(laststep == FALSE){
        # Find minimum following date on which we have data, using it to reset the
        # clusterdate so that it comprises the length of the clustering time-window
        clusterdate <- filter(eventdata, mt_time(eventdata) > clusterdate) %>%
          mt_time() %>%
          min() + days(clusterwindow)
        
        logger.trace(paste0("     Skipping to next date with data to cluster: ", as.Date(clusterdate)))  
      }
      
      next
    }
    
    
    #'  -------------------------------------------------------------
    ### 2.2. Perform Clustering                        -----------
    
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
    
    #' COMMENT BC: `as.dist(sf::st_distance(clusterpoints))` would produce the
    #' same result, but it's much slower. Using `dist` does however require
    #' coordinates to be in UTM, so it returns euclidean distances in meters to
    #' allow the correct application of cut-off `d`, below. This UTM check is
    #' done in the input validation section, above.
    mdist <- dist(cbind(clusterpoints$X, clusterpoints$Y))
    hc <- hclust(as.dist(mdist), method = "complete")
    
    # Add cluster information to data:
    clusterpoints$clust <- cutree(hc, h = d)
    genclustertable <- clusterpoints
    
    # Filter out clusters with only one location: (increasing to 3)
    cid <- as.vector(which(table(clusterpoints$clust) < 3))
    genclustertable %<>% filter(clust %!in% cid)
    
    if(nrow(genclustertable) == 0){
      logger.trace(paste0(as.Date(clusterdate), ":     Dismissing all clusters detected in current window as none contains more than 3 points."))
      next
    } 
    
    #' -----------------------------------------------------------------
    #### Append cluster data ---------------------------------------------
    
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
      mutate(clust_id = paste0("up", clust)) %>%
      dplyr::select(-clust)
    
    # Output location data:
    xydata <- clusterpoints %>% ungroup() %>%
      mutate(clust_id = paste0("up", clust)) %>%
      dplyr::select(-clust) %>%
      mutate(clust_id = if_else(clust_id %in% unique(clusts$clust_id), clust_id, "upNA"))
    
    
    #' ----------------------------------------------------
    ### 2.3. Identify Matching Clusters ----------------------
  
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
        arrange(clust_id) 
      

      # Generate distance matrix:
      # Rows are existing clusts, columns are new clusts
      dists <- st_distance(existingclust, newclust) |> units::set_units("m") |> units::drop_units() # ensuring dist in meters given `match_thresh` cut-off applied below
      rownames(dists) <- existingclust$clust_id
      colnames(dists) <- newclust$clust_id
      
      # Identify clusters within proximity threshold and arrange in table:
      closeClusterIndices  <- try(
        as_tibble(which(dists < match_thresh, arr.ind = T, useNames = T))
      )
      
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
    ###  2.4 Merging Clusters     ----------
    #' ----------------------------------------------
    
    #### Case 1: multiple-to-one ----------------------
    #'
    #' Multiple old clusters -> 1 new cluster
    #' In this case, we want to re-allocate the constituent locations
    #' of the new to the nearest 'old' cluster
    
    # Check to see if this is true:
    if (length(unique(matchingclustermap$updID)) != nrow(matchingclustermap)) {
      
      # browser()
      
      # Retrieve their IDs
      ids <- which(table(matchingclustermap$updID) > 1)
      for (u in 1:length(ids)) {
        
        # 1: Get merge, location, and old cluster data
        tempdat <- filter(matchingclustermap, updID == names(ids)[u])
        allocpoints <- filter(xydata, clust_id %in% tempdat$updID)
        existingclustpoints <- filter(clusterDataDwnld, clust_id %in% tempdat$existID)
        
        # Update cluster names:
        tempdat$updID <- paste0(tempdat$updID, ".", tempdat$existID)
        
        # Generate distances between new cluster locations and old cluster medians:
        pdists <- st_distance(allocpoints, existingclustpoints) %>% 
          units::drop_units() %>%
          as.data.frame() %>%
          mutate(
            # Get index of nearest cluster
            nearind = apply(., 1, which.min),
            allocclust = existingclustpoints$clust_id[nearind]
          )
        
        # Match to nearest:
        allocpoints$clust_id <- tempdat$updID[match(pdists$allocclust, tempdat$existID)]
        
        # update tracking data
        updclust_event_idx <- which(xydata$clust_id %in% names(ids)[u])
        xydata$clust_id[updclust_event_idx] <- allocpoints$clust_id
        
        
        # Update clustermap:
        matchingclustermap <- matchingclustermap %>% 
          filter(updID != names(ids)[u]) %>%
          bind_rows(tempdat)
        }
      }
    
    
    #### Case 2: one-to-multiple ----------------------
    #'
    #' 1 old cluster -> multiple new clusters
    #' In this case, we just want to keep the old cluster
    #' and we reallocate all new clusters to the same ID
    
    # Check to see if this is true:
    if (length(unique(matchingclustermap$existID)) != nrow(matchingclustermap)) {
      
      #browser()
      
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
        newclust <- filter(newclust, clust_id %!in% tempdat$updID)
        
        # 3: Deal with tracking data by assigning all to new cluster
        xydata$clust_id[xydata$clust_id %in% tempdat$updID] <- newname
      }
    }
    
    logger.trace(paste0(as.Date(clusterdate), ":       ", nrow(matchingclustermap), "  matched to existing clusters"))
    
    
    
    ####  Case 3: pre-existing clusters ----------------------
    #'
    #' Pre-existing clusters with no merges - i.e. those whose centroids are
    #' >175m away from centroids of clusters detected in current window, meaning
    #' no locations attributable to these clusters in this window

    # Retrieve the existing IDs not being merged:
    existIDs <- existingclust$clust_id[existingclust$clust_id %!in% matchingclustermap$existID]
    
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
    
    
    
    #### Case 4: new clusters ----------------------
    #'
    #' New clusters with no previous match (didn't previously exist)
    
    # Retrieve the IDs of these clusters:
    newIDs <- newclust$clust_id[newclust$clust_id %!in% matchingclustermap$updID]
    
    # We find the lowest 'free' value at which we can start their cluster IDs:
    if (length(newIDs) != 0) {
      if (!is.null(existingclust)) {
        startid <- max(clusterDataDwnld$clust_id) + 1
      } else {startid <- 1}
      
      # Generate table with these new IDs:
      newclustermap <- data.frame(existID = startid:(startid - 1 + length(newIDs)), 
                                  updID = newIDs)
    } else {
      newclustermap <- data.frame(existID = NULL, upID = NULL)
    }
    logger.trace(paste0(as.Date(clusterdate), ":       ", nrow(newclustermap), "  new clusters have no match"))
    
    
    #' ----------------------------------------------
    ### 2.4. Perform Matching ----------------------------
    
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
      rename(updID = clust_id) %>%
      left_join(updatedclustermap) %>% # Match clusters
      rename(clust_id = existID) %>%
      mutate(clust_id = as.numeric(clust_id)) %>%
      dplyr::select(-updID) %>% # Drop update column
      suppressMessages()

    # Secondly, cluster data:
    newclust <- newclust %>%
      rename(updID = clust_id) %>%
      left_join(updatedclustermap) %>%
      rename(clust_id = existID) %>%
      mutate(clust_id = as.numeric(clust_id)) %>%
      dplyr::select(-updID) %>%
      suppressMessages()

    # Combine previous and updated clusterdata:
    if (is.null(existingclust) | is.null(newclust)) {
      updatedClusters <- newclust # nothing to add in this case
    } else {
      updatedClusters <- rbind(existingclust[, c("clust_id", "geometry")], 
                               newclust) %>%
        group_by(clust_id) %>%
        summarise(geometry = st_union(geometry)) %>% # Join geometries together into multipoint
        st_centroid() # Find mean location of two centroids (we should change this later)
      
      #' COMEBCK BC: the `st_centroid` step doesn't seem to be relevant in the
      #' current clustering approach. It should probably be removed to avoid
      #' confusion with the centroid calculation via `calcGMedianSF()` performed 
      #' below
    }
    
    #browser()
    
    # Add cluster data into main location data:
    data %<>%
      filter(index %!in% xytagdata$index) %>% # remove rows whose clusters need to be updated
      #bind_rows(., xytagdata) %>%
      mt_stack(., xytagdata, .track_combine = "merge") %>% # add them back in with the new update
      arrange(mt_track_id(.), mt_time(.)) # reorder

    
    #' ----------------------------------------------------
    ### 2.5. Create TEMPORARY clustertable       ----------
    #'
    # This will contain only the essential information for
    # cluster-updating: location, time, clustID
    
    logger.trace(paste0(as.Date(clusterdate), ":     Generating temporary clustertable"))
    
    # Perform this step only if there are updates to perform on the tag data:
    if (nrow(filter(data, clust_id %in% updatedClusters$clust_id)) != 0) {
      
      # Bind necessary data on updated clusters:
      tempclusts <- data %>%
        filter(clust_id %in% updatedClusters$clust_id) %>% 
        mutate(datetime = mt_time(.))
      
      tempclusts %<>%
        group_by(clust_id) %>%
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
    ### 2.6. Return other clusters to the dataset             ---------------
    logger.trace(paste0(as.Date(clusterdate), ":     Merging back into clusterDataDwnld"))
    
    if (!is.null(clusterDataDwnld)) {
      
      clusterDataDwnld <- bind_rows(clustertable_update,
                                    filter(clusterDataDwnld, 
                                           clust_id %!in% clustertable_update$clust_id)) %>%
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
    

    
    #' ------------------------------------------------------------------
    ### 2.7. Generate output & move on rolling window        -----------
    
    # Log progress:
    logger.trace(paste0(as.Date(clusterdate), ":     COMPLETE. Total of ", nrow(clusterDataDwnld), " cumulative clusters generated"))
    if (laststep == FALSE) {
      logger.trace(paste0(as.Date(clusterdate), ": Increasing clusterdate to ", clusterdate + days(clusterstep)))
      clusterdate <- clusterdate + days(clusterstep)
    }
    
    # End of clustering loop
  }
  
  ## 3. Output processing ------------------
  
  if("clust_id" %!in% names(data)){
    logger.warn(paste0(
      "No clusters detected. Generating empty column 'clust_id', ",
      "which is a dependency of downstream cluster-related Apps"))
    
    data$clust_id <- NA
    
  }else{
    
    # identify 1-location clusters
    rem <- dplyr::count(data, clust_id) |> 
      filter(n == 1) |> 
      pull(clust_id)
    
    data <- data |> 
      mutate(
        # drop cluster id for 1 location clusters
        clust_id = ifelse(clust_id %in% rem, NA, clust_id),
        # concatenate user-specified code as the prefix of cluster_id
        clust_id = ifelse(!is.na(clust_id), paste0(clustercode, clust_id), NA)
        )  |> 
      # drop local auxiliary columns
      dplyr::select(-c("ID", "X", "Y"))
  }
  
  
  #' back-transforming coords to original projection, if they've
  #' were temporarily re-projected to UTM for clustering purposes
  if(transf_to_utm){
    logger.info(
      paste0("Reverting location coordinates to original projection of input data."))
    data <- data |> sf::st_transform(orig_crs)
  }

  
  logger.info("That's it - Finished Clustering!")
  
  return(data)
}



# Helper Functions ====================================================================

#' //////////////////////////////////////////////////////////////////////////////
# wee helpers
not_null <- Negate(is.null)
`%!in%` <- Negate(`%in%`)

#' //////////////////////////////////////////////////////////////////////////////
#' input validator for positive integer
is_positive_integer <- function(x){
  if(is.null(x)){
    stop(paste0("App input `", deparse(substitute(x)), "` specified as NULL. `", 
                deparse(substitute(x)), "` must be a positive integer."), 
         call. = FALSE)
    
  } else if(x < 1){
    stop(paste0("App input `", deparse(substitute(x)), "` must be a positive integer."), call. = FALSE)
    
  } else if(!is.integer(x)){
    stop(paste0("App input `", deparse(substitute(x)), "` must be a positive integer."), call. = FALSE)
    
  }
}


#' //////////////////////////////////////////////////////////////////////////////
#' input validator for clustering start/end timestamp being within time period covered by data
#' 
#' @param cluster_bound POSIXct, an end/start timestamp within which clustering is to be performed
#' @param min_time POSIXct, the minimum/earliest timestamp of the data under clustering
#' @param max_time POSIXct, the maximum/latest timestamp of the data under clustering
#' 
is_clusterbound_within_dt_period <- function(cluster_bound, min_time, max_time){
  
  if(length(cluster_bound) != 1) stop("cluster_bound must be of length 1", call. = FALSE)
  if(!is.POSIXct(cluster_bound)) stop("cluster_bound must be a POSIXct object", call. = FALSE)
  if(!is.POSIXct(min_time)) stop("min_time must be a POSIXct object", call. = FALSE)
  if(!is.POSIXct(max_time)) stop("max_time must be a POSIXct object", call. = FALSE)
  
  if(!dplyr::between(cluster_bound, min_time, max_time)){
    msg <- paste0(
      "App input `", deparse(substitute(cluster_bound)) ,"` is outside the range of timepoints ",
      "covered by the input data. Please provide timepoint whithin the range ['", 
      min_time, "', '", max_time, "'].")
    
    logger.fatal(msg)
    stop(msg, call. = FALSE)
  }  
  
}



#' //////////////////////////////////////////////////////////////////////////////
#' Converter of date-time string to a POSIXct object based on expected format
#' under MoveApps. Input validation also included
parse_mvapps_datetime <- function(x){
  
  msg <- paste0(
    "App input `", deparse(substitute(x)), "` must be a date-time string in ISO 8601 ",
    "format with timezone UTC, e.g. '2017-12-30 23:59:59'. Please provide a valid string."
  )
  
  if(is.POSIXct(x)){
    dt_tm <- lubridate::with_tz(x, tzone = "UTC")
    return(dt_tm)
  }
  
  
  if(!is.character(x)){
    logger.fatal(msg)
    stop(msg, call. = FALSE)
  }
  
  # parse string as POSIXct, respecting ISO 8601 format, as 
  # specified in https://docs.moveapps.org/#/appspec/current/settings/timestamp
  dt_tm <- tryCatch(
    as.POSIXct(
      x, 
      tz = "UTC", 
      tryFormats = c("%Y-%m-%dT%H:%M:%OSZ", 
                     "%Y-%m-%d %H:%M:%OS")),
    error = \(cond){
      logger.fatal(msg)
      stop(msg, call. = FALSE)  
    }
  )
  
  return(dt_tm)
}



#' //////////////////////////////////////////////////////////////////////////////
# Geometric medians
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

