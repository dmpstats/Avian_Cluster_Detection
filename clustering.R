# ~~~~~~~~~~~~~~~~~~
# clustering.R
# ~~~~~~~~~~~~~~~~~~
# Jul 2023

# ~~~~~~~~~~~~~
# Reworking for MoveApps / sf
# ~~~~~~~~~~~~~

clustering <- function(datmodsub, clusterstartdate, clusterenddate, clusterstep = 1, clusterwindow = 7, clustexpiration = 14) {
  
  
  # Filtering & initial clustering ------------------------------------------------------------
  
  
  # We initially have no previous cluster data
  clusterDataDwnld <- NULL
  fullclustertable <- NULL
  
  # Filter to data being clustered
  eventdata <- datmodsub %>%
    ungroup() %>%
    filter(between(mt_time(datmodsub), clusterstartdate - days(clusterwindow), clusterenddate + days(1)))
  
  
  # Set initial day to cluster:
  clusterdate <- ceiling_date(clusterstartdate, unit = "days")

  
  # Loop over the days within our clustering period:
  while(clusterdate <= ceiling_date(clusterenddate, unit = "days")) {
    
    # Unless we're on the first step, download cluster data
    if(!is.null(clusterDataDwnld)) {
      clusterDataDwnld <- fullclustertable
    }
    
    logger.debug(paste0(as.Date(clusterdate), ": Beginning clustering"))
    
    # Sample to relevant 7 days
    clusteringData <- eventdata %>%
      filter(between(mt_time(eventdata), clusterdate - days(clusterwindow), clusterdate))
    
    # Check data present:
    if(nrow(clusteringData) < 3) {
      logger.warn(paste0(as.Date(clusterdate), ": Clustering complete - not enough data within clustering period"))
      
      # Find minimum following date on which we have data
      clusterdate <- filter(eventdata, mt_time(eventdata) > clusterdate) %>%
        mt_time() %>%
        min() + days(clusterwindow)
      logger.warn(paste0("Skipping to next date with data to cluster: ", as.Date(clusterdate)))
      next
    }
    
    # We need to check the dataset isn't entirely travelling, as there wouldn't
    # be enough points to cluster:
    
    if(
      sum(clusteringData$behav != 'STravelling') < 2
    ) {
      logger.warn(paste0(as.Date(clusterdate), ": Not enough non-travelling behaviour for clustering (need >2 tracks). Trying next possible clusterdate"))
      clusterdate <- filter(eventdata, mt_time(eventdata) > clusterdate) %>%
        mt_time() %>%
        min() + days(clusterwindow)
      logger.warn(paste0("Skipping to clusterdate ", as.Date(clusterdate)))
      next
    }
    

    # Create new clusters
    logger.trace(paste0(as.Date(clusterdate), ": Creating new clusters"))
    clusters_new <- makeEventClusters(clusteringData, d = 500)
    
    

    
    # Matching new to old clusters ----------------------------------------------------
    
    logger.trace(paste0(as.Date(clusterdate), ": Creating matchingclustermap"))
    # Retrieve old clusters
    
    
    logger.trace(paste0(as.Date(clusterdate), ": Updating tagdata"))
    
    
    
    if (is.null(clusterDataDwnld)) {
      existingclust <- NULL # no previous clusters
    } else {
      existingclust <- clusterDataDwnld %>%
        filter(lastdatetime > clusterdate - days(clustexpiration)) %>%
        dplyr::select(x, y, xy.clust)
      logger.trace(paste0(as.Date(clusterdate), 
                          ": ", 
                          nrow(existingclust),
                          " existing clusters found in filtering period"))
    }
    
    # Retrieve new clusters
    newclust <- clusters_new$clusts 
    logger.trace(paste0(as.Date(clusterdate), 
                        ": ", 
                        nrow(newclust),
                        " clusters generated"))
    
    if (is.null(existingclust) | nrow(newclust) == 0) {
      matchingclustermap <- data.frame(existID = NULL, updID = NULL)
    } else {
      
      
      # find distances between existing and new clusters
      clustdist <- as.matrix(dist(rbind(existingclust[,1:2], newclust[,1:2]))) # x and y are utms
      diag(clustdist) <- NA
      # rows are existing clusters, columns are new clusters
      clustdist <- as.matrix(clustdist[1:nrow(existingclust), (nrow(existingclust)+1): (nrow(existingclust) + nrow(newclust))])
      
      #find near ones (200m), i.e. new clusters that match with an exsiting one. 
      closeClusterIndices <- try(as_tibble(which(clustdist<(175), arr.ind = TRUE, useNames = TRUE)) %>%
                                   rename(loc1_rowIndex = row, loc2_rowIndex = col), silent = TRUE)

      # map updIDs in newclusts to xy.clust ID in existing
      # MATCHES
      matchingclustermap <- data.frame(existID = existingclust$xy.clust[closeClusterIndices$loc1_rowIndex], 
                                       updID = newclust$xy.clust[closeClusterIndices$loc2_rowIndex])
      # Is this ^ the right way round?
      # CC: I've found some cases where the indexes seem to be swapped. Unsure what's causing it 
      
      # Short term fix
      # Removing any matched clusters with 'NA' as an index (can't find the cause right now) 
      if (nrow(matchingclustermap) != 0) {
        for (i in 1:nrow(matchingclustermap)) {
          if (is.na(matchingclustermap$existID[[i]]) | is.na(matchingclustermap$updID)[[i]])
          {matchingclustermap <- matchingclustermap[-i,]}
        }
      }
      
      logger.trace(paste0(as.Date(clusterdate), ": ", nrow(matchingclustermap), " clusters matched"))
      
       
    }
    
    # Merge clusters ---------------------------------------------------------------------
    # We want to keep existing clusters and allocate new points to new clusters
    
    
    
    
    logger.trace(paste0(as.Date(clusterdate), ": Merging with matched clusters"))
    if (length(unique(matchingclustermap$updID)) != dim(matchingclustermap)[1]) {
      
      # Find which tags are used more than once
      ids <- which(table(matchingclustermap$updID) > 1)
      
      for (u in 1:length(ids)) {
        
        tempdat <- filter(matchingclustermap, updID == names(ids)[u]) # selects IDs of clusters to be merged
        allocpoints <- filter(clusters_new$tagdata, xy.clust %in% tempdat$updID) # filters to tags associated with these points in NEW data
        existingclustpoints <- filter(clusterDataDwnld, xy.clust %in% tempdat$existID) # gets data for OLD clusters to be merged
        
        tempdat$updID <- paste0(tempdat$updID, ".", tempdat$existID) # tags these with previous cluster name
        
        # Allocate points to the nearest of existing IDs
        existclustlocs <- existingclustpoints %>% as.data.frame() %>% dplyr::select(x, y)
        allocpointsloc <- allocpoints %>% as.data.frame() %>% dplyr::select(x, y)
        pdists <- as.matrix(dist(rbind(existclustlocs, allocpointsloc)))
        pdists <- pdists[(nrow(existingclustpoints) + 1):nrow(pdists), 1:nrow(existingclustpoints)] 
        mindists <- apply(pdists, 1, min)
        
        for (d in 1:length(mindists)) {
          cid <- which(pdists[d,] == mindists[d])
          allocpoints$xy.clust[d] <- tempdat$updID[cid]
        }
        
        matchingclustermap <- filter(matchingclustermap, updID != names(ids)[u]) %>%
          bind_rows(., tempdat)
        clusters_new$tagdata <- filter(clusters_new$tagdata, xy.clust %!in% names(ids)[u]) %>%
          mt_stack(., allocpoints, .track_combine = "merge", .track_id_repair = "unique") # error here
        #  bind_rows(., allocpoints) # error here
        newclust <- filter(newclust, xy.clust != names(ids)[u])
        
      }
    }
    
    # When two new clusters match one existing cluster
    # edit to retain only the existing cluster
    
    if(length(unique(matchingclustermap$existID)) != dim(matchingclustermap)[1]) {
      
      ids <- which(table(matchingclustermap$existID) > 1)
      
      for (u in 1:length(ids)) {
        
        # find points allocated to updID
        tempdat <- filter(matchingclustermap, existID == names(ids)[u])
        
        allocpoints <- filter(clusters_new$tagdata, xy.clust %in% tempdat$updID)
        existclustpoints <- filter(clusterDataDwnld, xy.clust %in% tempdat$existID)
        
        newname <- paste0(tempdat$updID[1], ".", tempdat$existID)
        
        newclust <- filter(newclust, xy.clust %!in% tempdat$updID)
        
        allocpoints <- allocpoints %>% mutate(xy.clust = newname[1])
        clusters_new$tagdata <- filter(clusters_new$tagdata, xy.clust %!in% tempdat$updID) %>%
          mt_stack(., allocpoints, .track_combine = "merge") 
        
        tempdat$updID <- newname
        matchingclustermap <- filter(matchingclustermap, existID != names(ids)[u]) %>%
          bind_rows(., tempdat %>% slice(1))
        
        }
      
    }
    
    
    
    
    # EXISTING clusters
    existIDs <- existingclust$xy.clust[existingclust$xy.clust %!in% matchingclustermap$existID]
    
    if (length(existIDs) != 0) {
      existclustermap = data.frame(existID = existIDs,
                                   updID = NA)
    } else {
      existclustermap <- data.frame(existID = NULL, upID = NULL)
    }
    
    # NEW clusters
    newIDs <- newclust$xy.clust[newclust$xy.clust %!in% matchingclustermap$updID]
    
    # Find ID to start new clusters at
    if (length(newIDs) != 0) {
      if(!is.null(existingclust)) {
        startid <- max(clusterDataDwnld$xy.clust) + 1
      } else {startid <- 1}
      
      newclustermap <- data.frame(existID = startid:(startid - 1 + length(newIDs)), 
                                 updID = newIDs)
    } else {
      newclustermap <- data.frame(existID = NULL, upID = NULL)
      # newclustermap <- NULL
    }
    
    logger.trace(paste0(as.Date(clusterdate), ": Creating updatedclustermap"))
    updatedclustermap <- rbind(matchingclustermap, newclustermap, existclustermap, c(NA, "upNA"))
    colnames(updatedclustermap) <- c("existID", "updID") # fix for 0-cluster situation
    
    
    
    
    
    # Update the xytagdata (7-day data) ----------------------------------------------

    
    
    #f(!is.null(newclustermap) | !is.null(existclustermap)) {
     
      xytagdata <- clusters_new$tagdata %>% rename(updID = xy.clust) %>%
        left_join(., updatedclustermap) %>%
        rename(xy.clust = existID) %>%
        mutate(xy.clust = as.numeric(xy.clust)) %>%
        dplyr::select(-updID) %>%
        suppressMessages()
      
      # Change the new cluster IDs
      newclust <- newclust %>% rename(updID = xy.clust) %>%
        left_join(., updatedclustermap) %>%
        rename(xy.clust = existID) %>%
        mutate(xy.clust = as.numeric(xy.clust)) %>%
        dplyr::select(-updID) %>%
        suppressMessages()
      
      
      # Add back to 14-day clustertable
      if (is.null(existingclust) | is.null(newclust)) {
        updatedClusters <- newclust
      } else {
        updatedClusters <- full_join(existingclust, newclust) %>%
          group_by(xy.clust) %>%
          summarise(x = mean(x), y = mean(y)) %>%
          filter(xy.clust %in% newclust$xy.clust) %>%
          suppressMessages()
        
      } 
      
    #} else {
    #  
    #  # If no new clusters are created nor exist already, simply NA the full tagdata clusters
    #  xytagdata <- clusters_new$tagdata %>% 
    #    mutate(xy.clust = NA)
    #  updatedClusters <- NULL
    #  
    #}

    
    # Update main data
    datmodsub <- datmodsub %>%
      filter(index %!in% xytagdata$index) %>%
      bind_rows(., xytagdata) %>%
      arrange(mt_track_id(.), mt_time(.))
    
    
    
    # Create cluster table -------------------------------------------------------
    
    logger.trace(paste0(as.Date(clusterdate), ": Generating clustertable"))
    
    # We don't want to run this if there are no clusters present at this time
    if (nrow(filter(datmodsub, xy.clust %in% updatedClusters$xy.clust)) != 0) {
      
      clustertable_update <- makeclustertable(filter(datmodsub, xy.clust %in% updatedClusters$xy.clust),
                                              updatedClusters) #%>% 
      #      as.data.frame()
    } else {
      clustertable_update <- NULL
      }
    
    # Don't run this next section if there are no new clusters
  #  if(nrow(updatedClusters) != 0 & !is.null(updatedClusters)) {
    if(!is.null(updatedClusters)) {
      
      # Calculate revisits ---------------------------------------------------------
      
      
      
      logger.trace(paste0(as.Date(clusterdate), ": Calculating revisits"))
      cls = updatedClusters$xy.clust
      tempdat = NULL
      for(i in cls){
        tempdat <- datmodsub %>% 
          ungroup() %>%
          filter(between(mt_time(.), 
                         clustertable_update$firstdatetime[clustertable_update$xy.clust == i], 
                         clustertable_update$lastdatetime[clustertable_update$xy.clust == i]),
                 mt_track_id(.) %in% c(strsplit(clustertable_update$birds[clustertable_update$xy.clust == i], split=", ")[[1]])) %>%
          mutate(xy.clust = ifelse(xy.clust == i, xy.clust, 0),
                 xy.clust = ifelse(is.na(xy.clust), 0, xy.clust),
                 incluster = ifelse(xy.clust==i, 1, 0),
                 indaycluster = ifelse(xy.clust==i & between(hour(mt_time(.)), 10, 15), 1,0)) %>%
          group_by(mt_track_id(.)) %>%
          as.data.frame() %>%
          summarise(xy.clust = i,
                    visitsinevent = sum(rle(incluster)$values),
                    dayvisits = sum(rle(indaycluster)$values)) %>%
          summarise(xy.clust = i,
                    visitsinevent_tot = sum(visitsinevent),
                    visitsinevent_mean = mean(visitsinevent),
                    dayvisits_tot = sum(dayvisits),
                    dayvisits_mean = mean(dayvisits)) %>%
          bind_rows(tempdat, .) %>%
          as.data.frame()# %>%
        #dplyr::select(-geometry)
      }
      
      if(!is.null(tempdat)) {
        clustertable_update <- merge(clustertable_update, tempdat, by = "xy.clust")
      }
      
      
      
      # Calculate distance to night points ---------------------------------------
      # Select date range plus one day either side of cluster
      
      logger.trace(paste0(as.Date(clusterdate), ": calculating distance to night points"))
      tempdat = NULL
      for(i in cls){
        nightpts <- datmodsub %>% 
          ungroup() %>%
          filter(between(mt_time(.), 
                         clustertable_update$firstdatetime[clustertable_update$xy.clust == i] - days(1),
                         clustertable_update$lastdatetime[clustertable_update$xy.clust == i] + days(1)),
                 mt_track_id(.) %in% c(strsplit(clustertable_update$birds[clustertable_update$xy.clust == i], split=", ")[[1]]),
                 hour(mt_time(.)) > 21 | hour(mt_time(.)) < 3) %>%
          as.data.frame() %>%
          dplyr::select(x,y)
        if(nrow(nightpts)>0){
          clpts <- datmodsub %>%
            ungroup() %>%
            filter(between(mt_time(.),
                           clustertable_update$firstdatetime[clustertable_update$xy.clust == i], 
                           clustertable_update$lastdatetime[clustertable_update$xy.clust == i]),
                   mt_track_id(.) %in% c(strsplit(clustertable_update$birds[clustertable_update$xy.clust == i], split=", ")[[1]])) %>%
            as.data.frame() %>%
            dplyr::select(x,y)
          
          d <- as.matrix(dist(rbind(nightpts, clpts)))[1:nrow(nightpts),(nrow(nightpts)+1):(nrow(nightpts)+nrow(clpts))] # problem line
          if(!is.null(dim(d))){
            d <- apply(d, 2, min)
          }
          t <- data.frame(xy.clust = i, nightdist_mean = mean(d), nightdist_median = median(d), nightdist_sd = sd(d))
        }else{
          t <- data.frame(xy.clust = i, nightdist_mean = NA, nightdist_median = NA, nightdist_sd = NA)
        }
        tempdat <- bind_rows(tempdat, t) %>%
          as.data.frame() 
        
      }
      
      if(!is.null(tempdat)) {
        clustertable_update <- merge(clustertable_update, tempdat, by = "xy.clust")
      }
    
    }
    
    # Return other clusters --------------------------------------------------
    logger.trace(paste0(as.Date(clusterdate), ": Merging back into fullclustertable"))
    
    if(!is.null(clusterDataDwnld)) {
      fullclustertable <- bind_rows(clustertable_update,
                                    filter(clusterDataDwnld, xy.clust %!in% clustertable_update$xy.clust)) %>%
        arrange(desc(lastdatetime))
    } else {
      
      if(!is.null(clustertable_update)) {
        fullclustertable <- clustertable_update %>%
          arrange(desc(lastdatetime))
      } else {
        fullclustertable <- NULL
      }
      
    }
    
    
    
  
    logger.trace(paste0(as.Date(clusterdate), ": Increasing clusterdate to ", clusterdate + days(clusterstep)))
    # Increase date and repeat
    if(!is.null(fullclustertable)) {clusterDataDwnld <- fullclustertable} else {clusterDataDwnld <- NULL} # ifelse doesn't wor with null statements
    # clusterDataDwnld <- ifelse(!is.null(fullclustertable), fullclustertable, NULL) # return it to move2 here
    logger.debug(paste0(as.Date(clusterdate), ": COMPLETE. ", nrow(fullclustertable), " total clusters generated"))
    clusterdate <- clusterdate + days(clusterstep)
    
  }
  

  return(list(clustereventdata = datmodsub, clustereventtable = fullclustertable))
}

