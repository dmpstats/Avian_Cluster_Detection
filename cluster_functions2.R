# ~~~~~~~~~~~~~~~~~~
# CLUSTERING FUNCTIONS
# ~~~~~~~~~~~~~~~~~~
# Jan 2023
# Vulture Project
# Tanzania vultures

# ~~~~~~~~~~~~~
# Setup
# ~~~~~~~~~~~~~


# Shortened not in
`%!in%` <- Negate(`%in%`)

# ~~~~~~~~~~~~~
# Function for performing clustering
# ~~~~~~~~~~~~~


clustering <- function (datmodsub, clusterstartdate, clusterenddate, clusterDataDwnld = NULL, oldclassification = TRUE) {
  
  
  # Filter to the data being clustered:
  eventdata <- datmodsub %>%
    ungroup() %>%
    filter(between(datetime,  clusterstartdate - days(7) , clusterenddate + days(1)))
  
  # Set initial day to cluster:
  clusterdate <- ceiling_date(clusterstartdate)
  
  # Loop over the days within our clustering period:
  while(clusterdate <= ceiling_date(clusterenddate)) {
    
    # Unless we're on the first step, download cluster data
    if(!is.null(clusterDataDwnld)) {
      clusterDataDwnld <- fullclustertable
    }
    
    #if(clusterdate == as.POSIXct("2015-11-22 07:00:00", tz = "UTC")){browser()}
    
    #print(paste0("Clustering data leading up to ", toString(clusterdate)))
    
    # Sample to relevant 7 days
    clusteringData <-  eventdata %>%
      filter(between(datetime, clusterdate - days(7), clusterdate))
    
    if(nrow(clusteringData)==0){
      print(paste("Complete - no data", as.Date(clusterdate)))
      #browser()
      # find minimum date that is after current date (i.e. where we have data)
      clusterdate <- filter(eventdata, datetime > clusterdate) %>% 
        pull(datetime) %>% 
        min() + days(6)
      next
    }
    
    # Create new clusters
    clusters_new <- makeEventClusters(clusteringData, d=500, oldbehav = oldclassification)
    
    
    
    
    
    
    
    
    # match to existing clusters (from last 14 days)
    # create existing and new clust objects
    if(is.null(clusterDataDwnld)){
      existingclust <- NULL
    }else{
      existingclust <- clusterDataDwnld %>% 
        filter(lastdatetime > clusterdate - days(14)) %>%
        select(x, y, xy.clust)
    }
    newclust <- clusters_new$clusts
    
    
    if(is.null(existingclust) | is.null(newclust)){
      matchingclustermap <- data.frame(existID = NULL, updID = NULL)
    }else{
      # find distances between existing and new clusters
      clustdist <- as.matrix(dist(rbind(existingclust[,1:2], newclust[,1:2]))) # x and y are utms
      diag(clustdist) <- NA
      # rows are existing clusters, columns are new clusters
      clustdist <- as.matrix(clustdist[1:nrow(existingclust), (nrow(existingclust)+1): (nrow(existingclust) + nrow(newclust))])
      
      #find near ones (200m), i.e. new clusters that match with an exsiting one. 
      closeClusterIndices <- try(as_tibble(which(clustdist<(175), arr.ind = TRUE, useNames = TRUE)) %>%
                                   rename(loc1_rowIndex = row, loc2_rowIndex = col), silent = TRUE)
      if(class(closeClusterIndices)[1] == "try-error"){browser()}														 
      
      # map updIDs in newclusts to xy.clust ID in existing
      # MATCHES
      matchingclustermap <- data.frame(existID = existingclust$xy.clust[closeClusterIndices$loc1_rowIndex], 
                                       updID = newclust$xy.clust[closeClusterIndices$loc2_rowIndex])
    }
    
    
    # new cluster combo of existing clusters
    # keep existing and allocate new cluster points to one of new ones. 
    if(length(unique(matchingclustermap$updID)) != dim(matchingclustermap)[1]){
      ids <- which(table(matchingclustermap$updID)>1)
      for(u in 1:length(ids)){
        # find points allocated to upid
        tempdat <- filter(matchingclustermap, updID == names(ids)[u])
        allocpoints <- filter(clusters_new$tagdata, xy.clust %in% tempdat$updID) 
        existclustpoints <- filter(clusterDataDwnld, xy.clust %in% tempdat$existID) %>% select(x,y)
        
        tempdat$updID <- paste0(tempdat$updID, ".", tempdat$existID)
        # allocate points to nearest of existing ids
        pdists <- as.matrix(dist(rbind(existclustpoints, allocpoints[,c("x", "y")])))
        pdists <- pdists[(nrow(existclustpoints)+1):nrow(pdists), 1:nrow(existclustpoints)]
        mindists <- apply(pdists, 1, min)
        for(d in 1:length(mindists)){
          cid <- which(pdists[d,]==mindists[d])
          allocpoints$xy.clust[d] <- tempdat$updID[cid]
        }
        
        matchingclustermap <- filter(matchingclustermap, updID != names(ids)[u]) %>%
          bind_rows(., tempdat)
        clusters_new$tagdata <- filter(clusters_new$tagdata, xy.clust %!in% names(ids)[u]) %>%
          bind_rows(., allocpoints)
        newclust <- filter(newclust, xy.clust != names(ids)[u])
      }
    }
    
    # two new clusters matching one existing cluster
    # edit to retain existing cluster only. 
    if(length(unique(matchingclustermap$existID)) != dim(matchingclustermap)[1]){
      ids <- which(table(matchingclustermap$existID)>1)
      for(u in 1:length(ids)){
        # find points allocated to upid
        tempdat <- filter(matchingclustermap, existID == names(ids)[u])
        
        allocpoints <- filter(clusters_new$tagdata, xy.clust %in% tempdat$updID)
        existclustpoints <- filter(clusterDataDwnld, xy.clust %in% tempdat$existID)
        
        newname <- paste0(tempdat$updID[1], ".", tempdat$existID)
        
        newclust <- filter(newclust, xy.clust %!in% tempdat$updID)
        
        allocpoints <- allocpoints %>% mutate(xy.clust = newname[1])
        clusters_new$tagdata <- filter(clusters_new$tagdata, xy.clust %!in% tempdat$updID) %>%
          bind_rows(., allocpoints)
        
        tempdat$updID <- newname
        matchingclustermap <- filter(matchingclustermap, existID != names(ids)[u]) %>%
          bind_rows(., tempdat %>% slice(1))
      }
    }
    
    
    
    
    # EXISTING clusters
    existIDs <- existingclust$xy.clust[existingclust$xy.clust %!in% matchingclustermap$existID]
    
    if(length(existIDs)!=0){
      existclustermap = data.frame(existID = existIDs ,
                                   updID = NA)
    }else{existclustermap = NULL}
    
    # NEW clusters
    newIDs <- newclust$xy.clust[newclust$xy.clust %!in% matchingclustermap$updID]
    
    
    # find the ID to start the new clusters at
    if(length(newIDs)!=0){
      if(!is.null(existingclust)){
        startid <- max(clusterDataDwnld$xy.clust) + 1
      }else{startid <- 1}
      
      newclustermap = data.frame(existID = startid:(startid - 1 + length(newIDs)),
                                 updID = newIDs)
    }else{newclustermap=NULL}
    
    updatedclustermap <- rbind(matchingclustermap, newclustermap, existclustermap, c(NA, "upNA"))
    
    
    # change the xytagdata (7 days data)
    xytagdata <- clusters_new$tagdata %>% rename(updID = xy.clust) %>%
      left_join(., updatedclustermap) %>%
      rename(xy.clust = existID) %>%
      mutate(xy.clust = as.numeric(xy.clust)) %>%
      select(-updID)%>%
      suppressMessages()
    
    # change the newclust ids)
    newclust <- newclust %>% rename(updID = xy.clust) %>%
      left_join(., updatedclustermap) %>% 
      rename(xy.clust = existID) %>%
      mutate(xy.clust = as.numeric(xy.clust)) %>%
      select(-updID) %>%
      suppressMessages()
    
    
    # add back to 14 day cluster table by joining with existing clust
    if(is.null(existingclust) | is.null(newclust)){
      updatedClusters <- newclust
    }else{
      updatedClusters <- full_join(existingclust, newclust) %>%
        group_by(xy.clust) %>%
        summarise(x=mean(x), y=mean(y)) %>%
        filter(xy.clust %in% newclust$xy.clust) %>%
        suppressMessages()
    }
    
    
    # Add cluster information to datmodsub
    datmodsub <- datmodsub %>% 
      filter(index %!in% xytagdata$index) %>%
      bind_rows(., xytagdata) %>%
      arrange(Study, tag, datetime)
    
    
    
    
    
    
    
    
    
    # get table of cluster information
    clustertable_update <- makeclustertable(filter(datmodsub, xy.clust %in% updatedClusters$xy.clust),
                                            updatedClusters, oldbehav = oldclassification)
    
    #if(max(clustertable_update$nbirds)>1){browser()}
    
    ## calculate revisits
    cls = updatedClusters$xy.clust
    tempdat = NULL
    for(i in cls){
      tempdat <- datmodsub %>% 
        ungroup() %>%
        filter(between(datetime, 
                       clustertable_update$firstdatetime[clustertable_update$xy.clust == i], 
                       clustertable_update$lastdatetime[clustertable_update$xy.clust == i]),
               LocalID %in% c(strsplit(clustertable_update$birds[clustertable_update$xy.clust == i], split=", ")[[1]])) %>%
        mutate(xy.clust = ifelse(xy.clust == i, xy.clust, 0),
               xy.clust = ifelse(is.na(xy.clust), 0, xy.clust),
               incluster = ifelse(xy.clust==i, 1, 0),
               indaycluster = ifelse(xy.clust==i & between(hour, 10, 15), 1,0)) %>%
        group_by(LocalID) %>%
        summarise(xy.clust = i,
                  visitsinevent = sum(rle(incluster)$values),
                  dayvisits = sum(rle(indaycluster)$values)) %>%
        summarise(xy.clust = i,
                  visitsinevent_tot = sum(visitsinevent),
                  visitsinevent_mean = mean(visitsinevent),
                  dayvisits_tot = sum(dayvisits),
                  dayvisits_mean = mean(dayvisits)) %>%
        bind_rows(tempdat, .)
    }
      
   clustertable_update <- left_join(clustertable_update, tempdat) 
   
   ## distance to night points
   # select date range plus one day either side of cluster
  
  # browser()
   tempdat = NULL
   for(i in cls){
     nightpts <- datmodsub %>% 
       ungroup() %>%
       filter(between(datetime, 
                      clustertable_update$firstdatetime[clustertable_update$xy.clust == i] - days(1),
                      clustertable_update$lastdatetime[clustertable_update$xy.clust == i] + days(1)),
              LocalID %in% c(strsplit(clustertable_update$birds[clustertable_update$xy.clust == i], split=", ")[[1]]),
              hour > 21 | hour < 3) %>%
       select(x,y)
     if(nrow(nightpts)>0){
       clpts <- datmodsub %>%
         ungroup() %>%
         filter(between(datetime,
                        clustertable_update$firstdatetime[clustertable_update$xy.clust == i], 
                        clustertable_update$lastdatetime[clustertable_update$xy.clust == i]),
                LocalID %in% c(strsplit(clustertable_update$birds[clustertable_update$xy.clust == i], split=", ")[[1]])) %>%
         select(x,y)
       d <- as.matrix(dist(rbind(nightpts, clpts)))[1:nrow(nightpts),(nrow(nightpts)+1):(nrow(nightpts)+nrow(clpts))]
       if(!is.null(dim(d))){
         d <- apply(d, 2, min)
       }
       t <- data.frame(xy.clust = i, nightdist_mean = mean(d), nightdist_median = median(d), nightdist_sd = sd(d))
     }else{
       t <- data.frame(xy.clust = i, nightdist_mean = NA, nightdist_median = NA, nightdist_sd = NA)
     }
      tempdat <- bind_rows(tempdat, t)
   }
   
   clustertable_update <- left_join(clustertable_update, tempdat)
   
   ## classify importance
    
    # add back in values for clusters not in 14 days
    if(!is.null(clusterDataDwnld)){
      fullclustertable <- bind_rows(clustertable_update, 
                                    filter(clusterDataDwnld, xy.clust %!in% clustertable_update$xy.clust)) %>%
        arrange(desc(lastdatetime))
      
    }else{
      fullclustertable <- clustertable_update %>%
        arrange(desc(lastdatetime))
      
    }
    
    
    # Increase date and repeat
    clusterdate <- clusterdate + days(5) # move along 5 days
    clusterDataDwnld <- fullclustertable
    print(paste("Complete", as.Date(clusterdate)))
    
  }
  
  return(list(clustereventdata = datmodsub, clustereventtable = fullclustertable))
  
}



# ~~~~~~~~~~~~~
# makeEventClusters function
# ~~~~~~~~~~~~~


makeEventClusters <- function(data, d = 500, oldbehav = TRUE) {
  tags <- unique(data$tag)
  
  # filter on all tags, non-travelling behaviour and last 7 days
  if (oldbehav == TRUE) {
    tagdat <- filter(data, tag %in% tags, 
                     behavold %!in% c("STravelling", "Unknown"))
  } else {
    tagdat <- filter(data, tag %in% tags, 
                     behav %!in% c("XTravelling", "XSlowTravel", "XLeaving"))
  }
  

  # perform clustering
  xy <- sp::SpatialPointsDataFrame(tagdat[, c("x", "y")], data.frame(ID = seq(1:nrow(tagdat))))
  mdist <- dist(cbind(tagdat$x, tagdat$y))
  hc <- hclust(as.dist(mdist), method = "complete")
  #par(cex=1,font=3)
  #plot(hc, xlab = NA, main = "Dendrogram of Tag Clustering", labels = FALSE, ylim = c(0, 1000))
  #abline(h = 500, col = "red", lty = 'dashed')
  xy$clust <- cutree(hc, h=d)

  # find id of clusters with only 1 point and which rows in tag data
  cid <- as.vector(which(table(xy$clust)==1))
  rid <- which(xy$clust %in% cid == TRUE)
  
  # find the centroids of clusters
  cent <- matrix(ncol=2, nrow=max(xy$clust))
  for (i in 1:max(xy$clust)){
    # gCentroid from the rgeos package
    cent[i,] <- rgeos::gCentroid(subset(xy, clust == i))@coords
  }
    
  
  clusts <- as_tibble(cent) %>% 
    mutate(xy.clust = row_number()) %>%
    filter(xy.clust %!in% cid) %>%
    mutate(xy.clust = paste0("up", xy.clust)) %>%
    rename(x=V1, y=V2)
  
  xydata <- tagdat %>% ungroup () %>%
    mutate(xy.clust = xy$clust) %>%
    mutate(xy.clust = paste0("up", xy.clust)) %>%
    mutate(xy.clust = if_else(xy.clust %in% unique(clusts$xy.clust), xy.clust, "upNA"))

  
  return(list(tagdata = xydata, clusts = clusts))
}



# ~~~~~~~~~~~~~
# makeclustertable function
# ~~~~~~~~~~~~~


makeclustertable <- function(xytagdata, updatedClusters, oldbehav = TRUE){
  
  
  # calculate proportion of points feeding, roosting, resting
  if(oldbehav == TRUE) {
    # For old behaviour classification
    clust.behavold <- xytagdata %>%
      group_by(xy.clust) %>%
      count(behavold) %>%
      arrange(xy.clust)  %>%
      pivot_wider(values_from = n, names_from = behavold, values_fill=0) %>%
      mutate(Total = sum(across(starts_with ("S")), na.rm = T),
             propfeed = ifelse("SFeeding" %in% colnames(.), SFeeding/Total, 0),
             proprest = ifelse("SResting" %in% colnames(.), SResting/Total, 0),
             proproost = ifelse("SRoosting" %in% colnames(.), SRoosting/Total, 0))
  } else {
    # For new behaviour classification
    clust.behav <- xytagdata %>%
      group_by(xy.clust) %>%
      count(behav) %>%
      arrange(xy.clust)  %>%
      pivot_wider(values_from = n, names_from = behav, values_fill=0) %>%
      mutate(Total = sum(across(starts_with ("X")), na.rm = T),
             XDayRoost = ifelse("XDayRoost" %in% colnames(.), XDayRoost, 0),
             XNightRoost = ifelse("XNightRoost" %in% colnames(.), XNightRoost, 0),
             propfeed = ifelse("XFeeding" %in% colnames(.), XFeeding/Total, 0),
             proproost = (XDayRoost + XNightRoost)/Total)
  }

 
  
  # calculate number of days, birds and dates of each cluster
  clust.days <- xytagdata %>%
    group_by(xy.clust) %>%
    summarise(days = length(unique(yearmonthday)),
              nbirds = length(unique(tag)),
              birds = paste(unique(LocalID), collapse=", "),
              firstdatetime = min(datetime),
              lastdatetime = max(datetime),
              Studys = paste(unique(Study), collapse=", ")
    )
  
  if(oldbehav == TRUE) {
    clust.table <- left_join(clust.behavold, clust.days)%>%
      suppressMessages()
    # calculate the time spent at each cluster (only counting non-travelling points)
    clust.time <- xytagdata %>%
      group_by(xy.clust) %>%
      summarise(TimeTotal = sum(timediff_hrs),
                TimeDay = sum(timediff_hrs[hour>5 & hour < 17]),
                TimeFeed = sum(timediff_hrs[behavold=="SFeeding"]),
                MedianHourFeed = median(hour[behavold=="SFeeding"]),
                MedianHourDay = median(hour[hour>5 & hour < 17]),
                DistMedian = median(dist(cbind(x,y))),
                DistSD = sd(dist(cbind(x,y))),
                medloc =  Gmedian::Weiszfeld(data.frame(x=x, y=y))$median
      ) %>%
      mutate(MedianHourFeed = ifelse(is.na(MedianHourFeed), 0, MedianHourFeed),
             DistSD = ifelse(is.na(DistSD), 0, DistSD),
             x.med = medloc[,1],
             y.med = medloc[,2]) %>%
      select(-medloc)
  } else {
    clust.table <- left_join(clust.behav, clust.days)%>%
      suppressMessages()

    # calculate the time spent at each cluster (only counting non-travelling points)
    clust.time <- xytagdata %>%
      group_by(xy.clust) %>%
      summarise(TimeTotal = sum(timediff_hrs),
                TimeDay = sum(timediff_hrs[hour>5 & hour < 17]),
                TimeFeed = sum(timediff_hrs[behav=="XFeeding"]),
                MedianHourFeed = median(hour[behav=="XFeeding"]),
                MedianHourDay = median(hour[hour>5 & hour < 17]),
                DistMedian = median(dist(cbind(x,y))),
                DistSD = sd(dist(cbind(x,y))),
                medloc =  Gmedian::Weiszfeld(data.frame(x=x, y=y))$median
      ) %>%
      mutate(MedianHourFeed = ifelse(is.na(MedianHourFeed), 0, MedianHourFeed),
             DistSD = ifelse(is.na(DistSD), 0, DistSD),
             x.med = medloc[,1],
             y.med = medloc[,2]) %>%
      select(-medloc)
  }
  
  
  clust.table <- left_join(clust.table, clust.time) %>% ungroup()%>%
    suppressMessages()
  clust.table <- left_join(clust.table, updatedClusters)%>%
    suppressMessages()

  
  # add lats and longs
  xy <- data.frame(ID = 1:nrow(clust.table), x = clust.table$x, y = clust.table$y)
  coordinates(xy) <- c("x", "y")
  proj4string(xy) <- CRS("+init=epsg:32733")
  res <- as.data.frame(spTransform(xy, CRS("+proj=longlat +datum=WGS84")))
  
  xy <- data.frame(ID = 1:nrow(clust.table), x = clust.table$x.med, y = clust.table$y.med)
  coordinates(xy) <- c("x", "y")
  proj4string(xy) <- CRS("+init=epsg:32733")
  res.med <- as.data.frame(spTransform(xy, CRS("+proj=longlat +datum=WGS84")))
  
  clust.table %<>% mutate(lon=res$x, 
                          lat=res$y,
                          lon.med=res.med$x, 
                          lat.med=res.med$y)
  
  return(clust.table)
  
}





calcimportance <- function(clust.table){
  # calculate importance:
  clust.table <- clust.table %>% 
    mutate(importance = (propfeed*100) * days * TimeFeed * (nbirds*10)) %>%
    arrange(desc(importance),
            desc(firstdatetime))
  
  
  # band importance score
  #bandthresh <- quantile(clust.table$importance, probs=c(0, 0.25, 0.5, 0.7, 0.8, 0.9, 1))
  bandthresh <- c(0, 325, 1800, 2250, 10315, 60000, Inf)
  clust.table %<>%
    mutate(impBand = if_else(importance == 0, 0, as.numeric(cut(importance, breaks=c(bandthresh)))))
  
  
  
  risknames <- data.frame(impBand = c(0:6), impBandchr = c("No-feeding", "Quitelow", "Quitelow", "Quitelow", "Medium", "Quitehigh", "High"))
  
  clust.table <- left_join(clust.table, risknames)%>%
    suppressMessages()
  
  return(clust.table)
  
}
