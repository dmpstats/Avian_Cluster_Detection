# ~~~~~~~~~~~~~~~~~~
# clustering.R
# ~~~~~~~~~~~~~~~~~~
# Jul 2023

# ~~~~~~~~~~~~~
# Reworking for MoveApps / sf
# ~~~~~~~~~~~~~

makeclustertable <- function(xytagdata, updatedClusters, behavsystem = TRUE) {
  
  # Necessary columns
  xytagdata %<>%
    mutate(
      tag = mt_track_id(.),
      datetime = mt_time(.)
    )
  
  if(behavsystem == TRUE) {
    # Calculate proportion of feeding, roosting, and resting points
    clust <- xytagdata %>%
      as.data.frame() %>% # need to convert and re-convert for this method
      group_by(xy.clust) %>%
      count(behav) %>%
      arrange(xy.clust) %>%
      tidyr::pivot_wider(values_from = n, names_from = behav, values_fill = 0)%>%
      mutate(Total = sum(across(starts_with("S")), na.rm = T),
             propfeed = ifelse("SFeeding" %in% colnames(.), SFeeding/Total, 0), 
             proprest = ifelse("SResting" %in% colnames(.), SResting/Total, 0),
             proproost = ifelse("SRoosting" %in% colnames(.), SRoosting/Total, 0))
  } else {
    clust <- xytagdata %>%
      as.data.frame() %>% # need to convert and re-convert for this method
      group_by(xy.clust) %>%
      count(tag) %>%
      arrange(xy.clust) %>%
      tidyr::pivot_wider(values_from = n, names_from = tag, values_fill = 0)
  }
  
  # Calculate number of days, birds, and dates of each cluster
  clust.days <- xytagdata %>%
    group_by(xy.clust) %>%
    summarise(days = length(unique(date(datetime))),
              nbirds = length(unique(tag)),
              birds = paste(unique(tag), collapse = ", "),
              firstdatetime = min(datetime),
              lastdatetime = max(datetime)
              )
  
  clust.table <- left_join(clust, clust.days) %>%
    suppressMessages()
  
  # Calculate stationary time for each cluster (excluding travelling points)
  
  if(behavsystem == TRUE) {
    clust.time <- xytagdata %>%
      group_by(xy.clust) %>%
      summarise(TimeTotal = sum(timediff_hrs),
                TimeDay = sum(timediff_hrs[hour > 5 & hour < 17]),
                TimeFeed = sum(timediff_hrs[behav == "SFeeding"]),
                MedianHourFeed = median(hour[behav == "SFeeding"]),
                MedianHourDay = median(hour[hour > 5 & hour < 17]),
                DistMedian = median(dist(cbind(x, y))),
                DistSD = sd(dist(cbind(x, y))),
                medloc = Gmedian::Weiszfeld(data.frame(x = x, y = y))$median
      ) %>%
      
      mutate(MedianHourFeed = ifelse(is.na(MedianHourFeed), 0, MedianHourFeed),
             DistSD = ifelse(is.na(DistSD), 0, DistSD),
             x.med = medloc[, 1],
             y.med = medloc[, 2]) %>%
      dplyr::select(-medloc)
    
    clust.table <- left_join(clust.table, clust.time) %>% ungroup() %>%
      suppressMessages()
    clust.table <- left_join(clust.table, updatedClusters) %>%
      suppressMessages()
  } else {
    clust.time <- xytagdata %>%
      group_by(xy.clust) %>%
      summarise(TimeTotal = sum(timediff_hrs),
                TimeDay = sum(timediff_hrs[hour > 5 & hour < 17]),
                MedianHourDay = median(hour[hour > 5 & hour < 17]),
                DistMedian = median(dist(cbind(x, y))),
                DistSD = sd(dist(cbind(x, y))),
                medloc = Gmedian::Weiszfeld(data.frame(x = x, y = y))$median
      ) %>%
      
      mutate(DistSD = ifelse(is.na(DistSD), 0, DistSD),
             x.med = medloc[, 1],
             y.med = medloc[, 2]) %>%
      dplyr::select(-medloc)
    
    clust.table <- left_join(clust.table, clust.time) %>% ungroup() %>%
      suppressMessages()
    clust.table <- left_join(clust.table, updatedClusters) %>%
      suppressMessages()
  }

  
  
  # # Return to move2 object
  # clust.table <- move2::mt_as_move2(
  #   clust.table,
  #   coords = c("x", "y"),
  #   time_column = "firstdatetime",
  #   track_id_column = "xy.clust"
  # )
  
  # CC: Should we transform to lat/lon? I've used UTMs as standard output so far - leaving it this way for now

  
  return(clust.table)
  
}


