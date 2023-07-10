# ~~~~~~~~~~~~~~~~~~
# makeClusterFunctions.R
# ~~~~~~~~~~~~~~~~~~
# Jul 2023

# ~~~~~~~~~~~~~
# Reworking for MoveApps / sf
# ~~~~~~~~~~~~~

# Packages
library(move2)

# Shortened not in
`%!in%` <- Negate(`%in%`)


makeEventClusters <- function(data, d = 500) {
  
  # Get animal IDs
  tags <- unique(mt_track_id(data))

  # Filter all tags to non-travelling behaviour 
  tagdat <- filter(data, behav %!in% c("STravelling", "Unknown"))
  
  
  # Perform clustering on UTM coordinate data
  xy <- tagdat[, c("x", "y")] %>%
    #sf::st_as_sf(coords = c("x", "y"), crs = "+init=epsg:32733") %>%
    cbind(ID = seq(1:nrow(tagdat))) # Append ID
  
  mdist <- dist(cbind(tagdat$x, tagdat$y))
  hc <- hclust(as.dist(mdist), method = "complete")
  
  xy$clust <- cutree(hc, h=d)
  
  # find id of clusters with only 1 point and which rows in tag data
  cid <- as.vector(which(table(xy$clust)==1))
  rid <- which(xy$clust %in% cid == TRUE)
  
  # find the centroids of clusters
  cent <- matrix(ncol=2, nrow=max(xy$clust))
  #browser()
  
  for (i in 1:max(xy$clust)){
    # gCentroid from the rgeos package
    multi <- subset(xy, clust == i) %>%
      sf::st_coordinates() %>%
      sf::st_multipoint()
    cent[i,] <- sf::st_centroid(multi)
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

