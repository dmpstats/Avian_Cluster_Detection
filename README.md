# Avian Cluster Detection

MoveApps

Github repository: *https://github.com/dmpstats/Avian_Cluster_Detection*

## Description

This MoveApp utilizes a rolling time-window to spatially cluster location events, employing an iterative approach for cluster detection and cluster updating aimed at simulating real-time identification.

## Documentation

The App employs a spatial agglomerative clustering process on location events, covering either the entire input dataset or a specified time-period within it. This process operates over a rolling time-window that 'closes' clusters after a period of inactivity. The time-window then advances by the chosen time-step, iteratively reapplying the clustering process until reaching the last timestamp entry of the data under clustering. Clustered location events are annotated with a cluster ID code, which is appended to the input data as a column named `clust_id`. 

Below are the key steps of the clustering process:

1. Select the data for clustering by filtering the input data based on the chosen start (`clusterstart`) and end (`clusterend`) timepoints.

2. Set the initial clustering time-window using the start timepoint and the chosen duration for the rolling time-window (`clusterwindow`).

3. Apply hierarchical clustering [`stats::hclust`](https://rdrr.io/r/stats/hclust.html) to the distance matrix of locations within the current time-window. Form clusters by employing the chosen cut-height (`d`) to the resulting dendrogram, i.e. grouping locations so that within-group distances are smaller than the specified cut-height value.

4. Get the centroid of each detected cluster, calculated as the of median of cluster locations (via [`Gmedian::Weiszfeld`](https://rdrr.io/cran/Gmedian/man/Gmedian.html)).

5. Advance the current time-window by the specified time-step (`clusterstep`).

6. Repeat steps 3 and 4 to generate clusters in current iteration.

7. Match current clusters with those detected in the previous iteration based on a proximity threshold (`match_thresh`) between cluster centroids. Update clusters using the following actions for each "previous-to-current" case:

     - "One-to-One"" matches represent either expanding or stable clusters. Expanding clusters are updated with newest locations.
     
     - "None-to-One" matches are treated as newly formed clusters. Locations comprised in each new cluster are annotated with a freshly issued cluster ID.
     
     - "Multiple-to-One" matches are taken as evolving neighboring clusters. New locations are assigned to the nearest previous clusters.
     
     - "One-to-Multiple" matches are treated as previously existent clusters that are expanding sparsely. Locations on multiple current clusters are assigned to the previous cluster.
     
8. Re-calculate the centroids of updated clusters.

9. Tag stable clusters exceeding the chosen inactive expiration period (`clusterexpiration`) so that they are not considered in step 7 of next iteration.

10. Repeat steps 5 - 9 iteratively until the time-window's endpoint exceeds the final timestamp in the clustering data.


#### Some considerations on clustering parameters

- **Clustering Window** plays a pivotal role in this clustering process by determining the time period over which locations are eligible for grouping into clusters via the hierarchical clustering algorithm. Essentially, it sets the maximum allowed time gap between consecutive locations within a cluster.

- If the specified **Clustering Window** exceeds the time span covered by the data being clustered, cluster detection occurs only once on a fixed time window, i.e. the iterative steps are skipped. In such instances, the parameter **Maximum Inactive Cluster Duration** becomes irrelevant, and the maximum time gap between consecutive cluster element locations is dictated by the time span of the clustering data.

- The **Clustering Time-step** parameter presents a trade-off between the App runtime performance and the accuracy of cluster matching (step 7, above). While considering the temporal resolution of the input data, opting for smaller time-step values is typically advantageous, granted that computational time constraints permit their use.

- Specifying values other than NULL for either the **Clustering Start** or **End Timepoints** will result in outputs containing partially clustered location events.

- The Weiszfeld geometric median provides a more robust calculation of cluster centroids, which makes them less susceptible to the influence of outlier cluster locations. This stands in contrast to standard centroid calculation methods, such as `sf::st_centroid()`, which relies on the mean of cluster locations.


### MoveApps Worflow Dependencies

Choosing the **Use Behavioural Attribute** option requires the prior deployment of the ['Behavioural Classification for Vultures App'](https://www.moveapps.org/apps/browser/44bb2ffa-7d40-4fad-bff5-1269995ba1a2) ([Readme](https://github.com/dmpstats/Behavioural_Classification_for_Vultures)) in the workflow.


### Input data

A `move2::move2_loc` object.

### Output data

A `move2::move2_loc` object, with appended column `clust_id` (which cluster, if any, a location is associated with).

### Artefacts

None

### Settings

**Clustering Start Timepoint** (`clusterstart`) : Start date-time of period within input data to apply clustering. If `NULL` (default), the start timepoint is set to the earliest timestamp in input data.

**Clustering End Timepoint** (`clusterend`): Final date-time of period within input data to apply clustering. If `NULL` (default), the end timepoint is set to the latest timestamp in input data.

**Clustering Window** (`clusterwindow`) : Duration, in number of days, of the rolling time-window within which clustering is conducted. For example, setting this value to 7 (default) means clustering is applied iteratively to a 7-day rolling window.

**Clustering Time-step** (`clusterstep`): Number of days to advance between each iterative step of the rolling time-window. Defaults to a 1-day incremental step, which may be computationally demanding for large datasets. 

**Agglomerative Clustering Cut-Height** (`d`): The height at which to cut the dendrogram generated by the hierarchical clustering process (unit: meters; default: 500m).

**Cluster Matching Proximity Threshold** (`match_thresh`): The maximum centroid distance allowed for linking clusters across iterative clustering steps (units: meters; default: 175m). Consecutive clusters with centroids distanced beyond this value are considered separate entities.

**Use Behavioural Attribute** (`behavsystem`): Choose whether to consider animal behaviour in the analysis (default: `TRUE`). If selected, clustering is restricted to stationary-type classified locations (e.g. feeding, resting or roosting). Otherwise, clustering includes all available data. Note: this option requires prior deployment of the 'Behavioural Classification for Vultures' App in the workflow.

**Maximum Inactive Cluster Duration** (`clustexpiration`): Number of days of inactivity after which a cluster is considered 'closed' (default: 14 days). That is, if no new location has been added to a given cluster for this duration, no further locations will be allocated to that cluster.

**Cluster ID Prefix** (`clustercode`): Optional prefix to add to automatically generated cluster ID codes. For instance, setting it to 'A' (default) results in clusters labelled 'A.1', 'A.2', etc. Useful for post-processing, e.g. merging outputs from different datasets or studies.



### Most common errors

The app will halt processing an throw an error under the following conditions:

- Selecting **Use Behavioural Attribute** when column `behav` is not present in the input dataset. Column `behav` is created in the App 'Behavioural Classification for Vultures', which must be deployed earlier in the workflow.

- Specifying a **Clustering Start** or **End Timepoint** outside the time-period covered by the input data.

- Selecting a **Clustering Time-step** larger than **Clustering Window**, as that would lead to data being skipped between consecutive steps of the clustering process.

In addition, if the input data has significantly large or small gaps between events, adjusting the **Clustering Time-step** and **Clustering Window** is essential. Otherwise, outputs can either miss out potential clusters, or perform too few iterations.


### Null or error handling

- **Clustering Start Timepoint**: If `NULL`, start date will be set to earliest timestamp of the input data.

- **Clustering End Timepoint**: If `NULL`, end date will be set to the latest timestamp of the input data.

- The following parameters are mandatory, i.e. if any of these is set to `NULL`, the app will terminate abruptly and throw an error:
  - **Clustering Window**
  - **Clustering Time-step**
  - **Agglomerative Clustering Cut-Height**
  - **Cluster Matching Proximity Threshold** 
  - **Maximum Inactive Cluster Duration** 


