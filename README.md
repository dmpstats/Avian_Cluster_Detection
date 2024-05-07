# Avian Cluster Detection

MoveApps

Github repository: *https://github.com/dmpstats/Avian_Cluster_Detection*

## Description

This MoveApp utilizes a rolling time-window to spatially cluster location events, employing an iterative approach for cluster detection and updating aimed at simulating real-time identification.

## Documentation

The App employs a spatial agglomerative clustering process on location events, covering either the entire input dataset or a specified time-period within it. This process operates over a rolling time-window that 'closes' clusters after a period of inactivity. The time-window then advances by the chosen time-step, iterativelly reapplying the clustering process until reaching the last timestamp entry or the end date-time of the specified time-period. Cluster detection and updating are performed using the *median location* of clustered events, calculated via the [`Gmedian::Weiszfeld`](https://rdrr.io/cran/Gmedian/man/Gmedian.html) function. Clustered location events are annotated with a cluster ID code, which is appended to the input data as a column named `xy.clust`. 

Note: To use this process with no rolling window (i.e. perform a single clustering step of all input data):

  1. Select option **Clustering on Entire Data**.
  2. Specify **Clustering Window** to a number of days longer than the duration of your dataset (in days).



### MoveApps Worflow Dependencies

Choosing the **Use Behavioural Attribute** option requires the prior deployment of the ['Behavioural Classification for Vultures App'](https://www.moveapps.org/apps/browser/44bb2ffa-7d40-4fad-bff5-1269995ba1a2) ([Readme](https://github.com/dmpstats/Behavioural_Classification_for_Vultures)) in the workflow.



### Input data

A `move2::move2_loc` object.

### Output data

A `move2::move2_loc` object, with appended column `xy.clust` (which cluster, if any, a location is associated with).

### Artefacts

None

### Settings

**Clustering Start Timepoint** (`clusterstart`) : Start date-time of period within input data to apply clustering. If `NULL` (default), the start timepoint is set to the earliest timestamp in input data. Note: This input is only considered if **Clustering on Entire Data** is not selected. 

**Clustering End Timepoint** (`clusterend`): Final date-time of period within input data to apply clustering. If `NULL` (default), the end timepoint is set to the latest timestamp in input data. Note: This input is only considered if **Clustering on Entire Data** is not selected.

**Clustering Window** (`clusterwindow`) : Duration, in number of days, of the rolling time-window within which clustering is conducted. For example, setting this value to 7 (default) means clustering is applied iteratively to a 7-day rolling window.

**Clustering Time-step** (`clusterstep`): Number of days to advance between each iterative step of the rolling time-window. Defaults to a 1-day incremental step, which may be computationally demanding for large datasets. 

**Maximum Inactive Cluster Duration** (`clusterexpiration`): Number of days of inactivity after which a cluster is considered 'closed'. That is, if no new location has been added to a given cluster for this duration, no further locations will be allocated to that cluster.

**Use Behavioural Attribute** (`behavsystem`): Choose whether to consider animal behaviour in the analysis. If selected, clustering is restricted to stationary-type classified locations (e.g. feeding, resting or roosting). Otherwise, clustering includes all available data. Note: this option requires prior deployment of the 'Behavioural Classification for Vultures' App in the workflow.

**Agglomerative Clustering Cut-Height** (`d`): The height at which to cut the dendrogram generated by the hierarchical clustering process (unit: meters).

**Cluster ID Prefix** (`clustercode`): Optional prefix to add to automatically generated cluster ID codes. For instance, setting it to 'A' (default) results in clusters labelled 'A.1', 'A.2', etc. Useful for post-processing, e.g. merging outputs from different datasets or studies.



### Most common errors

The app will halt processing and return an error under the following conditions:

- Selecting **Use Behavioural Attribute** when column `behav` is not present in the input dataset. Column `behav` is created in the App 'Behavioural Classification for Vultures', which must be deployed earlier in the workflow.

- Specifying a **Clustering Start** or **Clustering End Timepoint** outside the time-period covered by the input data.

- Selecting a **Clustering Time-step** larger than **Clustering Window**, as that would lead to data being skipped between consecutive steps of the iterative clustering process.

In addition, if the input data has significantly large or small gaps between events, adjusting the **Clustering Time-step** and **Clustering Window** is essential. Otherwise, outputs can either miss out potential clusters, or perform too few iterations.


### Null or error handling

- **Clustering Start Timepoint**: If `NULL`, start date will be set to two weeks prior to the latest timestamp of the input data.

- **Clustering End Timepoint**: If `NULL`, end date will be set to the latest timestamp of the input data.
