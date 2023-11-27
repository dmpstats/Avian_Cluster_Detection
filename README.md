# Avian Cluster Detection

MoveApps

Github repository: *github.com/callumjclarke/Avian_Cluster_Detection*

## Description

This MoveApp spatially clusters data over a rolling window, generating data on the cluster's properties (hours spent, number of revisits, number of animals, and so on). Clusters are updated with each time step to simulate real-time identification.

## Documentation

This MoveApp is intended for use on the output of **Behavioural Classification for Vultures,** and its functionality in other workflows is not guaranteed.

A spatial agglomerative clustering process is applied to all locations within the specified timeframe, over a rolling-window that 'closes' clusters after a period of inactivity. The rolling window then moves forward by the selected time-step, and repeats the clustering, iteratively until the end of the specified timeframe is reached.

For each location associated with a cluster by this process, the cluster's name is appended to the tracking data as a column named *xy.clust*. At the end of the process, a cluster-table (detailing various cluster attributes) is also generated and released as a *move2* object.

To use this process with no rolling window (i.e. perform a single clustering step of all input data):

-   Set the **cluster start date** to the earliest possible date

-   Leave the **cluster end date** null: this will default to the final date containing data

-   Set the **cluster window** to a number of days longer than the duration of your dataset (in days)

Please note: the geometry of the outputted clusters is determined by their **median location,** not mean location. This is generated using the GMedian `Weiszfeld` function.

### Input data

*Move2* Location Data

### Output data

*Move2* Location Data, with appended column `xy.clust` (which cluster, if any, a location is associated with)

### Artefacts

`clustertable.rds`: A *Move2* object table listing certain attributes of each cluster (number of points, number of IDs, duration, etc.)

### Settings

`Cluster Start Date`: The earliest date on which to consider data for clustering (provided as an instant)

`Cluster End Date`: The final date on which to consider data for clustering (provided as an instant)

`Cluster Step`: The number of days to increase by between each 'rolling window' step (integer). The larger this is, the faster the app will run. If the chosen `Cluster Step` is larger than the `Cluster Window`, clusters cannot be merged.

`Cluster Window`: How many days should be considered in each 'rolling window' step (integer). For example, `Cluster Window = 7` means that the data will be clustered in a rolling 7-day window.

`Maximum inactive cluster event duration`: How many days of inactivity should pass before a cluster should be considered 'closed'. After a cluster is closed, no further data can be incorporated to this cluster.

`Behavioural Classification Clustering`: If the data has been processed by **Behavioural Classification for Vultures,** select this. If the data has not, do not select this - the clustering will be performed on all available data without behavioural classification.

`Agglomerative Clustering Cut-Height`: The height at which to cut the dendrogram after each iterative clustering. `d = 500` is default.

`Cluster ID Code`: An option string to append to the ID of clusters. For example, a clustercode of `abc` means that output clusters will be identified as `abc.1`, `abc.2`, and so on. This is useful if separate studies are being merged. Defaults to no clustercode.

### Most common errors

-   Selecting a `Cluster Step` larger than the selected `Cluster Window` means data could be skipped between rolling windows. For example, `Cluster Step = 7` but `Cluster Window = 7` means 5 days will be clustered and then the window will be increased by 7 days, skipping 2 days of data

-   If the data provided has significantly large or small gaps between events, adjusting the `Cluster Step` and `Cluster Window` is essential. Otherwise, outputs can either miss out potential clusters, or perform too few iterations

### Null or error handling

**Setting** `Behavioural Classification Clustering:` Selecting this without a column named `behav` will return the data with a warning

**Setting** `Cluster Start Date`: If an invalid date is provided, the default will be two weeks prior to the final point of data

**Setting** `Cluster End Date`: If an invalid date is provided, the default will be the timestamp of the final point of data
