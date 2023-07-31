# Avian Cluster Detection

MoveApps

Github repository: *github.com/callumjclarke/Avian_Cluster_Detection*

## Description

This MoveApp geographically clusters behaviourally-classified data over a rolling window, generating data on the cluster's properties (hours spent, number of revisits, number of animals, and so on). Clusters are updated with each time step to simulate real-time identification.

## Documentation

Input data to this MoveApp needs to have been processed by **Standardise Formats and Calculate Basic Statistics,** or a similar MoveApp, to append UTM data (i.e. columns named `x` and `y` containing UTM coordinate data).

If applying this to bird-tracking data, we recommend also using **Behavioural Classification for Vultures,** or an alternative MoveApp, which appends a column of predicted behaviours to each event. Currently, the only supported behaviours are `SFeeding`, `SResting`, `SRoosting` and `STravelling`. This will be incorporated in the clustering (and only stationary behaviours are clustered) - otherwise, all points will be used.

The data is filtered to points in the relevant time window (and, if selected, the relevant behaviours) and an agglomerative clustering process is applied. The rolling window then moves forward by the selected time-step and repeats the process. Any clusters that are geographically (centroids within 200m of one another) and temporally (both have points within the same `Cluster Window` duration) similar are merged into one renamed cluster; similarly, any clusters that have surpassed the given cluster expiration window without any further points being added are 'closed'.

If any given event is associated with a cluster, the cluster's name is appended to the tracking data. A cluster-table is generated at each step and merged with all previous tables. A final 'master' cluster-table is outputted as an artifact; both objects are in *Move2* format.

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

`clustertable.csv`: A *Move2* object table listing certain attributes of each cluster (number of points, number of IDs, duration, etc.)

### Settings

`Cluster Start Date`: The earliest date on which to consider data for clustering (provided as an instant)

`Cluster End Date`: The final date on which to consider data for clustering (provided as an instant)

`Cluster Step`: The number of days to increase by between each 'rolling window' step (integer). The larger this is, the faster the app will run. If the chosen `Cluster Step` is larger than the `Cluster Window`, clusters cannot be merged.

`Cluster Window`: How many days should be considered in each 'rolling window' step (integer). For example, `Cluster Window = 7` means that the data will be clustered in a rolling 7-day window.

`Maximum inactive cluster event duration`: How many days of inactivity should pass before a cluster should be considered 'closed'. After a cluster is closed, no further data can be incorporated to this cluster.

`Behavioural Classification Clustering`: If the data has been processed by **Behavioural Classification for Vultures,** select this. If the data has not, do not select this - the clustering will be performed on all available data without behavioural classification.

### Most common errors

-   Selecting a `Cluster Step` larger than the selected `Cluster Window` means data could be skipped between rolling-windows. For example, `Cluster Step = 7` but `Cluster Window = 7` means 5 days will be clustered and then the window will be increased by 7 days, skipping 2 days of data

-   If the data provided has significantly large or small gaps between events, adjusting the `Cluster Step` and `Cluster Window` is essential. Otherwise, outputs can either miss out potential clusters, or perform too few iterations

### Null or error handling

**Setting** `Behavioural Classification Clustering:` Selecting this without a column named `behav` will return the data with a warning

**Setting** `Cluster Start Date`: If an invalid date is provided, the default will be two weeks prior to the final point of data

**Setting** `Cluster End Date`: If an invalid date is provided, the default will be the timestamp of the final point of data
