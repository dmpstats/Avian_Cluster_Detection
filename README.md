# Avian Cluster Detection

MoveApps

Github repository: *github.com/yourAccount/Name-of-App* *(the link to the repository where the code of the app can be found must be provided)*

## Description

This MoveApp geographically clusters behaviourally-classified data over a rolling window, generating data on the cluster's properties (hours spent, number of revisits, number of animals, and so on). Clusters are updated with each time step to simulate real-time identification.

## Documentation

Input data to this MoveApp needs to have been processed by

1.  **Basic Data Processing for Merging Studies,** or a similar MoveApp, to append UTM data (i.e. columns named `x` and `y` containing UTM coordinate data)
2.  **Behavioural Classification for Avian Species,** or an alternative MoveApp, which appends a column of predicted behaviours to each event. Currently, the only supported behaviours are `SFeeding`, `SResting`, `SRoosting` and `STravelling`.

### Input data

*Indicate which type of input data the App requires. Currently only R objects of class `MoveStack` can be used. This will be extend in the future.*

*Example*: MoveStack in Movebank format

### Output data

*Indicate which type of output data the App produces to be passed on to subsequent apps. Currently only R objects of class `MoveStack` can be used. This will be extend in the future. In case the App does not pass on any data (e.g. a shiny visualization app), it can be also indicated here that no output is produced to be used in subsequent apps.*

*Example:* MoveStack in Movebank format

### Artefacts

*If the App creates artefacts (e.g. csv, pdf, jpeg, shapefiles, etc), please list them here and describe each.*

*Example:* `rest_overview.csv`: csv-file with Table of all rest site properties

### Settings

*Please list and define all settings/parameters that the App requires to be set by the App user, if necessary including their unit.*

*Example:* `Radius of resting site` (radius): Defined radius the animal has to stay in for a given duration of time for it to be considered resting site. Unit: `metres`.

### Most common errors

*Please describe shortly what most common errors of the App can be, how they occur and best ways of solving them.*

### Null or error handling

*Please indicate for each setting as well as the input data which behaviour the App is supposed to show in case of errors or NULL values/input. Please also add notes of possible errors that can happen if settings/parameters are improperly set and any other important information that you find the user should be aware of.*

*Example:* **Setting `radius`:** If no radius AND no duration are given, the input data set is returned with a warning. If no radius is given (NULL), but a duration is defined then a default radius of 1000m = 1km is set.
