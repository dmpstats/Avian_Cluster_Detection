# ------------------------- #
#         Preamble
# ------------------------- #

library(move2)
library(httr2)
library(purrr)
library(readr)
library(sf)

# Helpers
source("tests/app-testing-helpers.r")

# get App secret key for decrypting test dataset
app_key <- get_app_key()
 
# Read (encrypted) input datasets for testing
test_dt <- httr2::secret_read_rds("data/raw/vult_test_data.rds", key = I(app_key))
#map(test_dt, nrow)


set_interactive_app_testing()



# ---------------------------------------- #
# ----    Automated Unit testing        ----
# ---------------------------------------- #

testthat::test_file("tests/testthat/test_RFunction.R")


# ---------------------------------------- #
# ----   Interactive RFunction testing  ----
# ---------------------------------------- #


out_dt_nam <- rFunction(
  data = test_dt$nam, 
  clusterstart = "2024-03-05 18:00:09",
  clusterend = "2024-03-12 18:00:09",
  clusterwindow = 10L, 
  clusterstep = 2L)



out_dt_nam <- rFunction(
  data = test_dt$nam, 
  clusterstart = min(test_dt$nam$timestamp),
  clusterend = "2024-03-12 18:00:09",
  clusterwindow = as.integer(ceiling(diff(range(test_dt$nam$timestamp), units = "days"))), 
  clusterstep = 2L)


#' --------------------------
#' Example of input data with locations not in UTM

out_dt_nam <- rFunction(
  data = test_dt$nam, 
  clusterwindow = 10L, 
  clusterstep = 2L)

out_dt_nam_latlon <- rFunction(
  data = sf::st_transform(test_dt$nam, crs = 4326), 
  clusterwindow = 10L, 
  clusterstep = 2L)

identical(out_dt_nam_latlon$clust_id, out_dt_nam$clust_id)



#' --------------------------
#' example of case 1 of iterative cluster merging: Multiple old clusters -> 1 new cluster
#' 
#' Note: Uncomment `browser()` in associated section
#' 
out_dt_savahn <- rFunction(
  data = test_dt$savahn, 
  clusterwindow = 5L, 
  clusterstep = 1L)



#' --------------------------
#' example of case 2 of iterative cluster merging: 1 old clusters -> multiple new clusters
#' Note: Uncomment `browser()` in associated section
#' 
out <- rFunction(
  data = test_dt$ken_tnz, 
  clusterwindow = 4L, 
  clusterstep = 1L)





#' --------------------------
#' Testing no rolling window
#' 
out_no_rolling <- rFunction(
  data = test_dt$savahn,
  clusterstep = 2L, 
  clusterwindow = as.integer(ceiling(diff(range(test_dt$savahn$timestamp), units = "days"))), 
  clustexpiration = 8L, 
  behavsystem = TRUE)

#' Note how cluster expiration has no effect on no rolling window, as maximum
#' time lag between consecutive elements of some of the clusters is larger than
#' the specified 8 days
out_no_rolling |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  filter(max_lag > 8) |> 
  print(n = 50)





#' --------------------------
#' Demonstrating that expiration cut-off shouldn't be less than length of window
#' 
out_exp_lessthan_wdw <- rFunction(
  data = test_dt$savahn,
  clusterstep = 2L, 
  clusterwindow = 10L, 
  clustexpiration = 5L, 
  behavsystem = TRUE)


out_exp_lessthan_wdw |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  filter(max_lag > 5) |> 
  print(n = 50)



#' --------------------------
#' Demonstrating that, even if expiration is larger than window, we can still
#' end up with clusters whose point elements are separated by more than the
#' expiration cut-off
#' 
out_exp_morethan_wdw <- rFunction(
  data = test_dt$savahn,
  clusterstep = 2L, 
  clusterwindow = 10L, 
  clustexpiration = 12L, 
  behavsystem = TRUE)


out_exp_morethan_wdw |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  filter(max_lag > 12) 



#' --------------------------
#' ... But, if the rolling window time-step is increased, then output respects
#' the expiration cut-off!
#' 
out_exp_morethan_wdw_lrgstep <- rFunction(
  data = test_dt$savahn,
  clusterstep = 5L, 
  clusterwindow = 10L, 
  clustexpiration = 12L, 
  behavsystem = TRUE)


out_exp_morethan_wdw_lrgstep |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  filter(max_lag > 12)




#' --------------------------------------------------
#' test marginal effect of increasing clustering window

# 4 days
out_dt_wndw4 <- rFunction(
  data = test_dt$wcs,
  clusterstep = 2L, 
  clusterwindow = 4L, 
  clustercode = "A")

# nr clusters
length(unique(out_dt_wndw4$clust_id))

# max time lag between point elements of each cluster
out_dt_wndw4 |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)
            
out_dt_wndw4 |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)



# 12 days
out_dt_wndw12 <- rFunction(
  data = test_dt$wcs,
  clusterstep = 2L, 
  clusterwindow = 12L, 
  clustercode = "A")

length(unique(out_dt_wndw12$clust_id))

out_dt_wndw12 |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)

out_dt_wndw12 |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)


# CONCLUSION: perhaps potentially larger lags between point elements of clusters


#' --------------------------------------------------
#' test marginal effect of increasing cluster step - faster processing, but at what cost?

out_dt_stp1 <- rFunction(
  data = test_dt$ken_tnz,
  clusterstep = 1L, 
  clustercode = "A")

length(unique(out_dt_stp1$clust_id))

# 1-day time-step
out_dt_stp1 |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)

out_dt_stp1 |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)


# 5-day time-step
out_dt_stp5 <- rFunction(
  data = test_dt$ken_tnz,
  clusterstep = 5L, 
  clustercode = "A")

length(unique(out_dt_stp5$clust_id))

out_dt_stp5 |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)

out_dt_stp5 |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)


# CONCLUSION: No substancial difference



#' --------------------------------------------------
#' test marginal effect of increasing parameter d
#' 

out_dt_d500 <- rFunction(
  data = test_dt$nam,
  d = 500L,
  clustercode = "A")


length(unique(out_dt_d500$clust_id))

out_dt_d500 |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)

out_dt_d500 |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)



out_dt_d100 <- rFunction(
  data = test_dt$nam,
  d = 100L,
  clustercode = "A")


length(unique(out_dt_d100$clust_id))

out_dt_d100 |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)

out_dt_d100 |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)





#' -------------------------------------------------------
#' test marginal effect of matching proximity threshold
#' 

#' 175 meters (default)
out_1 <- rFunction(
  data = test_dt$nam,
  match_thresh = 175L,
  clustercode = "A")

table(out_1$clust_id)

#' 100 meters (default)
out_2 <- rFunction(
  data = test_dt$nam,
  match_thresh = 100L,
  clustercode = "A")

table(out_2$clust_id)





#' --------------------------------------------------
#'  Varied datasets with assorted reasonable settings
#'  
  
#' - nam
out_dt_nam <- rFunction(
  data = test_dt$nam,
  clusterstart = max(mt_time(test_dt$nam)) - days(5),
  clusterend = max(mt_time(test_dt$nam)) - days(1),
  clusterstep = 6L, 
  clusterwindow = 10L, 
  clustexpiration = 14L, 
  behavsystem = TRUE, 
  clustercode = "A")


out_dt_nam |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)


out_dt_nam |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)




#' - gaia
out_dt_gaia <- rFunction(
  data = test_dt$gaia,
  clusterstep = 1L, 
  clusterwindow = 8L, 
  clustexpiration = 14L, 
  behavsystem = TRUE, 
  d = 400L,
  clustercode = "gaia")

out_dt_gaia |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)


out_dt_gaia |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)



#' - sa_vfa
out_dt_sa_vfa <- rFunction(
  data = test_dt$sa_vfa,
  clusterstep = 3L, 
  clusterwindow = 10L, 
  clustexpiration = 14L, 
  behavsystem = TRUE, 
  clustercode = "sa_vfa")

out_dt_sa_vfa |> 
  tidyr::drop_na(clust_id) |> 
  group_split(clust_id) |> 
  map(~max(st_distance(.x))) |> 
  purrr::list_c() |> 
  sort(decreasing = TRUE)


out_dt_sa_vfa |> 
  tidyr::drop_na(clust_id) |> 
  group_by(clust_id) |> 
  arrange(timestamp) |>  
  mutate( time_lag = as.numeric(lead(timestamp) - timestamp, units = "days")) |>
  summarise(max_lag = max(time_lag, na.rm = TRUE)) |> 
  arrange(desc(max_lag)) |> 
  print(n = 20)



#' - savanha
out_dt_savah <- rFunction(
  data = test_dt$savahn,
  clusterstep = 3L, 
  clusterwindow = 10L, 
  clustexpiration = 14L, 
  behavsystem = TRUE, 
  clustercode = "savah")




# ---------------------------------------- #
# ----    MoveApps SDK testing          ----
# ---------------------------------------- #

# standard dataset with default inputs
run_sdk(
  data = test_dt$sa_vfa,
  clusterstep = 1L, 
  clusterwindow = 7L, 
  clustexpiration = 14L, 
  behavsystem = TRUE, 
  d = 500L,
  match_thresh = 100,
  clustercode = "A")
  
output <- readRDS("data/output/output.rds"); output
table(output$clust_id)



# no need to input strict integers 
run_sdk(
  data = test_dt$gaia,
  clusterstep = 1, 
  clusterwindow = 8, 
  clustexpiration = 14, 
  behavsystem = TRUE, 
  d = 500L,
  clustercode = "A")

output <- readRDS("data/output/output.rds"); output
table(output$clust_id)



# pass on strings to clusterstart/clusterend
run_sdk(
  data = test_dt$gaia,
  clusterstart = "2024-03-02 23:00:00",
  clusterstep = 1, 
  clusterwindow = 8, 
  clustexpiration = 14, 
  behavsystem = TRUE, 
  d = 500L,
  clustercode = "A")




