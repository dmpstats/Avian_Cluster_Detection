library(move2)
library(withr)
library(dplyr)
library(rlang)
library(httr2)
library(readr)
library(lubridate)

if(rlang::is_interactive()){
  library(testthat)
  source("tests/app-testing-helpers.r")
  set_interactive_app_testing()
  app_key <- get_app_key()
}

test_sets <- httr2::secret_read_rds(test_path("data/vult_unit_test_data.rds"), key = I(app_key))
input3 <- read_rds(test_path("data/input3_move2.rds"))


test_that("input validation is doing it's job", {
  
  # missing 'behav' column
  expect_error(
    rFunction(input3), 
    "Classified behaviour column 'behav' is not contained in the input data"
  )
  
  # invalid date-time format on `clusterstart` and `clusterend`
  expect_error(
    rFunction(input3, 
              wholedata = FALSE, 
              clusterstart = "INVALID_DATETIME", 
              clusterend = "2014-01-01 00:00:00"),
    "App input `clusterstart` must be a date-time string in ISO 8601 format"
  )
  
  expect_error(
    rFunction(input3, 
              wholedata = FALSE, 
              clusterstart = 1331, 
              clusterend = "2014-01-01 00:00:00"),
    "App input `clusterstart` must be a date-time string in ISO 8601 format"
  )
  
  expect_error(
    rFunction(input3, 
              wholedata = FALSE, 
              clusterstart = "2006-12-22 22:07:04", 
              clusterend = "INVALID_DATETIME"),
    "App input `clusterend` must be a date-time string in ISO 8601 format"
  )
  
  expect_error(
    rFunction(input3, 
              wholedata = FALSE, 
              clusterstart = "2006-12-22 22:07:04", 
              clusterend = 44363),
    "App input `clusterend` must be a date-time string in ISO 8601 format"
  )
  
  # `clusterstart` and `clusterend` outside data range
  expect_error(
    rFunction(input3, 
              wholedata = FALSE, 
              clusterstart = "2016-12-22 22:07:04", 
              clusterend = "2006-12-22 22:07:04"),
    "App input `clusterstart` is outside the range of timepoints covered by the input data."
  )
  
  expect_error(
    rFunction(input3, 
              wholedata = FALSE, 
              clusterstart = "2006-12-22 22:07:04", 
              clusterend = "2016-12-22 22:07:04"),
    "App input `clusterend` is outside the range of timepoints covered by the input data."
  )
  
  # `clusterstart` and `clusterend` ignored when `wholedata` is TRUE
  expect_no_error(
    rFunction(test_sets$wcs |> dplyr::slice(1:10), 
              wholedata = TRUE, 
              clusterstart = "INVALID_DATETIME", 
              clusterend = "9999-12-31 22:07:04")
  )
  
  # `clusterstart` exceeds `clusterend`
  expect_error(
    rFunction(test_sets$wcs, 
              wholedata = FALSE, 
              clusterstart = max(test_sets$wcs$timestamp) - days(1), 
              clusterend = max(test_sets$wcs$timestamp) - days(2)), 
    "Clustering Start Date-time \\(`clusterstart`\\) '2024-03-13 22:00:00' is equal, or exceeds, the Clustering End Date-time "
  )
  
  # other inputs
  expect_error(
    rFunction(input3,clusterstep = -1),
    "App input `clusterstep` must be a positive integer"
  )
  
  expect_error(
    rFunction(input3, clusterwindow = 55.3),
    "App input `clusterwindow` must be a positive integer"
  )
  
  expect_error(
    rFunction(input3, clustexpiration = "NOW"),
    "App input `clustexpiration` must be a positive integer"
  )
  
  expect_error(
    rFunction(input3, d = FALSE),
    "App input `d` must be a positive integer"
  )
  
})



test_that("output is a valid move2 object", {
  actual <- rFunction(data = test_sets$wcs |> dplyr::slice(1:100), wholedata = TRUE)
  # passses {move2} check
  expect_true(move2::mt_is_move2(actual))
  # check if 1st class is "move2"
  expect_true(class(actual)[1] == "move2")
})



test_that("output always have column 'xy.clust', even when no clusters found", {
  actual <- rFunction(data = test_sets$wcs |> dplyr::slice(1:2), wholedata = TRUE)
  expect_true("xy.clust" %in% names(actual))
})




test_that("clustering on subset of data produces partially clustered output", {
  
  tmstmp_start <- floor_date(max(test_sets$nam$timestamp), unit = "10 hour") - days(8)
  tmstmp_end <- max(test_sets$nam$timestamp) - days(2) + hours(2) + minutes(10) - seconds(27)
  
  output <- rFunction(
    data = test_sets$nam, 
    wholedata = FALSE, 
    clusterstart = tmstmp_start,
    clusterend = tmstmp_end,
    clusterwindow = 4L, 
    clusterstep = 2L)
  
  # no clusters on data earlier than clusterstart
  expect_true(
    actual |> 
      filter(timestamp <  tmstmp_start) |> 
      pull(xy.clust) |> 
      is.na() |> 
      all()
    )
  
  # no clusters on data earlier than clusterend
  expect_true(
    actual |> 
      filter(timestamp >  tmstmp_end) |> 
      pull(xy.clust) |> 
      is.na() |> 
      all()
  )
  
  #  clusters on subset of data between clusterstart & clusterend
  expect_false(
    actual |> 
      filter(between(timestamp, tmstmp_start, tmstmp_end)) |> 
      pull(xy.clust) |> 
      is.na() |> 
      all()
  )
  
})





test_that("Expected clustering outcome has not changed", {
  
  testthat::local_edition(3)
  
  # tanzania
  expect_snapshot_value(
     rFunction(
      test_sets$ken_tnz |> dplyr::filter(timestamp > max(timestamp) - lubridate::days(2)),
      wholedata = TRUE
    ) |>
      pull(xy.clust) |>
      table(),
    style = "json2"
  )

  # Namibia
  expect_snapshot_value(
    rFunction(
      test_sets$nam |> dplyr::filter(timestamp > max(timestamp) - lubridate::days(10)),
      wholedata = TRUE,
      clusterwindow = 7L
    ) |>
      pull(xy.clust) |>
      table(),
    style = "json2"
  )
   
  
  # savahn
  expect_snapshot_value(
    rFunction(
      test_sets$savahn |> dplyr::filter(timestamp > max(timestamp) - lubridate::days(10))
    ) |>
      pull(xy.clust) |>
      table(),
    style = "json2"
  )


  # WCS
  expect_snapshot_value(
    rFunction(
      test_sets$wcs,
      clusterstep = 3L
    ) |> 
      pull(xy.clust) |> 
      table(),
    style = "json2"
  )


})