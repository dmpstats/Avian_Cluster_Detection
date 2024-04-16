# ------------------------- #
#         Preamble
# ------------------------- #

library(move2)
library(httr2)
library(purrr)
library(readr)

# Helpers
source("tests/app-testing-helpers.r")

# # get App secret key for decrypting test dataset
# app_key <- get_app_key()
# 
# 
# # Read (encrypted) input datasets for testing
# test_dt <- secret_read_rds("data/raw/vult_test_data.rds", key = I(app_key))
# 
# 
# # thinning gaia and nam dataset for 5mins gap for faster testing.
# test_dt$gaia <- mt_filter_per_interval(test_dt$gaia, unit = "2 min")
# test_dt$nam_sop <- mt_filter_per_interval(test_dt$nam_sop, unit = "2 min")
# test_dt$metadata


dt_sa <- read_rds("data/raw/SA_VfA_MPIAB_processed_acc_merged_classified.rds")
dt_gaia <- read_rds("data/raw/gaia_processed_acc_merged_classified.rds")



# ---------------------------------------- #
# ----   Interactive RFunction testing  ----
# ---------------------------------------- #
set_interactive_app_testing()

out_dt_gaia <- rFunction(
  data = dt_gaia,
  clusterstep = 2, 
  clusterwindow = 7, 
  clustexpiration = 14, 
  behavsystem = TRUE, 
  clustercode = "A")


out_dt_gaia



out_dt_sa <- rFunction(
  data = dt_sa,
  clusterstep = 2, 
  clusterwindow = 7, 
  clustexpiration = 14, 
  behavsystem = TRUE, 
  clustercode = "A")


out_dt_$xy.clust





mt_track_id_column(out_dt_sa)
mt_time_column(
  out_dt_sa)

attributes(out_dt_sa)



cluster_tbl <- read_rds("data/output/clustertable.rds")
cluster_tbl

mt_track_id_column(cluster_tbl)
mt_time_column(cluster_tbl)




attr(out_dt_sa, "cluster_tbl") <- cluster_tbl

attributes(out_dt_sa)
attr(out_dt_sa, "cluster_tbl")
