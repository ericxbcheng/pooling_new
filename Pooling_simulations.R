library(tidyverse)
library(knitr)
library(kableExtra)
library(raster)
library(primers)
library(mc2d)
source("STD_simulation.R")
source("1D_pooling_simulation.R")
source("2D_pooling_simulation.R")

# n = 48
# q_val = 7
# k_val = 4
# n_pos = 3
# 
# n_pos_vals = 1:47
# thresh = 20
# n_iter = 1000
# 
# # STD
# STD = read.csv("STD_48_7_4.csv", header = TRUE, row.name = 1, stringsAsFactors = FALSE) %>%
#   as.matrix()
# STD_mat = read.csv("STD_48_7_4 matrix.csv", header = TRUE, row.name = 1, stringsAsFactors = FALSE) %>%
#   as.matrix()
# 
# set.seed(123)
# STD_results = tune_n_pos(n_pos_vals = n_pos_vals, n_iter = n_iter, n = n, thresh = thresh, STD_mat = STD, q = q_val)
# 
# saveRDS(object = STD_results, file = "STD_48_7_4_1000.rds")
# 
# # 1D pooling
# # Pool by rows
# set.seed(123)
# pool_row = tune_1d_n_pos(n_pos_vals = n_pos_vals, n_iter = n_iter, n = n, thresh = thresh, n_pool = 6)
# saveRDS(object = pool_row, file = "Pool_row_48_1000.rds")
# 
# # Pool by columns
# set.seed(123)
# pool_col = tune_1d_n_pos(n_pos_vals = n_pos_vals, n_iter = n_iter, n = n, thresh = thresh, n_pool = 8)
# saveRDS(object = pool_col, file = "Pool_col_48_1000.rds")
# 
# # 2D Pooling
# set.seed(123)
# pool_2d = tune_2d_n_pos(n_pos_vals = n_pos_vals, n_iter = n_iter, n = n, thresh = thresh, n_row = 6, n_col = 8)
# saveRDS(object = pool_2d, file = "Pool_2d_48_1000.rds")

######################### 2/16/2019
## Note: Use modified PERT(0, 0.7, 19.99, 80) + Gamma(2, mode = 40000)

thresh = 20
n_iter = 10000

STD_48 = read.csv("STD_48_7_2_1pos.csv", header = TRUE, row.name = 1, stringsAsFactors = FALSE)
STD_48_mat = read.csv("STD_48_7_2_1pos_mat.csv", header = TRUE, row.name = 1, stringsAsFactors = FALSE) %>%
  as.matrix()

STD_96 = read.csv("STD_96_5_3_1pos.csv", header = TRUE, row.name = 1, stringsAsFactors = FALSE) 
STD_96_mat = read.csv("STD_96_5_3_1pos_mat.csv", header = TRUE, row.name = 1, stringsAsFactors = FALSE) %>%
  as.matrix()

set.seed(123)
results_48 = tune_n_pos(n_pos_vals = 1:(48-1), n_iter = n_iter, n = 48, thresh = thresh, STD_mat = STD_48)
saveRDS(object = results_48, file = "STD_48_7_2_1000.rds")
rm(results_48)

set.seed(123)
results_96 = tune_n_pos(n_pos_vals = 1:(96-1), n_iter = n_iter, n = 96, thresh = thresh, STD_mat = STD_96)
saveRDS(object = results_96, file = "STD_96_5_3_1000.rds")
rm(results_96)

### 48-well
# Pool by rows
set.seed(123)
pool_row_48 = tune_1d_n_pos(n_pos_vals = 1:(48-1), n_iter = n_iter, n = 48, thresh = thresh, by = "row")
saveRDS(object = pool_row_48, file = "1D_row_48_1000.rds")
rm(pool_row_48)

# Pool by columns
set.seed(123)
pool_col_48 = tune_1d_n_pos(n_pos_vals = 1:(48-1), n_iter = n_iter, n = 48, thresh = thresh, by = "column")
saveRDS(object = pool_col_48, file = "1D_col_48_1000.rds")
rm(pool_col_48)

### 96-well
# Pool by rows
set.seed(123)
pool_row_96 = tune_1d_n_pos(n_pos_vals = 1:(96-1), n_iter = n_iter, n = 96, thresh = thresh, by = "row")
saveRDS(object = pool_row_96, file = "1D_row_96_1000.rds")
rm(pool_row_96)

# Pool by columns
set.seed(123)
pool_col_96 = tune_1d_n_pos(n_pos_vals = 1:(96-1), n_iter = n_iter, n = 96, thresh = thresh, by = "column")
saveRDS(object = pool_col_96, file = "1D_col_96_1000.rds")
rm(pool_col_96)

# 2D pooling
# sample size = 48
set.seed(123)
twoD_48 = tune_2d_n_pos(n_pos_vals = 1:(48-1), n_iter = n_iter, n = 48, thresh = thresh)
saveRDS(object = twoD_48, file = "2D_48_1000.rds")
rm(twoD_48)

# sample size = 96
set.seed(123)
twoD_96 = tune_2d_n_pos(n_pos_vals = 1:(96-1), n_iter = n_iter, n = 96, thresh = thresh)
saveRDS(object = twoD_96, file = "2D_96_1000.rds")
rm(twoD_96)