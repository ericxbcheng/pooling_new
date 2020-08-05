library(tidyverse)
source("Pooling_covid_functions.R")


n_iter = 10
Args_default = list(n = 20, n_pos = 2, pool_size = 5, thresh = 0)

sim_result = sim_iterate_covid(n_iter = n_iter, Args = Args_default)
sim_result_cleaned = clean_covid(list = sim_result)

test1 = tune_param_covid(Args = Args_default, n_iter = n_iter, param = "n_pos", val = 5)
test2 = tune_param_n_covid(vals = c(1,3,5), Args = Args_default, n_iter = 10, var_prim = "n_pos")
test3 = tune_param_sec_covid(Args = Args_default, var_prim = "n_pos", vals_prim = c(1,3,5), 
                     var_sec = "pool_size", n_iter = 10, vals_sec = c(2,4,6))

test4 = metrics_sec_covid(data = test3, vals_prim = c(1,3,5), vals_sec = c(2,4,6), n_iter = 10)


