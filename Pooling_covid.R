library(tidyverse)
source("Pooling_covid_functions.R")

n_iter = 10
n = 20000
n_pos = 1
pool_size = 5
thresh = 0

# Primary tuning param: n_pos
prev = seq(from = 0.001, to = 0.01, by = 0.001)
vals_prim = prev * n

# Secondary tuning param: pool_size
vals_sec = 1:10

# Simulation
Args_default = list(n = n, n_pos = n_pos, pool_size = pool_size, thresh = thresh)  

sim_data = tune_param_sec_covid(Args = Args_default, var_prim = "n_pos", vals_prim = vals_prim, 
                                var_sec = "pool_size", n_iter = n_iter, vals_sec = vals_sec)
saveRDS(object = sim_data, file = sprintf(fmt = "covid_sim_%d_0.001to0.01_1to10_%d.rds", n, n_iter))


sim_data_cleaned = metrics_sec_covid(data = sim_data, vals_prim = vals_prim, vals_sec = vals_sec, n_iter = n_iter)  
plot_tune2_ribbon_covid(data = sim_data_cleaned, xlab = "Prevalence (%)", legend_lab = "Type", xtick = prev * 100)


######################## Debug #####################

sim_result = sim_iterate_covid(n_iter = n_iter, Args = Args_default)
sim_result_cleaned = clean_covid(list = sim_result)

test1 = tune_param_covid(Args = Args_default, n_iter = n_iter, param = "n_pos", val = 5)
test2 = tune_param_n_covid(vals = c(1,3,5), Args = Args_default, n_iter = 10, var_prim = "n_pos")
test3 = tune_param_sec_covid(Args = Args_default, var_prim = "n_pos", vals_prim = c(1,3,5), 
                             var_sec = "pool_size", n_iter = 10, vals_sec = c(2,4,6))

test4 = metrics_sec_covid(data = test3, vals_prim = c(1,3,5), vals_sec = c(2,4,6), n_iter = 10)


