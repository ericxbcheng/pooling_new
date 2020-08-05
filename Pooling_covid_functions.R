make_samples = function(n, n_pos){
  
  stopifnot(n >= n_pos)
  
  a = rep.int(x = 1, times = n_pos)
  b = vector(mode = "numeric", length = n - n_pos)
  
  return(sample(x = c(a,b), size = n, replace = FALSE))
}

make_pool_idx = function(n, pool_size, for2D = FALSE){
  
  quotient = n %/% pool_size
  modulo = n %% pool_size
  
  if(modulo == 0){
    a = rep.int(x = 1:quotient, times = pool_size)
    return(a)
    
  } else {
    
    if(for2D == FALSE){
      a = rep.int(x = 1:quotient, times = pool_size)
    } else {
      a = rep(x = 1:quotient, each = pool_size)
    }
    
    b = seq.int(from = quotient + 1, to = quotient + modulo, by = 1)
    return(c(a,b))
  }
}

find_pool_pos_covid = function(conc, pool_idx, thresh){
  
  # Calculate pooled concentration
  pools = split(x = conc, f = pool_idx) %>%
    map_dbl(.x = ., .f = mean)
  
  # Determine positive or negative pool
  return(pools > thresh)
}


find_pool_pos_2d_covid = function(conc, n, pool_size, thresh, message = FALSE){
  
  #get number of sets to compare
  n_sets <- n %/% pool_size^2
  n_indiv = n %% pool_size^2
  
  if(message == TRUE){
    if(n %% pool_size^2 != 0) print("need to add results in for remainder")
  }
  
  
}

compare_2D_pools_ec = function(n, pool_size){
  
  #get number of sets to compare
  n_sets <- n %/% pool_size^2
  if(n %% pool_size^2 != 0) print("need to add results in for remainder")
  #compare by sets
  status<-n
  setID = rep(x = 1:n_sets, each = pool_size)
  for( i in 1:unique(setID) ) {
    d <- row_pools[i*(1:pool_size)] %o% col_pools[i*(1:pool_size)]
    status[(1:(pool_size)^2)*i] <- as.vector(d)
  }
  return(status)
}

calc_1d_metrics_covid = function(n, conc, pool_pos, pool_idx){
  
  class_true = {conc == 1} %>%
    as.numeric() %>%
    factor(x = ., levels = c(0, 1))
  
  class_pred = ifelse(test = pool_idx %in% which(pool_pos), yes = TRUE, no = FALSE) %>%
    as.numeric() %>%
    factor(x = ., levels = c(0, 1))
  
  # Make a contingency table and calculate sensitivity and specificity
  cont_table = table(class_true, class_pred)
  
  sensi = cont_table[2,2] / (cont_table[2,2] + cont_table[2,1])
  speci = cont_table[1,1] / (cont_table[1,1] + cont_table[1,2])
  
  out = c(sensi, speci)
  
  return(out)
}

sim_1d_outcome_covid = function(n, n_pos, pool_size, thresh = 0){
  
  stopifnot(n > n_pos | n > 0 | n_pos >= 0)
  
  conc = make_samples(n = n, n_pos = n_pos)
  
  pool_idx = make_pool_idx(n = n, pool_size = pool_size)
  
  pool_pos = find_pool_pos_covid(conc = conc, pool_idx = pool_idx, thresh = thresh)
  
  # Calculate sensitivity and specificity
  out = calc_1d_metrics_covid(n = n, conc = conc, pool_pos = pool_pos, pool_idx = pool_idx)
  
  return(out)
}

# Iterate once
gen_sim_1d_outcome_covid = function(n, n_pos, thresh, pool_size){
  function(...){
    sim_1d_outcome_covid(n = n, n_pos = n_pos, thresh = thresh, pool_size = pool_size)
  }
}

# First layer of iteration
sim_iterate_covid = function(n_iter, Args){
  
  # Check point: Is n_iter >= 1?
  stopifnot(n_iter >= 1)
  
  # Generate a sim_outcome_new() with loaded arguments
  a = do.call(what = gen_sim_1d_outcome_covid, args = Args)
  
  # Iterate that manufactured function for n_iter times
  b = map(.x = 1:n_iter, .f = a)
  
  return(b)
}

clean_covid = function(list){
  a = unlist(list)
  
  # Find the indices of each type of metric
  sens_ind = seq(from = 1, by = 2, length.out = length(a)/2)
  spec_ind = sens_ind + 1
  
  #Rearrange the results
  b = list("sensi" = a[sens_ind], "speci" = a[spec_ind])
  
  return(b)
}

# This function makes a list of arguments from the default arguments
update_arg_covid = function(Args, param, val){
  
  # Checkpoint
  stopifnot(param %in% names(Args))
  
  Args[[param]] = val
  return(Args)
}

# Create a function that passes one tuning parameter value to sim_iterate2()
tune_param_covid = function(Args, n_iter, param, val, ...){
  
  # Get tuning parameters
  a = update_arg_covid(Args = Args, param = param, val = val)
  
  b = sim_iterate_covid(n_iter = n_iter, Args = a) %>%
    clean_covid(list = .)
  
  # Add a vector of tuning parameter value
  b$param = rep.int(x = val, times = n_iter)
  
  return(b)
}

# A function that tunes over many values of the primary tuning parameter
tune_param_n_covid = function(vals, Args, n_iter, var_prim){
  return(map(.x = vals, .f = tune_param_covid, Args = Args, n_iter = n_iter, param = var_prim))
}


# A function that tunes over a crossed combination of the primary and secondary tuning parameter
tune_param_sec_covid = function(Args, var_prim, vals_prim, var_sec, n_iter, vals_sec){
  
  # Create a list of argument lists, each argument list corresponding to one secondary tuning value
  Args_sec = map(.x = vals_sec, .f = update_arg_covid, Args = Args, param = var_sec)
  
  # For each argument list in Args_sec, do iterate_tune1()
  sim_data = list()
  for(i in 1:length(vals_sec)){
    
    sim_data[[i]] = map(.x = vals_prim, .f = tune_param_covid, 
                        Args = Args_sec[[i]], n_iter = n_iter,param = var_prim)
    
  }
  return(sim_data)
}

metrics_sec_covid = function(data, vals_prim, vals_sec, n_iter){
  
  a = map(.x = data, .f = bind_rows) 
  b = bind_rows(a)
  param2 = rep(x = vals_sec, each = length(vals_prim) * n_iter)
  
  return(data.frame(b, param2))
}

# # Calculate the cost of reagents and pipetting
# calc_1d_cost = function(data, n, by){
#   
#   stopifnot(by %in% c("row", "column"))
#   stopifnot(n %in% c(48, 96))
#   
#   if(n == 48){
#     n_row = 6
#     n_col = 8
#   } else if(n == 96){
#     n_row = 8
#     n_col = 12
#   } else {
#     stop("Undefined n. Choose 48 or 96.")
#   }
#   
#   if(by == "row"){
#     a = do.call(what = rbind.data.frame, args = data) %>%
#       mutate(putative_pos = sensi*n_pos + (1-speci)*(n-n_pos),
#              n_test_total = n_row + putative_pos,
#              n_pipette = n + n_row + putative_pos)
#     return(a)
#     
#   } else if(by == "column"){
#     b = do.call(what = rbind.data.frame, args = data) %>%
#       mutate(putative_pos = sensi*n_pos + (1-speci)*(n-n_pos),
#              n_test_total = n_col + putative_pos,
#              n_pipette = n + n_col + putative_pos)
#     return(b)
#     
#   }
# }

plot_tune2_ribbon_covid = function(data, xlab, legend_lab, xtick){
  
  # Summarise the data
  a = data %>%
    gather(data = ., key = "metric", value = "value", -c(param, param2)) %>%
    group_by(param2, param, metric) %>%
    summarise(lb = quantile(x = value, probs = 0.025), 
              med = median(x = value),
              ub = quantile(x = value, probs = 0.975)) 
  
  # Visualize
  b = ggplot(data = a) +
    geom_ribbon(aes(x = param, ymin = lb, ymax = ub, group = as.factor(metric), fill = as.factor(metric)), alpha = 0.3) +
    geom_line(aes(x = param, y = med, color = as.factor(metric))) +
    geom_point(aes(x = param, y = med, color = as.factor(metric))) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1)) +
    scale_fill_discrete(name = legend_lab) +
    scale_color_discrete(name = legend_lab) +
    scale_x_continuous(breaks = unique(a$param), labels = xtick) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = xlab, y = "Metrics (2.5th - 97.5th percentile)") +
    theme_bw() +
    theme(legend.position = "top")+
    facet_wrap(~param2)
  
  return(b)
}
