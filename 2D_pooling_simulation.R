find_indi_pos = function(rows, cols){
  
  # Add the matrix for row pools and matrix for column pools
  # Convert it to a vector. We don't need to transpose the matrix because it should be consistent with how we generate conc
  # The index of each element indicates the sample ID (we count from left to right, from top to bottom on a 48-well plate)
  a = rows + cols %>%
    as.vector()
  
  b = ifelse(test = a == 2, yes = 1, no = 0) %>%
    factor(x = ., levels = c(0, 1))
  
  return(b)
}

calc_2d_metrics = function(conc, pool_r, pool_c, n, thresh){
  
  # Check the lengths of row pools and column pools
  stopifnot(length(pool_r) * length(pool_c) ==  n)
  
  # Remove the 0 concentration place holder
  a = conc[-1]
  
  # True AF status
  class_true = a >= thresh 
  class_true = class_true %>%
    as.numeric() %>%
    factor(x = ., levels = c(0, 1))
  
  # Generate the matrix of row pools and column pools
  mat_pool_r = find_pool_pos(c_pool = pool_r, n = n, by = "row", thresh = thresh)
  mat_pool_c = find_pool_pos(c_pool = pool_c, n = n, by = "column", thresh = thresh)
  
  # Predicted AF status
  class_pred = find_indi_pos(rows = mat_pool_r, cols = mat_pool_c)
  
  # Make a contingency table and calculate sensitivity and specificity
  cont_table = table(class_true, class_pred)
  
  sensi = cont_table[2,2] / (cont_table[2,2] + cont_table[2,1])
  speci = cont_table[1,1] / (cont_table[1,1] + cont_table[1,2])
  
  out = c(sensi, speci)
  
  return(out)
}

sim_2d_outcome = function(n, n_pos, thresh, lb, mode, tox){
  # Create a sequence of aflatoxin concentrations
  if(tox == "AF"){
    c_all = gen_elisa_af(n = n, n_pos = n_pos, lb = lb, mode = mode)
  } else if(tox == "FM"){
    c_all = gen_elisa_fm(n = n, n_pos = n_pos)
  }
  
  # Calculate the pooled concentrations
  pool_r = gen_1d_pool(conc = c_all, n = n, by = "row") 
  pool_c = gen_1d_pool(conc = c_all, n = n, by = "column")
  
  # Calculate sensitivity and specificity
  out = calc_2d_metrics(conc = c_all, pool_r = pool_r, pool_c = pool_c, n = n, thresh = thresh)
  out2 = c(out, n_pos)
  
  return(out2)
}

# Iterate once
gen_sim_2d_outcome = function(n, n_pos, thresh, lb, mode, tox){
  function(...){
    sim_2d_outcome(n = n, n_pos = n_pos, thresh = thresh, lb = lb, mode = mode, tox = tox)
  }
}

sim_2d_iterate = function(n_iter, ...){
  
  f_outcome = gen_sim_2d_outcome(...)
  a = map(.x = 1:n_iter, .f = f_outcome)
  b = clean(a)
  
  return(b)
}

tune_2d_n_pos = function(n_pos_vals, n_iter, n, thresh, lb, mode, tox){
  map(.x = n_pos_vals, .f = sim_2d_iterate, n_iter = n_iter, n = n, thresh = thresh, lb = lb, mode = mode, tox = tox)
}

calc_2d_cost = function(data, n){
  
  stopifnot(n %in% c(48, 96))
  
  if(n == 48){
    nrow = 6
    ncol = 8
  } else {
    nrow = 8
    ncol = 12
  }
  
  a = do.call(what = rbind.data.frame, args = data) %>%
    mutate(putative_pos = sensi*n_pos + (1-speci)*(n-n_pos),
           n_test_total = nrow + ncol + putative_pos,
           n_pipette = n * 2 + nrow + ncol + putative_pos)
  return(a)
}
