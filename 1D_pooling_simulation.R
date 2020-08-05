
## For simulated data
# gen_1d_pool = function(conc, n, n_pool){
#   # Remove the 0 concentration place holder
#   a = conc[-1]
# 
#   if(n_pool == 6){
#     b = matrix(data = a, nrow = n_pool, ncol = n/n_pool, byrow = FALSE) %>%
#       as.data.frame()
#     c = rowSums(b)/ncol(b)
#     return(c)
# 
#   }else if(n_pool == 8){
#     d = matrix(data = a, nrow = n/n_pool, ncol = n_pool, byrow = FALSE) %>%
#       as.data.frame()
#     e = colSums(d)/nrow(d)
#     return(e)
# 
#   } else {
#     stop("Unknown pooling method. Try either '6' or '8'.")
#   }
# }

gen_1d_pool = function(conc, n, by){
  
  stopifnot(by %in% c("row", "column"))
  stopifnot(n %in% c(48, 96))
  
  if(n == 48){
    n_row = 6
    n_col = 8
  } else if(n == 96){
    n_row = 8
    n_col = 12
  } else {
    stop("Undefined n. Choose 48 or 96.")
  }
  
  # Remove the 0 concentration place holder
  a = conc[-1]
  
  b = matrix(data = a, nrow = n_row, ncol = n_col, byrow = FALSE) %>%
    as.data.frame()
  
  if(by == "row"){
    c = rowMeans(b)
    return(c)
  } else if(by == "column"){
    c = colMeans(b)
    return(c)
  }
}

# find_pool_pos = function(c_pool, n, thresh) {
#   
#   a = matrix(data = 0, nrow = 6, ncol = 8)
#   n_pool = length(c_pool)
#   
#   # Pool by rows: n_pool = 6; Pool by columns: n_pool = 8
#   if (n_pool == 6) {
#     a[c_pool >= (thresh / (n / n_pool)),] = 1
#     
#   } else if (n_pool == 8) {
#     a[, c_pool >= (thresh / (n / n_pool))] = 1
#     
#   } else {
#     stop("Unknown pooling method. Please check the length of c_pool.")
#   }
#   
#   return(a)
# }

find_pool_pos = function(c_pool, n, by, thresh) {
  
  stopifnot(by %in% c("row", "column"))
  stopifnot(n %in% c(48, 96))
  
  if(n == 48){
    n_row = 6
    n_col = 8
  } else if(n == 96){
    n_row = 8
    n_col = 12
  } else {
    stop("Undefined n. Choose 48 or 96.")
  }
  
  a = matrix(data = 0, nrow = n_row, ncol = n_col)
  n_pool = length(c_pool)
  pos_pool = c_pool >= (thresh / (n / n_pool))
  
  if (by == "row") {
    a[pos_pool, ] = 1
    
  } else if (by == "column") {
    a[ ,pos_pool] = 1
    
  } else {
    stop("Unknown pooling method. Choose 'row' or 'column'.")
  }
  
  return(a)
}


calc_1d_metrics = function(n, thresh, conc, c_pool, by){
  # Remove the 0 concentration place holder
  a = conc[-1]
  
  # Find the kernels >= threshold
  class_true = a >= thresh 
  class_true = class_true %>%
    as.numeric() %>%
    factor(x = ., levels = c(0, 1))
  
  # Find positive pools
  b = find_pool_pos(c_pool = c_pool, n = n, by = by, thresh = thresh)
  
  # Find the predicted class
  class_pred = as.vector(b) %>%
    as.numeric() %>%
    factor(x = ., levels = c(0,1))
  
  # Make a contingency table and calculate sensitivity and specificity
  cont_table = table(class_true, class_pred)
  
  sensi = cont_table[2,2] / (cont_table[2,2] + cont_table[2,1])
  speci = cont_table[1,1] / (cont_table[1,1] + cont_table[1,2])
  
  out = c(sensi, speci)
  
  return(out)
}

sim_1d_outcome = function(n, n_pos, thresh, by, lb, mode, tox){
  # Create a sequence of aflatoxin concentrations
  if(tox == "AF"){
    c_all = gen_elisa_af(n = n, n_pos = n_pos, lb = lb, mode = mode)
  } else if (tox == "FM"){
    c_all = gen_elisa_fm(n = n, n_pos = n_pos)
  }
  
  # Calculate the pooled concentrations
  c_pool = gen_1d_pool(conc = c_all, n = n, by = by)
  
  # Calculate sensitivity and specificity
  out = calc_1d_metrics(n = n, thresh = thresh, conc = c_all, c_pool = c_pool, by =  by)
  out2 = c(out, n_pos)
  
  return(out2)
}

# Iterate once
gen_sim_1d_outcome = function(n, n_pos, thresh, by, lb, mode, tox){
  function(...){
    sim_1d_outcome(n = n, n_pos = n_pos, thresh = thresh, by = by, lb = lb, mode = mode, tox = tox)
  }
}

sim_1d_iterate = function(n_iter, ...){
  
  f_outcome = gen_sim_1d_outcome(...)
  a = map(.x = 1:n_iter, .f = f_outcome)
  b = clean(a)
  
  return(b)
}

tune_1d_n_pos = function(n_pos_vals, n_iter, n, thresh, by, lb, mode, tox){
  map(.x = n_pos_vals, .f = sim_1d_iterate, n_iter = n_iter, n = n, thresh = thresh, by = by, lb = lb, mode = mode, tox = tox)
}

# Calculate the cost of reagents and pipetting
calc_1d_cost = function(data, n, by){
  
  stopifnot(by %in% c("row", "column"))
  stopifnot(n %in% c(48, 96))
  
  if(n == 48){
    n_row = 6
    n_col = 8
  } else if(n == 96){
    n_row = 8
    n_col = 12
  } else {
    stop("Undefined n. Choose 48 or 96.")
  }
  
  if(by == "row"){
    a = do.call(what = rbind.data.frame, args = data) %>%
      mutate(putative_pos = sensi*n_pos + (1-speci)*(n-n_pos),
             n_test_total = n_row + putative_pos,
             n_pipette = n + n_row + putative_pos)
    return(a)
    
  } else if(by == "column"){
    b = do.call(what = rbind.data.frame, args = data) %>%
      mutate(putative_pos = sensi*n_pos + (1-speci)*(n-n_pos),
             n_test_total = n_col + putative_pos,
             n_pipette = n + n_col + putative_pos)
    return(b)
    
  }
}

