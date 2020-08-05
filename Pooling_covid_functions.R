make_samples = function(n, pos){
  
  stopifnot(n >= pos)
  
  a = rep.int(x = 1, times = pos)
  b = vector(mode = "numeric", length = n - pos)
  
  return(sample(x = c(a,b), size = n, replace = FALSE))
}

make_pool_idx = function(n, pool_size){
  
  quotient = n %/% pool_size
  modulo = n %% pool_size
  
  if(modulo == 0){
    a = rep.int(x = 1:quotient, times = pool_size)
    return(a)
    
  } else {
    a = rep.int(x = 1:quotient, times = pool_size)
    b = seq.int(from = quotient + 1, to = quotient + modulo, by = 1)
    return(c(a,b))
  }
}

find_pool_pos = function(conc, pool_idx, thresh){
  
  # Calculate pooled concentration
  pools = split(x = conc, f = pool_idx) %>%
    map_dbl(.x = ., .f = mean)
  
  # Determine positive or negative pool
  return(pools > thresh)
}

calc_1d_metrics_covid = function(n, conc, pool_pos, pool_idx){
  
  class_true = {conc == 1}
  class_pred = ifelse(test = pool_idx %in% which(pool_pos), yes = TRUE, no = FALSE)
  
  # Make a contingency table and calculate sensitivity and specificity
  cont_table = table(class_true, class_pred)
  
  sensi = cont_table[2,2] / (cont_table[2,2] + cont_table[2,1])
  speci = cont_table[1,1] / (cont_table[1,1] + cont_table[1,2])
  
  out = c(sensi, speci)
  
  return(out)
}
