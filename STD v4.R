## Update: wrap things into functions. 

## Create functions
sigmaq <- function(vector, j, i, gamma, q){
  #Calculate s_ij
  if (j < q) {
    s_ij <- 0
    for (c in 0:gamma) {
      s_ij <- s_ij + (j ^ c) * floor(i / (q ^ c))
    }
  } else if (j == q) {
    s_ij <- floor(i / (q ^ gamma))
  }
  #Calculate the number of places to shift vector
  move <- s_ij %% q
  binhf::shift(vector, move, dir = "right")
}

layer = function(vector, k, gamma, q, n){
  sapply(X = 0:(n-1), FUN = function(x){
    sigmaq(vector = vector, j = k, i = x, gamma = gamma, q = q)
  })
}

get_param = function(n, q_list, k_list, gamma_list, criteria){
  data.frame(q = q_list[criteria],
             gamma = gamma_list[criteria],
             k = k_list[criteria]) %>%
    mutate(t1 = q * k,
           t2 = q * (q + 1),
           t3 = q ^ 2 + ceiling((n - 1)/(q ^ gamma)) + 1)
}

#Inequality tests for determining which construction to use
get_construction = function(param, n, m){
  
  a = {param$k < param$q + 1} & 
    {{n/param$q < m} | {ceiling(n/param$q) == m}}
  b = {{param$k == param$q + 1} & 
    {n == param$q^(param$gamma+1) - 1}} & 
    {{param$q^param$gamma < m} | {param$q^param$gamma == m}}
  c = {{param$k == param$q + 1} & 
    {ceiling((n-1)/param$q^param$gamma) < param$q-1} & 
    {{param$q^param$gamma < m} | {param$q^param$gamma == m}}}
  
  return(list(a = a, b = b, c = c))
}

show_construction = function(param, criteria){
  
  cstr1 = param[criteria$a, ] %>% dplyr::select(q, gamma, k, t1) %>% rename(t = t1)
  cstr2 = param[criteria$b, ] %>% dplyr::select(q, gamma, k, t2) %>% rename(t = t2)
  cstr3 = param[criteria$c, ] %>% dplyr::select(q, gamma, k, t3) %>% rename(t = t3)
  
  return(list(cstr1 = cstr1, cstr2 = cstr2, cstr3 = cstr3))
}

get_best_cstr = function(list){
  a = bind_rows(list) %>%
    dplyr::filter(t == min(t)) %>%
    unlist(use.names = TRUE)
  
  return(a)
}

#best_cstr = get_best_cstr(list = cstr_list)

get_matrix = function(cstr, n){
  q_val = cstr[["q"]]
  k_val = cstr[["k"]]
  gamma_val = cstr[["gamma"]]
  
  C00 = c(1, rep(0, q_val - 1))
  
  layer_list = lapply(X = 0:(k_val -1), FUN = function(x){layer(C00, x, gamma_val, q_val, n)})
  M = as.matrix(do.call(rbind.data.frame, layer_list))
  
  colnames(M) = map2_chr(.x = rep("Sample ", times = ncol(M)), .y = seq(1:ncol(M)), .f = paste0)
  rownames(M) = map2_chr(.x = rep("Pool ", times = nrow(M)), .y = seq(1:nrow(M)), .f = paste0)
  
  return(M)
}

draw_STD = function(data, n){
  
  cstr = data$construction
  matrix = data$matrix
  
  q_val = cstr[["q"]]
  k_val = cstr[["k"]]
  gamma_val = cstr[["gamma"]]
  
  STD_scheme = raster(x = matrix, xmn = 0, xmx = n, ymn = 0, ymx = q_val*k_val)
  raster::plot(STD_scheme, xlab = "Sample", ylab = "Pool", 
       main = paste("STD (n = ", n, ", q = ",q_val, ", k = ", k_val, ") Pooling Scheme", sep = ""))
  layer_line = vector("numeric", length = k_val-1)
  for(i in 1:(k_val -1)){
    layer_line[i] = q_val*i
  }
  abline(h = layer_line, lty = 2)
  grid(nx = n, ny = nrow(matrix))
}

get_pools = function(matrix, cstr, n){
  
  q_val = cstr[["q"]]
  k_val = cstr[["k"]]
  gamma_val = cstr[["gamma"]]
  
  #Create a matrix of pooled samples
  pools = sapply(
    X = 1:nrow(matrix),
    FUN = function(x) {
      which(matrix[x, ] == 1)}, simplify = FALSE
  )
  
  # Initialize the final data frame
  result = matrix(0, nrow = q_val*k_val, ncol = ceiling(n/q_val))
  for(i in 1:length(pools)){
    result[i, 1:length(pools[[i]])] = pools[[i]]
  }
  
  #Add column names and row names
  colnames(result) = map2_chr(.x = rep("Sample ", times = ncol(result)), .y = seq(1:ncol(result)), .f = paste0)
  rownames(result) = map2_chr(.x = rep("Pool ", times = nrow(result)), .y = seq(1:nrow(result)), .f = paste0)
  
  return(result)
}

STD_generator = function(n, d, E, m){
  
  #Find all possible q and k values
  q_list = generate_primes(min = 1, max = n)
  gamma_list = ceiling((log(n) / log(q_list)) - 1)
  k_list = d * gamma_list + 2 * E + 1
  
  #Which (q, k) combinations pass the inequality test?
  pass = (k_list < q_list + 1) | (k_list == q_list + 1)
  param = get_param(n = n, q_list = q_list, k_list = k_list, gamma_list = gamma_list, criteria = pass)
  
  # Find the best construction
  cstr_criteria = get_construction(param = param, n = n, m = m)
  cstr_list = show_construction(param = param, criteria = cstr_criteria)
  best_cstr = get_best_cstr(list = cstr_list)
  
  # Generate the matrix M
  M = get_matrix(cstr = best_cstr, n = n)
  
  # Generate the pooling table
  final = get_pools(matrix = M, cstr = best_cstr, n = n)
  
  return(list(matrix = M, table = final, construction = best_cstr))
}



