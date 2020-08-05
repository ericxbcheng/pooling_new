make_samples = function(n, pos){
  
  stopifnot(n >= pos)
  
  a = rep.int(x = 1, times = pos)
  b = vector(mode = "numeric", length = n - pos)
  
  return(sample(x = c(a,b), size = n, replace = FALSE))
}

make_pools = function(n, pool_size){
  
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
