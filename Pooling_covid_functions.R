make_samples = function(n, pos){
  
  stopifnot(n >= pos)
  
  a = rep.int(x = 1, times = pos)
  b = vector(mode = "numeric", length = n - pos)
  return(c(a,b))
}
