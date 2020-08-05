#function space
make_samples = function(n, pos){
  
  stopifnot(n >= pos)
  
  a = rep.int(x = 1, times = pos)
  b = vector(mode = "numeric", length = n - pos)
  
  return(sample(x = c(a,b), size = n, replace = FALSE))
}

make_pools = function(n, pool_size, for2D=F){
  
  quotient = n %/% pool_size
  modulo = n %% pool_size
  
  if(modulo == 0){
    a = rep.int(x = 1:quotient, times = pool_size)
    return(a)
    
  } else {
    #makes the sequential pools
    a = rep.int(x = 1:quotient, times = pool_size)
    #makes the perpendicular pools for 2D
    if(for2D) a = rep(x = 1:quotient, each = pool_size)
    #adds in the remainder
    b = seq.int(from = quotient + 1, to = quotient + modulo, by = 1)
    return(c(a,b))
  }
}

compare_2D_pools = function(row_pools, col_pools, n, pool_size){
  
  #get number of sets to compare
  n_sets <- length(row_pools) %/% pool_size
  if(n > n_sets*pool_size^2) print("need to add results in for remainder")
  #compare by sets
  status<-n
  setID = rep(x = 1:n_sets, each = pool_size)
  for( i in 1:unique(setID) ) {
    d <- row_pools[i*(1:pool_size)] %o% col_pools[i*(1:pool_size)]
    status[(1:(pool_size)^2)*i] <- as.vector(d)
  }
  return(status)
}

#test space####
N = 10
POS = 2
POOL_SIZE=3
make_samples(n=N, pos=POS)
make_pools(n=N, pool_size = POOL_SIZE)
make_pools(n=N, pool_size = POOL_SIZE, for2D=F)
make_pools(n=N, pool_size = POOL_SIZE, for2D=T)

row_pools <- c(0,1,0)
col_pools <- c(1,0,0)

compare_2D_pools(row_pools, col_pools, n=9, pool_size = 3)
compare_2D_pools(row_pools=c(1,1,1), col_pools=c(1,1,0), n=9, pool_size = 3)
compare_2D_pools(row_pools=c(1,1,1), col_pools=c(1,1,0), n=10, pool_size = 3)
