# Create a function that generates the number of positive rows, the number of positive row combinations 
# Assumption: the scenario does not exist where all the kernels in a pool are negative, but the pool concentration exceeds the pool limit
gen_1d_data_theo = function(n = 48, n_pool){
  
  stopifnot(n_pool %in% c(6,8))
  
  # Make vectors for the possible number of contaminated kernels and the possible number of positive pools
  n_pos = 1:(n-1)
  pool_pos = 1:n_pool
  
  # n_pos = number of pos kernels, pool_pos = number of positive pools, pool_combn = number of combinations of positive pool, n_retest = number of re-tests needed
  df1 = expand.grid(n_pos, pool_pos) %>%
    rename(.data = ., n_pos = Var1, pool_pos = Var2) %>%
    arrange(.data = ., n_pos) %>%
    mutate(pool_combn = choose(n = n_pool, k = pool_pos), 
           n_retest = (n / n_pool)*pool_pos)
  
  # Only keep the obs where number of positive kernels >= number of positive pools
  df2 = df1 %>%
    dplyr::filter(.data = ., n_pos >= pool_pos)
  
  # layer = number of pools that can be filled up by 8 kernels. %/% means integer division. 
  # extra = number of kernels that are left out if we try to fill up the pools first. %% means modulo
  # layer = 0 means the current number of positive kernels < 8 and they can all appear in 1 pool; 
  # layer = 1 means there's at least 2 positive pools, because the number of positive kernels >= 9, etc.
  # layer = 1, extra = 0 means there are 1*8 + 0 = 8 pos kernels. They can perfectly fit into 1 pool
  # layer = 2, extra = 3 means there are 2*8 + 3 = 19 pos kernels. They can fill up 2 pools
  # This is useful for calculating pool_pos
  df3 = df2 %>%
    mutate(layer = n_pos %/% (n/n_pool),
           extra = n_pos %% (n/n_pool))
  
  # Filter the obs to remove impossible pool_pos at high n_pos
  # logi1 means find obs where after filling up all the pools with kernels, there are still extra kernels
  # logi2 means find obs where there's just enough of kernels that can fill up multiple pools without extra kernels
  logi1 = (df3$pool_pos >= df3$layer + 1) & (df3$extra >= 1)
  logi2 = (df3$pool_pos >= df3$layer) & (df3$extra == 0)
  
  df4 = df3 %>%
    dplyr::filter(.data = ., logi1|logi2)
  
  # Calculate the probability of different combinations of positive pools
  ## Calculate the total number of combinations of positive pools per n_pos
  df5 = df4 %>%
    dplyr::select(n_pos, pool_combn) %>%
    group_by(n_pos) %>%
    summarise(n_combn = sum(pool_combn))
  
  ## Join df5 and df4 by n_pos, then calculate the probability of combinations of positive pools, sensitivity, specificity 
  df6 = left_join(x = df4, y = df5, by = "n_pos") %>%
    mutate(P_pool_combn = pool_combn/n_combn, 
           sensi = 1,
           speci = (n - n_pos)/((n - n_pos) + (n_retest - n_pos)))
  
  return(df6)
}

calc_metrics_1d_theo = function(df){
  # Calculate the expected sensitivity and expected specificity
  df2 = df %>%
    mutate(frac_sensi = P_pool_combn * sensi,
           frac_speci = P_pool_combn * speci) %>%
    group_by(n_pos) %>%
    summarise(E_sensi = sum(frac_sensi),
              E_speci = sum(frac_speci))
  
  # Clean df2 by converting it to a long format data frame
  df3 = gather(data = df2, key = "Type", value = "Value", -n_pos)
  
  return(df3)
}

calc_cost_1d_theo = function(df, n_pool){
  
  # Calculate the expected number of retests for each n_pos
  df2 = df %>%
    group_by(n_pos) %>%
    mutate(frac_n_retest = n_retest*P_pool_combn) %>%
    summarise(.data = ., E_n_retest = ceiling(sum(frac_n_retest)))
  
  # Calculate the expected total number of tests (n_pool pools + # retests)
  df3 = df2 %>%
    mutate(E_total_tests = n_pool + E_n_retest)
  
  return(df3)
}

draw_metrics_1d_theo = function(data, n){
  ggplot(data = data) +
    geom_line(aes(x = n_pos, y = Value, color = Type)) +
    geom_point(aes(x = n_pos, y = Value, color = Type)) +
    scale_color_manual(labels = c("Sensitivity", "Specificity"), values = c("coral", "dodgerblue")) +
    labs(x = "Number of positive kernels", y = "Evaluation metrics") +
    coord_cartesian(xlim = c(0, n), y = c(0, 1)) +
    theme_bw()
}

draw_cost_1d_theo = function(data, n){
  ggplot(data = data) +
    geom_line(aes(x = n_pos, y = E_total_tests)) +
    geom_point(aes(x = n_pos, y = E_total_tests)) +
    geom_hline(yintercept = n, lty = 2) +
    labs(x = "Number of positive kernels", y = "Total number of tests needed") +
    theme_bw()
}