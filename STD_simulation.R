
## Generate concentration levels
rgamma_lim = function(n, alpha = 2, mode, lb, ub){
  # Alpha = 2 by default, theta is calculated by mode
  # Include lower bound and upper bound
  a = lb + rgamma(n = n, shape = alpha, scale = (mode-lb)/(alpha-1))
  # When a number is > upper bound, replace it with the mode
  if(!is.null(ub)){
    a[a >= ub] = mode
  } 
  return(a)
}

# From STD_identify.Rmd
classify = function(threshold_ind, conc, scheme){
  
  threshold_pool = threshold_ind / ncol(scheme)
  pos_pool_index = which(conc >= threshold_pool)
  
  pool_neg = scheme[-pos_pool_index, ]
  sample_neg = unique(unlist(pool_neg))
  
  pool_pos = scheme[pos_pool_index, ]
  pos1 = unique(unlist(pool_pos))
  sample_pos = pos1[!pos1 %in% sample_neg]
  
  result = list(sample_neg, sample_pos)
  names(result) = c("sample_neg", "sample_pos")
  return(result)
}

gen_elisa_af = function(n, n_pos){
  
  # Calculate the number of negative kernels
  n_neg = n - n_pos
  
  # By default, positive kernels follow a Gamma-like dist with a mode of 40000 ppb aflatoxin and a lower bound of 20 ppb
  # negative kernels follow a PERT distribution with min = 0, mode = 0.7, max < 20.
  c_pos = rgamma_lim(n = n_pos, alpha = 2, mode = 40000, lb = 20, ub = NULL)
  c_neg = rpert(n = n_neg, min = 0, mode = 0.7, max = 19.99, shape = 80)
  
  # Generate a random sequence 
  ind = sample(x = 1:n, replace = FALSE, size = n)
  
  # Randomize the positive kernels and negative kernels with the sequence
  c_all = c(c_pos, c_neg)[ind]
  names(c_all) = 1:n
  
  # Add an element 0, which is useful for STD pooling
  c_all = c("0" = 0, c_all)
  
  return(c_all)
}

gen_elisa_fm = function(n, n_pos){
  
  # Calculate the number of negative kernels
  n_neg = n - n_pos
  
  # c_pos ~ truncated normal(min = 0, max = 0.99, meanlog = -2.75, sdlog = 1.42)
  # c_neg ~ truncated normal(min = 1, max = Inf, meanlog = 3.62, sdlog = 1.74)
  c_neg = rlnormTrunc(n = n_neg, meanlog = -2.75, sdlog = 1.42, min = 0, max = 0.99)
  c_pos = rlnormTrunc(n = n_pos, meanlog = 3.62, sdlog = 1.74, min = 1, max = Inf)
  
  # Generate a random sequence 
  ind = sample(x = 1:n, replace = FALSE, size = n)
  
  # Randomize the positive kernels and negative kernels with the sequence
  c_all = c(c_pos, c_neg)[ind]
  names(c_all) = 1:n
  
  # Add an element 0, which is useful for STD pooling
  c_all = c("0" = 0, c_all)
  
  return(c_all)
}

# Pooling
gen_pool_af = function(STD_mat, conc){

  # Turn the matrix into a vector
  a = unlist(x = STD_mat, use.names = FALSE)
  
  # Match the concentrations from the "conc" to the corresponding location in the vector a
  b = data.frame("Kernel" = names(conc), "AF" = conc, stringsAsFactors = FALSE) 
  c = data.frame("Kernel" = as.character(a), stringsAsFactors = FALSE)
  d = left_join(x = c, y = b, by = "Kernel") 
  
  # Transform the concentration vector "d" back into an STD scheme
  e = matrix(data = d[["AF"]], nrow = nrow(STD_mat), ncol = ncol(STD_mat), dimnames = list(rownames(STD_mat), colnames(STD_mat))) %>%
    as.data.frame()
  
  # Calculate the pooled concentration
  c_pool = rowMeans(e)
  
  return(c_pool)
}

# From STD v3.R
calc_metrics = function(thresh, conc, n, result){
  
  # Make a contingency table
  putative_class = vector("numeric", length = n)
  putative_class[result$sample_pos] = 1
  true_class = ifelse(conc[-1] >= thresh, yes = 1, no = 0)
  
  # Manually convert the two vectors into factors
  putative_class = factor(x = putative_class, levels = c(0, 1))
  true_class = factor(x = true_class, levels = c(0, 1))
  cont_table = table(true_class, putative_class)
  
  # Calculate sensitivity and specificity
  sensi = cont_table[2,2] / (cont_table[2,2] + cont_table[2,1])
  speci = cont_table[1,1] / (cont_table[1,1] + cont_table[1,2])
  
  out = c(sensi, speci)
  
  return(out)
}

sim_outcome = function(n, n_pos, thresh, STD_mat,...){
  # Create a sequence of aflatoxin concentrations
  c_all = gen_elisa_af(n = n, n_pos = n_pos)
  
  # Calculate the pooled concentrations
  c_pool = gen_pool_af(STD_mat = STD_mat, conc = c_all)
  
  # Determine pos and neg kernels
  result = classify(threshold_ind = thresh, conc = c_pool, scheme = STD_mat)
  
  # Calculate sensitivity and specificity
  out = calc_metrics(thresh = thresh, conc = c_all, n = n, result = result)
  out2 = c(out, n_pos)
  
  return(out2)
}

# Iterate once
gen_sim_outcome = function(n, n_pos, thresh, STD_mat){
  function(...){
    sim_outcome(n = n, n_pos = n_pos, thresh = thresh, STD_mat = STD_mat)
  }
}

clean = function(list){
  a = unlist(list)
  
  # Find the indices of each type of metric
  sens_ind = seq(from = 1, by = 3, length.out = length(a)/3)
  spec_ind = sens_ind + 1
  n_pos_ind = sens_ind + 2
  
  #Rearrange the results
  b = list("sensi" = a[sens_ind], "speci" = a[spec_ind], "n_pos" = a[n_pos_ind])
  
  return(b)
}

# Iterate n times
sim_iterate = function(n_iter, ...){
  
  f_outcome = gen_sim_outcome(...)
  a = map(.x = 1:n_iter, .f = f_outcome)
  b = clean(a)
  
  return(b)
}

# Tune n_pos, each with n iterations
tune_n_pos = function(n_pos_vals, n_iter, n, thresh, STD_mat){
  map(.x = n_pos_vals, .f = sim_iterate, n_iter = n_iter, n = n, thresh = thresh, STD_mat = STD_mat)
}

## Visualization
draw_STD_2 = function(STD_mat, n, q, k){
  a = raster(x = STD_mat, xmn = 0, xmx = n, ymn = 0, ymx = q*k)
  plot(a, xlab = "Sample", ylab = "Pool", 
       main = paste0("STD(n = ", n, ", q =", q, ", k =", k, ") Pooling Scheme"))
  abline(h = seq(from = q, by = q, length.out = k-1), lty = 2)
  grid(nx = n, ny = nrow(STD_mat))
}

draw_metrics = function(df, method, n){
  if(method == "boxplot"){
    ggplot(data = df) +
      geom_boxplot(aes(x = n_pos, y = Value, group = interaction(n_pos, Type), color = Type)) +
      scale_color_manual(labels = c("Sensitivity", "Specificity"), values = c("coral", "dodgerblue")) +
      labs(x = "Number of positive kernels", y = "Evaluation metrics") +
      coord_cartesian(xlim = c(0, n), y = c(0, 1)) +
      theme_bw() +
      theme(legend.position = "top")
  } else if (method == "median_hilow"){
    ggplot(data = df, aes(x = n_pos, y = Value, color = Type)) +
      stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.95),
                   geom = "pointrange") +
      scale_color_manual(labels = c("Sensitivity", "Specificity"), values = c("coral", "dodgerblue")) +
      labs(x = "Number of positive kernels", y = "Evaluation metrics") +
      coord_cartesian(xlim = c(0, n), y = c(0, 1)) +
      theme_bw() +
      theme(legend.position = "top")
  } else {
    stop("Unknown method. Choose 'boxplot' or 'median_hilow'.")
  }
}

draw_cost = function(df, n, var, ylab = "The total number of tests needed"){
  ggplot(data = df, aes(x = n_pos, y = !!enexpr(var))) +
    stat_summary(geom = "pointrange", fun.data = median_hilow, fun.args = list(conf.int = 0.95)) +
    geom_hline(yintercept = n, lty = 2) +
    labs(x = "Number of positive kernels", y = ylab) +
    theme_bw()
}

# Contains facet_wrap()
# draw_metrics_all = function(df, n){
#   
#   ggplot(data = df, aes(x = n_pos, y = Value, color = Type)) +
#     stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.95),
#                  geom = "pointrange") +
#     scale_color_manual(labels = c("Sensitivity", "Specificity"), values = c("coral", "dodgerblue")) +
#     labs(x = "Number of positive kernels", y = "Evaluation metrics") +
#     coord_cartesian(xlim = c(0, n), y = c(0, 1)) +
#     theme_bw() +
#     theme(legend.position = "top") +
#     facet_wrap(~ Pooling)
# }

# draw_metrics_all = function(df, n){
#   
#   temp = df %>%
#     group_by(Pooling, Type, n_pos) %>%
#     summarise(med_val = median(Value),
#               q97.5 = stats::quantile(x = Value, probs = 0.975),
#               q2.5 = stats::quantile(x = Value, probs = 0.025))
#   
#   ggplot(data = temp) +
#     geom_ribbon(aes(x = n_pos, ymin = q2.5, ymax = q97.5, fill = Type), alpha = 0.2) +
#     geom_line(aes(x = n_pos, y = med_val, color = Type)) +
#     geom_point(aes(x = n_pos, y = med_val, color = Type)) +
#     facet_wrap( ~ Pooling) +
#     scale_color_manual(name = "Median",labels = c("Sensitivity", "Specificity"), values = c("coral", "dodgerblue")) +
#     scale_fill_manual(name = "Percentile (2.5th - 97.5th)",labels = c("Sensitivity", "Specificity"), values = c("coral", "dodgerblue")) +
#     labs(x = "Number of positive kernels", y = "Metrics") +
#     theme_bw() +
#     theme(legend.position = "top") 
# }

draw_metrics_all = function(df, n, ylab = "Metrics"){

  temp = df %>%
    group_by(Pooling, Type, n_pos) %>%
    summarise(med_val = median(Value),
              q97.5 = stats::quantile(x = Value, probs = 0.975),
              q2.5 = stats::quantile(x = Value, probs = 0.025))

  ggplot(data = temp) +
    geom_ribbon(aes(x = n_pos, ymin = q2.5, ymax = q97.5, group = Type), color = "grey", alpha = 0.2) +
    geom_line(aes(x = n_pos, y = med_val, group = Type)) +
    geom_point(aes(x = n_pos, y = med_val, shape = Type)) +
    facet_wrap( ~ Pooling) +
    scale_shape_manual(name = "Metric type",labels = c("Sensitivity", "Specificity"), values = c(16,17)) +
    labs(x = "Number of positive kernels", y = ylab) +
    theme_bw() +
    theme(legend.position = "top")
}


# draw_cost_all = function(df, n, var, ylab = "The total number of tests needed"){
#   ggplot(data = df, aes(x = n_pos, y = !!enexpr(var))) +
#     stat_summary(geom = "pointrange", fun.data = median_hilow, fun.args = list(conf.int = 0.95)) +
#     geom_hline(yintercept = n, lty = 2) +
#     labs(x = "The number of positive kernels", y = ylab) +
#     theme_bw() +
#     facet_wrap(~ Pooling)
# }

draw_cost_all = function(df, n, var, ylab){
  temp = df %>%
    group_by(Pooling, n_pos) %>%
    summarise(med_test = median(!!enexpr(var)),
              q2.5 = quantile(x = !!enexpr(var), probs = 0.025),
              q97.5 = quantile(x = !!enexpr(var), probs = 0.975))
  
  ggplot(data = temp) +
    geom_ribbon(aes(x = n_pos, ymin = q2.5, ymax = q97.5), color = "grey", alpha = 0.2) +
    geom_line(aes(x = n_pos, y = med_test)) +
    geom_point(aes(x = n_pos, y = med_test)) +
    facet_wrap( ~ Pooling) +
    geom_hline(yintercept = n, lty = 2) +
    labs(x = "Number of positive kernels", y = ylab) +
    theme_bw() +
    theme(legend.position = "top") 
}

calc_STD_cost = function(data, n, scheme){
  do.call(what = rbind.data.frame, args = data) %>%
    mutate(putative_pos = sensi*n_pos + (1-speci)*(n-n_pos),
           n_test_total = nrow(scheme) + putative_pos,
           n_pipette = nrow(scheme) * ncol(scheme) + nrow(scheme) + putative_pos)
}