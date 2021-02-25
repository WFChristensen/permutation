# gendat() generates data for a two independent samples, X and Y, from the
# lognormal distribution with sd_x = sd_y = 1 so that 
# skewness = (exp(sd_x ^ 2) + 2) * sqrt(exp(sd_x ^ 2) - 1) = 6.18
# n_x is the sample size of group X, n_y is the sample size of group Y, and 
# delta is the desired difference in the means of the two groups.
# gendat() returns the data generated for both samples (vals_x and vals_y)
gendat <- function(n_x, n_y, delta) {
  
  sd_x <- 1
  sd_y <- 1
  
  if (delta < 0) {
    
    mean_y <- 0
    tau <- log(exp(mean_y + (sd_y^2 / 2)) + delta) - mean_y - (sd_y^2 / 2)
    mean_x <- mean_y + tau
    
  } else if (delta > 0) {
    
    mean_x <- 0
    tau <- log(exp(mean_x + (sd_x^2 / 2)) - delta) - mean_x - (sd_x^2 / 2)
    mean_y <- mean_x + tau
    
  } else {
    
    mean_x <- 0
    mean_y <- 0
    
  }
  
  vals_x <- exp(rnorm(n = n_x, mean = mean_x, sd = sd_x))
  vals_y <- exp(rnorm(n = n_y, mean = mean_y, sd = sd_y)) 
  
  list(vals_x, vals_y)
  
}


# generate the permutation distribution ----------------------------------------
twosamp_permdist <- function(vals_x, vals_y, num_samp, center = "mean") {
  
  vals <- c(vals_x, vals_y)
  perm_dist_new <- rep(NA, num_samp)
  
  for (i in 1:num_samp) {
    
    permuted_vals <- sample(1:length(vals), 
                            size = length(vals_x), 
                            replace = FALSE)
    
    if (center == "mean") {
      
      perm_dist_new[i] <- mean(vals[permuted_vals]) - 
        mean(vals[-permuted_vals])
      
    } else {
      
      perm_dist_new[i] <- median(vals[permuted_vals]) -
        median(vals[-permuted_vals])
      
    }
    
  }
  
  perm_dist_new
  
}


# calculate both one-sided p-values --------------------------------------------
onesided_pvals <- function(perm_dist, obs_diff, num_samp) {
  # upper tail p-value
  pval_up <- (sum(perm_dist >= obs_diff) + 1) / (num_samp + 1)
  # lower tail p-value
  pval_low <- (sum(perm_dist <= obs_diff) + 1) / (num_samp + 1)
  return(c(pval_up, pval_low))
  
}


# calculate the traditional two-sided p-value ----------------------------------
twosided_pval_trad <- function(perm_dist, obs_diff, num_samp) {
  
  pval2_trad <- (sum(abs(perm_dist) >= abs(obs_diff)) + 1) / (num_samp + 1)
  pval2_trad
  
}


# calculate the alternative two-sided p-value (by doubling the smallest 
# one-sided p-value) -----------------------------------------------------------
twosided_pval_double <- function(pval_up, pval_low) {
  
  pval2_double <- min(pval_up, pval_low) * 2
  pval2_double
  
}


# Box Cox transform (1) a vector of values from two populations (treatments) or 
# (2) a permutation distribution
# for (1): vals = vals, trt = trt, obs_diff = NA
# for (2): vals = vals, trt = 1, obs_diff = obs_diff ---------------------------
boxcox_trans <- function(vals_x, vals_y, trt = FALSE, obs_diff = NA, 
                         center = FALSE) {
  
  if (center == TRUE) {
    
    vals_centered <- c(vals_x - mean(vals_x), vals_y - mean(vals_y))
    
  } else {
    
    vals_centered <- c(vals_x, vals_y)
    
  }
  
  # shift values up to only have positive values for boxcox
  if (min(vals_centered) < 0) { 
    
    vals_shifted <- vals_centered - min(vals_centered) + 0.01
    obs_diff_shifted <- obs_diff - min(vals_centered) + 0.01
    
    if (!is.na(obs_diff_shifted) & obs_diff_shifted < 0) {
      
      obs_diff_shifted <- 0
      
    }
    
  } else if (min(vals_centered) == 0) { 
    
    vals_shifted <- vals_centered + 0.01
    obs_diff_shifted <- obs_diff + 0.01
    
  } else {
    
    vals_shifted <- vals_centered
    obs_diff_shifted <- obs_diff
    
  }
  
  # boxcox transform vals: get optimal lambda value
  if (is.na(obs_diff)) {  # trans vals
    
    bc <- MASS::boxcox(vals_shifted ~trt,
                       lambda = seq(-6, 6, 0.01),
                       plotit = FALSE)  
    
  } else {  # trans perm dist
    
    bc <- MASS::boxcox(vals_shifted ~1,
                       lambda = seq(-6, 6, 0.01),
                       plotit = FALSE) 
    
  }
  
  bc_df <- data.frame(bc$x, bc$y)  
  bc_max <- bc_df[which.max(bc_df$bc.y), ]
  lambda <- bc_max$bc.x
  
  # apply bc lambda value transform to original, uncentered data
  vals <- c(vals_x, vals_y)
  
  # shift values up to only have positive values for boxcox
  if (min(vals) < 0) { 
    
    vals_orig <- vals - min(vals) + 0.01
    obs_diff_orig <- obs_diff - min(vals) + 0.01
    if (!is.na(obs_diff_orig) & obs_diff_orig < 0) {
      obs_diff_orig <- 0
      
    }
    
  } else if (min(vals) == 0) {  
    
    vals_orig <- vals + 0.01
    obs_diff_orig <- obs_diff + 0.01
    
  } else {
    
    vals_orig <- vals
    obs_diff_orig <- obs_diff
    
  }
  
  # using optimal lambda value, transform vals
  if (lambda == 0) {
    
    vals_trans <- log(vals_orig)
    obs_diff_trans <- log(obs_diff_shifted)
    
  } else {
    
    vals_trans <- (vals_orig^lambda - 1) / lambda
    obs_diff_trans <- (obs_diff_shifted^lambda - 1) / lambda
    
  }
  
  list(vals_trans, obs_diff_trans)
  
} 


# returns the results of a permutation test for the difference in means/medians 
# between two samples  ---------------------------------------------------------
twosamp_perm <- function(vals_x, vals_y, num_samp, center = "mean", 
                         transform_perm_dist = FALSE) {
  
  if (!transform_perm_dist) {
    
    obs_diff <- ifelse(center == "mean", 
                       mean(vals_x) - mean(vals_y),
                       median(vals_x) - median(vals_y))
    
    perm_dist <- twosamp_permdist(vals_x, 
                                  vals_y, 
                                  num_samp = num_samp,
                                  center = center)
    
    pvals <- onesided_pvals(perm_dist, obs_diff, num_samp)
    pval_up <- pvals[1]
    pval_low <- pvals[2]
    pvals2_trad <- twosided_pval_trad(perm_dist, obs_diff, num_samp)
    pvals2_double <- twosided_pval_double(pval_up, pval_low)
    
  } else {  # Box Cox transform permutation dist then get two-sided p-value
    
    obs_diff <- mean(vals_x) - mean(vals_y)
    
    perm_dist <- twosamp_permdist(vals_x, 
                                  vals_y, 
                                  num_samp = num_samp,
                                  center = center)
    
    bc_res <- boxcox_trans(perm_dist, trt = 1, obs_diff = obs_diff)
    perm_dist_trans <- bc_res[[1]]
    obs_diff_trans <- bc_res[[2]]
    
    pvals <- onesided_pvals(perm_dist_trans, 
                            obs_diff_trans, 
                            num_samp)
    pval_up <- pvals[1]
    pval_low <- pvals[2]
    pvals2_trad <- twosided_pval_trad(perm_dist_trans, 
                                      obs_diff_trans, 
                                      num_samp)
    pvals2_double <- twosided_pval_double(pval_up, pval_low)
    
  }
  
  return(c(pvals2_trad, pvals2_double))
  
}


# Code provided by Dr. Marco Marozzi
# Computes the p-value from the tri-aspect test
# M. Marozzi (2004), A Bi-Aspect Nonparametric Test for the Two-Sample Location 
# Problem, Computational Statistics and Data Analysis, 44, p. 639-648.
# M. Marozzi (2007), Multivariate Tri-Aspect Non-Parametric Testing, Journal of 
# Nonparametric Statistics, 19, 6, p. 269-282.
# Please notify any usage to marco.marozzi@unive.it ----------------------------
Tabc = function(x1, x2, alt, B = 1000) {
  x = c(x1, x2)
  n1 = length(x1)
  n2 = length(x2)
  n = n1 + n2
  mediana = median(x) 
  ranghi = rank(x) 
  # changed below from sum() to mean() so code works with unbalanced samples:
  ta.ob = mean(x1) - mean(x2)
  tb.ob = length(x1[x1 >= mediana]) - length(x2[x2 >= mediana]) 
  tc.ob = sum(ranghi[1:n1]) - sum(ranghi[(n1 + 1):n])
  
  ta.perm = vector(, B)
  tb.perm = vector(, B)
  tc.perm = vector(, B)
  tabc.perm = vector(, B)
  
  for (b in 1:B) {
    x.perm = sample(x)
    x1.perm = x.perm[1:n1]  
    x2.perm = x.perm[(n1 + 1):(n1 + n2)] 
    ranghi.perm = rank(x.perm)
    ta.perm[b] = mean(x1.perm) - mean(x2.perm)
    tb.perm[b] = length(x1.perm[x1.perm >= mediana]) - 
      length(x2.perm[x2.perm >= mediana])
    tc.perm[b] = sum(ranghi.perm[1:n1]) - sum(ranghi.perm[(n1 + 1):n])
  }
  
  
  if (alt == "greater") {
    pv.ta.ob = length(ta.perm[ta.perm >= ta.ob]) / B
    pv.tb.ob = length(tb.perm[tb.perm >= tb.ob]) / B
    pv.tc.ob = length(tc.perm[tc.perm >= tc.ob]) / B
  }
  if (alt == "less") {   
    pv.ta.ob = length(ta.perm[ta.perm <= ta.ob]) / B
    pv.tb.ob = length(tb.perm[tb.perm <= tb.ob]) / B
    pv.tc.ob = length(tc.perm[tc.perm <= tc.ob]) / B
  }
  if (alt == "two.sided") {
    pv.ta.ob = length(abs(ta.perm)[abs(ta.perm) >= abs(ta.ob)]) / B
    pv.tb.ob = length(abs(tb.perm)[abs(tb.perm) >= abs(tb.ob)]) / B
    pv.tc.ob = length(abs(tc.perm)[abs(tc.perm) >= abs(tc.ob)]) / B
  }
  
  tabc.ob = min(pv.ta.ob, pv.tb.ob, pv.tc.ob)
  
  for (b in 1:B) {
    if (alt == "greater") {
      pv.ta.perm = length(ta.perm[ta.perm >= ta.perm[b]]) / B
      pv.tb.perm = length(tb.perm[tb.perm >= tb.perm[b]]) / B
      pv.tc.perm = length(tc.perm[tc.perm >= tc.perm[b]]) / B
    } 
    if (alt == "less") {
      pv.ta.perm = length(ta.perm[ta.perm <= ta.perm[b]]) / B
      pv.tb.perm = length(tb.perm[tb.perm <= tb.perm[b]]) / B
      pv.tc.perm = length(tc.perm[tc.perm <= tc.perm[b]]) / B
    } 
    if (alt == "two.sided") {
      pv.ta.perm = length(abs(ta.perm)[abs(ta.perm) >= abs(ta.perm[b])]) / B
      pv.tb.perm = length(abs(tb.perm)[abs(tb.perm) >= abs(tb.perm[b])]) / B
      pv.tc.perm = length(abs(tc.perm)[abs(tc.perm) >= abs(tc.perm[b])]) / B
    }
    
    tabc.perm[b] = min(pv.ta.perm, pv.tb.perm, pv.tc.perm)
  }
  
  pv.tabc = length(tabc.perm[tabc.perm <= tabc.ob]) / B
  # print(pv.tabc)
  return(pv.tabc)
}