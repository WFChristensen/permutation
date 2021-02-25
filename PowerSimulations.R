# Main function for obtaining two-sided p-values from different methods 
# comparing the centers of two populations

source("SimulationFunctions.R")

# set parameters ---------------------------------------------------------------
# n_x: sample size of group X
# epsilon = n_y/n_x: dictates the sample size of group Y
# delta = mu.x - mu-y: true difference in means between X and Y
# num_perm_dist_sampled: number of samples to draw when computing the 
#                        permutation distribution
# num_replications: desired number of simulated data sets and results
n_x <- 15  
epsilon <- 3 
delta <- -1
num_perm_dist_sampled <- 1000  
num_replications <- 10

# create initial vectors to hold results ---------------------------------------
t_test_student_pval <- rep(NA, num_replications)
t_test_student_boxcox_pval <- rep(NA, num_replications)
t_test_welch_pval <- rep(NA, num_replications)
t_test_welch_boxcox_pval <- rep(NA, num_replications)
wilcoxon_pval <- rep(NA, num_replications)
perm_mean_pval2_trad <- rep(NA, num_replications)
perm_mean_pval2_double <- rep(NA, num_replications)
perm_median_pval2_trad <- rep(NA, num_replications)
perm_median_pval2_double <- rep(NA, num_replications)
perm_mean_boxcox_pval <- rep(NA, num_replications)
tri_aspect_pval <- rep(NA, num_replications)
tri_aspect_boxcox_pval <- rep(NA, num_replications)


for (j in 1:num_replications) {
  
  # generate data --------------------------------------------------------------
  n_y <- n_x * epsilon
  sim_dat <- gendat(n_x, n_y, delta)
  vals_x <- sim_dat[[1]]
  vals_y <- sim_dat[[2]]
  
  # Box-Cox transform data -----------------------------------------------------
  group <- c(rep("grp1", n_x), rep("grp2", n_y))
  bc_res <- boxcox_trans(vals_x, 
                         vals_y, 
                         trt = group, 
                         obs_diff = NA,
                         center = TRUE)
  vals_trans <- bc_res[[1]]
  vals_x_trans <- vals_trans[1:sum(group == "grp1")]
  vals_y_trans <- vals_trans[(sum(group == "grp1") + 1):length(group)]
  
  # two-sample Students t-test [Students t(x)] ---------------------------------
  t_test_student_res <- t.test(vals_x, 
                               vals_y, 
                               alternative = "two.sided",
                               var.equal = TRUE)
  t_test_student_pval[j] <- t_test_student_res$p.value
  
  # two-sample Students t-test transformed [boxcox(Students t(x))] -------------
  t_test_student_trans_res <- t.test(vals_x_trans,
                                     vals_y_trans,
                                     alternative = "two.sided",
                                     var.equal = TRUE)
  t_test_student_boxcox_pval[j] <- t_test_student_trans_res$p.value
  
  # two-sample Welch t-test [(Welch t(x))] -------------------------------------
  t_test_welch_res <- t.test(vals_x,
                             vals_y,
                             alternative = "two.sided",
                             var.equal = FALSE)
  t_test_welch_pval[j] <- t_test_welch_res$p.value
  
  # two-sample Welch t-test transformed [boxcox(Welch t(x))] -------------------
  t_test_welch_trans_res <- t.test(vals_x_trans,
                                   vals_y_trans,
                                   alternative = "two.sided",
                                   var.equal = FALSE)
  t_test_welch_boxcox_pval[j] <- t_test_welch_trans_res$p.value
  
  # Wilcoxon Rank-Sum Test -----------------------------------------------------
  wilcox_res <- wilcox.test(c(vals_x, vals_y) ~ group)
  wilcoxon_pval[j] <- wilcox_res$p.value
  
  # permutation: mean ----------------------------------------------------------
  pval2_perm_mean <- twosamp_perm(vals_x, 
                                  vals_y, 
                                  num_perm_dist_sampled,
                                  center = "mean")
  perm_mean_pval2_trad[j] <- pval2_perm_mean[1]
  perm_mean_pval2_double[j] <- pval2_perm_mean[2]
  
  # permutation: median --------------------------------------------------------
  pval2_perm_median <- twosamp_perm(vals_x, 
                                    vals_y, 
                                    num_perm_dist_sampled,
                                    center = "median")
  perm_median_pval2_trad[j] <- pval2_perm_median[1]
  perm_median_pval2_double[j] <- pval2_perm_median[2]
  
  # permutation: mean, boxcox transformed data ---------------------------------
  pval2_perm_boxcox_dat <- twosamp_perm(vals_x_trans,
                                        vals_y_trans,
                                        num_perm_dist_sampled,
                                        center = "mean")
  perm_mean_boxcox_pval[j] <- pval2_perm_boxcox_dat[1]
  
  # tri-aspect test (Marozzi, 2007) --------------------------------------------
  tri_aspect_pval[j] <- Tabc(vals_x, vals_y, alt = "two.sided", 
                             B = num_perm_dist_sampled)
  
  # tri-aspect test (Marozzi, 2007) boxcox transformed data --------------------
  tri_aspect_boxcox_pval[j] <- Tabc(vals_x_trans, vals_y_trans, 
                                    alt = "two.sided", 
                                    B = num_perm_dist_sampled)
  
}

# summarize results ------------------------------------------------------------
results <- as.data.frame(cbind(t_test_student_pval, 
                               t_test_student_boxcox_pval,
                               t_test_welch_pval,
                               t_test_welch_boxcox_pval,
                               wilcoxon_pval,
                               perm_mean_pval2_trad,
                               perm_mean_pval2_double,
                               perm_median_pval2_trad,
                               perm_median_pval2_double,
                               perm_mean_boxcox_pval,
                               tri_aspect_pval, 
                               tri_aspect_boxcox_pval))

names(results) <- c("tStudent",
                    "tStudentBoxCox",
                    "tWelch",
                    "tWelchBoxCox",
                    "Wilcoxon",
                    "PermMean",   
                    "PermMeanAlt",
                    "PermMedian",  
                    "PermMedianAlt",
                    "PermMeanBoxCox",    
                    "Tri-Aspect",
                    "Tri-AspectBoxCox")

results_summarized <- data.frame("method" = names(results), 
                                 "prop" = colMeans(results < 0.01),
                                 row.names = 1:length(names(results)))
results_summarized

# plot summarized results ------------------------------------------------------
dotchart(results_summarized$prop, 
         labels = results_summarized$method,
         xlab = "Proportion of p-values < 0.01")

