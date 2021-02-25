# Code for reproducing the examples in Section 5 of "When your permutation test is doomed to fail"  ----
library(e1071)
source("SimulationFunctions.R")

nickwkend <- c(0.002660, 0.000434, 0.002150, 0.002510, 0.000434, 0.001370, 0.000434, 0.000434, 0.000434)
nickwkday <- c(0.010600, 0.002390, 0.003190, 0.002480, 0.017200, 0.008840, 0.002980, 0.002470, 0.001490, 
               0.000930, 0.044800, 0.003060, 0.007150, 0.084800, 0.000434, 0.000434, 0.000434, 0.005460, 
               0.008720, 0.000434)

KennedyD <- c(3.4522440, 0.0000000, 1.2489590, 1.7293560, 0.0000000, 0.2198285, 0.8149959, 0.7189073)
KennedyM <- c(0.0000000, 0.9592326, 0.3311258, 0.3508156, 0.0000000, 1.7930010, 0.4157140, 0.0000000, 
              0.0000000, 0.1365188, 0.0000000, 0.4039996, 0.0966931, 0.7992895, 0.0000000, 0.4101723, 
              0.0000000, 0.1635590, 0.3138873, 0.0000000, 0.1708088, 0.0000000, 0.0000000, 1.5616130,
              0.5188067, 0.3278689)

exampledata <- "IR"    # Can be "nickel" or "IR"

# Create x and y data vectors... -----------------------------------------------------------------------

if (exampledata=="nickel"){
  vals_x <- nickwkend
  vals_y <- nickwkday
}
if (exampledata=="IR"){
  vals_x <- KennedyD
  vals_y <- KennedyM
}

n_x <- length(vals_x)
n_y <- length(vals_y)


# Rules of Thumb calculations (from Section 4 of the manuscript) ---------------------------------------
ST <- function(nx){1.13*nx^.33} # Calculates the "skewness threshold" from equation (4) 
skewness(c(vals_x,vals_y))
ST(n_x)

PA  <- function(skew,nx,ny)  # Calculates the "power asymmetry" from equation (5) 
{
  exp(40*skew^2.4*(ny-nx)^.3*nx^(-3.1)) 
}
PA(skewness(c(vals_x,vals_y)), length(vals_x), length(vals_y))


# Set seed and select the number of permutations to be used in permutation tests -----------------------
set.seed(8675309)
num_perm_dist_sampled <- 50000

# Box-Cox transform data -------------------------------------------------------------------------------
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
t_test_student_pval <- t_test_student_res$p.value

# two-sample Students t-test transformed [boxcox(Students t(x))] -------------
t_test_student_trans_res <- t.test(vals_x_trans,
                                   vals_y_trans,
                                   alternative = "two.sided",
                                   var.equal = TRUE)
t_test_student_boxcox_pval <- t_test_student_trans_res$p.value

# two-sample Welch t-test [(Welch t(x))] -------------------------------------
t_test_welch_res <- t.test(vals_x,
                           vals_y,
                           alternative = "two.sided",
                           var.equal = FALSE)
t_test_welch_pval <- t_test_welch_res$p.value

# two-sample Welch t-test transformed [boxcox(Welch t(x))] -------------------
t_test_welch_trans_res <- t.test(vals_x_trans,
                                 vals_y_trans,
                                 alternative = "two.sided",
                                 var.equal = FALSE)
t_test_welch_boxcox_pval <- t_test_welch_trans_res$p.value

# Wilcoxon Rank-Sum Test -----------------------------------------------------
wilcox_res <- wilcox.test(c(vals_x, vals_y) ~ group)
wilcoxon_pval <- wilcox_res$p.value

# permutation: mean ----------------------------------------------------------
pval2_perm_mean <- twosamp_perm(vals_x, 
                                vals_y, 
                                num_perm_dist_sampled,
                                center = "mean")
perm_mean_pval2_trad <- pval2_perm_mean[1]
perm_mean_pval2_double <- pval2_perm_mean[2]

# permutation: median --------------------------------------------------------
pval2_perm_median <- twosamp_perm(vals_x, 
                                  vals_y, 
                                  num_perm_dist_sampled,
                                  center = "median")
perm_median_pval2_trad <- pval2_perm_median[1]
perm_median_pval2_double <- pval2_perm_median[2]

# permutation: mean, boxcox transformed data ---------------------------------
pval2_perm_boxcox_dat <- twosamp_perm(vals_x_trans,
                                      vals_y_trans,
                                      num_perm_dist_sampled,
                                      center = "mean")
perm_mean_boxcox_pval <- pval2_perm_boxcox_dat[1]

# tri-aspect test (Marozzi, 2007) --------------------------------------------
tri_aspect_pval <- Tabc(vals_x, vals_y, alt = "two.sided", 
                           B = num_perm_dist_sampled)

# tri-aspect test (Marozzi, 2007) boxcox transformed data --------------------
tri_aspect_boxcox_pval <- Tabc(vals_x_trans, vals_y_trans, 
                                  alt = "two.sided", 
                                  B = num_perm_dist_sampled)


# summarize results ------------------------------------------------------------
results.example <- as.data.frame(cbind(t_test_student_pval, 
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

names(results.example) <- c("tStudent",
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

results.example

