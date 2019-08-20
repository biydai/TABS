# load functions
source("functions.R")

# Example 1 (analyze a real data set with the two stage framework)
data(lung) # built-in data set in survival package
lung <- na.omit(lung) # for simplicity, I only took the complete cases
survY <- Surv(lung$time, lung$status)
X <- lung[,c("age","sex","ph.ecog","ph.karno","meal.cal","wt.loss")]
two_stage_permutation(survY, X, itrB = 100, S =0.5)
two_stage_bootstrap(survY, X, itrB = 100, seed = 666, S =0.5)
two_stage_bootstrap(survY, X, itrB = 100, seed = 666, S =0.5, null_C = 0.6)
