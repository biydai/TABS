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

# Example 2 (simulation under the null, with permutation)
two_stage_sim(nn = 300, # sample size
              rate_l = 0.025, rate_c = 0.025,
              beta = rep(0,10),mean = 0, sd = log(1.5),
              nfold = 10,
              itrB = 50, #recommend >= 200
              S = 0.5,
              alpha_1 = 0.25, alpha_2 = 0.2,
              null_C = 0.5,
              method = "permutation",
              rep = 5, #recommend >= 500
              seedID = 999
)

# Example 3 (simulation under the null, with bootstrap)
two_stage_sim(nn = 300, # sample size
              rate_l = 0.025, rate_c = 0.025,
              beta = rep(0,10),mean = 0, sd = log(1.5),
              nfold = 10,
              itrB = 50, #recommend >= 200
              S = 0.5,
              alpha_1 = 0.25, alpha_2 = 0.2,
              null_C = 0.5,
              method = "bootstrap",
              rep = 5, #recommend >= 500
              seedID = 999
)



# Example 4: Simulation under the alternative
# step 1: tune the beta coefficients to hit targeted C index
# here, we target at 0.65
tuning_C(nn = 10000,
         mean = 0,
         sd = log(2),
         rate_l = 0.025,
         p = 7,
         by = 0.01,
         start = 0.5,
         end = 1,
         targetC = 0.65,
         seed = 111)

# step 2: tune the rate_c to hit targeted censoring percentage
tuning_censoring(nn = 10000,
                 mean = 0,
                 sd = log(2), 
                 rate_l = 0.025,
                 start = 0.001,
                 end = 0.025*0.3*(1-0.3) + 0.05,
                 by = 0.0001,
                 beta = c(rep(0.1,5),0.84),
                 target = 0.3,
                 seed = 111)

# step 3: use tuned parameters in the simulation function

two_stage_sim(nn = 500, # sample size
              rate_l = 0.025, rate_c = 0.0102,
              beta = c(rep(0.1,5),0.84),mean = 0, sd = log(1.5),
              nfold = 10,
              itrB = 50, #recommend >= 200
              S = 0.5,
              alpha_1 = 0.25, alpha_2 = 0.2,
              null_C = 0.5,
              method = "permutation",
              rep = 5, #recommend >= 500
              seedID = 999
)

two_stage_sim(nn = 500, # sample size
              rate_l = 0.025, rate_c = 0.0097,
              beta = c(rep(0.1,5),0.84),mean = 0, sd = log(1.5),
              nfold = 10,
              itrB = 50, #recommend >= 200
              S = 0.5,
              alpha_1 = 0.25, alpha_2 = 0.2,
              null_C = 0.5,
              method = "bootstrap",
              rep = 5, #recommend >= 500
              seedID = 999
)

