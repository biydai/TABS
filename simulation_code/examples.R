# load functions
source("functions.R")

# Example 1 (simulation under the null, with permutation)
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

# Example 2 (simulation under the null, with bootstrap)
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



# Example 3: Simulation under the alternative
# step 1: tune the beta coefficients to hit targeted C index
# here, we target at 0.6

tuned_C <- tuning_C(nn = 10000,
                    mean = 0,
                    sd = log(1.5),
                    rate_l = 0.025,
                    p = 5,
                    by = 0.01,
                    start = 0.5,
                    end = 1.1,
                    targetC = 0.6,
                    seed = 123)
tuned_C
tuned_C$beta
tuned_C$beta_j


#### Step 2: Tune Censoring Distribution to Hit Targeted Censoring Percentage

tuned_censor<- tuning_censoring(nn = 10000,
                                mean = 0,
                                sd = log(1.5), 
                                rate_l = 0.025,
                                start = 0.001,
                                end = 0.025*0.3*(1-0.3) + 0.05,
                                by = 0.0001,
                                beta = c(tuned_C$beta,tuned_C$beta_j),
                                target = 0.3,
                                seed = 111)
tuned_censor
tuned_censor$rate_c

#### Step 3: Input the parameters obtained from Step 1 and 2 into the simulation function:
# Example of Permutation under the alternative
result_P <- two_stage_sim(nn = 300, # sample size
                          rate_l = 0.025, rate_c = tuned_censor$rate_c,
                          beta = c(tuned_C$beta,tuned_C$beta_j),mean = 0, sd = log(1.5),
                          nfold = 10,
                          itrB = 50, #recommend >= 200
                          S = 0.5,
                          alpha_1 = 0.25, alpha_2 = 0.2,
                          method = "permutation",
                          rep = 10, #recommend >= 500
                          seedID = 999
)
head(result_P)
# power:
mean(result_P[,"final_validation"])

# probability of Early Stopping:
mean(result_P[,"early_stop"])

## Bootstrap (under Null: C = 0.6):
  
result_B <- two_stage_sim(nn = 300, # sample size
                          rate_l = 0.025, rate_c =  tuned_censor$rate_c,
                          beta = c(tuned_C$beta,tuned_C$beta_j),mean = 0, sd = log(1.5),
                          nfold = 10,
                          itrB = 50, #recommend >= 200
                          S = 0.5,
                          alpha_1 = 0.25, alpha_2 = 0.2,
                          null_C = 0.6,
                          method = "bootstrap",
                          rep = 20, #recommend >= 500
                          seedID = 123
)
head(result_B)
# type one error rate:
mean(result_B[,"final_validation"])

# probability of Early Stopping:
mean(result_B[,"early_stop"])

