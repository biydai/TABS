# R code for conducting two stage adaptive design

library(survival)
library(survcomp)
library(caret)


# function that creates folds
.cvFolds <- function(Y, V){
  # Create CV folds (stratify by outcome)   
  Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
  Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}     
  return(folds)
}

# function that computes cross validated linear predictor
cvcox <- function(survY, nn ,k = nfold, X){
  #id <- caret::createFolds(survY[,2], k)
  id <- .cvFolds(survY[,2],k)
  sig <- numeric(nn)
  for(i in 1:k){
    test_id <- id[[i]]
    df <- data.frame(survY[-test_id,],X[-test_id,])
    colnames(df)[1] <- "surv"
    #training
    fit <- survival::coxph(surv ~. , data = df)
    #fitting
    sig[test_id] <- predict(fit, newdata = X[test_id,]) 
  }
  sig
}


# function that runs analysis 

two_stage_permutation <- function(survY, #output from surv function in survival package 
                                  X, # matrix or data frame of variables of interest, without missing value
                                  itrB, # number of iteration used in permutation
                                  seed = NULL, # seed number for reproducing results
                                  nfold = 10, # cross-validation fold
                                  S = 0.5, # proportion used in first stage
                                  alpha_1 = 0.25, # split of type one error
                                  alpha_2 = 0.2,
                                  printmodel = FALSE
  
){
  if(!is.null(seed)){set.seed(seed)}
  # split data into stage1 and stage2 (stratified by censoring status)
  status <- survY[,2]
  ID <- caret::createDataPartition(status,p = S)
  nn_S <- length(ID[[1]])
  survY_1 <- survY[ID[[1]],]
  survY_2 <- survY[-ID[[1]],]
  
  X_1 <- X[ID[[1]],]
  X_2 <- X[-ID[[1]],]
  
  # fit model and get signature (with CV)
  sig <- cvcox(survY = survY_1, nn = nn_S ,k = nfold, X = X_1)
  
  # Estimate C
  survC <- survival::concordance(survY_1~sig, reverse = TRUE)$concordance
  
  # Permutation
  B <- itrB
  C_p<- numeric(B)
  resamp_id <- numeric(nn_S)
  for(i in 1:B){
    resamp_id <- sample(1:nn_S,replace = FALSE)
    survY_r <- survY_1[resamp_id,]
    psig <- cvcox(survY = survY_r, nn = nn_S ,k = nfold, X = X_1)
    C_p[i] <- survival::concordance(survY_r~psig,reverse = TRUE)$concordance
  }
  
  threshold <- quantile(C_p,1 - alpha_1)
  
  early_stop <- 0
  
  #stage 1
  if(survC < threshold){
    early_stop <- 1
    return(c(survC = survC,
             threshold =  threshold,
             early_stop = early_stop,
             C_2 = NA,
             se_2 = NA,
             final_reject = 0,
             cens = 1 - mean(status),
             cens1 = 1 - mean(survY_1[,2]),
             cens2 = 1 - mean(survY_2[,2]))
    )
  }else{
    #stage 2
    df <- data.frame(survY_1,X_1)
    fit1 <- coxph(survY_1 ~ ., data = df)
    sig2 <- predict(fit1, newdata = X_2)
    #test on stage 2
    test2 <- concordance.index(sig2, survY_2[,1], survY_2[,2], method = "noether")
    C_2 <- survival::concordance(survY_2~sig2,reverse = TRUE)$concordance
    se_2 <- test2$se
    t <- (C_2 - 0.5)/se_2
    final_reject <- t > abs(qnorm(alpha_2))
    
    if(printmodel == FALSE){
    return(c(survC = survC,
             threshold = threshold,
             early_stop = early_stop,
             C_2 = C_2,
             se_2 = se_2,
             final_reject = final_reject,
             cens = 1 - mean(status),
             cens1 = 1 - mean(survY_1[,2]),
             cens2 = 1 - mean(survY_2[,2])))
    }
    if(printmodel == TRUE){
      return(list(c(survC = survC,
               threshold = threshold,
               early_stop = early_stop,
               C_2 = C_2,
               se_2 = se_2,
               final_reject = final_reject,
               cens = 1 - mean(status),
               cens1 = 1 - mean(survY_1[,2]),
               cens2 = 1 - mean(survY_2[,2])),
               final_model = fit1))
    }
  }
}

two_stage_bootstrap <- function(survY, #output from surv function in survival package 
                                  X, # matrix or data frame of variables of interest, without missing value
                                  itrB, # number of iteration used in boostrap
                                  seed = NULL, # seed number for reproducing results
                                  nfold = 10, # cross-validation fold
                                  S = 0.5, # proportion used in first stage
                                  alpha_1 = 0.25, # split of type one error
                                  alpha_2 = 0.2,
                                  null_C = 0.5, # the null hypothesis of interest
                                  printmodel = FALSE
                                  
){
  if(!is.null(seed)){set.seed(seed)}
  # split data into stage1 and stage2 (stratified by censoring status)
  status <- survY[,2]
  ID <- caret::createDataPartition(status,p = S)
  nn_S <- length(ID[[1]])
  survY_1 <- survY[ID[[1]],]
  survY_2 <- survY[-ID[[1]],]
  
  X_1 <- X[ID[[1]],]
  X_2 <- X[-ID[[1]],]
  
  # fit model and get signature (with CV)
  sig <- cvcox(survY = survY_1, nn = nn_S ,k = nfold, X = X_1)
  
  # Estimate C
  survC <- survival::concordance(survY_1~sig, reverse = TRUE)$concordance
  
  # smoothed boostrap
  B <- itrB
  resamp_id <- numeric(nn_S)
  #varY <- var(survY_1[,1])
  h <- 0.001
  m <- mean(survY_1[,1])
  C_p <- rep(NA,B) 
  
  for(b in 1:itrB){
    resamp_id <- sample(1:nn_S,replace = TRUE)
    survY_r <- survY_1[resamp_id,]
    X_r <- X_1[resamp_id,]
    psig <- cvcox(survY = survY_r, nn = nn_S ,k = nfold, X = X_r)
    C_p[b] <- try(survival::concordance(survY_r~psig, reverse = TRUE)$concordance)
  }
  
  threshold <- quantile(C_p,alpha_1,na.rm = TRUE)
  
  early_stop <- 0
  
  #stage 1
  if(null_C > threshold){
    early_stop <- 1
    return(c(survC = survC,
             threshold =  threshold,
             early_stop = early_stop,
             C_2 = NA,
             se_2 = NA,
             final_reject = 0,
             cens = 1 - mean(status),
             cens1 = 1 - mean(survY_1[,2]),
             cens2 = 1 - mean(survY_2[,2]))
    )
  }else{
    #stage 2
    df <- data.frame(survY_1,X_1)
    fit1 <- coxph(survY_1 ~ ., data = df)
    sig2 <- predict(fit1, newdata = X_2)
    #test on stage 2
    test2 <- concordance.index(sig2, survY_2[,1], survY_2[,2], method = "noether")
    C_2 <- survival::concordance(survY_2~sig2, reverse = TRUE)$concordance
    se_2 <- test2$se
    t <- (C_2 - null_C)/se_2
    final_reject <- t > abs(qnorm(alpha_2))
    
    if(printmodel == FALSE){
      return(c(survC = survC,
               threshold = threshold,
               early_stop = early_stop,
               C_2 = C_2,
               se_2 = se_2,
               final_reject = final_reject,
               cens = 1 - mean(status),
               cens1 = 1 - mean(survY_1[,2]),
               cens2 = 1 - mean(survY_2[,2])))
    }
    if(printmodel == TRUE){
      return(list(c(survC = survC,
                  threshold = threshold,
                  early_stop = early_stop,
                  C_2 = C_2,
                  se_2 = se_2,
                  final_reject = final_reject,
                  cens = 1 - mean(status),
                  cens1 = 1 - mean(survY_1[,2]),
                  cens2 = 1 - mean(survY_2[,2])),
                  final_model = fit1))
    }
  }
}
 