# this file contains wrapper functions for performing multinomial logit regression within ZAP

## nnet
zap_multinom_nnet = function(Hmat, X_tilde){
  k <- ncol(X_tilde)
  if (k == 1){
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ 1, trace = F)
  }else{
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ X_tilde[, -1], trace = F)
  }
  param_multinom <-  c(t(summary( multinom_mod )$coeff))
  return(param_multinom)
}

## glmnet
zap_multinom_glmnet = function(Hmat, X_tilde){
  k <- ncol(X_tilde)
  if (k < 2){
    # throw an error
    stop("glmnet cannot take less than 2 covariates")
  }else{
    multinom_mod<-glmnet::glmnet(
                        x = X_tilde[, -1],
                        y = Hmat[, c(2,1,3)]  ,
                        family = "multinomial",
                        lambda= c(0),
                        intercept = TRUE)
  }
  beta_matrix =
    cbind(multinom_mod$beta[[2]] - multinom_mod$beta[[1]] ,
          multinom_mod$beta[[3]] - multinom_mod$beta[[1]] )

  beta_matrix = rbind(
    multinom_mod$a0[-1]- multinom_mod$a0[1],
    beta_matrix
  )

  param_multinom <-  c(as.matrix(beta_matrix))
  return(param_multinom)
}



