# this gets a reasonable initialization for the EM algorithm in the spirit of Lei and Fithian

truncate = function(x, a, b){

  return(  pmin(b, pmax( x, a)))
}



get_piMat_init = function( unmask_set, mask_set, X, s_l, s_r,  l_set, r_set ){


  # first do the binomial regression for data on the "right"
  m <- nrow(X)
  J_tilde_indicator <- numeric(m)
  J_tilde_indicator[mask_set] = 1
  data_r <- data.frame(X = X[r_set, ], response = J_tilde_indicator[r_set])
  data_l <- data.frame(X = X[l_set, ], response = J_tilde_indicator[l_set])
  logit_r <- glm(response ~ . , family = binomial, data = data_r)
  logit_l <- glm(response ~ . , family = binomial, data = data_l)

  # get conditional pi_r
  pi_J_tilde_r <-  predict(logit_r, newdata = data.frame(X = X), type = "response")  # this gives predicted probabilities
  J_tilde_r_null_val <-   1 - 0.5/(0.5 - 2*( 1- s_r))
  pi_r_cond_est<- (1  - pi_J_tilde_r )*J_tilde_r_null_val + pi_J_tilde_r
  pi_r_cond_est <- truncate(pi_r_cond_est, 0, 1)


  # get conditional pi_l
  pi_J_tilde_l <-  predict(logit_l, newdata = data.frame(X = X), type = "response")
  J_tilde_l_null_val <-   1 - 0.5/(0.5 - 2*s_l)
  pi_l_cond_est<- (1  - pi_J_tilde_l )*J_tilde_l_null_val + pi_J_tilde_l
  pi_l_cond_est <- truncate(pi_l_cond_est, 0, 1)


  # Now we try to obtain the fitted probabilities for  U_i to be > 0.5
  D <- numeric(m)
  D[r_set] <- 1
  logit_right_side <- glm(response ~ . , family = binomial, data = data.frame(response = D, X = X))
  prob_to_right <- predict(logit_right_side, newdata = data.frame(X = X), type = "response")
  ## final multinom regression with full data
  # pi_r_est <- pi_r_cond_est*length(r_set)/m
  # pi_l_est <- pi_l_cond_est*length(l_set)/m
  pi_r_est <- pi_r_cond_est*prob_to_right
  pi_l_est <- pi_l_cond_est*(1 - prob_to_right)
  pi_null_est <- 1 - pi_l_est - pi_r_est
  full_dat <- data.frame(left = pi_l_est, null = pi_null_est,  right = pi_r_est, X)

  # multinom_mod  <- nnet::multinom( as.matrix(full_dat[, c(2, 1, 3)]) ~ X1 + X2  , data=full_dat)
  #prediction <- predict(multinom_mod, type = "probs")
  # theta_l <- summary(multinom_mod)$coeff[1, ]
  # theta_r <- summary(multinom_mod)$coeff[2, ]
  # return(c(theta_l, theta_r))
  return(as.matrix(full_dat[, 1:3]))

  #
  # a <-  exp(cbind(1, X)%*%t(summary(mlogit_mod)$coeff))
  # head( (rowSums(cbind(1, a)))^{-1}*a)
  #
  # H_null <- prediction[, 1]
  # H_l <- prediction[, 2]
  # H_r  <- prediction[, 3]
  #

}




get_beta_init = function(U,  X_tilde,  l_set, r_set,extraParam, beta0){

  k = ncol(X_tilde)
  U_l <- U[l_set]
  U_r <- U[r_set]
  gamma_l <- extraParam[1]
  gamma_r <- extraParam[2]
  beta_l <- optim(par = beta0[1:k],
                  fn = neg_LogLike_LeftBeta,
                  gr = neg_LogLike_LeftBeta_grad,
                  U = U_l ,
                  X_tilde =  X_tilde[l_set, , drop = F] ,
                  gamma_l =  gamma_l, method = "BFGS")$par

  beta_r <- optim(par = beta0[(k+1):(2*k)],
                  fn = neg_LogLike_RightBeta,
                  gr = neg_LogLike_RightBeta_grad,
                  U = U_r ,
                  X_tilde =  X_tilde[r_set, , drop = F] ,
                  gamma_r =  gamma_r, method = "BFGS")$par

  return(c(beta_l, beta_r))

}



