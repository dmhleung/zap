#' @title EM for finite-sample ZAP
#' @description an internal function used to perform the EM for zap_finite
#' @examples TBD


# this will be my EM algorithm for the finite version

# the output will be the fitted parameters for the current adaptive update


EM_finite= function(U, U_mirror, tol,  unmask_set,  mask_set,
                    pprob = NULL,  maxit, beta0, piMat, extraParam, X_tilde){

  k <- ncol(X_tilde)


    Hmat_unmask <- form_Hmat_unmask_init(U= U[unmask_set],
                                         beta0 = beta0,
                                         piMat= piMat[unmask_set, , drop = F],
                                         X_tilde = X_tilde[unmask_set, , drop = F],
                                         extraParam = extraParam)



    Hmat_Ymat_list_mask <- form_Hmat_Yab_mask_init(U_mask  = U[mask_set],
                                              U_mirror_mask = U_mirror[mask_set],
                                              beta0 = beta0,
                                              piMat= piMat[mask_set, , drop = F],
                                              X_tilde = X_tilde[mask_set, , drop = F],
                                              extraParam = extraParam)

    Hmat <- overall_Hmat( Hmat_mask = Hmat_Ymat_list_mask$Hmat, Hmat_unmask = Hmat_unmask,
                          mask_set = mask_set, unmask_set = unmask_set )


    Ymat <- overall_Ymat( Ymat_mask = Hmat_Ymat_list_mask$Ymat, U_unmask = U[unmask_set],
                          mask_set = mask_set, unmask_set = unmask_set )


  ### this part begins the breaking down the two optimization problems

    param_beta <- optim(par = beta0,
                        gr = neg_Q_fn_finite_beta_part_grad_cpp,
                        fn = neg_Q_fn_finite_beta_part,
                        Ymat = Ymat, Hmat = Hmat,
                        extraParam = extraParam, X_tilde = X_tilde,
                        method = "BFGS")$par

  # now do the multinom regression

  if (k ==1){
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ 1, trace = F)
  } else {
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ X_tilde[, -1], trace = F)
  }
  param_multinom <-  c(t(summary( multinom_mod )$coeff))
  param <- c(param_beta, param_multinom)

  squarem_obj <- SQUAREM::squarem(par = param, fixptfn=EM_finite_fix_point_fn,
                  # objfn = EM_finite_obj_fn,
                  control=list(tol=tol, maxiter= maxit),
          U = U, U_mirror = U_mirror, unmask_set = unmask_set,
          mask_set = mask_set, extraParam = extraParam, X_tilde = X_tilde )

  return_list = list(param = squarem_obj$par)
  return(return_list)

}

EM_finite_obj_fn = function(param, U, U_mirror, unmask_set,  mask_set, extraParam, X_tilde ){
  paraMat <- createParaMat(param, X_tilde)
  loglike = LogLike_finite(U = U,  U_mirror = U_mirror,
                 paraMat = paraMat, extraParam = extraParam,
                 mask_set = mask_set, unmask_set = unmask_set)
  return(-loglike)
}



## define the fix point function for applying squarem() in SQUAREM package

EM_finite_fix_point_fn = function(param, U, U_mirror, unmask_set,  mask_set, extraParam, X_tilde ){
  k <- ncol(X_tilde)
  paraMat <- createParaMat(param, X_tilde)
  Hmat_unmask <- form_Hmat_unmask(U= U[unmask_set],
                                  paraMat = paraMat[unmask_set, , drop= F],
                                  extraParam = extraParam)

  Hmat_Ymat_list_mask <- form_Hmat_Yab_mask(U_mask  = U[mask_set],
                                            U_mirror_mask = U_mirror[mask_set],
                                            paraMat = paraMat[mask_set, , drop= F],
                                            extraParam = extraParam)

  Hmat <- overall_Hmat( Hmat_mask = Hmat_Ymat_list_mask$Hmat, Hmat_unmask = Hmat_unmask,
                        mask_set = mask_set, unmask_set = unmask_set )

  Ymat <- overall_Ymat( Ymat_mask = Hmat_Ymat_list_mask$Ymat, U_unmask = U[unmask_set],
                        mask_set = mask_set, unmask_set = unmask_set )

  param_beta <- optim(par = param[1:(2*k)],
                      gr = neg_Q_fn_finite_beta_part_grad_cpp,
                      fn = neg_Q_fn_finite_beta_part,
                      Ymat = Ymat, Hmat = Hmat,
                      extraParam = extraParam, X_tilde = X_tilde,
                      method = "BFGS")$par


  if ( k ==1){
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ 1, trace = F)
  } else {
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ X_tilde[, -1], trace = F)
  }
  param_multinom <-  c(t(summary( multinom_mod )$coeff))
  param <- c(param_beta, param_multinom)

  return(param)
}




### the function below may not be useful anymore (has been replaced by cpp implementation)

# neg_Q_fn_finite_beta_part_grad_fun = function(param_beta, Hmat, Ymat, extraParam,
#                     X_tilde){
#   m  = nrow(X_tilde)
#   ### form
#   H_l = Hmat[, 1]
#   H_r = Hmat[, 3]
#   Yl_a = Ymat[, 1]
#   Yl_b = Ymat[, 2]
#   Yr_a = Ymat[, 3]
#   Yr_b = Ymat[, 4]
#   k = ncol(X_tilde)
#
#   exp_neg_X_beta_l <- as.numeric(exp( - X_tilde%*%param_beta[1:k] ))  # a vector of length m
#   exp_neg_X_beta_r <- as.numeric(exp( - X_tilde%*%param_beta[(k+1):(2*k)] )) # a vector of length m
#   k_l <- 1/(1 + exp_neg_X_beta_l) # a vector of length m
#   k_r <- 1/(1 + exp_neg_X_beta_r) # a vector of length m
#   Df_Dk_l <- H_l*(Yl_a - digamma(k_l) + digamma (k_l + extraParam[1]) ) #df/dk_l, this s a long vector now
#   Df_Dk_r <- H_r*(Yr_b - digamma(k_r) + digamma (k_r + extraParam[2]) ) #df/dk_r, this s a long vector now
#
#   Dk_l_Dbeta_l <- X_tilde*(k_l^2)* exp_neg_X_beta_l   # d k_l / d beta_l , an m-by-k matrix
#   Dk_r_Dbeta_r <- X_tilde*(k_r^2)* exp_neg_X_beta_r   # d k_r / d beta_r , an m-by-k matrix
#
#   grad_mat = cbind( Dk_l_Dbeta_l*Df_Dk_l , Dk_r_Dbeta_r*Df_Dk_r ) # matrix of dimension m-by-2k
#
#   return(-colSums( grad_mat ))
#
# }
