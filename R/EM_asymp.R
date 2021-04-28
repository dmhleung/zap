

# @title EM for asymptotic ZAP
# @description an internal function used to perform the EM for zap_asymp
# @examples TBD



EM_asymp= function(z, param,  tol= 0.0001 , pprob = NULL,
                   maxit, extraParam, X_tilde, uvals.cutoff){

  k <- ncol(X_tilde)
  U<- pnorm(z)
  U[U <= uvals.cutoff] <- uvals.cutoff
  U[U >= 1-  uvals.cutoff] <- 1 -  uvals.cutoff
  squarem_obj <- SQUAREM::squarem(par = param, fixptfn=EM_asymp_fix_point_fn,
                        # objfn = EM_asymp_obj_fn,
                        control=list(tol=tol, maxiter= maxit),
                        U = U,  extraParam = extraParam, X_tilde = X_tilde )

  return_list = list(param = squarem_obj$par,
                     loglike = -EM_asymp_obj_fn(squarem_obj$par, U,  extraParam, X_tilde))
}




# @title EM asymp fix point function for SQUAREM
# @description an internal function used to perform the EM for zap_asymp using the  SQUAREM package
# @examples TBD
EM_asymp_fix_point_fn = function(param, U, extraParam, X_tilde ){
  k <- ncol(X_tilde)
  paraMat <- createParaMat(param, X_tilde)
  Hmat<- form_Hmat_unmask(U= U, paraMat = paraMat, extraParam = extraParam)
  Yl_a = trunc_log_arma(U)
  Yl_b = trunc_log_arma(1 - U)
   Umat = cbind( Yl_a ,  Yl_b, Yl_a  , Yl_b)   # Yr_a  = Yl_a, Yr_b  = Yl_b
  param_beta <- optim(par = param[1:(2*k)],
                      gr = neg_Q_fn_finite_beta_part_grad_cpp,
                      fn = neg_Q_fn_finite_beta_part,
                      Ymat =  Umat , Hmat = Hmat,
                      extraParam = extraParam, X_tilde = X_tilde,
                      method = "BFGS")$par

  if (k == 1){
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ 1, trace = F)
  }else{
    multinom_mod <- nnet::multinom(formula = Hmat[, c(2,1,3)] ~ X_tilde[, -1], trace = F)
  }
  param_multinom <-  c(t(summary( multinom_mod )$coeff))
  param <- c(param_beta, param_multinom)

  return(param)
}


# @title Negative log-likelihood of the 3-component beta mixture
# @description an internal function used to return the objective function value for the EM of zap_asymp
# @examples TBD
EM_asymp_obj_fn = function(param, U,  extraParam, X_tilde ){
  paraMat <- createParaMat(param, X_tilde)
  loglike = LogLike_asymp(U = U, paraMat = paraMat, extraParam = extraParam)
  return(-loglike)
}


