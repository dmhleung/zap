

# @title EM for asymptotic ZAP
# @description an internal function used to perform the EM for zap_asymp
# @examples TBD



EM_asymp= function(z,
                   param,
                   tol= 0.0001,
                   pprob = NULL,
                   maxit,
                   extraParam,
                   X_tilde,
                   uvals.cutoff,
                   workhorse){

  k <- ncol(X_tilde)
  U<- stats::pnorm(z)
  U[U <= uvals.cutoff] <- uvals.cutoff
  U[U >= 1-  uvals.cutoff] <- 1 -  uvals.cutoff

  squarem_obj <- try(SQUAREM::squarem(par = param, fixptfn=EM_asymp_fix_point_fn,
                             # objfn = EM_asymp_obj_fn,
                             control=list(tol=tol, maxiter= maxit),
                             U = U,  extraParam = extraParam, X_tilde = X_tilde, workhorse = workhorse ))

  if (inherits(squarem_obj, "try-error")){
    attributes(squarem_obj)$condition
    print("squarem has an error; we perform a usual EM loop now without the help of squarem!")

    # We will have to perform an usual EM procedure now without the help of squarem
    new_param = EM_asymp_fix_point_fn(param, U, extraParam, X_tilde)
    eval_count = 1
    norm_diff_param = sqrt(sum(abs(new_param - param)^2))
    while(norm_diff_param > tol & eval_count <= maxit){
      param = new_param
      new_param = EM_asymp_fix_point_fn(param, U, extraParam, X_tilde)
      norm_diff_param = sqrt(sum(abs(new_param - param)^2))
      eval_count = eval_count +1
    }
    param = new_param
  } else {
    param = squarem_obj$par
  }
  return_list = list(param = param,
                     loglike = -EM_asymp_obj_fn(param, U,  extraParam, X_tilde))
  return(return_list)
}




# @title EM asymp fix point function for SQUAREM
# @description an internal function used to perform the EM for zap_asymp using the  SQUAREM package
# @examples TBD
EM_asymp_fix_point_fn = function(param, U, extraParam, X_tilde, workhorse ){
  k <- ncol(X_tilde)
  paraMat <- createParaMat(param, X_tilde)
  Hmat<- form_Hmat_unmask(U= U, paraMat = paraMat, extraParam = extraParam)
  Yl_a = trunc_log_arma(U)
  Yl_b = trunc_log_arma(1 - U)
  Umat = cbind( Yl_a ,  Yl_b, Yl_a  , Yl_b)   # Yr_a  = Yl_a, Yr_b  = Yl_b
  param_beta <- stats::optim(par = param[1:(2*k)],
                      gr = neg_Q_fn_finite_beta_part_grad_cpp,
                      fn = neg_Q_fn_finite_beta_part,
                      Ymat =  Umat , Hmat = Hmat,
                      extraParam = extraParam, X_tilde = X_tilde,
                      method = "BFGS")$par

  if (workhorse == "nnet"){
    param_multinom <- zap_multinom_nnet(Hmat, X_tilde)
  }

  if (workhorse == "glmnet"){
    param_multinom <- zap_multinom_glmnet(Hmat, X_tilde)
  }

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


