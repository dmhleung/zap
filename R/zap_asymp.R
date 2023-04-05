#' @title Asymptotic ZAP
#'
#' @description Asymptotic version of ZAP based on a 3-component beta mixture
#'   working model.
#'
#' @param z a vector of z-values of length n
#' @param param_init a numeric vector of the form \eqn{(\beta_l, \beta_r, \theta_l, \theta_r)} that initialized he values for \eqn{\beta_l}, \eqn{\beta_r}, \eqn{\theta_l} and \eqn{\theta_r}, in that order; see the paper.  By default, if there is no input from the user ( \code{param_init = NULL} ), the intercept coefficients in \eqn{\beta_l}, \eqn{\beta_r}, \eqn{\theta_l} and \eqn{\theta_r} are set to be \eqn{log(5)}, \eqn{log(5)}, \eqn{log(0.1)} and \eqn{log(0.1)} respectively; all the other coefficients are set to be zero.
#' @param X a covariate matrix of dimension n by p, where  n is the length of z and p is
#'   the number of covariates; p must at least be 1, i.e. there must be at least one covariate, intercept excluded. Must be a matrix (even if it has just one column).
#' @param gamma_l the choice of  \eqn{\gamma_l}, the second shape parameter for the left-leaning beta density. Default to 4.
#' @param gamma_r the choice of  \eqn{\gamma_r}, the second shape parameter for the right-leaning beta density. Default to 4.
#' @param gen_size number of random uniformly distributed realizations generated
#'   to estimate the mirror statistics. Default to 50000.
#' @param maxit maximum number of iterations performed for the EM algorithm. Default to 100.
#' @param tol tolerance level to terminate the EM algorithm if \deqn{\|(\beta_l, \beta_r, \theta_l, \theta_r)^{(t)} - (\beta_l, \beta_r, \theta_l, \theta_r)^{(t+1)}\|_2 \leq tol,} where \eqn{\|\cdot\|_2} is the Euclidean norm. Default to 0.0001.
#' @param alpha  targeted FDR level. Default to 0.05.
#' @param ep additive factor in the numerator of the FDP estimates in ZAP;
#'   default to 1. Changing it may result in an inflated FDR of the procedure.
#' @param uvals.cutoff an ``epsilon-small" number. The u-values will be cut off
#'   to be within the interval \deqn{[uvals.cutoff, \quad 1 - uvals.cutoff].} This
#'   stabilizes the beta mixture density estimation. Default is 1e-15.
#' @param workhorse the workhorse for running the multinomial logit regression for the class probabilities under the hood; can choose between "glmnet" and "nnet", which correspond to the functions \code{glmnet::glmnet()} and \code{nnet::multinom()} respectively. Default to "glmnet", which is recommended for p >= 2, i.e. when the covariate matrix \code{X} has two or more columns. If p = 1, it is automatically fixed to be "nnet" (Since \code{glmnet::glmnet()} can only take two or more covariates).
#' @details TBD
#' @return a list with containing the followings: \itemize{ \item rej_index: a
#'   vector of indices corresponding to the rejected hypotheses. Can be of
#'   length 0 if no hypotheses are rejected. \item stat: a vector of resulting
#'   test statistics \item mirror_stat: a vector of "mirror" test statistics
#'   \item param: a vector of estimated parameters concatenating  \eqn{\beta_l}, \eqn{\beta_r}, \eqn{\theta_l} and
#'   \eqn{\theta_r} in that order }
#' @export
#' @examples #TBD


zap_asymp = function(z,
                     X,
                     param_init = NULL,
                     gamma_l = 4, gamma_r = 4,
                     gen_size = 50000, alpha = 0.05,
                     ep =1 , uvals.cutoff = 1e-15,
                     tol = 0.0001,
                     maxit = 100,
                     workhorse = "glmnet"){

  if (!is.matrix(X)) stop("X must be a matrix!")  # make sure X is a matrix
  if (ncol(X) < 2) workhorse = "nnet"

  # initialize param_init if no input from the user
  if (is.null(param_init)) param_init <- c(log(5), rep(0, ncol(X)), log(5), rep(0, ncol(X)),
                 log(.1), rep(0, ncol(X)), log(.1), rep(0, ncol(X)))

  extraParam <- c(gamma_l, gamma_r)

  EM_fit <- try(
    R.utils::withTimeout(
      EM_asymp(z, param = param_init ,
               tol = tol,
               maxit = maxit,
               extraParam = extraParam,
               X_tilde = cbind(1, X),
               uvals.cutoff = uvals.cutoff,
               workhorse = workhorse),
      timeout = 60)
  )

  if (inherits(EM_fit, "try-error")){
    EM_fit = list(param = param_init)
  }

  blfdr_parts <-create_blfdr_parts(u = stats::pnorm(z),
                                   param = EM_fit$param,
                                   extraParam = extraParam,
                                   X = X,
                                   gen_size = gen_size)

  blfdr_vec <- blfdr_parts$blfdr_vec
  blfdr_flip_vec <- blfdr_parts$blfdr_flip_vec

  bottom <- sapply(X = blfdr_vec,
                   FUN = function(cutpt, vec){max(1, sum( vec <= cutpt ))},
                   vec = blfdr_vec)

  top_BC <-  ep +  sapply(X = blfdr_vec,
                          FUN = function(cutpt, vec){ sum( vec <= cutpt )},
                          vec = blfdr_flip_vec)

  FDR_est <- ( top_BC /bottom)
  cutpt_max_BC <- max( blfdr_vec[FDR_est<= alpha])  # this can be -Inf
  rej_index_blfdr_BC  <- which(  blfdr_vec <=    cutpt_max_BC ) # this can be empty

  output = list()
  output$rej_index = rej_index_blfdr_BC
  output$stat = blfdr_vec
  output$mirror_stat = blfdr_flip_vec
  output$param = EM_fit$param
  output$FDR = FDR_est


  return(output)
}
