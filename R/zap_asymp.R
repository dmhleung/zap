
## something



#' @title Asymptotic ZAP
#'
#' @description Asymptotic version of ZAP based on a 3-component beta mixture
#'   working model.
#'
#' @param z a vector of z-values
#' @param X a matrix of dimension n by p, where  n is the length of z and p is
#'   the number of covariates. p should at least be 1.
#' @param beta_incpt_init a numeric vector of length 2 concatenating initial values
#'   for the INTERCEPT coefficients of the vectors \eqn{\beta_l} and
#'   \eqn{\beta_r} in that order. Default to c(log(5),  log(5)).
#' @param theta_incpt_init a numeric vector of length 2 concatenating initial values
#'   for the INTERCEPT coefficients of the vectors \eqn{\theta_l} and
#'   \eqn{\theta_r} in that order. Default to c(log(0.1),  log(0.1)).
#' @param gamma a numeric vector of length 2 concatenating the second shape parameters
#'   \eqn{\gamma_l}, \eqn{\gamma_r}  for the left-leaning and right-leaning beta
#'   densities in that order. Default to c(4, 4).
#' @param gen_size number of random uniformly distributed realizations generated
#'   to estimate the mirror statistics. Default to 50000.
#' @param alpha  targeted FDR level. Default to 0.05.
#' @param ep additive factor in the numerator of the FDP estimates in ZAP.
#'   Default is 1. Changing it may result in an inflated FDR of the procedure.
#' @param uvals.cutoff an ``epsilon-small" number. The u-values will be cut off
#'   to be within the interval \deqn{[uvals.cutoff, 1 - uvals.cutoff].} This
#'   stabilizes the beta mixture density estimation. Default is 1e-15.
#' @details TBD
#' @return a list with containing the followings: \itemize{ \item rej_index: a
#'   vector of indices corresponding to the rejected hypotheses. Can be of
#'   length 0 if no hypotheses are rejected. \item stat: a vector of resulting
#'   test statistics \item mirror_stat: a vector of "mirror" test statistics
#'   \item param: a vector of estimated parameters concatenating  \eqn{\beta_l}, \eqn{\beta_r}, \eqn{\theta_l} and
#'   \eqn{\theta_r} in that order }
#' @export
#' @examples


zap_asymp = function(z,X,
                     beta_incpt_init = c(log(5),  log(5)),
                     theta_incpt_init =  c(log(.1), log(.1)),
                     gamma = c(4, 4),
                     alpha = 0.05,
                     gen_size= 50000,
                     maxit = 100,
                     tol = 0.0001,
                     ep = 1,
                     uvals.cutoff = 1e-15){
  p <- ncol(X)
  beta_init <- c(beta_incpt_init[1], rep( 0, p) , beta_incpt_init[2],  rep( 0, p))
  theta_init <- c(theta_incpt_init[1], rep( 0, p) , theta_incpt_init[2],  rep( 0, p))


  EM_fit <- EM_asymp(z, param = c(beta_init, theta_init) ,
                       tol = tol,
                       maxit = maxit,
                       extraParam = gamma,
                       X_tilde = cbind(1, X),
                     uvals.cutoff = uvals.cutoff)

  blfdr_parts <-create_blfdr_parts(u = pnorm(z),
                                   param = EM_fit$param,
                                   extraParam = gamma,
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
  return(output)
}
