## this gives the finite sample version of zap
## this will direct give the rejection index

# starting value for s_r should be 0.8< s_r < 1
# starting value for s_l should be 0 <  s_l <  0.2


# x should be a dataset

#' @title Finite-sample ZAP
#' @description Finite-sample version of ZAP based on a 3-component beta mixture
#'   working model.
#' @param z a vector of z-values
#' @param X a matrix of dimension n by p, where  n is the length of z and p is
#'   the number of covariates. p should at least be 1.
#' @param beta_incpt_init a numeric vector of length 2 concatenating initial
#'   values for the INTERCEPT coefficients of the vectors \eqn{\beta_l} and
#'   \eqn{\beta_r} in that order. Default to c(log(5),  log(5)).
#' @param gamma a numeric vector of length 2 concatenating the second shape parameters
#'   \eqn{\gamma_l}, \eqn{\gamma_r}  for the left-leaning and right-leaning beta
#'   densities in that order. Default to c(4, 4).
#' @param s_l,s_r left and right initial thresholding functions, each should be
#'   a vector of the same length as z. If equal to NULL, s_l will default to a
#'   vector with all components equal to 0.2, and s_r will default to a vector
#'   with all components equal to 0.8.
#' @param alpha targeted FDR level.  Default to 0.05.
#' @param maxit maximum number of iterations for every EM update
#' @param tol tolerance
#' @param nfits maximum number of EM updates during the ZAP procedure. Default to 100.
#' @param uvals.cutoff an ``epsilon-small" number. The u-values will be cut off
#'   to be within the interval \deqn{[uvals.cutoff, 1 - uvals.cutoff].} This
#'   stabilizes the beta mixture density estimation. Default is 1e-15.
#' @param ms indicator of whether model selection should be performed. Default
#'   to FALSE. Under construction, please don't touch this option for now.
#' @param x Under construction, please don't touch this option for now.
#' @param df Under construction, please don't touch this option for now.
#' @details TBD
#' @return a vector of indices of the rejected hypotheses
#' @export
#' @examples #TBD


zap_finite = function(z, X,
                      beta_incpt_init = c(log(5),  log(5)),
                      gamma = c(4, 4),
                      s_l = NULL, s_r =  NULL,
                      alpha = 0.05,
                      maxit = 10,
                      tol = 0.0001,
                      nfits = 100,
                      uvals.cutoff = 1e-15,
                      ms = F, x= NULL, df = 6
                      ){

  m = length(z)
  p = ncol(X)
  if (is.null(s_l)) s_l = rep(0.2, m)
  if (is.null(s_r)) s_r = rep(0.8, m)

  U = stats::pnorm(z)
  U[U <= uvals.cutoff] <- uvals.cutoff
  U[U >= 1-  uvals.cutoff] <- 1 -  uvals.cutoff

  r_set = which(U > 0.5) # this indicate where the u value is to the right
  l_set = setdiff( 1:m, r_set)
  # form the mirror image of U that depends on D
  U_mirror = numeric(m)
  U_mirror[r_set] = 1.5 - U[r_set]
  U_mirror[l_set] = 0.5 - U[l_set]
  U_outer = numeric(m)
  U_outer [r_set] = pmax( U[r_set] , U_mirror[r_set] )
  U_outer [l_set] = pmin( U[l_set] , U_mirror[l_set] )
  unmask_set= which((s_l <= U & U<=  (0.5 - s_l))| ((1.5 -   s_r) <= U & U <= s_r ))
  mask_set = setdiff(1:m, unmask_set)
  # index for masked u-values

  # initiate FDP
  rej_index = which( (U <  s_l |   U > s_r) )
  R = length(rej_index)
  A = sum(  U>  (0.5 - s_l) &   U  < (1.5 -   s_r) )
  FDPest = (1 + A)/max(R, 1)

  #### Do model selection step here
  if (ms == F){
    X_tilde <- cbind(1, X)
    beta_init <- c(beta_incpt_init[1], rep(0, p) , beta_incpt_init[2], rep(0, p) )
  }else{
     X <- ns_mod_select_cv(U = U, U_mirror = U_mirror, U_outer = U_outer,
                           extraParam = gamma,
                          l_set = l_set, r_set = r_set,
                          unmask_set = unmask_set,
                          mask_set = mask_set,
                          x = x, df = df, s_l = s_l, s_r = s_r,
                          beta0_incpt_l = beta_incpt_init[1],
                          beta0_incpt_r = beta_incpt_init[2],
                          maxit = maxit, tol = tol, K= 10)
     X_tilde <- cbind(1, X)
     beta_init <- c(beta_incpt_init[1], rep(0, p) , beta_incpt_init[2], rep(0, p) )
  }

  ## begin the while loop if FDPest > alpha
  count = 0
  while(FDPest > alpha & length(mask_set)!= 0 ) {
    if (count%%(m%/%nfits) == 0){
      print(paste("FDP_est = ", FDPest))

      piMat <-  get_piMat_init(unmask_set, mask_set, X, s_l, s_r,  l_set, r_set )

      beta_for_EM <- get_beta_init(U_outer,
                                   X_tilde,  l_set, r_set, gamma,
                                   beta0 =  beta_init)

      print("beta submitted into EM is ")
      print(beta_for_EM)
      param <- EM_finite(U = U , U_mirror = U_mirror,
                        unmask_set = unmask_set,  mask_set = mask_set,
                        maxit = maxit, tol= tol, beta0 = beta_for_EM, piMat = piMat,
                        extraParam = gamma, X_tilde = X_tilde)$param


      linear_predictors_mask <- X_tilde[mask_set, , drop = F]%*%matrix(param, nrow = 1 + p)


      blfdr_mask_outer <-   parallel::mcmapply(FUN = blfdr_stat, U = U_outer[mask_set],
                                               eta1 = linear_predictors_mask[, 1],
                                               eta2 = linear_predictors_mask[, 2],
                                               eta3 = linear_predictors_mask[, 3],
                                               eta4 = linear_predictors_mask[, 4],
                                               MoreArgs = list(extraParam = gamma))


    }


     # print(paste0("mask set length is " , length(blfdr_mask_outer)))
    max_position <- which.max(blfdr_mask_outer)
    new_s <- U_outer[mask_set[max_position] ]
    if (new_s <= 0.5){
      s_l[mask_set[max_position]] = new_s
    }else{
      s_r[mask_set[max_position]] = new_s
    }
    unmask_set= which((s_l <= U & U<=  (0.5 - s_l))| ((1.5 -   s_r) <= U & U <= s_r ))
    mask_set = setdiff(1:m, unmask_set)
    blfdr_mask_outer <- blfdr_mask_outer[- max_position]
    rej_index = which( (U <  s_l |   U > s_r) )
    R = length(rej_index)
    A = sum(  U>  (0.5 - s_l) &   U  < (1.5 -   s_r) )
    FDPest = (1 + A)/max(R, 1)
    count = count +1

  }

  return(rej_index)

}








# inner_bdy[r_set] <-  s_r[r_set]
# inner_bdy[l_set] <-  s_l[l_set]
# blfdr_mask_inner_bdy <-   parallel::mcmapply(FUN = blfdr_stat,
#                                    u_main = inner_bdy[mask_set],
#                                    eta1 = linear_predictors_mask[, 1],
#                                    eta2 = linear_predictors_mask[, 2],
#                                    eta3 = linear_predictors_mask[, 3],
#                                    eta4 = linear_predictors_mask[, 4],
#                                    MoreArgs = list(extraParam = extraParam)
#                                    )


# if (mask_set[max_position] %in% r_set)
#    search_interval = c(U_outer[mask_set[max_position] ] , 1)
# else
#   search_interval = c(0 , U_outer[mask_set[max_position]]  )


# new_s <- uniroot(f = blfdr_stat_shift,
#         eta1 = linear_predictors_mask[max_position, 1],
#         eta2 = linear_predictors_mask[max_position, 2],
#         eta3 = linear_predictors_mask[max_position, 3],
#         eta4 = linear_predictors_mask[max_position, 4],
#         extraParam = extraParam,
#         shift = blfdr_mask_max - 0.000000001,
#         interval = search_interval) $root
#
# blfdr_stat_shift( u_main =U_outer[mask_set[max_position]] ,
#                   eta1 = linear_predictors_mask[max_position, 1],
#                   eta2 = linear_predictors_mask[max_position, 2],
#                   eta3 = linear_predictors_mask[max_position, 3],
#                   eta4 = linear_predictors_mask[max_position, 4],
#                   extraParam = extraParam,
#                   shift = blfdr_mask_max )
#
#
#
# # we can now see which one needs solving for a root
# uniroot_mask_indexm <- which ( )
