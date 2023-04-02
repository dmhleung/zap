#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]



using namespace Rcpp;



//[[Rcpp::export]]
double neg_LogLike_LeftBeta(arma::vec beta_l, arma::vec U, arma::mat X_tilde, double gamma_l){
  arma::vec k_l = 1/(1 + exp(-X_tilde*beta_l));
  arma::vec log_B_l = lgamma(k_l) + lgamma(gamma_l ) - lgamma(k_l + gamma_l );
  arma::vec loglike_vec = (k_l -1)%arma::trunc_log(U )+ (gamma_l - 1)*arma::trunc_log(1 - U) - log_B_l ;
  return(-accu(loglike_vec)) ;
}



// [[Rcpp::export]]
arma::vec neg_LogLike_LeftBeta_grad(arma::vec beta_l, arma::vec U, arma::mat X_tilde, double gamma_l) {

  arma::vec exp_neg_X_beta_l  =  exp(-X_tilde*beta_l)   ;

  arma::vec k_l = 1/(1 + exp_neg_X_beta_l);
  NumericVector k_l_rcpp = wrap( k_l   );
  arma::vec digamma_k_l  = digamma(k_l_rcpp);
  arma::vec digamma_k_l_gamma_l  =  digamma(k_l_rcpp + gamma_l);

  arma::vec Df_Dk_l = trunc_log(U) - digamma_k_l + digamma_k_l_gamma_l  ;
  arma::mat Dk_l_Dbeta_l = X_tilde.each_col() % (square(k_l)%exp_neg_X_beta_l);

  arma::vec gradient = trans(Dk_l_Dbeta_l)*Df_Dk_l;
  return(-gradient) ;
}



//[[Rcpp::export]]
double neg_LogLike_RightBeta(arma::vec beta_r, arma::vec U, arma::mat X_tilde, double gamma_r){
  arma::vec k_r = 1/(1 + exp(-X_tilde*beta_r));
  arma::vec log_B_r = lgamma(k_r) + lgamma(gamma_r ) - lgamma(k_r + gamma_r );
  arma::vec loglike_vec = (gamma_r - 1)*arma::trunc_log(U)+   (k_r -1 )%arma::trunc_log(1 - U) - log_B_r ;
  return(-accu(loglike_vec)) ;
}



// negative right beta log-likelihood gradient
//
// This is a function that evaluate the right leaning beta negative log likelihood gradient for a vector of u-values
//param beta_r The first shape parameter
//param U A vector of u-values to evaluate the beta density on
//param X_tilde A matrix of the covariates augmented with 1's as intercepts. Should have the same number of rows as the length of U
//param gamma_r The second shape parameter, should be greater than 2
//return the evaluated negative log likelihood gradient, a vector of dimension p

// [[Rcpp::export]]
arma::vec neg_LogLike_RightBeta_grad(arma::vec beta_r, arma::vec U, arma::mat X_tilde, double gamma_r) {

  arma::vec exp_neg_X_beta_r  =  exp(-X_tilde*beta_r)   ;

  arma::vec k_r = 1/(1 + exp_neg_X_beta_r);
  NumericVector k_r_rcpp = wrap( k_r   );
  arma::vec digamma_k_r  = digamma(k_r_rcpp);
  arma::vec digamma_k_r_gamma_r  =  digamma(k_r_rcpp + gamma_r);

  arma::vec Df_Dk_r = trunc_log(1 - U) - digamma_k_r + digamma_k_r_gamma_r  ;
  arma::mat Dk_r_Dbeta_r = X_tilde.each_col() % (square(k_r)%exp_neg_X_beta_r);

  arma::vec gradient = trans(Dk_r_Dbeta_r)*Df_Dk_r;
  return(-gradient) ;
}



// log likelihood of the 3-component beta mixture evaluated at U (vector of u-values)
// [[Rcpp::export]]
double LogLike_asymp(arma::vec U, arma::mat paraMat, arma::vec extraParam){

  arma::vec k_l = paraMat.col(0);
  arma::vec k_r = paraMat.col(1);
  arma::vec exp_Xtheta_l = paraMat.col(2);
  arma::vec exp_Xtheta_r = paraMat.col(3);
  arma::vec base =  1 + exp_Xtheta_l + exp_Xtheta_r ;
  arma::vec B_l = exp(lgamma(k_l) + lgamma(extraParam(0) ) - lgamma(k_l + extraParam(0) ));
  arma::vec B_r = exp(lgamma(k_r) + lgamma(extraParam(1)) - lgamma(k_r + extraParam(1)));

  arma::vec like_vec =  1 +
    exp_Xtheta_l%exp( (k_l -1)%trunc_log(U)+ (extraParam(0) - 1)*trunc_log(1 - U))/B_l +
    exp_Xtheta_r%exp((extraParam(1) - 1)*trunc_log(U)+ (k_r -1 )%trunc_log(1 - U))/B_r;

  return(accu(trunc_log(like_vec) - trunc_log(base)));
}


// negative Q function for the beta part when running zap (finite)
// [[Rcpp::export]]
double neg_Q_fn_finite_beta_part(arma::vec param_beta, arma::mat Ymat, arma::mat Hmat,
                                 arma::vec extraParam, arma::mat X_tilde) {


  arma::vec Yl_a = Ymat.col(0);
  arma::vec Yl_b = Ymat.col(1);
  arma::vec Yr_a = Ymat.col(2);
  arma::vec Yr_b = Ymat.col(3);

  int k = X_tilde.n_cols;
  arma::vec k_l = 1/(1 + exp(-X_tilde*param_beta.subvec(0, k - 1)));
  arma::vec k_r = 1/(1 + exp(-X_tilde*param_beta.subvec(k, 2*k-1)));
  // vec log_B_l = lgamma(k_l) + lgamma(extraParam(0)) - lgamma(k_l + extraParam(0) );
  // vec log_B_r = lgamma(k_r) + lgamma(extraParam(1)) - lgamma(k_r + extraParam(1));
  arma::vec log_B_l = lgamma(k_l) - lgamma(k_l + extraParam(0) );
  arma::vec log_B_r = lgamma(k_r) - lgamma(k_r + extraParam(1));


  /*
   vec loglike_vec =
   Hmat.col(0)%( (k_l -1)%Yl_a + (extraParam(0) - 1)*Yl_b - log_B_l ) +
   Hmat.col(2)%((extraParam(1) - 1)*Yr_a + (k_r - 1)%Yr_b - log_B_r);
   */

  arma::vec loglike_vec =
    Hmat.col(0)%( (k_l -1)%Yl_a  - log_B_l ) +
    Hmat.col(2)%((k_r - 1)%Yr_b - log_B_r);

  return(-accu(loglike_vec)) ;
}



// gradient of negative Q function for the beta part when running zap (finite)
// [[Rcpp::export]]
arma::vec neg_Q_fn_finite_beta_part_grad_cpp(arma::vec param_beta, arma::mat Ymat, arma::mat Hmat,
                                             arma::vec extraParam, arma::mat X_tilde) {


  arma::vec Yl_a = Ymat.col(0);
  arma::vec Yl_b = Ymat.col(1);
  arma::vec Yr_a = Ymat.col(2);
  arma::vec Yr_b = Ymat.col(3);
  int k = X_tilde.n_cols;

  arma::vec exp_neg_X_beta_l  =  exp(-X_tilde*param_beta.subvec(0, k - 1))   ;
  arma::vec exp_neg_X_beta_r  =   exp(-X_tilde*param_beta.subvec(k, 2*k-1))   ;
  arma::vec k_l = 1/(1 + exp_neg_X_beta_l);
  arma::vec k_r = 1/(1 + exp_neg_X_beta_r) ;
  NumericVector k_l_rcpp = wrap( k_l   );
  NumericVector k_r_rcpp =   wrap( k_r    );

  arma::vec digamma_k_l  = digamma(k_l_rcpp);
  arma::vec digamma_k_l_gamma_l  =  digamma(k_l_rcpp + extraParam(0));
  arma::vec digamma_k_r  = digamma(k_r_rcpp );
  arma::vec digamma_k_r_gamma_r  =  digamma(k_r_rcpp  + extraParam(1));


  arma::vec Df_Dk_l = Hmat.col(0)%(Yl_a - digamma_k_l + digamma_k_l_gamma_l ) ;
  arma::vec Df_Dk_r = Hmat.col(2)%(Yr_b - digamma_k_r + digamma_k_r_gamma_r ) ;

  arma::mat Dk_l_Dbeta_l = X_tilde.each_col() % (square(k_l)%exp_neg_X_beta_l);

  arma::mat Dk_r_Dbeta_r =  X_tilde.each_col() %   (square(k_r)%exp_neg_X_beta_r);

  arma::vec gradient(2*k);

  gradient.subvec(0, k-1)  = trans(Dk_l_Dbeta_l)*Df_Dk_l ;
  gradient.subvec(k, 2*k-1)  = trans(Dk_r_Dbeta_r )*Df_Dk_r  ;

  return(-gradient) ;
}



// [[Rcpp::export]]
double LogLike_finite(arma::vec U,  arma::vec U_mirror,  arma::mat paraMat, arma::vec extraParam, arma::uvec mask_set, arma::uvec unmask_set){

  arma::vec U_unmask = U(unmask_set - 1);
  arma::vec U_mask = U(mask_set - 1);
  arma::vec U_mirror_mask = U_mirror(mask_set - 1);

  arma::vec k_l = paraMat.col(0);
  arma::vec k_r = paraMat.col(1);
  arma::vec exp_Xtheta_l = paraMat.col(2);
  arma::vec exp_Xtheta_r = paraMat.col(3);
  arma::vec base =  1 + exp_Xtheta_l + exp_Xtheta_r ;
  arma::vec B_l = exp(lgamma(k_l) + lgamma(extraParam(0) ) - lgamma(k_l + extraParam(0) ));
  arma::vec B_r = exp(lgamma(k_r) + lgamma(extraParam(1)) - lgamma(k_r + extraParam(1)));

  arma::vec k_l_unmask = k_l(unmask_set - 1);
  arma::vec k_r_unmask = k_r(unmask_set - 1);
  arma::vec exp_Xtheta_l_unmask = exp_Xtheta_l(unmask_set - 1);
  arma::vec exp_Xtheta_r_unmask = exp_Xtheta_r(unmask_set - 1);
  arma::vec B_l_unmask = B_l(unmask_set - 1);
  arma::vec B_r_unmask = B_r(unmask_set - 1);

  arma::vec k_l_mask = k_l(mask_set - 1);
  arma::vec k_r_mask = k_r(mask_set - 1);
  arma::vec exp_Xtheta_l_mask = exp_Xtheta_l(mask_set - 1);
  arma::vec exp_Xtheta_r_mask = exp_Xtheta_r(mask_set - 1);
  arma::vec B_l_mask = B_l(mask_set - 1);
  arma::vec B_r_mask = B_r(mask_set - 1);

  arma::vec like_vec_unmask =  1 +
    exp_Xtheta_l_unmask%exp( (k_l_unmask -1)%trunc_log(U_unmask)+ (extraParam(0) - 1)*trunc_log(1 - U_unmask))/B_l_unmask +
    exp_Xtheta_r_unmask%exp((extraParam(1) - 1)*trunc_log(U_unmask)+ (k_r_unmask -1 )%trunc_log(1 - U_unmask))/B_r_unmask;

  arma::vec like_vec_mask =  2 +
    exp_Xtheta_l_mask%exp( (k_l_mask -1)%trunc_log(U_mask)+ (extraParam(0) - 1)*trunc_log(1 - U_mask))/B_l_mask +
    exp_Xtheta_l_mask%exp( (k_l_mask -1)%trunc_log(U_mirror_mask)+ (extraParam(0) - 1)*trunc_log(1 - U_mirror_mask))/B_l_mask +
    exp_Xtheta_r_mask%exp((extraParam(1) - 1)*trunc_log(U_mask)+ (k_r_mask  -1 )%trunc_log(1 - U_mask))/B_r_mask +
    exp_Xtheta_r_mask%exp((extraParam(1) - 1)*trunc_log(U_mirror_mask)+ (k_r_mask -1 )%trunc_log(1 - U_mirror_mask))/B_r_mask;

  return(accu(trunc_log(like_vec_unmask)) + accu(trunc_log(like_vec_mask)) - accu(trunc_log(base))  );
}

