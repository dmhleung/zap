#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


//[[Rcpp::export]]

arma::vec blfdr_fn(double u_main, arma::vec unif_dat,
             double eta1, double eta2,
             double eta3, double eta4,
             arma::vec extraParam){

  double k_l = 1/(1 +exp(- eta1));
  double k_r = 1/(1 + exp(-eta2));
  double exp_Xtheta_l = exp(eta3);
  double exp_Xtheta_r = exp(eta4);
  double B_l = exp(std::lgamma(k_l) + std::lgamma(extraParam(0) ) - std::lgamma(k_l + extraParam(0) ));
  double B_r = exp(std::lgamma(k_r) + std::lgamma(extraParam(1)) - std::lgamma(k_r + extraParam(1)));

  arma::vec whole_blfdr_dat_vec =
    1/( 1+
    exp_Xtheta_l*exp( (k_l -1)*arma::trunc_log(unif_dat )+ (extraParam(0) - 1)*arma::trunc_log(1 - unif_dat))/B_l +
    exp_Xtheta_r*exp((extraParam(1) - 1)*arma::trunc_log(unif_dat)+ (k_r -1 ) *arma::trunc_log(1 - unif_dat))/B_r
    );

  double blfdr_main_stat =   1/( 1+
                                 exp_Xtheta_l*exp( (k_l -1)*arma::trunc_log(u_main )+
                                 (extraParam(0) - 1)*arma::trunc_log(1 - u_main))/B_l +
                                 exp_Xtheta_r*exp((extraParam(1) - 1)*arma::trunc_log(u_main)+
                                 (k_r -1 ) *arma::trunc_log(1 - u_main))/B_r
  );

  double sum = accu(whole_blfdr_dat_vec < blfdr_main_stat);
  double p = 1 - sum/(unif_dat.n_elem);
  arma::vec P(1);
  P.fill(p);
  //the probability value

  arma::vec  blfdr_flip = arma::quantile(whole_blfdr_dat_vec, P );

  arma::vec result(2) ;
  result(0) = blfdr_main_stat;
  result(1) = blfdr_flip(0);
  return(result);
}




