#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

//[[Rcpp::export]]

double blfdr_stat(double U,
                  double eta1, double eta2,
                  double eta3, double eta4,
                  arma::vec extraParam){

  double k_l = 1/(1 +exp(- eta1));
  double k_r = 1/(1 + exp(-eta2));
  double exp_Xtheta_l = exp(eta3);
  double exp_Xtheta_r = exp(eta4);
  double B_l = exp(lgamma(k_l) + lgamma(extraParam(0) ) - lgamma(k_l + extraParam(0) ));
  double B_r = exp(lgamma(k_r) + lgamma(extraParam(1)) - lgamma(k_r + extraParam(1)));

  double blfdr_stat =   1/( 1+
                            exp_Xtheta_l*exp( (k_l -1)*arma::trunc_log(U )+
                            (extraParam(0) - 1)*arma::trunc_log(1 - U))/B_l +
                            exp_Xtheta_r*exp((extraParam(1) - 1)*arma::trunc_log(U)+
                            (k_r -1 ) *arma::trunc_log(1 - U))/B_r
  );

  return(blfdr_stat);
}
