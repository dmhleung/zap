#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;




// [[Rcpp::export]]
arma::vec trunc_log_arma(arma::vec U) {
  return(arma::trunc_log(U)) ;
}



// [[Rcpp::export]]
arma::mat createParaMat(arma::vec param, arma::mat X_tilde){
  int k = X_tilde.n_cols;
  int N = X_tilde.n_rows;
  arma::mat M = arma::zeros(N, 4);
  M.col(0) =  1/(1 + exp(-X_tilde*param.subvec(0, k-1)));
  M.col(1) =  1/(1 + exp(-X_tilde*param.subvec(k, 2*k-1)));
  M.col(2) =  exp(X_tilde*param.subvec(2*k, 3*k-1));
  M.col(3) =  exp(X_tilde*param.subvec(3*k, 4*k-1));
  return(M);
}



// [[Rcpp::export]]
arma::mat form_Hmat_unmask(arma::vec U, arma::mat paraMat, arma::vec extraParam){

  arma::vec k_l = paraMat.col(0);
  arma::vec k_r = paraMat.col(1);
  arma::vec exp_Xtheta_l = paraMat.col(2);
  arma::vec exp_Xtheta_r = paraMat.col(3);
  arma::vec B_l = exp(arma::lgamma(k_l) + lgamma(extraParam(0) ) - arma::lgamma(k_l + extraParam(0) ));
  arma::vec B_r = exp(arma::lgamma(k_r ) + lgamma(extraParam(1)) - arma::lgamma(k_r  + extraParam(1)));


  arma::mat pprob = arma::zeros<arma::mat>(U.n_elem, 3);
  pprob.col(0) = exp_Xtheta_l%exp( (k_l -1)%arma::trunc_log(U)+ (extraParam(0) - 1)*arma::trunc_log(1 - U))/B_l ;
  pprob.col(1).ones() ;
  pprob.col(2) = exp_Xtheta_r%exp((extraParam(1) - 1)*arma::trunc_log(U)+ (k_r -1 )%arma::trunc_log(1 - U))/B_r;
  // Here we have omitted the base denominator because it doesn't matter
  pprob =  arma::normalise(pprob, 1, 1);
  return(pprob);
}





// [[Rcpp::export]]
arma::mat form_Hmat_unmask_init(arma::vec U, arma::vec beta0, arma::mat piMat,
                                arma::mat X_tilde, arma::vec extraParam){

  int k = X_tilde.n_cols ;
  arma::vec k_l = 1/(1 + exp(-X_tilde*beta0.subvec(0, k-1)));
  arma::vec k_r = 1/(1 + exp(-X_tilde*beta0.subvec(k, 2*k-1)));
  arma::vec B_l = exp(lgamma(k_l) + lgamma(extraParam(0) ) - lgamma(k_l + extraParam(0) ));
  arma::vec B_r = exp(lgamma(k_r ) + lgamma(extraParam(1)) - lgamma(k_r  + extraParam(1)));


  arma::mat pprob = arma::zeros<arma::mat>(U.n_elem, 3);
  pprob.col(0) = piMat.col(0)%exp( (k_l -1)%arma::trunc_log(U)+ (extraParam(0) - 1)*arma::trunc_log(1 - U))/B_l ;
  pprob.col(1) = piMat.col(1);
  pprob.col(2) = piMat.col(2)%exp((extraParam(1) - 1)*arma::trunc_log(U)+ (k_r -1 )%arma::trunc_log(1 - U))/B_r;
  // Here we have omitted the base denominator because it doesn't matter
  pprob =  normalise(pprob, 1, 1);
  return(pprob);
}



// [[Rcpp::export]]
Rcpp::List form_Hmat_Yab_mask_init(arma::vec U_mask, arma::vec U_mirror_mask, arma::vec beta0,
                                   arma::mat piMat, arma::mat X_tilde, arma::vec extraParam){

  arma::vec U= U_mask;
  arma::vec U_mirror= U_mirror_mask;
  //redefine the names for simplicity

  int k = X_tilde.n_cols ;
  arma::vec k_l = 1/(1 + exp(-X_tilde*beta0.subvec(0, k-1)));
  arma::vec k_r = 1/(1 + exp(-X_tilde*beta0.subvec(k, 2*k-1)));
  arma::vec B_l = exp(lgamma(k_l) + lgamma(extraParam(0) ) - lgamma(k_l + extraParam(0) ));
  arma::vec B_r = exp(lgamma(k_r ) + lgamma(extraParam(1)) - lgamma(k_r  + extraParam(1)));

  arma::vec h_l = exp( (k_l -1)%arma::trunc_log(U)+ (extraParam(0) - 1)*arma::trunc_log(1 - U))/B_l;
  arma::vec h_l_mirror = exp( (k_l -1)%arma::trunc_log(U_mirror)+ (extraParam(0) - 1)*arma::trunc_log(1 - U_mirror))/B_l;
  arma::vec h_r = exp((extraParam(1) - 1)*arma::trunc_log(U)+ (k_r -1 )%arma::trunc_log(1 - U))/B_r;
  arma::vec h_r_mirror = exp((extraParam(1) - 1)*arma::trunc_log(U_mirror)+ (k_r -1 ) %arma::trunc_log(1 - U_mirror))/B_r;

  arma::mat pprob = arma::zeros<arma::mat>(U.n_elem, 3);
  pprob.col(0) = piMat.col(0)%(h_l + h_l_mirror) ;
  pprob.col(1) = 2*piMat.col(1) ;
  pprob.col(2) = piMat.col(2)%(h_r + h_r_mirror);
  pprob =  arma::normalise(pprob, 1, 1);

  arma::mat Ymat = arma::zeros<arma::mat>(U.n_elem, 4);
  arma::vec Yl_a = (h_l%arma::trunc_log(U) + h_l_mirror%arma::trunc_log(U_mirror))/(h_l + h_l_mirror);
  arma::vec Yl_b = (h_l%arma::trunc_log(1 - U) + h_l_mirror%arma::trunc_log(1 - U_mirror))/(h_l + h_l_mirror);
  arma::vec Yr_a = (h_r%arma::trunc_log(U) + h_r_mirror%arma::trunc_log(U_mirror))/(h_r + h_r_mirror);
  arma::vec Yr_b = (h_r%arma::trunc_log(1 - U) + h_r_mirror%arma::trunc_log(1 - U_mirror))/(h_r + h_r_mirror);
  Ymat.col(0) = Yl_a ;
  Ymat.col(1) = Yl_b ;
  Ymat.col(2) = Yr_a ;
  Ymat.col(3) = Yr_b ;

  return Rcpp::List::create(Rcpp::Named("Hmat")=pprob,
                            Rcpp::Named("Ymat")=Ymat);
}




// [[Rcpp::export]]
Rcpp::List form_Hmat_Yab_mask(arma::vec U_mask, arma::vec U_mirror_mask,
                              arma::mat paraMat, arma::vec extraParam){

  arma::vec U= U_mask;
  arma::vec U_mirror= U_mirror_mask;
  //redefine the names for simplicity

  arma::vec k_l = paraMat.col(0);
  arma::vec k_r = paraMat.col(1);
  arma::vec exp_Xtheta_l = paraMat.col(2);
  arma::vec exp_Xtheta_r = paraMat.col(3);
  arma::vec B_l = exp(arma::lgamma(k_l) + lgamma(extraParam(0) ) - lgamma(k_l + extraParam(0) ));
  arma::vec B_r = exp(arma::lgamma(k_r ) + lgamma(extraParam(1)) - lgamma(k_r  + extraParam(1)));

  arma::vec h_l = exp( (k_l -1)%arma::trunc_log(U)+ (extraParam(0) - 1)*arma::trunc_log(1 - U))/B_l;
  arma::vec h_l_mirror = exp( (k_l -1)%arma::trunc_log(U_mirror)+ (extraParam(0) - 1)*arma::trunc_log(1 - U_mirror))/B_l;
  arma::vec h_r = exp((extraParam(1) - 1)*arma::trunc_log(U)+ (k_r -1 ) %arma::trunc_log(1 - U))/B_r;
  arma::vec h_r_mirror = exp((extraParam(1) - 1)*arma::trunc_log(U_mirror)+ (k_r -1 ) %arma::trunc_log(1 - U_mirror))/B_r;

  arma::mat pprob = arma::zeros<arma::mat>(U.n_elem, 3);
  pprob.col(0) = exp_Xtheta_l%(h_l + h_l_mirror) ;
  pprob.col(1).fill(2) ;
  pprob.col(2) = exp_Xtheta_r%(h_r + h_r_mirror);
  pprob =  normalise(pprob, 1, 1);

  arma::mat Ymat = arma::zeros<arma::mat>(U.n_elem, 4);
  arma::vec Yl_a = (h_l%arma::trunc_log(U) + h_l_mirror%arma::trunc_log(U_mirror))/(h_l + h_l_mirror);
  arma::vec Yl_b = (h_l%arma::trunc_log(1 - U) + h_l_mirror%arma::trunc_log(1 - U_mirror))/(h_l + h_l_mirror);
  arma::vec Yr_a = (h_r%arma::trunc_log(U) + h_r_mirror%arma::trunc_log(U_mirror))/(h_r + h_r_mirror);
  arma::vec Yr_b = (h_r%arma::trunc_log(1 - U) + h_r_mirror%arma::trunc_log(1 - U_mirror))/(h_r + h_r_mirror);
  Ymat.col(0) = Yl_a ;
  Ymat.col(1) = Yl_b ;
  Ymat.col(2) = Yr_a ;
  Ymat.col(3) = Yr_b ;

  return Rcpp::List::create(Rcpp::Named("Hmat")=pprob,
                            Rcpp::Named("Ymat")=Ymat);
}



// [[Rcpp::export]]
arma::mat overall_Hmat(arma::mat Hmat_mask, arma::mat Hmat_unmask, arma::uvec mask_set, arma::uvec unmask_set){

  arma::mat Hmat =  arma::zeros<arma::mat>(Hmat_mask.n_rows + Hmat_unmask.n_rows, 3);
  Hmat.rows(unmask_set-1 ) = Hmat_unmask;
  Hmat.rows(mask_set-1 ) = Hmat_mask;
  return Hmat;
}


// [[Rcpp::export]]
arma::mat overall_Ymat(arma::mat Ymat_mask, arma::vec U_unmask, arma::uvec mask_set, arma::uvec unmask_set){
  int n = mask_set.n_elem+ unmask_set.n_elem;

  arma::mat Ymat =  arma::zeros<arma::mat>(n, 4);
  arma::vec Yl_a(n);
  arma::vec Yl_b(n);
  arma::vec Yr_a(n);
  arma::vec Yr_b(n);

  Yl_a(unmask_set-1) = arma::trunc_log(U_unmask);
  Yl_a(mask_set-1) = Ymat_mask.col(0);
  Yl_b(unmask_set-1) = arma::trunc_log(1 - U_unmask);
  Yl_b(mask_set-1) = Ymat_mask.col(1);


  Yr_a(unmask_set-1) = arma::trunc_log(U_unmask);
  Yr_a(mask_set-1) = Ymat_mask.col(2);
  Yr_b(unmask_set-1) = arma::trunc_log(1 - U_unmask);
  Yr_b(mask_set-1) = Ymat_mask.col(3);


  Ymat.col(0) = Yl_a;
  Ymat.col(1) = Yl_b;
  Ymat.col(2) = Yr_a;
  Ymat.col(3) = Yr_b;

  return Ymat;
}


