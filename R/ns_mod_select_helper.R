
  
  
  
  
ns_mod_select_helper <- function(U, U_mirror, U_outer, extraParam,
                               l_set, r_set,  unmask_set,  mask_set,
                               x, df, s_l, s_r, 
                               beta0_incpt_l = log(5),
                               beta0_incpt_r = log(5),
                               maxit, tol){
  
  X_ls <- ns_expand(x , df)
  beta0_ls <- list()
  for (i in 1:df){
   
    beta0_ls[[i]] <- c(beta0_incpt_l, rep(0, i*ncol(x)),  beta0_incpt_r, rep(0, i*ncol(x)))
  }
  X_tilde_ls<- lapply(X = X_ls, FUN  = function(X) {cbind(1, X)})


  piMat_ls <- parallel::mcmapply(
    FUN = get_piMat_init,
    X = X_ls,
    mc.cores = 8,
    SIMPLIFY = F, 
    MoreArgs =list(
    unmask_set = unmask_set, 
    mask_set = mask_set,
    s_l = s_l, s_r = s_r,  
    l_set = l_set, r_set = r_set)
    )
  
  beta0_ls <-  
    parallel::mcmapply(
      FUN = get_beta_init,
      X_tilde = X_tilde_ls, 
      beta0 = beta0_ls, 
      mc.cores = 8,
      SIMPLIFY = F, 
      MoreArgs = list(U_outer = U_outer, 
                      l_set = l_set, r_set =r_set,
                      extraParam = extraParam))

  param_ls <- 
    parallel::mcmapply(
      FUN = EM_finite,
      X_tilde = X_tilde_ls,
      piMat = piMat_ls,
      beta0 = beta0_ls,
      mc.cores = 8,
      MoreArgs = list(U = U, U_mirror = U_mirror,  
                      unmask_set = unmask_set, 
                      mask_set = mask_set,
                      maxit = maxit, tol= tol,
                      extraParam = extraParam))
  
  return_ls = list()
  return_ls$param_ls = param_ls

  return(return_ls)
  
  
}