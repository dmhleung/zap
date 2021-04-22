# A model selection function that returns the "best" B-spline basis matrix for a natural cubic spline

# Input: 
#1. univeriat covariate x (numeric vector)
#2. the u values and their mirrors (two numeric vectors)
#3. unmasked and masked sets


# Output: a B-spline basis matrix



# use caret::createFolds to create partitions


ns_mod_select_cv <- function(U, U_mirror, U_outer, extraParam,
                              l_set, r_set,  unmask_set,  mask_set,
                              x, df, s_l, s_r, 
                              beta0_incpt_l = log(5),
                              beta0_incpt_r = log(5),
                              maxit, tol, K= 10){
  print(paste0("Attempting model selection with ", K, "-fold-CV" ))
  # set.seed(130)
  m <- length(U)
  partition <- caret::createFolds(1:m, k = K)
  # LogLike_vec <- numeric(df)
  LogLike_mat <- matrix(0, nr = K, nc = df)
  
  LogLike_no_covar <- numeric(K)
  for (k in 1:K){
      print(paste0("Performing model selection (CV) for fold " , k, " of the data"  ))
      test_set_index <-  partition[[k]]  
      train_set_index <- setdiff(1:m, test_set_index)
      
      s_l_train <- s_l[train_set_index]
      s_r_train <- s_r[train_set_index]
      x_train <- x[train_set_index, , drop = F]  # this's still a dataframe
      x_test <- x[test_set_index, , drop = F]# this's still a dataframe
      l_set_train <- which(train_set_index %in% l_set)
      r_set_train <- which(train_set_index %in% r_set)
      unmask_set_train <- which(train_set_index %in% unmask_set)
      mask_set_train <- which(train_set_index %in% mask_set)
      U_train <- U[train_set_index]
      U_mirror_train <- U_mirror[train_set_index]
      U_outer_train <- U_outer[train_set_index]
      
      # each element should be a matrix for the following two lists
      X_ls_test <- ns_expand(x_test , df)  
      X_tilde_ls_test <- lapply(X = X_ls_test, FUN  = function(X) {cbind(1, X)})
 
      param_ls<- 
        ns_mod_select_helper(U = U_train, 
                             U_mirror = U_mirror_train, 
                             U_outer = U_outer_train, 
                             extraParam = extraParam,
                             l_set = l_set_train, r_set = r_set_train,  
                             s_l = s_l_train, s_r = s_r_train, 
                             unmask_set = unmask_set_train,  
                             mask_set = mask_set_train,
                              x = x_train, df = df, 
                             beta0_incpt_l =beta0_incpt_l ,
                              beta0_incpt_r =beta0_incpt_r,
                               maxit = maxit, tol = tol)$param_ls
    
      paraMat_ls_test <-parallel::mcmapply(
        FUN = createParaMat,
        param = param_ls,
        X_tilde =X_tilde_ls_test, 
        mc.cores = 8, SIMPLIFY = F)
      
      LogLike_mat[k, ] <- parallel::mcmapply(
          FUN = LogLike_finite,
          paraMat = paraMat_ls_test, 
          mc.cores = 8,
          MoreArgs = list(U = U[test_set_index], U_mirror = U_mirror[test_set_index],  
                          unmask_set = which(test_set_index  %in% unmask_set), 
                          mask_set = which(test_set_index %in% mask_set),
                          extraParam = extraParam))
 
  }
  
  LogLike_mean <- apply(X = LogLike_mat, MARGIN = 2, FUN = mean)

  LogLike_sd <-  apply(X = LogLike_mat, MARGIN = 2, FUN = sd)
  max_index <- which.max(LogLike_mean)
  lowest_thrshold <- LogLike_mean[max_index] - LogLike_sd[max_index]
  
  chosen_df <- min(which( LogLike_mean > lowest_thrshold))
  print(paste0("the chosen df is ", chosen_df))
  
   return(ns_expand(x , df )[[chosen_df]])
}


