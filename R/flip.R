


# a function with
# input: test statistics, list of "shell" functions, number of std nomrals generated
# shell_fn_ls has the same length as stat
# returns the "flip" quantile

# Along all the functions in this file, we will only use create_blfdr_parts() for our package, for now



create_blfdr_parts= function(u, param, extraParam,  X, gen_size){
  linear_predictors = cbind(1, X)%*%matrix(param, nrow = 1 + ncol(X))
  blfdr_mat =   parallel::mcmapply(FUN = blfdr_fn,
                                   u_main = u,
                                   eta1 = linear_predictors[, 1],
                                   eta2 = linear_predictors[, 2],
                                   eta3 = linear_predictors[, 3],
                                   eta4 = linear_predictors[, 4],
                                   MoreArgs = list(extraParam = extraParam,
                                                   unif_dat = stats::runif(gen_size)
                                                   ),
                                   mc.cores = 8)

  return(list( blfdr_vec =  blfdr_mat[1, ],
               blfdr_flip_vec = blfdr_mat[2, ]))
}



# ## this get the distribution probabilities for the BH variation
# blfdr_dp_sum_vec <- function(blfdr_vec,
#                             param,
#                             extraParam,
#                             X, gen_size){
#   N = nrow(X)
#   linear_predictors = cbind(1, X)%*%matrix(param, nrow = 1 + ncol(X))
#
#     result <- parallel::mcmapply(
#               FUN = blfdr_dp_vec,
#               eta1 = linear_predictors[, 1],
#               eta2 = linear_predictors[, 2],
#               eta3 = linear_predictors[, 3],
#               eta4 = linear_predictors[, 4],
#               MoreArgs = list(extraParam = extraParam,
#                               unif_dat = stats::runif(gen_size),
#                               cutoff_vec = blfdr_vec),
#               mc.cores = 8)
#
#   return(rowSums(result))
#
# }
#
# ### will not be used at the moment
#
#
# get_blfdr_stuff = function(u_main_vec,
#                            param,
#                            extraParam,
#                            X, gen_size){
#   N = nrow(X)
#   linear_predictors = cbind(1, X)%*%matrix(param, nrow = 1 + ncol(X))
#   unif_dat = stats::runif(gen_size)
#   blfdr_dp_vec = numeric(N)
#   blfdr_flip_vec = numeric(N)
#   blfdr_vec = parallel::mcmapply(FUN = blfdr_stat,
#                                  U = u_main_vec,
#                                  eta1 = linear_predictors[, 1],
#                                  eta2 = linear_predictors[, 2],
#                                  eta3 = linear_predictors[, 3],
#                                  eta4 = linear_predictors[, 4],
#                                  MoreArgs = list(extraParam = extraParam),
#                                  mc.cores = 8)
#   blfdr_vec = as.numeric(blfdr_vec)
#
#   for (j in 1:N){  # each loop corresponds to q info about one statistic
#     print(j)
#     # this is a vector of length gen_size
#     blfdr_null_dat_vec = generate_blfdr_null_dat_vec(
#       unif_dat = unif_dat,
#       eta1 = linear_predictors[j, 1],
#       eta2 = linear_predictors[j, 2],
#       eta3 = linear_predictors[j, 3],
#       eta4 = linear_predictors[j, 4],
#       extraParam = extraParam
#     )
#     # blfdr_flip_vec[j] <- blfdr_flip_stat(blfdr_vec[j], blfdr_null_dat_vec)
#     blfdr_dp_vec_update = parallel::mcmapply(
#       FUN = count_less_than_cutoff_prop,
#       cutoff = blfdr_vec,
#       MoreArgs = list(blfdr_null_dat_vec = blfdr_null_dat_vec),
#       mc.cores = 8)
#     blfdr_dp_vec = blfdr_dp_vec + blfdr_dp_vec_update
#   }
#   return_list  = list(blfdr_vec = blfdr_vec,
#                       # blfdr_flip_vec = blfdr_flip_vec,
#                       blfdr_dp_vec = blfdr_dp_vec)
#   return(return_list)
# }
#
#
#
#
#
# ########################################################
# ########################################################
# ########################################################
# ##create a list of functions
# # para = parameters(flexmix_obj2 )
# #  null_col = which.min(abs(para["coef.(Intercept)", ] ))
# # concom_para = attributes(attributes(flexmix_obj2 )$concomitant)$coef
# # prob_element = exp(cbind(1, X)%*%concom_para)
# # prob = prob_element*(1/rowSums( prob_element))
# create_plfdr_fn = function(prob_pair, para,  null_col){
#   plfdr_fn = function(z){
#     num = prob_pair[null_col]*stats::dnorm(z, mean = para[1, null_col], sd= para[2, null_col])
#     denom = num + prob_pair[-null_col]*stats::dnorm(z, mean = para[1, -null_col], sd= para[2, -null_col])
#     return(num/denom)
#
#   }
#   return(plfdr_fn)
#
# }
# # fn_ls = apply(prob, MARGIN = 1, FUN  = create_plfdr_fn, para = para, null_col = null_col)
# #
# flip = function(stat_vec , shell_fn_ls, norm_gen_size){
#   N = length(stat_vec)
#   if (N != length(shell_fn_ls )) stop("length of stat_vec must be same as shell_fn_ls")
#   norm.dat= stats::rnorm(n = norm_gen_size)
#   # use mclapply
#   dat_ls <- parallel::mclapply(shell_fn_ls, FUN = function(f) {f(norm.dat)}, mc.cores = 8)
#
#   obtain_flip = function(datvec, stat){
#     stats::quantile(datvec, 1 - sum(datvec < stat)/norm_gen_size )
#   }
#
#   stat_flip =
#     mapply(FUN = obtain_flip,
#          dat_ls,
#          stat_vec
#          )
#  return(stat_flip)
#
# }
#
#
# create_beta_lfdr_parts = function(z, param, extraParam,  X, gen_size){
#
#   create_beta_lfdr_fn= function(param, extraParam){
#
#     beta_lfdr_fn = function(z){
#       beta_lfdr =
#         (stats::dbeta(stats::pnorm(z), 1/(1 +exp(- param[1])), extraParam[1])*exp(param[3]) +
#            stats::dbeta(stats::pnorm(z), extraParam[2],  1/(1 + exp(-param[2])))*exp(param[4])   +
#            1)^(-1)
#       return(beta_lfdr)
#     }
#     return(beta_lfdr_fn)
#   }
#
#   linear_predictors = cbind(1, X)%*%matrix(param, ncol = 4)
#   beta_lfdr_fn_ls <- apply(X = linear_predictors,
#                            MARGIN  =1,
#                           FUN= create_beta_lfdr_fn,
#                            extraParam = extraParam)
#
#
#   beta_lfdr <- mapply(function(f, z){f(z)}, beta_lfdr_fn_ls, z)
#
#   # datavec= pnorm(z)
#   # datavec[which(datavec ==1)] = 0.9999999
#
#   beta_lfdr_flip <- flip(beta_lfdr, beta_lfdr_fn_ls, gen_size)
#
#   return(list(beta_lfdr = beta_lfdr, beta_lfdr_flip = beta_lfdr_flip))
#
# }



