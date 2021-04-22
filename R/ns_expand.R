# x should be a data frame

ns_expand = function(x , df){
  X_ls <- list()
  for (i in 1:df){
    if (i == 1){
      X_ls[[i]] <- as.matrix(x)
      colnames(X_ls[[i]]) = NULL
      rownames(X_ls[[i]]) = NULL
    }else{
      X_ls[[i]]  <-do.call(what = cbind,
                           args  = lapply(X = as.list(x),  FUN = ns, df = i))
      colnames(X_ls[[i]]) = NULL
      rownames(X_ls[[i]]) = NULL
    }
  }
  return(X_ls)
}
