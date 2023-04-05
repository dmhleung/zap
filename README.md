
<!-- README.md is generated from README.Rmd. Please edit that file -->

# zap

<!-- badges: start -->
<!-- badges: end -->

The goal of zap is to perform z-value based covariate-adaptive testing,
as described in [our zap
paper](https://academic.oup.com/jrsssb/article/84/5/1886/7072884).

## Installation

You can install the development version of zap from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dmhleung/zap")
```

## Example

This is a basic example which shows you how to use `zap_asymp()`, taken
from Setup 1 in Section 4.1 of the
[paper](https://academic.oup.com/jrsssb/article/84/5/1886/7072884), for
the simulation parameters
![\zeta = 0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Czeta%20%3D%200.5 "\zeta = 0.5"),
![\epsilon = 1.9](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%20%3D%201.9 "\epsilon = 1.9")
and
![\eta = -2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ceta%20%3D%20-2 "\eta = -2"):

``` r
set.seed(123)
RNGkind("L'Ecuyer-CMRG") # make sure we can reproduce the result even using parallel computing
library(zap)
## basic example code
dens = -2
inform = 0.5
esize = 1.9
symm = 1 # symm  = 1  # if setup is S1
altvar = 1
m = 5000
P=2
ep = 1
alpha = 0.05
X = matrix(rnorm(m*P, 0, sqrt(1/2)), m,P)
direction = sample(x = c(-1, 1),
                     size = m,
                     replace = TRUE,
                     prob = c(1 - symm,  symm)
)
coeff_reg = c(dens,rep(inform, P))  # change how asymmetrically X influence the probabilities
psi = as.numeric(crossprod(t(cbind(1, X)), coeff_reg ))
prior_mean = 2*esize*direction/(1 + exp(- inform*rowSums(X)))
wsuccess = 1/(1+exp(-psi))
gammatrue = rbinom(m,1,wsuccess)
thetatrue = prior_mean
thetatrue[!gammatrue] = 0
nbr_true_sig = sum(gammatrue )

sd = rep(1, m)
sd[gammatrue==1] = sqrt(altvar)
z = rnorm(m, mean = thetatrue, sd = sd)
pvals = 2*pnorm(abs(z), lower.tail  = FALSE)
sgn = sign(z) 

# apply the zap asymp function
zap_asymp_obj <- zap_asymp(z = z, X = X)
names(zap_asymp_obj)
#> [1] "rej_index"   "stat"        "mirror_stat" "param"       "FDR"
zap_asymp_obj$rej_index
#>   [1]   87  118  141  219  245  267  334  365  378  379  381  384  392  434  519
#>  [16]  540  602  624  629  631  655  671  675  680  713  877  902  947  969  985
#>  [31] 1038 1054 1091 1192 1202 1318 1328 1398 1454 1459 1493 1649 1672 1685 1729
#>  [46] 1770 1779 1817 1821 1846 1847 1911 1917 1941 1968 1992 2004 2029 2055 2103
#>  [61] 2151 2165 2185 2210 2280 2295 2375 2404 2427 2456 2462 2505 2572 2589 2692
#>  [76] 2695 2738 2809 2813 2870 2884 2885 2919 2952 2980 3005 3006 3020 3023 3059
#>  [91] 3083 3148 3166 3209 3218 3239 3246 3294 3352 3356 3373 3375 3440 3494 3501
#> [106] 3516 3535 3553 3625 3649 3655 3667 3686 3698 3768 3773 3868 3932 4043 4126
#> [121] 4292 4301 4319 4349 4354 4379 4419 4482 4496 4514 4527 4620 4626 4652 4655
#> [136] 4661 4686 4746 4765 4785 4823 4828 4847 4865 4882 4937 4946 4998
```

This package and its documentation is a work in progress. My apologies
for any inconvenience at the moment!
