
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

This is a basic example which shows you how to solve a common problem,
taken from Setup 1 in Section 4.1 of the
[paper](https://academic.oup.com/jrsssb/article/84/5/1886/7072884), for
the simulation parameters
![\zeta = 0.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Czeta%20%3D%200.5 "\zeta = 0.5"),
![\epsilon = 1.9](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%20%3D%201.9 "\epsilon = 1.9")
and
![\eta = -2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ceta%20%3D%20-2 "\eta = -2"):

``` r
set.seed(123)
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
#>   [1]   27   49  145  157  161  169  226  243  257  264  273  277  306  310  322
#>  [16]  351  363  419  452  466  471  492  517  556  584  600  630  693  732  749
#>  [31]  842  851  917  975 1008 1011 1021 1023 1070 1119 1126 1131 1139 1170 1198
#>  [46] 1244 1354 1451 1483 1488 1494 1510 1524 1536 1544 1555 1585 1661 1672 1696
#>  [61] 1745 1749 1771 1793 1838 1877 1886 1887 1897 1971 1985 1993 2039 2059 2099
#>  [76] 2120 2146 2164 2211 2232 2233 2307 2318 2414 2445 2449 2486 2527 2581 2663
#>  [91] 2668 2739 2754 2770 2786 2832 2848 2883 2947 2973 3014 3082 3100 3111 3146
#> [106] 3161 3185 3189 3197 3217 3253 3287 3313 3356 3381 3385 3423 3464 3512 3520
#> [121] 3521 3524 3547 3575 3583 3609 3657 3773 3813 3841 3861 3879 3889 3947 4041
#> [136] 4059 4064 4090 4094 4150 4170 4179 4208 4210 4275 4277 4302 4348 4397 4423
#> [151] 4456 4544 4549 4567 4649 4655 4730 4748 4784 4797 4828
```

This package and its documentation is a work in progress. My apologies
for any inconvenience at the moment!
