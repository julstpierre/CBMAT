
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CBMAT: Copula Based Multivariate Association Test for bivariate mixed phenotypes

<!-- badges: start -->

<!-- badges: end -->

## Installation

``` r
# development version from GitHub
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("julstpierre/CBMAT")
```

## Example

``` r
library(CBMAT)
attach(data_mixed)
mixed.fit <-CBMAT(y1=y.bin,
                  fam1="binomial(link=probit)",
                  y2=y.Gamma,
                  fam2="Gamma(link=log)",
                  x=x,
                  G=G,
                  copfit=c("Gaussian","Clayton","Franck","Gumbell"),
                  weight=FALSE,
                  weight.para1=1,
                  weight.para2=25,
                  pval.method="min"
)

mixed.fit
#> $p.value
#> [1] 0.4540164
#> 
#> $alpha
#> [1] 0.3429957
#> 
#> $tau
#> [1] 0.1463919
#> 
#> $gamma.y1
#> [1] -1.991215  1.481467  1.881062
#> 
#> $gamma.y2
#> [1] 0.2565104 1.0136100 0.7674400
#> 
#> $cop
#> [1] "Clayton"

cont.fit <-CBMAT(y1=y.gauss,
                 fam1="gaussian()",
                 y2=y.Gamma,
                 fam2="Gamma(link=log)",
                 x=x,
                 G=G,
                 copfit=c("Gaussian","Clayton","Franck","Gumbell"),
                 weight=FALSE,
                 weight.para1=1,
                 weight.para2=25,
                 pval.method="min"
)
cont.fit
#> $p.value
#> [1] 0.5554563
#> 
#> $alpha
#> [1] 0.2672714
#> 
#> $tau
#> [1] 0.172244
#> 
#> $gamma.y1
#> [1] 1.773083 1.600274 1.815735
#> 
#> $gamma.y2
#> [1] 0.2501279 1.0248267 0.7646350
#> 
#> $cop
#> [1] "Gaussian"
```
