#' Simulated mixed bivariate phenotype
#'
#' @description To simulate genotypes, we used the 1000 Genomes Project database.
#' Variants within 500kbs of the BRCA1 gene, for which several known mutations
#' are associated with a higher risk of developing breast, ovarian and prostate
#' cancers, were selected. To avoid any multicollinearity, we pruned variants
#' based on linkage disequilibrium r^2 > 0.7. Further, a total of 503 subjects
#' with a European genetic ancestry are selected in order to avoid any population
#' structure. One discrete and one continuous traits are simulated using a gaussian copula
#' to model the joint dependence.
#' 
#' @format This data frame has 503 rows and the following 35 columns:
#' \describe{
#'   \item{y1}{discrete trait simulated from a latent gaussian variable}
#'   \item{y2}{continous trait simulated from a Gamma distribution}
#'   \item{x1}{intercept}
#'   \item{x2}{discrete covariate}
#'   \item{x3}{continous covariate}
#'   \item{V1:V30}{genetic variants}
#' }
#' @source \url{https://www.internationalgenome.org/data/}
"df_mixed"