## This file contains checks of input function parameters that are
## used by several other functions.

check_copula <- function(copfit) {
  if (! all(copfit %in% c("Gaussian", "Gumbel", "Frank", "Clayton"))) {
    stop("Unkown copula model specified; the copfit parameter ",
         "should be one of: 'Gaussian', 'Gumbel', 'Frank'",
         ", or 'Clayton'" )
  }
}

check_fam1 <- function(fam) {
  if (! fam %in% c("binomial(link=probit)", "gaussian()", "Gamma(link=log)", "Student(link=)")) {
    stop("Unkown model specified; the family parameter ",
         "should be one of: 'binomial(link=probit)', 'gaussian()', 'Gamma(link=log)'",
         ", or 'Student(link=)'" )
  }
}

check_fam2 <- function(fam) {
  if (! fam %in% c("gaussian()", "Gamma(link=log)", "Student(link=)")) {
    stop("Unkown model specified; the family parameter ",
         "should be one of: 'gaussian()', 'Gamma(link=log)'",
         ", or 'Student(link=)'" )
  }
}

check_pheno <- function(y) {
  if (is.null(y)) {
    stop("Argument y is empty, you did not specify a phenotype and/or genotype")
  }
}

check_covariates <- function(covariates, y, label) {
  if (!is.null(covariates)) {
    if (!class(covariates)[1] == "matrix") {
      stop("Argument covariates is not a matrix")
    }
  }
  
  if (nrow(covariates) != length(y)) {
    msg <- paste0("The number of rows of the ",label," matrix (",
                  nrow(covariates),
                  ") is not equal to the number of elements in the",
                  " phenotype (",
                  length(y),
                  ")")
    stop(msg)
  }
}
