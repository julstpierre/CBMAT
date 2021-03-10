#' @title CBMAT main function
#' @description FUNCTION_DESCRIPTION
#' @param y1 vector designates the first phenotype (one row-entry per individual), of length \eqn{n}.
#' @param fam1 a description of the error distribution and link function to be used in the marginal model of the first phenotype. Can be one of the following:
#' \itemize{
#' \item \code{"binomial(link=probit)"}
#' \item \code{"gaussian()"}
#' \item \code{"Gamma(link=log)"}
#' \item \code{"Student(link=)"}
#' }
#' @param y2 vector designates the second phenotype (one row-entry per individual), of length \eqn{n}.
#' @param fam2 a description of the error distribution and link function to be used in the marginal model of the second phenotype. Can be one of the following:
#' \itemize{
#' \item \code{"gaussian()"}
#' \item \code{"Gamma(link=log)"}
#' \item \code{"Student(link=)"}
#' }
#' @param x matrix of covariates including intercept (dimension:\eqn{n \times k}, with \eqn{k} the number of covariates), Default: cbind(rep(1, length(y1)))
#' @param G matrix of SNPs (dimension:\eqn{n \times p}, with \eqn{p} the number of SNPs)
#' @param copfit character, selects the copula(s) to use for modelling phenotypes dependence. 
#' Can be any one of the following:
#' \itemize{
#' \item \code{"Gaussian"}, Gaussian copula
#' \item \code{"Clayton"}, Clayton copula 
#' \item \code{"Gumbel"}, Gumbel copula 
#' \item \code{"Frank"}, Frank copula
#' }, Default: c("Gaussian", "Clayton", "Franck", "Gumbell")
#' @param weight logical, should weights be used to increase power for rare variants, Default: FALSE
#' @param weight.para1 alpha parameter of beta distribution used to simulate weights, Default: 1
#' @param weight.para2 beta parameter of beta distribution used to simulate weights, Default: 25
#' @param pval.method character, which method to use to calculate p-value.
#' Can be one of the following:
#' \itemize{
#' \item \code{"min"} (default)
#' \item \code{"Fischer"}, Fischer
#' \item \code{"MFKM"}, MFKM
#' }
#' @details When "weight=TRUE", a weighted test is used. The weight for each SNP is a beta fuction of the corresponding minor allele frequency (MAF). 
#' There are two parameters, weight.para1 and weight.para2, for beta function. For example, when weight.para1=weight.para2=0.5, the corresponding weight is 1/sqrt(MAF*(1-MAF)); 
#' when weight.para=1 and weight.para=25, the corresponding weight is the one suggested by SKAT.
#' @return A list containing results of the association test. The output list contains the following results:
#'   \itemize{
#'   \item \code{p.value}: the p-value of the region-based association test
#'   \item \code{gamma.y1}: the intercept and covariates estimates from the marginal model of the first phenotype
#'   \item \code{gamma.y2}: the intercept and covariates estimates from the marginal model of the second phenotype
#'   \item \code{alpha}: the estimated phenotypic dependence parameter, as a parameter of the used copula
#'   \item \code{tau}: Kendall's tau, which is a function of the estimated phenotypic dependence parameter alpha
#'   \item \code{cop}: based on the AIC criterion calculated based on the copula-model 
#'   }
#' @rdname CBMAT
#' @author Julien St-Pierre and Karim Oualkacha
#' @export 
CBMAT <- function(y1=NULL,
				  fam1=NULL,
				  y2=NULL,
				  fam2=NULL,
				  x=cbind(rep(1,length(y1))),
				  G=NULL,
				  copfit=c("Gaussian","Clayton","Franck","Gumbell"),
				  weight=FALSE,
				  weight.para1=1,
				  weight.para2=25,
				  pval.method="min"
)
{

	G<-as.matrix(G)
	p<-dim(G)[2]
	maf <- apply(G, 2, mean)/2
	if (!weight) {
	W <- diag(p)
	} else {
	W <- diag(dbeta(maf, weight.para1, weight.para2)^2, nrow = p, 
				ncol = p)
	}
	
	#y1
	if (fam1=="Student(link=)"){
		if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
			stop("Package \"LaplacesDemon\" needed for this function to work. Please install it.",
			call. = FALSE)
		}
		model1 <- estim(y1,x,fam="Student(link=)")
		fam1 <- list(family=sub("\\(.*","",model1$fam),link=gsub(".*=|\\).*","",model1$fam))
		g1 <- gaussian()$linkfun
		g1.inv <- gaussian()$linkinv
		v1 <- gaussian()$variance
		gamma.y1.init<-model1$coefficients
		phi1.init<-model1$phi
		log.phi1.init<-log(phi1.init)
		mu01.init<-g1.inv(c(x%*%gamma.y1.init))
		df1.init <- model1$df
		k<-length(model1$coefficients)
		
	} else {
		fam1 <- eval(parse(text=fam1))
		g1<-fam1$linkfun
		g1.inv<-fam1$linkinv
		v1<-fam1$variance
		reg1<-glm(y1~-1+x,family=fam1)
		gamma.y1.init<-reg1$coefficients
		phi1.init<-summary(reg1)$dispersion
		log.phi1.init<-log(phi1.init)
		mu01.init<-g1.inv(c(x%*%gamma.y1.init))
		df1.init <- NULL
		k<-length(reg1$coefficients)
	}
	
	#y2
	if (fam2=="Student(link=)"){
		if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
			stop("Package \"LaplacesDemon\" needed for this function to work. Please install it.",
			call. = FALSE)
		}
		model2 <- estim(y2,x,fam="Student(link=)")
		fam2 <- list(family=sub("\\(.*","",model2$fam),link=gsub(".*=|\\).*","",model2$fam))
		g2 <- gaussian()$linkfun
		g2.inv <- gaussian()$linkinv
		v2 <- gaussian()$variance
		gamma.y2.init<-model2$coefficients
		phi2.init<-model2$phi
		log.phi2.init<-log(phi2.init)
		mu02.init<-g2.inv(c(x%*%gamma.y2.init))
		df2.init <- model2$df
		k<-length(model2$coefficients)
		
	} else {
		fam2 <- eval(parse(text=fam2))
		g2<-fam2$linkfun
		g2.inv<-fam2$linkinv
		v2<-fam2$variance
		reg2<-glm(y2~-1+x,family=fam2)
		gamma.y2.init<-reg2$coefficients
		phi2.init<-summary(reg2)$dispersion
		log.phi2.init<-log(phi2.init)
		mu02.init<-g2.inv(c(x%*%gamma.y2.init))
		df2.init <- NULL
	}
	
	#Create cdf and pdf functions for y1 and y2
	if (fam1$family=="Gamma"){
		F1<-function(mu,phi,df) pgamma(y1,shape=1/phi,scale=mu*phi)
	} else if (fam1$family=="gaussian"){
		F1<-function(mu,phi,df) pnorm(y1,mu,sqrt(phi))
	} else if (fam1$link=="probit"){
		F1<-function(mu,phi,df) mu
	} else if (fam1$family=="Student"){
		F1<-function(mu,phi,df) LaplacesDemon::pst(y1,mu=mu,sigma=sqrt(phi^2*(df-2)/df),nu=df)
	}
  
	if (fam2$family=="Gamma"){
		F2<-function(mu,phi,df) pgamma(y2,shape=1/phi,scale=mu*phi)
		f2<-function(mu,phi,df) dgamma(y2,shape=1/phi,scale=mu*phi)
	} else if (fam2$family=="gaussian"){
		F2<-function(mu,phi,df) pnorm(y2,mu,sqrt(phi))
		f2<-function(mu,phi,df) dnorm(y2,mu,sqrt(phi))
	} else if (fam2$family=="Student"){
		F2<-function(mu,phi,df) LaplacesDemon::pst(y2,mu=mu,sigma=sqrt(phi^2*(df-2)/df),nu=df)
		f2<-function(mu,phi,df) LaplacesDemon::dst(y2,mu=mu,sigma=sqrt(phi^2*(df-2)/df),nu=df)
	}
	
	#Initial estimate for copula parameter
	cop.family <- sapply(1:length(copfit),function(i) switch (copfit[i],
        "Gaussian" = 1,
        "Clayton"  = 3,
        "Franck"   = 5,
        "Gumbell"  = 4)
	)
	cop<- VineCopula::BiCopSelect(F1(mu01.init,phi1.init,df1.init),F2(mu02.init,phi2.init,df2.init),familyset=cop.family,rotations=FALSE)$family
	
	if (fam1$link=="probit"){
		cor.est <- cor(y1,y2,method = "kendall")
	} else {
		cor.est <- cor(F1(mu01.init,phi1.init,df1.init),F2(mu02.init,phi2.init,df2.init),method = "kendall")
	}
	#--- cop param domain for L-BFGS-G method
	  if (cop==1) {
		lower.par.cop = -.98
		upper.par.cop = .98
	  }
	  
	  if (cop==3) {
		lower.par.cop = .001
		upper.par.cop = 10
	  }
	  
	  if ( ( cor.est > 0 ) & (cop == 5) ) {
		lower.par.cop = .001
		upper.par.cop = 10
	  }
	  
	  if ( ( cor.est < 0 ) & (cop == 5) ) {
		lower.par.cop = -10
		upper.par.cop = -.001
	  }
	
	alpha.init <- VineCopula::BiCopTau2Par(cop,max(cor.est,VineCopula::BiCopPar2Tau(cop,lower.par.cop)))
	
	#Estimate nuisance parameters under null hypothesis of no association
	if (fam1$link=="probit"){
		r=2*k+2+as.numeric(!is.null(df2.init)) # number of nuisance parameters
		pars.init<-c(alpha.init,gamma.y1.init,gamma.y2.init,log.phi2.init,df2.init);	#log.phi1.init is removed
		outoptim<-optim(pars.init, negloglik,cop=cop,y1=y1,y2=y2,x=x,
						fam1=fam1,fam2=fam2,g1.inv=g1.inv,g2.inv=g2.inv,
						method = "L-BFGS-B", 
						lower=c(lower.par.cop,rep(-Inf,2*k),-1,2.5),
						upper=c(upper.par.cop,rep(Inf,2*k+2))
						)
		df2 <- if(is.null(df2.init)){NULL} else{outoptim$par[2*k+3]}
		pars.est=rbind(c(outoptim$par[1],
						outoptim$par[2:(k+1)],
						outoptim$par[(k+2):(2*k+1)],
						exp(outoptim$par[2*k+2]),
						df2
						)
					)					
		
	} else {
		r=2*k+3+as.numeric(!is.null(df1.init))+as.numeric(!is.null(df2.init)) # number of nuisance parameters
		pars.init<-c(alpha.init,gamma.y1.init,gamma.y2.init,log.phi1.init,log.phi2.init,df1.init,df2.init);	
		outoptim<-optim(pars.init, negloglik,cop=cop,y1=y1,y2=y2,x=x,
						fam1=fam1,fam2=fam2,g1.inv=g1.inv,g2.inv=g2.inv,
						method = "L-BFGS-B", 
						lower=c(lower.par.cop,rep(-Inf,2*k+1),-0.99,2.5,2.5),
						upper=c(upper.par.cop,rep(Inf,2*k+4))
						)
		df1 <- if(is.null(df1.init)){NULL} else{outoptim$par[2*k+4]}
		df2 <- if(is.null(df2.init)){NULL} else{outoptim$par[r]}
		pars.est=rbind(c(outoptim$par[1],
						outoptim$par[2:(k+1)],
						outoptim$par[(k+2):(2*k+1)],
						exp(outoptim$par[2*k+2]),
						exp(outoptim$par[2*k+3]),
						df1,
						df2
						)
					)
	
	}
  
  #Calculate score vector and Hessian
  dlog <- -numDeriv::jacobian(negloglik.hessian,c(pars.est,rep(0,2*p)),cop=cop,y1=y1,y2=y2,cov=x,G=G,fam1=fam1,fam2=fam2,g1.inv=g1.inv,g2.inv=g2.inv)[(r+1):(r+2*p)]
  hess <- numDeriv::hessian(negloglik.hessian,c(pars.est,rep(0,2*p)),cop=cop,y1=y1,y2=y2,cov=x,G=G,fam1=fam1,fam2=fam2,g1.inv=g1.inv,g2.inv=g2.inv,method.args = list(d=0.001))
  UUT <- hess
  
  d2log <- -hess[(r+1):(r+2*p),(r+1):(r+2*p)]
  I.theta.theta <- hess[1:r,1:r]
  I.beta.theta <- hess[(r+1):(r+2*p),1:r]
  C <- cbind(-I.beta.theta%*%solve(I.theta.theta),diag(2*p))
  d2log.est <- C%*%UUT%*%t(C)
  sqrt.var.est <- Re(expm::sqrtm(d2log.est))
  Y.tilde=solve(sqrt.var.est, dlog, tol = 1e-20)
  
  f.rho <- function(rho){
	  U1 <- 1/2*t(dlog)%*%matrix(rbind(cbind(W, rho*W),(cbind(rho*W,W))),nrow=2*p)%*%dlog
	  U2 <- 1/2*sum(diag(d2log%*%matrix(rbind(cbind(W, rho*W),(cbind(rho*W,W))),nrow=2*p)))
	  U=U1+U2
	  matK <- sqrt.var.est%*%matrix(rbind(cbind(W, rho*W),(cbind(rho*W,W))),nrow=2*p)%*%sqrt.var.est
	  lambda<-eigen(Re(expm::sqrtm(-d2log))%*%matrix(rbind(cbind(W, rho*W),(cbind(rho*W,W))),nrow=2*p)%*%Re(expm::sqrtm(-d2log)),symmetric = TRUE)$values
	  lambda.est<-eigen(matK,symmetric = TRUE)$values
	  values <- lambda.est / sum(lambda.est)
	  pval.rho<-CompQuadForm::davies(2*U+sum(lambda),lambda.est)$Qq
	  K <- matK / sum(lambda.est)
	  v.Q <- 2*mult(K,K,dim(matK)[1])
	  Q <- t(Y.tilde) %*% K %*% Y.tilde
	  return(list("pval.rho"=pval.rho,"lambda"=lambda,"lambda.est"=lambda.est,"U"=U,"copula"=cop,"AIC"=AIC,"v.Q"=v.Q,"matK"=K,"Q"=Q,"values"=values))
	  }
  
  f.rho.out<-sapply(seq(0,.9,.1),function(rho) f.rho(rho))
  
  #----------- p.min  ---------#
  if (pval.method=="min"){
  	pval.rho <- unlist(f.rho.out[1,])
  	p.min <- min(pval.rho)
      
  	V.Qs <- unlist(f.rho.out[7,])
  	matK<-f.rho.out[8,]
  	d=length(V.Qs)
  	bb <- lapply(1:(length(V.Qs)-1), function(j)
  			lapply(j:(length(V.Qs)-1), function(i, x) 
  				2*mult(x[[j]], x[[i+1]], dim(x[[i]])[1]), x = matK
  			)
  		  )
  	Cov.Gam <- NULL
  	Cov.Gam <- diag(V.Qs)
    
  	for (j in 1:(length(V.Qs)-1))  {
  	Cov.Gam[j,(1+j):length(V.Qs)] = unlist(bb[[j]])}
  
  	Cov.Gam = Cov.Gam + t(Cov.Gam) - diag(V.Qs)
  	Corr.Gam = cov2cor(Cov.Gam)
  	#Corr.Gam[Corr.Gam>1]<-1 #make sure correlation is not > 1
  	cop.R = copula::normalCopula(copula::P2p(Corr.Gam), dim = d, dispstr = "un")
  
  	pval.min <- 1 - copula::pCopula(c(rep(1-p.min,d)), cop.R)
  	p.value <- pval.min
  }
	#------------- p.MFKM -----------#
  if (pval.method=="MFKM"){
    if (!requireNamespace(c("rARPACK","matrixcalc"), quietly = TRUE)) {
      stop("Package \"rARPACK\" and \"matrixcalc\" needed for this method to work. Please install it.",
           call. = FALSE)
    }
  	Qs <- unlist(f.rho.out[9,])
  	Q.sum = sum(Qs)
  	K.sum = Reduce("+", matK) 
  	values.sum = rARPACK::eigs(K.sum,matrixcalc::matrix.rank(K.sum,method="qr"))$values #JStP02MAY20:Added method="qr" to remove warning from chol method. 
  	p.MFKM = CompQuadForm::davies(Q.sum,values.sum)$Qq
  	p.value <- p.MFKM
  }
	#----------- p.Q.Fisher  ---------#
  if (pval.method=="Fischer"){
  	values <- t(matrix(unlist(f.rho.out[10,]),nrow=(2*p)))
  	Q.values <- cbind(Qs,values)
  	Q.Fisher <- sum(apply(Q.values,1,function(x) -2*log(davies.SURV(x[1],x[-1]))))
  	p.Q.Fisher <- perm.Q.Fisher(Q.Fisher, matK, values, nb.perm=1000)$p.Q.Fisher.perm
  	p.value <- p.Q.Fisher
  }
  
  #------------ Output -----------------------#
  return(list("p.value"=p.value,
			  "alpha"=pars.est[1],
			  "tau"=VineCopula::BiCopPar2Tau(cop,pars.est[1]),
			  "gamma.y1"=pars.est[2:(k+1)],
			  "gamma.y2"=pars.est[(k+2):(2*k+1)],
			  "cop" = ifelse(cop==1,"Gaussian", 
								  ifelse(cop==3,"Clayton", 
								  ifelse(cop==5,"Franck",  
								  ifelse(cop==4,"Gumbell"))))
			 )
		)
}
