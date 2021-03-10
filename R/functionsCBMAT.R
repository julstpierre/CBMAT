#----------------------------------------------#
# Calcul de la trace d'un produit de matrices
#----------------------------------------------#

mult.i=function(i,A,B) {(A[i,]%*%B[,i])[1,1]}
mult = function(A,B,n) {sum(sapply(1:n,mult.i,A=A,B=B))}

#----------------------------------------------#
## Calcul de la fonction de survie des formes quadratiques par 
## l'approximation de davies en forme matricielle pour les integrales 
#----------------------------------------------#
davies.SURV.2 = function(i,q,values) {
  qqq = CompQuadForm::davies(q[i],lambda=values)$Qq
  if (qqq > 1) qqq = 1
  if (qqq < 0) qqq =0 
  return(qqq)
}
davies.SURV = function(q,values) {sapply(1:(length(q)),davies.SURV.2,q=q,values=values)}

#----------------------------------------------------#
## perm.Q.Fisher
#----------------------------------------------------#
perm.Q.Fisher <- function(Q.Fisher, matK, values, nb.perm=1000){
  
  MCMC.Q.Fisher = NULL
  for (s.MCMC in 1:nb.perm){
    MCMC.Z = rnorm(dim(matK[[1]])[1],0,1)
    Q.perm.all <- unlist(lapply(matK, function(x) {
      Q.perm <- t(MCMC.Z) %*% x %*% MCMC.Z
      Q.perm
    }
    ))
    Q.perm.all.values = cbind(Q.perm.all,values)
    Q.Fisher.perm = sum(apply(Q.perm.all.values,1,function(x) -2*log(davies.SURV(x[1],x[-1]))))
    MCMC.Q.Fisher = c(MCMC.Q.Fisher, Q.Fisher.perm)
  }
  p.Q.Fisher <- mean( MCMC.Q.Fisher > as.numeric(Q.Fisher) )
  
  m1=mean(MCMC.Q.Fisher)
  m2=mean((MCMC.Q.Fisher-m1)^2)
  m4=mean((MCMC.Q.Fisher-m1)^4)
  g=(m4/(m2^2))-3
  df=12/g
  stat=df+(Q.Fisher-m1)*sqrt(2*df/m2)
  corrected.p.value = 1-pchisq(stat,df=df)
  results<-list(p.Q.Fisher=p.Q.Fisher, p.Q.Fisher.perm = corrected.p.value)
  results
}

#--------------------------------------------------------------------#
#This functions computes the negative log-likehood for estimation of 
#parameters under H0
#--------------------------------------------------------------------#
negloglik <- function(pars,cop,y1,y2,x,fam1,fam2,g1.inv,g2.inv){
  
  alpha<- pars[1];
  k<-dim(x)[2]
  r<-length(pars)
  gamma.y1<-pars[2:(k+1)]
  gamma.y2<-pars[(k+2):(2*k+1)]
  mu1=g1.inv(c(x%*%gamma.y1))
  mu2=g2.inv(c(x%*%gamma.y2))
  
  if (fam1$link=="probit"){
	log.phi2<-pars[2*k+2]
	phi2<-exp(log.phi2)
	if (fam2$family=="Student"){
		df2=pars[r]
		if (df2<=2){return(Inf)}
	} else {
		df2=NA
	}
	LL<-joint_prob_mixed(y1,y2,fam1,fam2,mu1,mu2,phi2,df2,alpha,cop)
  } else {
	log.phi1<-pars[2*k+2]
	log.phi2<-pars[2*k+3]
	phi1<-exp(log.phi1)
	phi2<-exp(log.phi2)
	if (fam1$family=="Student"){
		df1=pars[2*k+4]
		if (df1<=2){return(Inf)}
	} else {
		df1=NA
	}
	if (fam2$family=="Student"){
		df2=pars[r]
		if (df2<=2){return(Inf)}
	} else {
		df2=NA
	}
	LL<-joint_prob(y1,y2,fam1,fam2,mu1,mu2,phi1,phi2,df1,df2,alpha,cop)
  }
  
  return(-sum(log(LL)))
}

#--------------------------------------------------------------------#
#This functions computes the negative log-likehood for calculation of 
#the score's statistics
#--------------------------------------------------------------------#
negloglik.hessian<- function(pars,cop,y1,y2,cov,G,fam1,fam2,g1.inv,g2.inv){
  
  k<-dim(cov)[2]
  p<-dim(G)[2]
  r<-length(pars)-2*p
  alpha<-pars[1]
  gamma.y1<-pars[2:(k+1)]
  gamma.y2<-pars[(k+2):(2*k+1)]
  
  if (fam1$link=="probit"){
	phi2<-pars[2*k+2]
	if (fam2$family=="Student"){df2=pars[r]} else{df2=NA}
	beta<-pars[(r+1):(r+2*p)]
	LL<-joint_prob_test_mixed(y1,y2,fam1,fam2,cov,G,gamma.y1,gamma.y2,beta,phi2,df2,alpha,cop,g1.inv,g2.inv)
  } else {
	phi1<-pars[2*k+2]
	phi2<-pars[2*k+3]
	if (fam1$family=="Student"){df1=pars[2*k+4]} else{df1=NA}
	if (fam2$family=="Student"){df2=pars[r]} else{df2=NA}
	beta<-pars[(r+1):(r+2*p)]
	LL<-joint_prob_test(y1,y2,fam1,fam2,cov,G,gamma.y1,gamma.y2,beta,phi1,phi2,df1,df2,alpha,cop,g1.inv,g2.inv)
  }

	return(-sum(log(LL)))   
}


#-----------------------------------------------------------------------#
# joint_prob is the model LL used to estimate parameters under H0
# when y1 and y2 are continuous
#-----------------------------------------------------------------------#
joint_prob<-function(y1,y2,fam1,fam2,mu1,mu2,phi1,phi2,df1,df2,alpha,cop){
  
  if (fam1$family=="Gamma"){
    F.y1<-pgamma(y1,shape=1/phi1,scale=mu1*phi1)
    f.y1<-dgamma(y1,shape=1/phi1,scale=mu1*phi1)
  } else if (fam1$family=="gaussian"){
    F.y1<-pnorm(y1,mu1,sqrt(phi1))
    f.y1<-dnorm(y1,mu1,sqrt(phi1))
  } else if (fam1$family=="Student"){
	F.y1<-LaplacesDemon::pst(y1,mu1,sigma=sqrt(phi1^2*(df1-2)/df1),nu=df1)
	f.y1<-LaplacesDemon::dst(y1,mu1,sigma=sqrt(phi1^2*(df1-2)/df1),nu=df1)
  }
  
  if (fam2$family=="Gamma"){
    F.y2<-pgamma(y2,shape=1/phi2,scale=mu2*phi2)
    f.y2<-dgamma(y2,shape=1/phi2,scale=mu2*phi2)
  } else if (fam2$family=="gaussian"){
    F.y2<-pnorm(y2,mu2,sqrt(phi2))
    f.y2<-dnorm(y2,mu2,sqrt(phi2))
  } else if (fam2$family=="Student"){
	F.y2<-LaplacesDemon::pst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
	f.y2<-LaplacesDemon::dst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
  }
  
  LL<-f.y1*f.y2*VineCopula::BiCopPDF(F.y1,F.y2,cop,alpha)
  return(LL)
}

#-----------------------------------------------------------------------#
# joint_prob_mixed is the model LL used to estimate parameters under H0
# when y1 is discrete and y2 is continuous
#-----------------------------------------------------------------------#
joint_prob_mixed<-function(y1,y2,fam1,fam2,mu1,mu2,phi2,df2,alpha,cop){
  
  if (fam1$link=="probit"){
    F.y1<-1-mu1
  }
  
  if (fam2$family=="Gamma"){
    F.y2<-pgamma(y2,shape=1/phi2,scale=mu2*phi2)
    f.y2<-dgamma(y2,shape=1/phi2,scale=mu2*phi2)
  } else if (fam2$family=="gaussian"){
    F.y2<-pnorm(y2,mu2,sqrt(phi2))
    f.y2<-dnorm(y2,mu2,sqrt(phi2))
  } else if (fam2$family=="Student"){
	F.y2<-LaplacesDemon::pst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
	f.y2<-LaplacesDemon::dst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
  }
  
  div.C.y2<-VineCopula::BiCopHfunc2(F.y1,F.y2,cop,alpha)
  LL<-f.y2*((y1==0)* div.C.y2 + (y1==1)*(1-div.C.y2))
  return(LL)
}

#-----------------------------------------------------------------------#
# joint_prob test is the model LL used to calculate the score's statistic
# when y1 and y2 are continuous
#-----------------------------------------------------------------------#
joint_prob_test<-function(y1,y2,fam1,fam2,cov,G,gamma.y1,gamma.y2,beta,phi1,phi2,df1,df2,alpha,cop,g1.inv,g2.inv){
  
  p<-dim(G)[2]
  mu1<-g1.inv(c(cov%*%gamma.y1)+c(G%*%beta[1:p]))
  mu2<-g2.inv(c(cov%*%gamma.y2)+c(G%*%beta[(p+1):(2*p)]))
  
  
  if (fam1$family=="Gamma"){
    F.y1<-pgamma(y1,shape=1/phi1,scale=mu1*phi1)
    f.y1<-dgamma(y1,shape=1/phi1,scale=mu1*phi1)
  } else if (fam1$family=="gaussian"){
    F.y1<-pnorm(y1,mu1,sqrt(phi1))
    f.y1<-dnorm(y1,mu1,sqrt(phi1))
  } else if (fam1$family=="Student"){
	F.y1<-LaplacesDemon::pst(y1,mu1,sigma=sqrt(phi1^2*(df1-2)/df1),nu=df1)
	f.y1<-LaplacesDemon::dst(y1,mu1,sigma=sqrt(phi1^2*(df1-2)/df1),nu=df1)
  }
  
  if (fam2$family=="Gamma"){
    F.y2<-pgamma(y2,shape=1/phi2,scale=mu2*phi2)
    f.y2<-dgamma(y2,shape=1/phi2,scale=mu2*phi2)
  } else if (fam2$family=="gaussian"){
    F.y2<-pnorm(y2,mu2,sqrt(phi2))
    f.y2<-dnorm(y2,mu2,sqrt(phi2))
  } else if (fam2$family=="Student"){
	F.y2<-LaplacesDemon::pst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
	f.y2<-LaplacesDemon::dst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
  }
  
  LL<-f.y1*f.y2*VineCopula::BiCopPDF(F.y1,F.y2,cop,alpha)
  return(LL)
}

#-----------------------------------------------------------------------#
# joint_prob test_mixed is the LL used to calculate the score's statistic
# when y1 is discrete and y2 continuous
#-----------------------------------------------------------------------#
joint_prob_test_mixed<-function(y1,y2,fam1,fam2,cov,G,gamma.y1,gamma.y2,beta,phi2,df2,alpha,cop,g1.inv,g2.inv){
  
  p<-dim(G)[2]
  mu1<-g1.inv(c(cov%*%gamma.y1)+c(G%*%beta[1:p]))
  mu2<-g2.inv(c(cov%*%gamma.y2)+c(G%*%beta[(p+1):(2*p)]))
  
  
  if (fam1$link=="probit"){
    F.y1<-1-mu1
  }
  
  if (fam2$family=="Gamma"){
    F.y2<-pgamma(y2,shape=1/phi2,scale=mu2*phi2)
    f.y2<-dgamma(y2,shape=1/phi2,scale=mu2*phi2)
  } else if (fam2$family=="gaussian"){
    F.y2<-pnorm(y2,mu2,sqrt(phi2))
    f.y2<-dnorm(y2,mu2,sqrt(phi2))
  } else if (fam2$family=="Student"){
	F.y2<-LaplacesDemon::pst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
	f.y2<-LaplacesDemon::dst(y2,mu2,sigma=sqrt(phi2^2*(df2-2)/df2),nu=df2)
  }
  
  div.C.y2<-VineCopula::BiCopHfunc2(F.y1,F.y2,cop,alpha,check.pars=FALSE) #JStP 11NOV2019 : added check.pars=FALSE to allow calculation of hessian when copula parameter is close to treshold
  LL<-f.y2*((y1==0)* div.C.y2 + (y1==1)*(1-div.C.y2))
  return(LL)
}

###########MLE for fitting gamma, student and normal distributions################
LL.student<-function(y,cov,pars){
  nu<-pars[1]
  gamma<-pars[2:length(pars)]
  mu=c(cov%*%gamma)
  sigma=sqrt((nu-2)/nu)
  LL<-LaplacesDemon::dst(y, mu=mu, nu=nu,sigma=sigma)
  return(-sum(log(LL)))
  
}

estim<-function(y,x,fam){
  
  if (fam=="gaussian()") {
    reg <- glm(y~-1+x)
    coefficients<-reg$coefficients
    phi<-summary(reg)$dispersion
    mu<-c(x%*%coefficients)
    df <- NA
    aic = AIC(reg)
  } else if (fam=="Gamma(link=log)" & all(y>0)) {
    reg<-glm(y~-1+x,family=Gamma(link=log))
    coefficients<-reg$coefficients
    phi<-summary(reg)$dispersion
    mu<-Gamma(link=log)$linkinv(c(x%*%coefficients))
    df <- NA
    aic = AIC(reg)
  } else if (fam=="Student(link=)") {
    reg<-glm(y~-1+x)
    gamma.init<-reg$coefficients
    optim.st<-try(optim(c(10,gamma.init),LL.student,y=y,cov=x,method = "L-BFGS-B",lower=c(2.5,rep(-Inf,length(gamma.init)))),silent=TRUE)
	if (class(optim.st)=="try-error"){
		log.LL <- NA
		aic = NA
		df <- NA
		coefficients<- NA
		phi <- NA
	} else{
		log.LL <- -optim.st$value
		aic = -2*log.LL+2*(length(optim.st$par)+1)
		df <-optim.st$par[1]
		coefficients<-optim.st$par[-1]
		phi <- sd(y)
	}
  }
  else {
    aic=NA
    coefficients=NA
    df=NA
    phi=NA
  }
  
  return (list(family=fam,aic=aic,coefficients=coefficients,df=df,phi=phi))
  
}

best.model <- function(y,x){
  
  model.norm <- estim(y,x,fam="gaussian()")
  model.gamma <- estim(y,x,fam="Gamma(link=log)")
  model.st <- estim(y,x,fam="Student(link=)")
  min.aic <- min(model.norm$aic
                ,model.gamma$aic
                ,model.st$aic
                ,na.rm=T)
  
  if (model.norm$aic==min.aic){
    return (list(family=model.norm$fam,aic=model.norm$aic,coefficients=model.norm$coefficients,phi=model.norm$phi))
  } 
  
  else if (!is.na(model.st$aic) & model.st$aic==min.aic){
    return (list(family=model.st$fam,aic=model.st$aic,coefficients=model.st$coefficients,phi=model.st$phi,df=model.st$df))
  } 
  
  else if (!is.na(model.gamma$aic) & model.gamma$aic==min.aic){
    return (list(family=model.gamma$fam,aic=model.gamma$aic,coefficients=model.gamma$coefficients,phi=model.gamma$phi))
  } 
  
}

