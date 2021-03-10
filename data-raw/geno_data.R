######################################
# R Source code file for creating simulated dataset to be included in the CBMAT package
# Author: Julien St-Pierre
# Created: March 04, 2021
# Updated: 
# Notes:
#####################################

if (!require("pacman")) install.packages("pacman") 
pacman::p_load(usethis)
pacman::p_load(mvtnorm)
pacman::p_load(VineCopula)

set.seed(1245)
#----------------------
# Initialize parameters
#----------------------
p=30              # number of snps on the gene
v=0.1             # fraction of snps that are non zero
tau=0.3           # variance component parameter
rho=0             # covariance parameter for beta
cop=1             # gaussian copula to model joint dependence
Kendall.tau = 0.2 # Kendall tau for joint dependence
phi1=1            # Dispersion parameter for y1
phi2=1            # Dispersion parameter for y2

#--------------------------
# Simulate genotype matrix
#--------------------------
G.all <- read.table("data-raw/G.012")
s <- sample(1:(dim(G.all)[2]-p+1),1)
G <- as.matrix(G.all[,s:(s+p-1)]);colnames(G) <- paste("V",1:p,sep="")

#-----------------------------------------------------
# Simulate coefficients for the genetic region effect
#-----------------------------------------------------
W <- diag(p)
diag(W)[sample(1:p,p*(1-v))] <- 0
beta <- t(mvtnorm::rmvnorm(1,sigma=tau*matrix(rbind(cbind(W, rho*W),(cbind(rho*W,W))),nrow=2*p)))

#------------------------------------------------------
# Simulate non-genetic covariates (including intercept)
#------------------------------------------------------
k=3 							
gamma.y1 <- c(runif(k,0,2))
gamma.y2 <- c(runif(k,0,2))
x <- cbind(1,rbinom(nrow(G),1,0.5),rnorm(nrow(G),0,1))

#-------------------------------------------------
# Simulate univariate cdfs using a gaussian copula
#-------------------------------------------------
alpha <- VineCopula::BiCopTau2Par(cop,Kendall.tau)
u <- VineCopula::BiCopSim(nrow(G),cop,alpha)

#-------------------------------------------------------
#Simulate discrete trait using gaussian latent variable
#-------------------------------------------------------
mu1.c <- c(x%*%gamma.y1+G%*%beta[1:p])
y1.c <- qnorm(u[,1],mu1.c,sqrt(phi1))
y1 <- as.numeric(y1.c>quantile(y1.c,0.70))

#-------------------------------------------------------
#Simulate a continuous Gamma trait
#-------------------------------------------------------
g2.inv <- Gamma(link="log")$linkinv
mu2 <- g2.inv(c(x%*%gamma.y2+G%*%beta[(p+1):(2*p)]))
y2 <- qgamma(u[,2],shape=1/phi2,scale=mu2*phi2)


#write data
data_mixed <- list("y.bin"=y1,"y.gauss"=y1.c,"y.Gamma"=y2,"x"=x,"G"=G)
usethis::use_data(data_mixed, overwrite = TRUE)
