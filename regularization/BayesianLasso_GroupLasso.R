# Data generation
n = 100
p = 100
beta <- c(2,2, -2, rep(0, p-3))
X <- matrix(rnorm(n*p), nrow = n)
eps <- rnorm(n , mean = 0, sd = sqrt(2))
y <- X%*%beta + eps 

library(mvtnorm)
library(invgamma)
library(statmod)

#+---------------Bayesian Group LASSO---------------------+
hat.beta <- solve(t(X)%*%X + diag(1000, p))%*%t(X)%*%y
hat.tau2 <-  rep(1,p)
hat.sig2 <- 1
hat.lambda2 <- 1


MCMC.size <- 5000
MCMC.beta <- matrix(0,MCMC.size, p)
for (j in 1: MCMC.size){
  
  par1 <- solve(t(X)%*%X + diag(hat.tau2^{-1}))
  mean_vec <- par1%*%t(X)%*%y 
  cov_mat <- hat.sig2*par1
  
  # beta
  hat.beta <-as.numeric(
    rmvnorm( 1, mean_vec, cov_mat)
  )
  MCMC.beta[j,] <- hat.beta 
  # sig2
  shape.par = (n+p)/2
  rate.par = 0.5*(sum((y-X%*%hat.beta)^2) +
                    t(hat.beta)%*%diag(hat.tau2^{-1})%*%hat.beta)
  hat.sig2 <- rinvgamma(1, shape = shape.par, rate = rate.par)
  
  # 1/tau2
  G <- p/5
  hat.tau2 <- numeric(p)
  for (g in 1:G){
    hat.beta.g <- hat.beta[(5*(g-1)+1) : (5*(g-1)+5)]
    ig.mean.par.g <- sqrt(hat.lambda2*hat.sig2/hat.beta.g^2)
    ig.dispersion.par.g <- rep(1/hat.lambda2, 5) # dispersion = 1/var
    hat.tau2.g <- 1/rinvgauss(1, ig.mean.par.g,dispersion = ig.dispersion.par.g) ## tau2 is a vector.
    hat.tau2[(5*(g-1)+1) : (5*(g-1)+5)] <- hat.tau2.g
  }
  
  # lambda2
  r = 1
  delta = 1.78
  gamma.par1 <- (p + G)/2 + r
  gamma.par2 <- sum(hat.tau2)/2/5 + delta # 5 for each group
  hat.lambda2 <- rgamma(1, gamma.par1,rate = gamma.par2)
  
  cat(j,",")
}
beta.bayes.group.lasso <- apply(MCMC.beta[2001:5000, ],2,mean)

#+------------------- Bayesian LASSO---------------------------+

# Initial values
## this beta initial value is not used in the loop
hat.beta <- solve(t(X)%*%X + diag(1000, p))%*%t(X)%*%y 
hat.tau2 <-  rep(1,p)
hat.sig2 <- 1
hat.lambda2 <- 1

MCMC.size <- 5000
MCMC.beta <- matrix(0,MCMC.size, p)
for (j in 1: MCMC.size){
  
  # library(mvtnorm)
  par1 <- solve(t(X)%*%X + diag(hat.tau2^{-1}))
  mean_vec <- par1%*%t(X)%*%y 
  cov_mat <- hat.sig2*par1
  
  # beta
  hat.beta <-as.numeric(
    rmvnorm( 1, mean_vec, cov_mat)
  )
  MCMC.beta[j,] <- hat.beta 
  # sig2
  # library(invgamma)
  shape.par = (n+p)/2
  rate.par = 0.5*(t(y-X%*%hat.beta)%*%(y-X%*%hat.beta) +
                    t(hat.beta)%*%diag(hat.tau2^{-1})%*%hat.beta)
  hat.sig2 <- rinvgamma(1, shape = shape.par, rate = rate.par)
  
  # 1/tau2
  # library(statmod)
  ig.mean.par <- sqrt(hat.lambda2*hat.sig2/hat.beta^2)
  ig.dispersion.par <- rep(1/hat.lambda2, p)
  hat.tau2 <- 1/rinvgauss(p, mean = ig.mean.par,
                          dispersion = ig.dispersion.par) ## tau2 is a vector.
  
  # lambda2
  r = 1
  delta = 1.78
  gamma.par1 <- p + r
  gamma.par2 <- 0.5*sum(hat.tau2)+ delta
  hat.lambda2 <- rgamma(1, gamma.par1,rate = gamma.par2)
  
  cat(j,",")
}
beta.bayes.lasso <- apply(MCMC.beta[2001:5000, ],2,mean)

# plot
plot(1:p, beta, pch = 15, cex = 0.4, col = 1, ylim = c(-3,3))
points(1:p, beta.bayes.group.lasso, pch = 15, cex = 0.4, col = 2)
points(1:p, beta.bayes.lasso, pch = 15,cex = 0.4, col = 4)
legend("topright", inset=.05, title="Method",
       c("True","Group LASSO","LASSO"), fill=c(1, 2, 4),cex = 0.7, horiz=TRUE)
