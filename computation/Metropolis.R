set.seed(123)
y <- c(0,0)
rho <- 0.8
# Metropolis Algorithm 
library(mvtnorm)
M <- 3000
theta <- matrix(0, M, 2)

theta[1,] <- c(2.5, -2.5)
V <- 0.1 * matrix(c(1, rho,rho, 1), 2, 2) # Cov for Proposal
Sig <- matrix(c(1, rho,rho, 1), 2, 2) # Cov matrix for the posterior

for(i in 2:M){
      theta.star <- rmvnorm(1, theta[(i-1),], V)
      log.R <- dmvnorm(theta.star, y, Sig, log=T) - dmvnorm(theta[(i-1),], y, Sig, log=T)
      if(log(runif(1))<=log.R){
            theta[i, ] <- theta.star
      }
      else{
            theta[i, ] <- theta[(i-1), ]
      }
}
head(theta)
plot(theta[,1], theta[,2], type = 'l', ylim=c(-4, 4), xlim=c(-4,4))
points(theta[,1], theta[,2], pch=15)
