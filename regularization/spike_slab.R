# Ex  4 spike and slab

# Data generation
n = 100
p = 200
beta <- c(2,2, -2, rep(0, p-3))
X <- matrix(rnorm(n*p), nrow = n)
eps <- rnorm(n , mean = 0, sd = sqrt(2))
y <- X%*%beta + eps 


plot(1:p, beta, pch = 20, col = 1, ylim = c(-3,3))
# spike-and-slab no bias
library(mvtnorm)
library(invgamma)


lambda1 <- 1000
lambda0 <- 0.001
hat.gamma <- rep(0,p)
hat.beta <- solve(t(X)%*%X + diag(1/lambda0, p))%*%t(X)%*%y
hat.sig2 <- 1
hat.w <- 0.5 # could also use a layer.

MC.size <- 1000
MCMC.beta <- matrix(0, MC.size, p)

for (k in 1:MC.size){
#update beta
Lambda.gamma <- hat.gamma*lambda1 + (1-hat.gamma)*lambda0
tXX <- solve(t(X)%*%X + diag(1/Lambda.gamma))
mean.beta <- tXX%*%t(X)%*%y
var.beta <- hat.sig2*tXX
hat.beta <- as.numeric(rmvnorm(1,mean.beta, var.beta))

# update sig2
a.sig2 <- (n+p+1)/2
b.sig2 <- sum((y-X%*%hat.beta)^2)/2 + sum(hat.beta^2/Lambda.gamma)/2+1/2
hat.sig2 <- rinvgamma(1, shape = a.sig2, rate = b.sig2)

# update gamma
prob.1 <- log(hat.w) + 
  dnorm(hat.beta, 0, sqrt(hat.sig2*lambda1), log = T) #log w + log P1
prob.0 <- log(1-hat.w) + 
  dnorm(hat.beta, 0, sqrt(hat.sig2*lambda0), log = T)
prob.10 <- 1/(1+exp(prob.0-prob.1))
hat.gamma <- rbinom(p,1, prob.10) # Bernoulli sample

#plot(hat.gamma, ylim = c(0,1), main = k)

MCMC.beta[k, ] <- hat.beta 
}
# called stochastic search var selection
hat.beta.ssvs <- apply(MCMC.beta[-(1:200),], 2, mean)

plot(1:p, beta, pch = 20, col = 2)
points(1:p, hat.beta, pch = 20, col = 3)
