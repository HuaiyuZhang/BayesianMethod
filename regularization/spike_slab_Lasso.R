# Exam question 3
# spike and slab and LASSO
library(mvtnorm)
library(invgamma)
library(statmod)
library(rmutil)


# Data generation
# n = 100
# p = 200
# beta <- c(2,2, -2, rep(0, p-3))
# X <- matrix(rnorm(n*p), nrow = n)
# eps <- rnorm(n , mean = 0, sd = sqrt(2))
# y <- X%*%beta + eps 
# plot(1:p, beta, pch = 20, col = 1, ylim = c(-3,3))

# data
library(ISLR)
data(Hitters)
Data <- na.omit(Hitters)
y <- Data[,19]
X <- Data[,-19]
X$League <- ifelse(X$League == 'A', 1, 0)
X$Division <- ifelse(X$Division =='E', 1, 0)
X$NewLeague <- ifelse(X$NewLeague=='A', 1, 0)
X <- as.matrix(scale(X)[,])
y <- (y-mean(y))/sd(y)

n <- nrow(X)
p <- ncol(X)





# Initial values
lambda1 <- 10
lambda0 <- 0.001
hat.gamma <- rep(0,p)
# hat.beta <- solve(t(X)%*%X + diag(1/lambda0, p))%*%t(X)%*%y
hat.sig2 <- 1
hat.w <- 0.5 # could also use a layer.
hat.tau2 <-  rep(1,p) # lasso part

MC.size <- 5000
MCMC.beta <- matrix(0, MC.size, p)
MCMC.gamma <- matrix(0, MC.size, p)
for (k in 1:MC.size){
  #update beta
  Lambda.gamma <- hat.gamma*hat.tau2 + (1-hat.gamma)*lambda0
  tXX <- solve(t(X)%*%X + diag(1/Lambda.gamma))
  mean.beta <- tXX%*%t(X)%*%y
  var.beta <- hat.sig2*tXX
  hat.beta <- as.numeric(rmvnorm(1,mean.beta, var.beta))
  
  # 1/tau2
  for (j in 1:p){
    if (hat.gamma[j] == 1){
      ig.mean.par <- sqrt(hat.sig2/hat.beta[j]^2/lambda1)
      ig.dispersion.par <- lambda1
      hat.tau2[j] <- 1/statmod::rinvgauss(1, ig.mean.par,dispersion = ig.dispersion.par) 
    }
    else{
      hat.tau2[j] <- rexp(1, rate = 1/(2*lambda1))
    }
  }
   
  # update sig2 
  a.sig2 <- (n+p)/2
  b.sig2 <- sum((y-X%*%hat.beta)^2)/2 + sum(hat.beta^2/Lambda.gamma)/2 
  hat.sig2 <- rinvgamma(1, shape = a.sig2, rate = b.sig2)
  
  
  # update gamma
  for (j in 1:p){
    prob.1 <- log(hat.w) + 
      dnorm(hat.beta[j], 0, sqrt(hat.sig2*hat.tau2[j]), log = T)
    prob.0 <- log(1-hat.w) + 
      dnorm(hat.beta[j], 0, sqrt(hat.sig2*lambda0), log = T)
    prob.10 <- 1/(1+exp(prob.0-prob.1))
    
    hat.gamma[j] <- rbinom(1,1, prob.10) # Bernoulli sample
  }
 
  MCMC.beta[k, ] <- hat.beta
  MCMC.gamma[k, ] <- hat.gamma
  cat(k, '~')
}


hat.beta <- apply(MCMC.beta[-(1:2000),], 2, mean)
res <- data.frame(name = colnames(X), hat.beta = hat.beta)

#  variable selection based on gamma.
post.indicator <- apply(MCMC.gamma[-(1:2000),], 2, 
                        function(x) mean(x)>0.5)
res[post.indicator,]
