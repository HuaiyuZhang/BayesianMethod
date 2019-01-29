library(LearnBayes)
set.seed(1234)
data(donner)
y <- donner$survival
n <- length(y)
Male <- donner$male
Age <- donner$age
X <- cbind(rep(1, n), Male, Age)
logit <- function(x, log = F){
      if(!log){
            1/(1+1/exp(x))
      }
      else{
            log(1/(1+1/exp(x)))
      }
}
# Laplace
neg.log.post <- function(b){
      Xb <- as.numeric(X%*%b)
      (-1)*sum(y*plogis(Xb, 0, 1, log = T) + (1-y)*plogis(-Xb, 0, 1, log = T))-
            sum(dnorm(b, 0, sqrt(10^5), log = T))
}
Asymp.post <- optim(c(0,0,0), neg.log.post, hessian = T)
Laplace.Post.Mean <- Asymp.post$par
Laplace.Post.SD <- sqrt(diag(solve(Asymp.post$hessian)))

# 