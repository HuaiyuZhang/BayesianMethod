# Example for MH algorithm: generate random numbers from Rayleigh distribution
set.seed(1234)
m <- 10000
sigma2 <- 4
# log of the target distribution
logf.ray <- function(x){
      log(x)-log(sigma2)-(x^2)/(2*sigma2)
}

x <- numeric(m)  # collector
x[1] <- rchisq(1, df = 1) # initial value
for(t in 2:m){
      x.t_1 <- x[t-1]
      x.star <- rchisq(1, df = x.t_1)
      log.num <- logf.ray(x.star) + dchisq(x.t_1, df = x.star, log = T)
      log.den <- logf.ray(x.t_1) + dchisq(x.star, df = x.t_1, log = T)
      log.R <- log.num - log.den
      if(log(runif(1)) <= log.R){
            x[t] <- x.star
      }
      else{
           x[t] <- x.t_1 
      }
}
plot(x, type = 'l')
