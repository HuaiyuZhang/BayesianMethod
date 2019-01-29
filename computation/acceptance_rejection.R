# A simulation for Beta(2,2) using acceptance-rejection algorithm

c = 1.5 # optimal
#c = 0.5
f = function(x) 6*x*(1-x) # dbeta(2,2)
B = 1000
y = NULL
i = 0 # number of iteration
k = 0 # number of acceptance
while (k < B){
  Z = runif(1)
  u = runif(1)
  if (u<= f(Z)/c*1){
    y_new = Z
    y = c(y, y_new)
    k = k+1 # update the number of acceptance
  }  
  i = i + 1
}
length(y)
k/i
# Draw Q-Q plot
Empirical.Q <- quantile(y, seq(0.1, 0.9, 0.05))
Q <- qbeta(seq(0.1,0.9,0.05), 2,2)
plot(Q, Empirical.Q)
abline(0,1)

