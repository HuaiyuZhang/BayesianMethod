# Bayesian SVM, kernel version.
# Created: "Tue Dec 11 09:23:31 2018"

#+-----------Kernel function---------------+
# Input: tau: smoothness parameter; X: n by p matrix, original space; 
#        Y: m by p matrix 
# Output: n-dim vector
gaussian_kernel <- function(tau, X, Y){
  n <- dim(X)[1]
  p <- dim(X)[2]
  f <- function(x,y)  exp(-tau/(n*p)^2* t(x-y)%*%(x-y))
  outer( 1:nrow(Y), 1:nrow(X), Vectorize( function(j,i) f(Y[j,], X[i,]) ) )
}

#+-----------Gibbs sampler---------------+
#Input: X: data matrix; y: labels; B: number of iterations
#       tau: kernel parameter; lambda: prior smoothness parameter
#Output: posterior samples
bayes_svm_gibbs <- function(X, y, B = 5000, lambda = 1, tau = 2){
  n <- dim(X)[1]
  #initialization
  beta0 <- rep(0, n)
  v0 <- rep(1, n)
  #
  beta <- beta0  # dim-n
  v <- v0 # dim-n
  library(statmod) # Inverse gaussian
  library(mvtnorm) # multi-normal
  # iteration tracker
  b <-  1
  BETA <- matrix(nrow = B, ncol = n)
  K <- gaussian_kernel(tau, X, X) #dim n by n
  while(T){
    Z <- y*K/sqrt(v) #dim n by n
    w <- (1+v)/sqrt(v) # dim-n
    cov_mat <- solve(t(Z)%*%Z+lambda*diag(1, nrow = n)) #dim n by n
    mean_vec <- cov_mat%*%t(Z)%*%w  # dim-n
    
    # update beta
    beta <- as.numeric(rmvnorm( 1, mean_vec, cov_mat))
    # update v
    invgauss_mean <-  1/abs(1-y*K%*%beta)
    v <- 1/rinvgauss(n, invgauss_mean,dispersion = 1)  # dim-n
    # Output
    BETA[b,] <- beta
    # Control flow
    b <- b + 1
    # if (b%%10 == 0) cat(b,'~')
    if (b>B) break
  }
  BETA[-(1:1000),] # Burning
}



#+----------------Data processing-------------+
library(SIS)
data("leukemia.train")
data("leukemia.test")

X <- as.matrix(leukemia.train[,1:7129])
y <- leukemia.train[,7130]
y <- ifelse(y==0, 1, -1 )

#+-----------------------------------------+
# Question (a)
hyper_mat <- matrix(c(1,2,1,0.5,10,2,10,0.5), nrow = 4, byrow = T)
post_mean <- NULL
post_var <- NULL
for (i in 1:4){
  BETA_post <- bayes_svm_gibbs(X, y, lambda = hyper_mat[i,1], tau = hyper_mat[i,2],  B = 5000)
  post_mean <- rbind(post_mean, apply(BETA_post,2, mean))
  post_var <-  rbind(post_var, apply(BETA_post,2, var) )
}
# xtable::xtable(cbind(t(post_mean),t(post_var)))




#+--------Question (b)-----------------+
X_test <- as.matrix(leukemia.test[,1:7129])
y_test_01 <- leukemia.test[,7130]
y_test <- ifelse(y_test_01==0, 1, -1 )

hyper_mat <- matrix(c(1,2,
                      1,0.5,
                      10,2,
                      10,0.5), nrow = 4, byrow = T)
y_pred <- NULL
for (i in 1:4){
  BETA_post <- bayes_svm_gibbs(X, y, lambda = hyper_mat[i,1], tau = hyper_mat[i,2],  B = 5000)
  K_new <- gaussian_kernel(hyper_mat[i,2], X, X_test)
  PRED_MAT <- BETA_post%*%t(K_new) 
  y_pred <- rbind(y_pred,
                  apply(PRED_MAT, 2, 
                        function(x) ifelse(mean(x>0)>0.5, 1, -1)))
}


for (i in 1:4){
  y_predict <- factor(y_pred[i,], levels = c('-1', '1'))
  y_true <- factor(y_test, levels = c('-1', '1'))
  print(table(y_predict, y_true))
}


#+-------------- Question (c)----------------+
hyper_mat <- matrix(c(0.1, 2,
                      0.1, 4,
                      0.05, 2,
                      0.05, 4),
                    nrow = 4, byrow = T)
y_pred <- NULL
for (i in 1:4){
  BETA_post <- bayes_svm_gibbs(X, y, lambda = hyper_mat[i,1], tau = hyper_mat[i,2],  B = 5000)
  K_new <- gaussian_kernel(hyper_mat[i,2], X, X_test)
  PRED_MAT <- BETA_post%*%t(K_new) 
  y_pred <- rbind(y_pred,
                  apply(PRED_MAT, 2, 
                        function(x) ifelse(mean(x>0)>0.5, 1, -1)))
}

for (i in 1:4){
  y_predict <- factor(y_pred[i,], levels = c('-1', '1'))
  y_true <- factor(y_test, levels = c('-1', '1'))
  print(table(y_predict, y_true))
}
# library(MLmetrics)
# Specificity(y_true, y_predict)
# Sensitivity(y_true, y_predict)
# F1_Score(y_true = y_true, y_pred = y_predict )

