# Generate Beta-binomial
n=10
a<-2
b<-2
N<-50000
NN<-100000

# Direct method
znew<-rbeta(N,a,b)
xnew<-rbinom(N,n,znew)
table(xnew)/N

# Gibbs Sampler
zx<-0.5
xznew<-NULL
zxnew<-NULL
for(i in 1:NN)
{
  xz<-rbinom(1,n,zx)
  xznew<-c(xznew,xz)
  zx<-rbeta(1,a+xz,n-xz+b)
  zxnew<-c(zxnew,zx)
}
x_gibbs <- xznew[50001:NN]
z_gibbs <- zxnew[50001:NN]
table(x_gibbs)/N


# Generate t-distribution

nu<-2
sig2<-2
mu<-1

# Direct method
znew<-rgamma(N,nu/2,nu/2)
xnew<-rnorm(N,mu,sqrt(sig2/z4new))
q.direct <- quantile(xnew,seq(0.05, 0.95, 0.05))

# Gibbs sampler
zx<-0.5
xznew<-NULL
zxnew<-NULL
for(i in 1:NN)
{
  xz<-rnorm(1,mu,sqrt(sig2/zx))
  xznew<-c(xznew,xz)
  zx<-rgamma(1,(niu+1)/2,((xz-mu)^2)/(2*sig2)+ nu/2)
  zxnew<-c(zxnew,zx)
}
x_gibbs<-xznew[N+1:NN]
z_gibbs<-zxnew[N+1:NN]
q.gibbs <- quantile(x_gibbs,seq(0.05, 0.95, 0.05))
# Compare through QQ plot
plot(q.direct, q.gibbs)
abline(0,1)


z <- rnorm(50000)
q.norm <- quantile(z, seq(0.05, 0.95, 0.05))
plot(q.direct, q.norm)
abline(0,1)
qqplot(z, xnew)
