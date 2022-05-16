# Demonstration of MAM

library(MatrixManifold) # available at https://github.com/linulysses/matrix-manifold
dyn.load('SBF_continuous_pointwise.dll')
source('mam.R')


set.seed(1)


d <- 3  # dimension of SPD
D <- d*(d+1)/2  # intrinsic dimension of SPD(d)
sig <- 0.05  # noise level
mu <- diag(d) # Frechet mean
p <- 3  # number of predictors
n <- 100  # sample size
metric <- 'LogCholesky'

# regression function
f <- function(x) Reduce('+', lapply(1:p, function(k) matrix(sin(2*k*pi*x[k])/k,d,d)))

# sample predictors
X <- matrix(runif(n*p),n,p)

# sample size
N <- matrix(sig*rnorm(n*D),n,D)

# set up manifold structure
mfd <- matrix.manifold('spd',metric,d)

# an orthonormal basis at mu for adding noise
B <- ortho.basis(d,mu,metric)

# generate noisy response
Y <- array(0,c(d,d,n))
for(i in 1:n)
{
    # noise
    E <- 0
    for(j in 1:D) E <- E + N[i,j]*B[,,j]
    
    Y[,,i] <- rie.exp(mfd,mu,f(X[i,])+E)
}

# manifold additive regression
res <- mam(train=list(X=X,Y=Y),X_test=X,mfd=mfd,cf.return=T)

# sample MSE
mean(sapply(1:n, function(i) geo.dist(mfd,res$yhat[,,i],Y[,,i])^2))

# anova to test whether all component functions are zero
pval <- sapply(1:n, function(k) anova(X,Y,res$yhat, X[k,], res$fhat[,,k,], res$muhat, res$h, metric))
pval.adj <- p.adjust(pval,method='fdr') # adjusted p values
