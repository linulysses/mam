################  Manifold Additive Regression #################

#' Manifold additive model
#' @param train a list of two members, X, a matrix of which each represents a sample, Y, the response represented by a d*d*n array where d is the matrix dimension and n is the sample size
#' @param valid optional validation dataset of the same structure of \code{train}
#' @param h.rng the range of bandwidth
#' @param h a optionally specified bandwidth; automatically selected if not provided
#' @param X_test observations of predictors at which the regression is estimated
#' @param cf.return whether component functions are returned
#' @param mfd manifold structure; default is SPD with LogCholesky metric
mam <- function(train,valid=NULL,h.rng=c(1e-4,1),h=NULL,X_test=NULL,cf.return=F,mfd=NULL)
{
    tmp <- transform(train$Y,mfd=mfd)
    muhat <- tmp$mu
    train$Y <- tmp$R
    if(!is.null(valid))
    {
        valid$Y <- transform(valid$Y,muhat)$R
    }
    
    if(is.null(X_test)) X_test <- train$X
    
    n <- dim(train$Y)[3]
    d <- dim(train$Y)[1]
    R <- array(0,c(d,d,nrow(X_test)))
    f <- NULL  # estimated f function
    
    sel.h <- array(0,c(d,d,dim(train$X)[2]))
    if(!is.null(h)) sel.h <- h
    
    # if sel.h is a scalar
    if(is.null(dim(sel.h))) sel.h <- array(sel.h,c(d,d,dim(train$X)[2]))
    
    for(i in 1:d)
        for(j in 1:i)
        {
            res <- estimate(train,i,j,valid,h.rng,sel.h[i,j,],X_test)
            
            # if bandwidth is too small, then try to increase it
            maxit <- 10
            it <- 1
            while(any(is.nan(res$yhat)) || any(is.na(res$yhat)) || any(is.infinite(res$yhat)))
            {
                h0 <- res$h * 1.5
                res <- estimate(train,i,j,valid,h.rng,h0,X_test)
                it <- it + 1
                if(it > maxit) stop('I tried my best ...')
            }

            R[i,j,] <- res$yhat
            R[j,i,] <- res$yhat
            
            sel.h[i,j,] <- res$h
            sel.h[j,i,] <- res$h
            
            if(is.null(f) && cf.return) 
            {
                f <- array(0,c(d,d,ncol(res$mhat), ncol(X_test)))
            }
            if(cf.return){
                f[i,j,,] <- t(res$mhat[,,1])
                f[j,i,,] <- t(res$mhat[,,1])
            } 
        }
    
    
    if(cf.return)
    {
        w <- array(NA,dim(f))
        
        for(k in 1:dim(f)[4])
        {
            for(j in 1:dim(f)[3])
            {
                tmp <- parallel.transport(mfd,muhat,diag(rep(1,d)),f[,,j,k])
                tmp <- rie.exp(mfd,diag(rep(1,d)),tmp) 
                if(any(is.infinite(tmp)) || any(is.na(tmp)) || any(is.nan(tmp)))
                {
                    w[,,j,k] <- matrix(NA,d,d)
                }
                else
                {
                    E <- eigen(tmp)$values
                    if(any(E <= 1e-10))
                    {
                        tmp <- tmp + diag(rep(max(1e-10,abs(min(E))),d))
                    }
                    w[,,j,k] <- tmp
                }
            }
            
        }
        
        return(list(yhat=transform.back(R,muhat,mfd),muhat=muhat,what=w,fhat=f,h=sel.h))
    }
    else transform.back(R,muhat,mfd)
}

# transform data from manifold to tangent spaces
transform <- function(dat,mu=NULL,mfd=NULL)
{
    n <- dim(dat)[3]
    d <- dim(dat)[1]
    if(is.null(mu)) mu <- frechet.mean(mfd,dat)
    R <- array(0,dim(dat))
    for(i in 1:n)
    {
        R[,,i] <- rie.log(mfd,mu,dat[,,i])
    }
    return(list(R=R,mu=mu))
}


# transform tangent vectors onto manifold
transform.back <- function(R,muhat,mfd)
{
    n <- dim(R)[3]
    d <- dim(R)[1]
    for(i in 1:n)
    {
        R[,,i] <- rie.exp(mfd,muhat,R[,,i])
        if(any(is.infinite(R[,,i])) || any(is.na(R[,,i])) || any(is.nan(R[,,i])))
        {
            R[,,i] <- matrix(0,d,d)
        }
        else
        {
            E <- eigen(R[,,i])$values
            if(any(E <= 1e-10))
            {
                R[,,i] <- R[,,i] + diag(rep(max(1e-10,abs(min(E))),d))
            }
        }
    }
    return(R)
}

# estimate component functions coordinate-by-coordinate
estimate <- function(train,i,j,valid=NULL,h.rng=c(0,1),h=NULL,X_test=NULL)
{
    Y_training <- train$Y[i,j,]
    X_training <- train$X
    if(!is.null(valid))
    {
        X_valid <- valid$X
        Y_valid <- valid$Y[i,j,]
    }
    else
    {
        X_valid <- NULL
        Y_valid <- NULL
    }
    
    d <- dim(train$Y[,,1])[1]
    
    D <- dim(X_training)[2]
    
    if(is.null(h) || any(h==0))
    {
        h <- select.bw(X_training,Y_training,X_valid,Y_valid,h.rng)
        h <- rep(h,D)
        
    }
    
    if(is.null(X_test))
        X_test <- X_training 
    
    result <- SBF(X_test=X_test,X_training=X_training,Y_training=Y_training,h=h)
    result$h <- h
    return(result)
    
}

# select bandwith by validation dataset or by cross-validation 
select.bw <- function(X_training,Y_training,X_valid,Y_valid,h.rng)
{
    
    rmse <- function(v1,v2) sqrt(mean((v1-v2)^2))
    
    if(!is.null(X_valid))
    {
        hs <- seq(h.rng[1],h.rng[2],length.out=20)
        err <- rep(0,length(hs))
        D <- dim(X_training)[2]
        if(!is.null(X_valid))
        {
            for(k in 1:length(hs))
            {
                h <- rep(hs[k],D)
                result <- SBF(X_test=X_valid,X_training=X_training,Y_training=Y_training,h=h)
                err[k] <- rmse(Y_valid,result$yhat)
            }
        }
    }
    else
    {
        hs <- seq(h.rng[1],h.rng[2],length.out=20)
        n <- nrow(X_training)
        K <- 5
        cv <- cv.partition(n,K)
        err <- matrix(0,length(hs),K)
        D <- dim(X_training)[2]
        for(j in 1:length(hs))
            for(k in 1:K)
            {
                tr.idx <- cv$training[[k]]
                te.idx <- cv$test[[k]]
                X <- X_training[tr.idx,]
                Y <- Y_training[tr.idx]
                h <- rep(hs[j],D)
                result <- SBF(X_test=X_training[te.idx,],X_training=X,Y_training=Y,h=h)
                err[j,k] <- rmse(Y_training[te.idx],result$yhat)
            }
        err <- apply(err,1,mean)
    }
    return(hs[which.min(err)])
}

# generate an orthonormal basis of the tangent space at mu
ortho.basis <- function(d,mu,metric)
{
    #raw <- gen.sym(d,n=d*(d+1)/2)
    sym.mfd <- matrix.manifold('sym','Frobenius',d)
    raw <- rmatrix(sym.mfd,n=d*(d+1)/2)
    mfd <- matrix.manifold('spd',metric,d)
    B <- array(0,dim(raw))
    n <- dim(raw)[3]
    B[,,1] <- raw[,,1] / sqrt(rie.metric(mfd,mu,raw[,,1],raw[,,1]))
    for(j in 2:n)
    {
        B[,,j] <- raw[,,j]
        for(k in 1:(j-1))
        {
            B[,,j] <- B[,,j] - rie.metric(mfd,mu,B[,,k],B[,,j]) * B[,,k]
        }
        B[,,j] <- B[,,j] / sqrt(rie.metric(mfd,mu,B[,,j],B[,,j]))
    }
    B
}


#' the function mimic the cvpartition in MATLAB
#' @param n sample size
#' @param K the number of cv folders
#' @return a list of the follow elements
#' \item{training}{a list of \code{K} elements. Each element is a vector of indices of samples for training}
#' \item{test}{a list of \code{K} elements. Each element is a vector of indices of samples for test}
#' \item{num.test.sets}{an integer signifying the number of test sets}
#' @export cv.partition
cv.partition <- function(n,K)
{
    if(K <= 1) stop('K must be larger than 1')
    
    I <- sample.int(n)
    a <- n %% K
    b <- (n-a)/K
    
    trunk.size <- c(rep(b+1,a),rep(b,K-a))
    start.idx <- c(1,1+cumsum(trunk.size[1:(K-1)]))
    end.idx <- start.idx + trunk.size - 1
    
    test <- lapply(1:K,function(k){
        I[start.idx[k]:end.idx[k]]
    })
    
    training <- lapply(1:K,function(k){
        I[-(start.idx[k]:end.idx[k])]
    })
    
    return(list(training=training,test=test,num.test.sets=length(test)))
}


################  ANOVA  ###################################


# this is the kernel used in 'SBF_continuous_pointwise.dll'
epan <- function(u) ifelse( (u>-1) & (u<1), 3*(1-u^2)/4, 0 )

#' nornalization constant
#' @param K the kernel function
#' @param v the location
#' @param h the bandwidth
nconst <- function(K,v,h)
{
    a <- max(-1,-v/h)
    b <- min(1,(1-v)/h)
    L <- 100
    r <- (b-a)/L
    
    sum(K((seq(a+r/2,b-r/2,length.out = L)))*r)
}

# kernel function in Jeon & Park (2020)
kern <- function(u,v,h)
{
    a <- nconst(epan,v,h)
    epan((u-v)/h)/ h / a
}

#' @param X the data matrix of predictors, scaled to [0,1]
#' @param Y d*d*n array of SPD responses
#' @param Yhat the fitted responses
#' @param h used bandwidth to fit MAM
#' @param newx the vector of new predictor
#' @param fhat the estimate component function (in tangent space at muhat) at newx
#' @param muhat the estimated Frechet mean
#' @param metric Riemannian metric
anova <- function(X,Y,Yhat,newx,fhat,muhat,h,metric)
{
    n <- nrow(X)
    q <- ncol(X)
    d <- dim(Y)[1]
    
    a <- n^{1/5}*h
    K2 <- 0.6
    
    px <- sapply(1:q,function(k){
        mean(kern(X[,k],newx[k],h[k]))
    })
    
    B <- ortho.basis(d,muhat,metric)
    mfd <- matrix.manifold('spd',metric,d)
    
    D <- dim(B)[3]
    
    # project data into tangent spaces and the selected orthonormal basis
    U <- matrix(0,n,D)
    for(i in 1:n)
    {
        Zhat <- rie.log(mfd,muhat,Yhat[,,i])
        Z <- rie.log(mfd,muhat,Y[,,i])
        E <- Z - Zhat
        
        for(j in 1:D)
            U[i,j] <- rie.metric(mfd,muhat,E,B[,,j])
    }
    
    
    # estimate covariance structure (in tangent space)
    Sig <- sum(1/(a*px)) * K2 * cov(U)
    
    
    
    V <- rep(0,D)
    for(k in 1:q)
        for(j in 1:D) V[j] <- V[j] + rie.metric(mfd,muhat,fhat[,,k],B[,,j])
    
    Z <- solve(expm::sqrtm(Sig)) %*% V
    chi <- n^{4/5} * sum(Z^2) # follow chi-square dist
    
    # return p value
    pchisq(chi,df=D,lower.tail = F)
}



################  Smooth Backfitting #######################


dyn.load('SBF_continuous_pointwise.dll')

# Function for obtaining additional grid points for numerical integration
# X_test: N times d matrix whose elements take values in [0,1]
# We use sort(X_test[,j],seq(0,1,length=add)) as the grid points for numerical integration for jth integral
integration.grid=function(X_test,add)
{
    add_vector=seq(0,1,length=add)
    added_X_test=matrix(,add+nrow(X_test),ncol(X_test))
    for(j in 1:ncol(X_test))
    {
        added_X_test[,j]=sort(c(add_vector,X_test[,j]))
    }
    added_X_test
}

# Function for fitting smooth backfitting
# X_test: N times d matrix whose elements take values in [0,1]
# X_training: n times d matrix whose elements take values in [0,1]
# Y_training: n-dimensional vector
# add: number of additional grid points for numerical integration
# max_iteration: maximum iteration number for smooth bacfkfitting algorithm
# epsilon: convergence criterion for smooth backfitting algorithm
# h: d-dimensional vector of bandwidths
SBF=function(X_test,X_training,Y_training,add=0,max_iteration=50,epsilon=10^-4,h)
{
    n=nrow(X_training)
    d=ncol(X_training)
    N=nrow(X_test)
    
    add=as.integer(add)
    max_iteration=as.integer(max_iteration)
    h=as.double(h)
    
    g=integration.grid(X_test,add)
    ngrid=nrow(g)
    actual_iteration=as.integer(0)
    T=as.integer(1)
    Y_training=matrix(Y_training,n,T)
    yhat=matrix(0,N,T)
    mhat=array(0,dim=c(d,ngrid,T))
    
    result=.Fortran('SBF_continuous_pointwise',n0=N,n=n,d=d,ngrid=ngrid,T=T,X0=X_test,X=X_training,Y=Y_training,g=g,h=h,
                    max_iteration=max_iteration,epsilon=epsilon,actual_iteration=actual_iteration,
                    yhat=yhat,mhat=mhat)
    
    return(list(yhat=result$yhat,mhat=result$mhat,actual_iteration=result$actual_iteration))
    if(any(result$actual_iteration==max_iteration)) print('Some smooth backfitting algorithm did not converge.')
}


