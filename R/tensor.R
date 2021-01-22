### Iterative Projection Estimation Functions of Tensor Factor Model

# list of functions:

# get x.hat from a tensor x and all loading matrics in list Q
# vector time series
# UP (Upfolding Procedure) one step for multi-dim tensor time series
# TIPUP one step for any Dim tensor time series
# TOPUP one step for any Dim tensor time series
# BIC estimators for determing the numbers of factors
# eigenvalue ratio estimators for determing the numbers of factors
# iterative versions of all three methods for any Dim tensor time series

#'functions for algorithm
#'
# get x.hat from a tensor x and all loading matrices in list Q by projection
#'@name get.hat
#'@rdname get.hat
#'@aliases get.hat
#'@export
#'@param x tensor of any dimension, x: \eqn{d1 * d2 * d3 * \cdots * d_k * n}
#'@param Q list of k matrices
#'@param d.seq 1:k
#'@return a tensor object x.hat
get.hat <- function(x,Q,d.seq){
  ans <- ttl(x,lapply(Q,t),d.seq)
  ans <- ttl(ans,Q,d.seq)
  ans
}

#'vector time series
#'
#'@name vts
#'@rdname vts
#'@aliases vts
#'@export
#'@param x a n*d matrix
#'@param h0 Pre-scribed parameter h
#'@param r First r eigenvectors
#'@return a list containing the following:\describe{
#'\item{\code{M}}
#'\item{\code{Q}}{The eigenvectors of matrix M}
#'}
vts <- function(x,h0,r){
  # for vector time series
  # x is a matrix: n * d
  n <- nrow(x)
  d <- ncol(x)
  M <- matrix(0,d,d)
  for(h in 1:h0){
    Omega <- t(x[1:(n-h),]) %*% x[(h+1):n,] / (n-h)
    M <- M + Omega %*% t(Omega)
  }
  Q <- eigen(M)$vectors[,1:r,drop=FALSE]
  list("M"=M,"Q"=Q)
}

#'UP one step for any Dim tensor time series
#'
#'@name up.init.tensor
#'@rdname up.init.tensor
#'@aliases up.init.tensor
#'@export
#'@param x tensor of any dimension : \eqn{d1 * d2 * d3 * \cdots * d_k * n}
#'@param r initial guess of # of factors
#'@param oneside.true  If oneside.true==TRUE, then only compute the one sided column space, not the other sides, this option is useful for the iterative method
#'@param norm.true If norm.true==TRUE, normalize the tensor
#'@return a list containing the following:\describe{
#'\item{\code{Q}}{Orthonormal matrix Q}
#'\item{\code{lambda}}{Singular values}
#'\item{\code{norm.percent}}{x after standardization}
#'\item{\code{x.hat}}{A tensor object x.hat}
#'}
up.init.tensor <- function(x,r,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # UP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_k * n
  # if oneside.true==TRUE, then only compute the one sided column space,
  # not the other sides, this option is useful for the iterative method
  dd <- dim(x)
  d <- length(dd)
  ans.Q <- NULL
  ans.lambda <- NULL
  if(oneside.true==TRUE){
    k=1
  } else{
    k=d-1
  }
  for(i in 1:k){
    xmat <- matrix(aperm(x,c(i,c(1:d)[-i])),nrow=dd[i])
    svd.xmat <- svd(xmat,nu=r[i],nv=0)
    ans.Q <- c(ans.Q,list(svd.xmat$u[,1:r[i],drop=FALSE]))
    ans.lambda <- c(ans.lambda,list(svd.xmat$d))
  }
  norm.percent <- NULL
  x.hat <- NULL
  if(norm.true==TRUE){
    x.tnsr <- as.tensor(x)
    x.hat <- get.hat(x.tnsr,ans.Q,1:k)
    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}

#'TIPUP one step for any Dim tensor time series
#'
#'@name tipup.init.tensor
#'@rdname tipup.init.tensor
#'@aliases tipup.init.tensor
#'@export
#'@param x tensor of any dimension \eqn{d1 * d2 * d3 * \cdots * d_d * n}
#'@param r initial guess of # of factors
#'@param h0 Pre-scribed parameter h
#'@param oneside.true  If oneside.true==TRUE, then only compute the one sided column space, not the other sides, this option is useful for the iterative method
#'@param norm.true If norm.true==TRUE, normalize the tensor
#'@return a list containing the following:\describe{
#'\item{\code{M}}{Estimator}
#'\item{\code{Q}}{Orthonormal matrix Q}
#'\item{\code{lambda}}{singular values}
#'\item{\code{norm.percent}}{x after standardization}
#'\item{\code{x.hat}}{A tensor object x.hat}
#'}
tipup.init.tensor <- function(x,r,h0=1,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # TIPUP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_d * n
  # if oneside.true==TRUE, then only compute the one sided column space,
  # not the other sides, this option is useful for the iterative method
  dd <- dim(x)
  d <- length(dd) # d >= 2
  n <- dd[d]
  dd.prod <- prod(dd) / n
  x.matrix <- matrix(x,ncol=n)
  if(oneside.true==TRUE){
    k=1
  } else{
    k=d-1
  }
  ans.M <- ans.Q <- ans.lambda <- NULL
  for(i in 1:k){
    M.temp <- matrix(0,dd[i],dd[i])
    for(h in 1:h0){
      x.left <- array(x.matrix[,1:(n-h)],c(dd[-d],n-h))
      x.right <- array(x.matrix[,(h+1):n],c(dd[-d],n-h))
      Omega <- tensor(x.left,x.right,c(1:d)[-i],c(1:d)[-i])/(n-h)
      M.temp <- M.temp + Omega %*% t(Omega)
    }
    #M.temp <- M.temp / dd.prod * dd[i]
    ans.eig <- eigen(M.temp)
    ans.M <- c(ans.M,list(M.temp))
    ans.Q <- c(ans.Q,list(ans.eig$vectors[,1:r[i],drop=FALSE]))
    ans.lambda <- c(ans.lambda,list(ans.eig$values))
  }
  norm.percent <- NULL
  x.hat <- NULL
  if(norm.true==TRUE){
    x.tnsr <- as.tensor(x)
    x.hat <- get.hat(x.tnsr,ans.Q,1:k)
    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("M"=ans.M,"Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}


#'TOPUP one step for any Dim tensor time series
#'
#'@name topup.init.tensor
#'@rdname topup.init.tensor
#'@aliases topup.init.tensor
#'@export
#'@param x tensor of any dimension
#'@param x \eqn{d1 * d2 * d3 * \cdots * d_d * n}
#'@param r initial guess of # of factors
#'@param h0 Pre-scribed parameter h
#'@param oneside.true  If oneside.true==TRUE, then only compute the one sided column space, not the other sides, this option is useful for the iterative method
#'@param norm.true If norm.true==TRUE, calculate the normalized residual of the tensor
#'@return a list containing the following:\describe{
#'\item{\code{M}}{Estimator}
#'\item{\code{Q}}{Orthonormal matrix Q}
#'\item{\code{lambda}}{singular values}
#'\item{\code{norm.percent}}{normalized residual}
#'\item{\code{x.hat}}{A tensor object x.hat}
#'}
topup.init.tensor <- function(x,r,h0=1,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # TOPUP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_d * n
  # if oneside.true==TRUE, then only compute the one sided column space,
  # not the other sides, this option is useful for the iterative method

  dd <- dim(x)
  d <- length(dd) # d >= 2
  n <- dd[d]
  dd.prod <- prod(dd) / n
  x.matrix <- matrix(x,ncol=n)
  if(oneside.true==TRUE){
    k=1
  } else{
    k=d-1
  }
  ans.M <- ans.Q <- ans.lambda <- NULL
  # allocate 0 for all k matrices with different sizes
  for(i in 1:k){
    M.temp <- matrix(0,dd[i],dd[i])
    ans.M <- c(ans.M,list(M.temp))
  }
  for(h in 1:h0){
    x.left <- array(x.matrix[,1:(n-h)],c(dd[-d],n-h))
    x.right <- array(x.matrix[,(h+1):n],c(dd[-d],n-h))
    Omega <- tensor(x.left,x.right,d,d)/(n-h)
    for(i in 1:k){
      ans.M[[i]] <- ans.M[[i]] + tensor(Omega,Omega,c(1:(2*(d-1)))[-i],c(1:(2*(d-1)))[-i])
    }
  }
  for(i in 1:k){
    ans.eig <- eigen(ans.M[[i]])
    ans.Q <- c(ans.Q,list(ans.eig$vectors[,1:r[i],drop=FALSE]))
    ans.lambda <- c(ans.lambda,list(ans.eig$values))
  }
  norm.percent <- NULL
  x.hat <- NULL
  if(norm.true==TRUE){
    x.tnsr <- as.tensor(x)
    x.hat <- get.hat(x.tnsr,ans.Q,1:(d-1))
    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("M"=ans.M,"Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}

#'iterative versions of all three methods for any Dim tensor time series
#'
#'@name iter.tensor.bic
#'@rdname iter.tensor.bic
#'@aliases iter.tensor.bic
#'@export
#'@param x tensor of any dimension, \eqn{x: d1 * d2 * d3 * \cdots * d_K * n}
#'@param r initial guess of # of factors
#'@param h0 Pre-scribed parameter h
#'@param method the method chosen among UP,TIPUP,TOPUP
#'@param tol level of error tolerance
#'@param niter number of iterations
#'@param tracetrue if TRUE, record the dis value for each iteration
#'@return a list containing the following:\describe{
#'\item{\code{Q}}
#'\item{\code{Qinit}}{initial value of Q}
#'\item{\code{Qfirst}}{value of Q at first iteration}
#'\item{\code{x.hat}}
#'\item{\code{x.hat.init}}{initial value of xhat}
#'\item{\code{x.hat.first}}{value of xhat at first iteration}
#'\item{\code{factor.num}}{number of factors}
#'\item{\code{timer}}{a timer}
#'\item{\code{norm.percent}}{x after standardization}
#'\item{\code{dis}}{difference of fnorm.resid for each iteration}
#'\item{\code{niter}}{number of interations}
#'\item{\code{fnorm.resid}}{normalized residual}
#'}
iter.tensor.bic <- function(x,r,h0=1,method,tol=1e-4,niter=100,tracetrue=FALSE){
  # UP iterative
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]
  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  timer <- rep(NA,3)
  time.now <- proc.time()[3]
  factor.num <- array(NA, c(d-1,5,niter))
  if(method=="UP"){
    ans.init <- up.init.tensor(x,r,norm.true=TRUE)
  }
  if(method=="TIPUP"){
    ans.init <- tipup.init.tensor(x,r,h0,norm.true=TRUE)
  }
  if(method=="TOPUP"){
    ans.init <- topup.init.tensor(x,r,h0,norm.true=TRUE)
  }
  ddd=dd[-d]
  for(i in 1:(d-1)){
    factor.num[i,,1]=tensor.bic(ans.init$lambda[[i]],h0,ddd[i],ddd[-i],n)
  }
  timer[1] <- proc.time()[3]-time.now
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,niter)
  x.hat <- get.hat(x.tnsr,ans.init$Q,d.seq)
  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q
  time.now <- proc.time()[3]
  while((dis > tol) & (iiter < niter)){             #while((dis > tol) & (iiter < niter)){
    for(i in 1:(d-1)){
      x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
      if(method=="UP"){
        ans.Q[[i]] <- up.init.tensor(x.new,c(r[i],r[-i]),oneside.true=TRUE,norm.true=FALSE)$Q[[1]]
      }
      if(method=="TIPUP"){
        ans.iter <- tipup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
        ans.Q[[i]] <- ans.iter$Q[[1]]
      }
      if(method=="TOPUP"){
        ans.iter <- topup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
        ans.Q[[i]] <- ans.iter$Q[[1]]
      }
      ddd=dd[-d]
      factor.num[i,,1+iiter]=tensor.bic(ans.iter$lambda[[1]],h0,ddd[i],ddd[-i],n)
      r[i]=factor.num[i,3,1+iiter]
    }

    x.hat <- get.hat(x.tnsr,ans.Q,d.seq)
    fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
    dis <- abs(fnorm.resid[iiter+1] - fnorm.resid[iiter])
    if(iiter==1){
      Qfirst <- ans.Q
      x.hat.first <- x.hat@data
      timer[2] <- proc.time()[3]-time.now
    }
    iiter <- iiter + 1
    if(tracetrue==TRUE){
      print(Q[[1]])
      cat("iiter=",iiter,", dis=",dis,"\n", sep="")
    }
    #    if(prod(factor.num[,3,iiter]==factor.num[,3,iiter-1])){
    #      break
    #    }
  }
  factor.num[,,niter]=factor.num[,,iiter-1]
  timer[3] <- proc.time()[3]-time.now
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]
  norm.percent <- c(ans.init$norm.percent,fnorm.resid[1],tail(fnorm.resid,1))
  list("Q"=ans.Q,"Qinit"=ans.init$Q, "Qfirst"=Qfirst,
       "x.hat"=x.hat@data,"x.hat.init" = ans.init$x.hat, "x.hat.first" =x.hat.first,
       "factor.num"=factor.num, ###"eigen.gap"=eigen.gap,
       "timer"=timer,
       "norm.percent"=norm.percent,
       "dis"=dis,"niter"=iiter,"fnorm.resid"=fnorm.resid)
}



#'iterative versions of all three methods for any Dim tensor time series
#'
#'@name iter.tensor.ratio
#'@rdname iter.tensor.ratio
#'@aliases iter.tensor.ratio
#'@export
#'@param x tensor of any dimension , \eqn{d1 * d2 * d3 * \cdots * d_K * n}
#'@param r initial guess of # of factors
#'@param h0 Pre-scribed parameter h
#'@param method the method chosen among UP,TIPUP,TOPUP
#'@param tol level of error tolerance
#'@param niter number of interations
#'@param tracetrue if TRUE, record the dis value for each iteration
#'@return a list containing the following:\describe{
#'\item{\code{Q}}
#'\item{\code{Qinit}}{initial value of Q}
#'\item{\code{Qfirst}}{value of Q at first iteration}
#'\item{\code{x.hat}}
#'\item{\code{x.hat.init}}{initial value of xhat}
#'\item{\code{x.hat.first}}{value of xhat at first iteration}
#'\item{\code{factor.num}}{number of factors}
#'\item{\code{timer}}{a timer}
#'\item{\code{norm.percent}}{x after standardization}
#'\item{\code{dis}}{difference of fnorm.resid for each iteration}
#'\item{\code{niter}}{number of interations}
#'\item{\code{fnorm.resid}}{normalized residual}
#'}
iter.tensor.ratio <- function(x,r,h0=1,method,tol=1e-4,niter=100,tracetrue=FALSE){
  # x: tensor of any dimension
  # UP iterative
  # x: d1 * d2 * d3 * ... * d_K * n
  # method: UP, TIPUP, or TOPUP
  # r: initial estimator of # of factors
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]
  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  timer <- rep(NA,3)
  time.now <- proc.time()[3]
  factor.num <- array(NA, c(d-1,5,niter))
  if(method=="UP"){
    ans.init <- up.init.tensor(x,r,norm.true=TRUE)
  }
  if(method=="TIPUP"){
    ans.init <- tipup.init.tensor(x,r,h0,norm.true=TRUE)
  }
  if(method=="TOPUP"){
    ans.init <- topup.init.tensor(x,r,h0,norm.true=TRUE)
  }
  ddd=dd[-d]
  for(i in 1:(d-1)){
    factor.num[i,,1]=tensor.ratio(ans.init$lambda[[i]],h0,ddd[i],ddd[-i],n)
  }
  timer[1] <- proc.time()[3]-time.now
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,niter)
  x.hat <- get.hat(x.tnsr,ans.init$Q,d.seq)
  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q
  time.now <- proc.time()[3]
  while((dis > tol) & (iiter < niter)){                 #while((dis > tol) & (iiter < niter)){
    for(i in 1:(d-1)){
      x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
      if(method=="UP"){
        ans.Q[[i]] <- up.init.tensor(x.new,c(r[i],r[-i]),oneside.true=TRUE,norm.true=FALSE)$Q[[1]]
      }
      if(method=="TIPUP"){
        ans.iter <- tipup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
        ans.Q[[i]] <- ans.iter$Q[[1]]
      }
      if(method=="TOPUP"){
        ans.iter <- topup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
        ans.Q[[i]] <- ans.iter$Q[[1]]
      }
      ddd=dd[-d]
      factor.num[i,,1+iiter]=tensor.ratio(ans.iter$lambda[[1]],h0,ddd[i],ddd[-i],n)
      r[i]=factor.num[i,1,1+iiter]
    }
    x.hat <- get.hat(x.tnsr,ans.Q,d.seq)
    fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
    dis <- abs(fnorm.resid[iiter+1] - fnorm.resid[iiter])
    if(iiter==1){
      Qfirst <- ans.Q
      x.hat.first <- x.hat@data
      timer[2] <- proc.time()[3]-time.now
    }
    iiter <- iiter + 1
    if(tracetrue==TRUE){
      print(Q[[1]])
      cat("iiter=",iiter,", dis=",dis,"\n", sep="")
    }
    #    if(prod(factor.num[,1,iiter]==factor.num[,1,iiter-1])){
    #      break
    #    }
  }
  factor.num[,,niter]=factor.num[,,iiter-1]
  timer[3] <- proc.time()[3]-time.now
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]
  norm.percent <- c(ans.init$norm.percent,fnorm.resid[1],tail(fnorm.resid,1))
  list("Q"=ans.Q,"Qinit"=ans.init$Q, "Qfirst"=Qfirst,
       "x.hat"=x.hat@data,"x.hat.init" = ans.init$x.hat, "x.hat.first" =x.hat.first,
       "factor.num"=factor.num, ###"eigen.gap"=eigen.gap,
       "timer"=timer,
       "norm.percent"=norm.percent,
       "dis"=dis,"niter"=iiter,"fnorm.resid"=fnorm.resid)
}

#'BIC estimators for determing the numbers of factors
#'@name tensor.bic
#'@rdname tensor.bic
#'@aliases tensor.bic
#'@export
#'@param reigen list of eigenvalues
#'@param h0 Pre-scribed parameter h
#'@param p1 p1
#'@param p2 p2
#'@param n n
#'@return factor.p1: Estimated number of factors
tensor.bic<-function(reigen,h0=1,p1,p2,n){
  #p1
  #p2
  #n
  if(length(p2)>1){
    p2=prod(p2)
  }
  factor.p1=numeric(8)
  p=p1*p2
  m1=ceiling(p1/3)

  lambda.p1<-reigen[p1:1]
  cumlambda.p1<-cumsum(lambda.p1)
  cumlambda.p1<-cumlambda.p1[(p1-1):1]

  #threshold
  ic=cumlambda.p1[1:m1]/p^2+(1:m1)*h0*(1/n)*log(p*n/(p+n))
  factor.p1[1]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2+(1:m1)*h0*(1/n)*log(min(p,n))
  factor.p1[2]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2+(1:m1)*h0*(1/n+1/p)*log(p*n/(p+n))
  factor.p1[3]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2+(1:m1)*h0*(1/n+1/p)*log(min(p,n))
  factor.p1[4]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2+(1:m1)*h0*(1/n+1/p)*log(min(p1,n))
  factor.p1[5]<-which.min(ic)

  factor.p1
}


#'Eigenvalue ratio estimators for determining the numbers of factors
#'@name tensor.ratio
#'@rdname tensor.ratio
#'@aliases tensor.ratio
#'@export
#'@param reigen list of eigenvalues
#'@param h0 Pre-scribed parameter h
#'@param p1 p1
#'@param p2 p2
#'@param n n
#'@return factor.p1: Estimated number of factors
#'
tensor.ratio<-function(reigen,h0=1,p1,p2,n){

  if(length(p2)>1){
    p2=prod(p2)
  }
  factor.p1=numeric(5)
  p=p1*p2
  m1=ceiling(p1/3)

  lambda.p1<-reigen[p1:1]
  cumlambda.p1<-cumsum(lambda.p1)
  cumlambda.p1<-cumlambda.p1[(p1-1):1]

  #ratio
  ratio<-(lambda.p1[(p1-1):(p1-m1)] +h0*0.01)/(lambda.p1[p1:(p1-m1+1)] +h0*0.01)
  factor.p1[1]<-which.min(ratio)
  ratio<-(lambda.p1[(p1-1):(p1-m1)] +p^2/n^2)/(lambda.p1[p1:(p1-m1+1)] +p^2/n^2)  #p^2*1/n
  factor.p1[2]<-which.min(ratio)
  ratio<-(lambda.p1[(p1-1):(p1-m1)] +p^2*(1/n^2/p1^2))/(lambda.p1[p1:(p1-m1+1)] +p^2*(1/n^2/p1^2))  #p^2*(1/n+1/p)
  factor.p1[3]<-which.min(ratio)
  ratio<-(lambda.p1[(p1-1):(p1-m1)] +p^2*(1/n^2/p1^2+1/n^2/p2^2))/(lambda.p1[p1:(p1-m1+1)] +p^2*(1/n^2/p1^2+1/n^2/p2^2)) #p^2*(1/n+1/p1)
  factor.p1[4]<-which.min(ratio)
  ratio<-(lambda.p1[(p1-1):(p1-m1)] +p^2*(1/n^2/p1+1/n^2/p2))/(lambda.p1[p1:(p1-m1+1)] +p^2*(1/n^2/p1+1/n^2/p2))  #p^2*(1/n)*log(p*n/(p+n))
  factor.p1[5]<-which.min(ratio)

  factor.p1
}



























#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
## the following functions are no longer algorithms, they are for simulations
# they are exactly the same as v5 and v6

# distance function with three types
dist.fn <- function(H1,H2,dist.type="all"){
  # H1 and H2 are orthonormal
  # compute three distances between them
  # type="qiwei","spectra","frob","all"
  error <- rep(NA,3)
  maxd <- max(ncol(H1),ncol(H2))
  error.matrix <- H1 %*% t(H1) - H2 %*% t(H2)
  if(dist.type=="qiwei" | dist.type=="all"){
    # (1-1/max(q_1,q_2) * tr(H_1H_1^TH_2H_2^T))^{1/2} from Qiwei's paper
    error[1] <- sqrt(max(1 - 1/maxd*sum(diag(H1 %*% t(H1) %*% H2 %*% t(H2))),0))
  }
  if(dist.type=="spectra" | dist.type=="all"){
    # \|H_1H_1^T-H_2H_2^T\|_spectra
    error[2] <- svd(error.matrix, nu=0,nv=0)$d[1]
  }
  if(dist.type=="frob" | dist.type=="all"){
    # \|H_1H_1^T-H_2H_2^T\|_F
    error[3] <- sqrt(sum(error.matrix^2))
  }
  error
}

# generate orthonormal matrix Q with size p*r
gen.Q <- function(p,r){
  # A: p*r, iid N(0,1)
  # return Q: p*r, qr(A)$Q
  A <- matrix(rnorm(p*r),p,r)
  Q <- qr.Q(qr(A))
}

# generate a matrix F with AR coefficients of length n
gen.F <- function(ar1.coef,n,innovsd){
  # ar1.coef: r1*r2
  # F: r1*r2*n
  r1r2 = dim(ar1.coef)
  r1 = r1r2[1]
  r2 = r1r2[2]
  bF = array(NA,dim=c(r1,r2,n))

  for(ir in 1:r1){
    for(jr in 1:r2){
      bF[ir,jr,] = arima.sim(n=n, model=list(ar=ar1.coef[ir,jr]),sd = innovsd)
    }
  }
  bF
}



gen.F.tensor3 <- function(ar1.coef,n,innovsd){
  # ar1.coef: r1*r2*r3
  # F: r1*r2*r3*n
  r = dim(ar1.coef)
  r1 = r[1]
  r2 = r[2]
  r3 = r[3]
  bF = array(NA,dim=c(r1,r2,r3,n))

  for(ir1 in 1:r1){
    for(ir2 in 1:r2){
      for(ir3 in 1:r3){
        bF[ir1,ir2,ir3,] = arima.sim(n=n, model=list(ar=ar1.coef[ir1,ir2,ir3]),sd = innovsd)
      }
    }
  }
  bF
}


# generate the error tensor E with given two covariance matrices
gen.E <- function(n,p1,p2,cov1.half,cov2.half){
  # cov1.half: p1 * p1
  # cov2.half: p2 * p2
  # return p1 * p2 * n
  E <- array(rnorm(n*p1*p2),c(p1,p2,n)) # p1 * p2 * n
  E <- tensor(E,cov1.half,1,1) # p2 * n * p1
  E <- tensor(E,cov2.half,1,1) # n * p1 * p2
  aperm(E,c(2,3,1)) #p1 * p2 * n
}

# generate the error tensor E with given three covariance matrices
gen.E.tensor3 <- function(n,p1,p2,p3,cov1.half,cov2.half,cov3.half){
  # cov1.half: p1 * p1
  # cov2.half: p2 * p2
  # cov3.half: p3 * p3
  # return p1 * p2 * p3 * n
  E <- array(rnorm(n*p1*p2*p3),c(p1,p2,p3,n)) # p1 * p2 * p3 * n
  E <- tensor(E,cov1.half,1,1) # p2 * p3 * n * p1
  E <- tensor(E,cov2.half,1,1) # p3 * n * p1 * p2
  E <- tensor(E,cov3.half,1,1) # n * p1 * p2 * p3
  aperm(E,c(2:4,1)) #p1 * p2 * n
}

# generate X = bF \times_1 Q1 \times_2 Q2
gen.X <- function(Q1,Q2,bF){
  # Q1: p1*r1 orthonormal
  # Q2: p2*r2 orthonormal
  # bF: r1*r2*n

  X <- tensor(bF, Q1, 1, 2) # r2*n*p1
  X <- tensor(X, Q2, 1, 2) # n*p1*p2
  X <- aperm(X,c(2,3,1)) # p1*p2*n
  X
}

# generate X = bF \times_1 Q1 \times_2 Q2 \times_3 Q3
gen.X.tensor3 <- function(Q1,Q2,Q3,bF){
  # Q1: p1*r1 orthonormal
  # Q2: p2*r2 orthonormal
  # Q3: p3*r3 orthonormal
  # bF: r1*r2*r3*n

  X <- tensor(bF, Q1, 1, 2) # r2*r3*n*p1
  X <- tensor(X, Q2, 1, 2) # r3*n*p1*p2
  X <- tensor(X, Q3, 1, 2) # n*p1*p2*p3
  X <- aperm(X,c(2:4,1)) # p1*p2*p3*n
  X
}

# Generate covariance matrix with diagonal = 1, off-diagonal == offdiag.val
gen.cov.eqaloffdiag <- function(p,offdiag.val){
  cov.matrix <- matrix(offdiag.val,p,p) + diag(rep(1-offdiag.val,p))
  cov.matrix
}

# find squart root of matrix x
sqrt.m <- function(x){
  x.eig <- eigen(x)
  values <- x.eig$values
  values.sqrt <- rep(0,nrow(x))
  values.sqrt[values > 0] <- sqrt(values[values > 0])
  x.eig$vectors %*% diag(values.sqrt) %*% t(x.eig$vectors)
}

# compute the errors for estimations of Q1, Q2, and Q for three types of distances
error.space.fn.tensor <- function(Qlist,Qhatlist,Q){
  # Q=kronecker(Q1, Q2, Q3,...)
  # Qlist, Qhatlist: list of orthonormal matrices
  d <- length(Qlist)
  error <- rep(NA, d+1)
  for(i in 1:d){
    error[i] <- dist.fn(Qhatlist[[i]],Qlist[[i]], dist.type = "spectra")[2]
  }
  Qhat <- kronecker_list(Qhatlist)
  error[d+1] <- dist.fn(Qhat,Q, dist.type = "spectra")[2]
  error
}

error.all.fn <- function(ans, Qlist, Q, x0){
  # return: error: 3 columns for ans$Qinit, ans$Qfirst and ans$Q respectively
  k <- length(ans$Q)
  error <- matrix(NA,k+4,3)
  error[1:(k+1),1] <- error.space.fn.tensor(Qlist,ans$Qinit,Q)
  error[1:(k+1),2] <- error.space.fn.tensor(Qlist,ans$Qfirst,Q)
  error[1:(k+1),3] <- error.space.fn.tensor(Qlist,ans$Q,Q)
  error[k+2,] <- ans$norm.percent
  error[k+3,1] <- sqrt(sum((ans$x.hat.init-x0)^2))/sqrt(sum((x0)^2))
  error[k+3,2] <- sqrt(sum((ans$x.hat.first-x0)^2))/sqrt(sum((x0)^2))
  error[k+3,3] <- sqrt(sum((ans$x.hat-x0)^2))/sqrt(sum((x0)^2))
  error[k+4,] <- ans$timer
  error
}
