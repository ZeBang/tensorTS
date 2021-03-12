#' Tensor factor estimation method
#'
#' Estimation function for tensor-valued time series factor model.
#'@name tenFM.est
#'@rdname tenFM.est
#'@aliases tenFM.est
#'@usage tenFM.est(x,r,method='TIPUP')
#'@export
#'@param x \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param r input rank of factor tensor
#'@param h0 lag number
#'@param dim dimension of coefficient matrices
#'@param method character string, specifying the type of the estimation method to be used. \describe{
#'  \item{\code{"TIPUP",}}{TIPUP method}
#'  \item{\code{"TOPUP",}}{TOPUP method}
#'}
#'@param iter boolean, specifying using an iterative approach or an non-iterative approach
#'@param vmax boolean, specifying using varimax rotation on the factor matrix or not
#'@param tol tolerance for each iteration
#'@param maxiter integer, the max number of iterations allowed
#'@return a list containing the following:\describe{
#'\item{\code{Ft}}{estimated factor matrices of dimension \eqn{T*d_1*d_2 *\cdots *d_K}}
#'\item{\code{Ft.all}}{Summation of factor matrices over time, dimension \eqn{d_1,d_2,\cdots,d_K} }
#'\item{\code{Q}}{loading matrices}
#'\item{\code{x.hat}}{fitted part}
#'\item{\code{niter}}{number of iterations}
#'\item{\code{fnorm.resid}}{Frobenius norm of residual}
#'}
#'@examples
#'Q1 = gen.Q(5,2)
#'Q2 = gen.Q(5,2)
#'ar_coef = array(c(0.3,0.5,0.4,0.6),c(2,2))
#'bF = gen.F(ar_coef,20,1)
#'x = aperm(gen.X(Q1,Q1,bF),c(3,1,2))
#'x = as.tensor(x) # Generate x and permute the axes
#'result <- tenFM.est(x,c(2,2),method='TOPUP')  # Estimate the factor and loadings
#'Ft <- result$Ft
tenFM.est=function(x,r,h0=1,method='TIPUP',iter=TRUE,vmax=FALSE,tol=1e-4,maxiter=100){
  x <- aperm(x@data,c(2:length(dim(x)),1))
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]

  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  eigen.gap <- array(NA, c(d-1,max(dd[-d]),maxiter))

  if(method=="TIPUP"){
    ans.init <- tipup.init.tensor(x,r,h0,norm.true=TRUE)
  }else if(method=="TOPUP"){
    ans.init <- topup.init.tensor(x,r,h0,norm.true=TRUE)
  }else{
    stop('Wrong method !')
  }
  ddd=dd[-d]
  for(i in 1:(d-1)){
    eigen.gap[i,1:dd[i],1]=ans.init$lambda[[i]]
  }
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter)
  x.hat <- get.hat(x.tnsr,ans.init$Q,d.seq)
  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q
  Ft <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
  Ft.all <- apply(Ft@data,c(1:(d-1)),sum)
  fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm

  if(iter==TRUE){
    while((dis > tol) & (iiter < maxiter)){
      for(i in 1:(d-1)){
        x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
        if(method=="TIPUP"){
          ans.iter <- tipup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else if(method=="TOPUP"){
          ans.iter <- topup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else{
          stop('Wrong method !')
        }
        ddd=dd[-d]
        eigen.gap[i,1:dd[i],1+iiter]=ans.iter$lambda[[1]]

        if(vmax==TRUE){
          for(j in 1:(d-1)){
            if(r[j] > 1){
              Qj <- ans.Q[[j]]
              Qj.ans <- varimax(Qj)
              Qjrot <- Qj %*% Qj.ans$rotmat
              for(k in 1:r[j]){
                Qjrot[,k]=Qjrot[,k]*sign(sum(Qjrot[,k]))
              }
              ans.Q[[j]] = Qjrot
            }
          }
        }

        x.hat <- get.hat(x.tnsr,ans.Q,d.seq)

        Ft <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
        Ft.all <- apply(Ft@data,c(1:(d-1)),sum)

        fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
        dis <- abs(fnorm.resid[iiter+1] - fnorm.resid[iiter])
        if(iiter==1){
          Qfirst <- ans.Q
          x.hat.first <- x.hat@data
        }
        iiter <- iiter + 1
      }
    }
  }else{
    iiter <- iiter + 1
  }
  eigen.gap[,,maxiter]=eigen.gap[,,iiter-1]
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]
  norm.percent <- c(ans.init$norm.percent,fnorm.resid[1],tail(fnorm.resid,1))
  return(list("Ft"=aperm(Ft@data,c(d,1:(d-1))),"Ft.all"=Ft.all,"Q"=ans.Q,"x.hat"=aperm(x.hat@data,c(d,1:(d-1))),"niter"=iiter,"fnorm.resid"=fnorm.resid[iiter]))
}


#' Tensor factor rank determination
#'
#' Function for rank determination of the factor
#'@name Rank.est
#'@rdname Rank.est
#'@aliases Rank.est
#'@usage Rank.est(x,r,rank='BIC',method='TIPUP')
#'@export
#'@param x \eqn{m_1 * \cdots * m_K * T} tensor-valued time series
#'@param r input rank of factor matrix
#'@param h0 lag number
#'@param rank character string, specifying the type of the rank determination method to be used. \describe{
#'  \item{\code{"BIC",}}{Bayesian information criterion}
#'  \item{\code{"ER",}}{Eigen ratio method}
#'}
#'@param method character string, specifying the type of the factor estimation method to be used. \describe{
#'  \item{\code{"TIPUP",}}{TIPUP method}
#'  \item{\code{"TOPUP",}}{TOPUP method}
#'}
#'@param inputr boolean, if TRUE, use input rank for each iteration; if FLASE, update the rank r in each iteration
#'@param iter boolean, specifying using an iterative approach or an non-iterative approach
#'@param penalty takes value in {1,2,3,4,5}, decide which penalty function to use
#'@param delta1 weakest factor strength, a tuning parameter used for BIC method only
#'@param tol tolerance for each iteration
#'@param maxiter integer, the max number of iterations allowed
#'@return a list containing the following:\describe{
#'\item{\code{factor.num}}{estimated number of factors \eqn{n_1,n_2,\cdots,n_K}}
#'}
#'@examples
#'Q1 = gen.Q(5,2)
#'Q2 = gen.Q(5,2)
#'ar_coef = array(c(0.3,0.5,0.4,0.6),c(2,2))
#'bF = gen.F(ar_coef,20,1)
#'x = aperm(gen.X(Q1,Q1,bF),c(3,1,2))
#'x = as.tensor(x) # Generate x and permute the axes
#'rank <- Rank.est(x,c(2,2),rank='BIC',method='TIPUP')  # Estimate the rank
Rank.est = function(x,r,h0=1,rank='BIC',method='TIPUP',inputr=FALSE,iter=TRUE,penalty=1,delta1=0,tol=1e-4,maxiter=100){
  x <- aperm(x@data,c(2:length(dim(x)),1))
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]
  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  factor.num <- array(NA, c(d-1,5,maxiter))
  eigen.gap <- array(NA, c(d-1,max(dd[-d]),maxiter ))

  if(method=="TIPUP"){
    ans.init <- tipup.init.tensor(x,r,h0,norm.true=TRUE)
  }else if(method=="TOPUP"){
    ans.init <- topup.init.tensor(x,r,h0,norm.true=TRUE)
  }
  else{
    stop('Wrong method !')
  }
  ddd=dd[-d]
  for(i in 1:(d-1)){
    if(rank=='BIC'){
      sigmas.hat=1
      factor.num[i,,1]=tenFM.bic(ans.init$lambda[[i]]/sigmas.hat,h0,ddd[i],ddd[-i],n,delta1)
    }else if(rank=='ER'){
      factor.num[i,,1]=tenFM.ratio(ans.init$lambda[[i]],ddd[i],ddd[-i],n)
    }else{
      stop('Wrong rank method !')  # should catch exception
    }
    eigen.gap[i,1:dd[i],1]=ans.init$lambda[[i]]
  }
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter)
  x.hat <- get.hat(x.tnsr,ans.init$Q,d.seq)
  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q

  if(iter==TRUE){
    while((dis > tol) & (iiter < maxiter)){
      for(i in 1:(d-1)){
        x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
        if(method=="TIPUP"){
          ans.iter <- tipup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else if(method=="TOPUP"){
          ans.iter <- topup.init.tensor(x.new,c(r[i],r[-i]),h0,oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else{
          stop('Wrong estimation method input !')  # should catch exception
        }
        ddd=dd[-d]

        if(rank=='BIC'){
          factor.num[i,,1+iiter]=tenFM.bic(ans.iter$lambda[[1]]/sigmas.hat,h0,ddd[i],ddd[-i],n,delta1)
        }else if(rank=='ER'){
          factor.num[i,,1+iiter]=tenFM.ratio(ans.iter$lambda[[1]],ddd[i],ddd[-i],n)
        }else{
          stop('Wrong rank method !')  # should be error
        }
        if(inputr==FALSE){
          r[i]=factor.num[i,penalty,1+iiter]
        }

        eigen.gap[i,1:dd[i],1+iiter]=ans.iter$lambda[[1]]
      }
      x.hat <- get.hat(x.tnsr,ans.Q,d.seq)
      fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
      dis <- abs(fnorm.resid[iiter+1] - fnorm.resid[iiter])
      if(iiter==1){
        Qfirst <- ans.Q
        x.hat.first <- x.hat@data
      }
      iiter <- iiter + 1
    }
  }else{
    iiter <- iiter + 1
  }

  factor.num[,,maxiter]=factor.num[,,iiter-1]
  eigen.gap[,,maxiter]=eigen.gap[,,iiter-1]
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]
  norm.percent <- c(ans.init$norm.percent,fnorm.resid[1],tail(fnorm.resid,1))
  return(factor.num[,penalty,maxiter])
}


#' BIC estimators for determing the numbers of factors
#'
#'@name tenFM.bic
#'@rdname tenFM.bic
#'@aliases tenFM.bic
#'@export
#'@param reigen list of eigenvalues
#'@param h0 lag number
#'@param p1 p1
#'@param p2 p2
#'@param n n
#'@return factor.p1: Estimated number of factors
tenFM.bic<-function(reigen,h0=1,p1,p2,n,delta1=0){
  #p1
  #p2
  #n
  delta=2*delta1
  if(length(p2)>1){
    p2=prod(p2)
  }
  factor.p1=numeric(5)
  p=p1*p2
  m1=ceiling(p1/3)

  lambda.p1<-reigen[p1:1]
  cumlambda.p1<-cumsum(lambda.p1)
  cumlambda.p1<-cumlambda.p1[(p1-1):1]

  #threshold
  ic=cumlambda.p1[1:m1]/p^2*p^delta+(1:m1)*h0*(1/n)*log(p*n/(p+n))
  factor.p1[1]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2*p^delta+(1:m1)*h0*(1/n)*log(min(p,n))
  factor.p1[2]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2*p^delta+(1:m1)*h0*(1/n+1/p)*log(p*n/(p+n))
  factor.p1[3]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2*p^delta+(1:m1)*h0*(1/n+1/p)*log(min(p,n))
  factor.p1[4]<-which.min(ic)
  ic=cumlambda.p1[1:m1]/p^2*p^delta+(1:m1)*h0*(1/n+1/p)*log(min(p1,n))
  factor.p1[5]<-which.min(ic)

  factor.p1
}


#' Eigenvalue ratio estimators for determining the numbers of factors
#'
#'@name tenFM.ratio
#'@rdname tenFM.ratio
#'@aliases tenFM.ratio
#'@export
#'@param reigen list of eigenvalues
#'@param h0 lag number
#'@param p1 p1
#'@param p2 p2
#'@param n n
#'@return factor.p1: Estimated number of factors
#'
tenFM.ratio<-function(reigen,p1,p2,n){
  #p1
  #p2
  #n
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
  ratio<-lambda.p1[(p1-1):(p1-m1)]/lambda.p1[p1:(p1-m1+1)]
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


#' get fitted x.hat(signal part)
#'
#' get x.hat from a tensor x and all loading matrices in list Q by projection
#'@name get.hat
#'@rdname get.hat
#'@aliases get.hat
#'@export
#'@param x tensor of any dimension, x: \eqn{d1 * d2 * d3 * \cdots * d_K * T}
#'@param Q list of k matrices
#'@param d.seq 1:k
#'@return a tensor object x.hat
get.hat <- function(x,Q,d.seq){
  ans <- ttl(x,lapply(Q,t),d.seq)
  ans <- ttl(ans,Q,d.seq)
  ans
}


#' get M matrix and first r eigenvectors of M
#'
#'@name vts
#'@rdname vts
#'@aliases vts
#'@export
#'@param x a n*d matrix
#'@param h0 lag number
#'@param r First r eigenvectors
#'@return a list containing the following:\describe{
#'\item{\code{M}}{.}
#'\item{\code{Q}}{The first r eigenvectors of matrix M}
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


#' TIPUP one step for any Dim tensor time series
#'
#'@name tipup.init.tensor
#'@rdname tipup.init.tensor
#'@aliases tipup.init.tensor
#'@export
#'@param x tensor of any dimension \eqn{d1 * d2 * d3 * \cdots * d_K * T}
#'@param r input rank of factor
#'@param h0 lag number
#'@param oneside.true  If TRUE, only compute the one sided column space, not the other sides, this option is useful for the iterative method
#'@param norm.true If TRUE, calculate the normalized residual of the tensor
#'@return a list containing the following:\describe{
#'\item{\code{M}}{Estimator}
#'\item{\code{Q}}{Orthonormal matrix Q}
#'\item{\code{lambda}}{singular values}
#'\item{\code{norm.percent}}{Frobenius norm of residual}
#'\item{\code{x.hat}}{Fitted signal part, a tensor object}
#'}
tipup.init.tensor <- function(x,r,h0=1,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # TIPUP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_K * n
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


#' TOPUP one step for any Dim tensor time series
#'
#'@name topup.init.tensor
#'@rdname topup.init.tensor
#'@aliases topup.init.tensor
#'@export
#'@param x tensor of any dimension \eqn{d1 * d2 * d3 * \cdots * d_K * T}
#'@param r input rank of factor
#'@param h0 lag number
#'@param oneside.true  If TRUE, only compute the one sided column space, not the other sides, this option is useful for the iterative method
#'@param norm.true If TRUE, calculate the normalized residual of the tensor
#'@return a list containing the following:\describe{
#'\item{\code{M}}{Estimator}
#'\item{\code{Q}}{Orthonormal matrix Q}
#'\item{\code{lambda}}{singular values}
#'\item{\code{norm.percent}}{Frobenius norm of residual}
#'\item{\code{x.hat}}{Fitted signal part, a tensor object}
#'}
topup.init.tensor <- function(x,r,h0=1,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # TOPUP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_K * T
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


#' generate orthonormal loading Q
#'
#' For simulation only, generate loading matrix Q with size p*r
#'@name gen.Q
#'@rdname gen.Q
#'@aliases gen.Q
#'@export
#'@param p row dimension
#'@param r column dimension, should not be larger than p
#'@return An orthonormal matrix Q with size p*r
gen.Q <- function(p,r){
  A <- matrix(rnorm(p*r),p,r)
  Q <- qr.Q(qr(A))
  Q
}

#' generate matrix \eqn{F_t}
#'
#' For simulation only, generate a matrix F with AR coefficients of length n
#'@name gen.F
#'@rdname gen.F
#'@aliases gen.F
#'@export
#'@param ar1.coef r1*r2
#'@param n length of time
#'@param innovsd sd of factor ar series
#'@return An array of size r1*r2*n
gen.F <- function(ar1.coef,n,innovsd){
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

#' generate tensor \eqn{F_t}
#'
#' For simulation only, generate a tensor F with AR coefficients of length n
#'@name gen.F.tensor3
#'@rdname gen.F.tensor3
#'@aliases gen.F.tensor3
#'@export
#'@param ar1.coef r1*r2*r3
#'@param n length of time
#'@param innovsd sd of factor ar series
#'@return An array of size r1*r2*r3*n
gen.F.tensor3 <- function(ar1.coef,n,innovsd){
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


#' Generate Error \eqn{E_t} for Matrix
#'
#' For simulation only, generate the error tensor E with given two covariance matrices
#'@name gen.E
#'@rdname gen.E
#'@aliases gen.E
#'@export
#'@param n length of time
#'@param p1 p1
#'@param p2 p2
#'@param cov1.half p1 * p1 covariance matrix
#'@param cov2.half p2 * p2 covariance matrix
#'@return A tensor of size p1 * p2 * n
gen.E <- function(n,p1,p2,cov1.half,cov2.half){
  E <- array(rnorm(n*p1*p2),c(p1,p2,n)) # p1 * p2 * n
  E <- tensor(E,cov1.half,1,1) # p2 * n * p1
  E <- tensor(E,cov2.half,1,1) # n * p1 * p2
  aperm(E,c(2,3,1)) #p1 * p2 * n
}

#' generate error \eqn{E_t} for tensor
#'
#' For simulation only, generate the error tensor E with given three covariance matrices
#'@name gen.E.tensor3
#'@rdname gen.E.tensor3
#'@aliases gen.E.tensor3
#'@export
#'@param n length of time
#'@param p1 p1
#'@param p2 p2
#'@param p3 p3
#'@param cov1.half p1 * p1 covariance matrix
#'@param cov2.half p2 * p2 covariance matrix
#'@param cov3.half p3 * p3 covariance matrix
#'@return A tensor of size p1 * p2 * p3 * n
gen.E.tensor3 <- function(n,p1,p2,p3,cov1.half,cov2.half,cov3.half){
  E <- array(rnorm(n*p1*p2*p3),c(p1,p2,p3,n)) # p1 * p2 * p3 * n
  E <- tensor(E,cov1.half,1,1) # p2 * p3 * n * p1
  E <- tensor(E,cov2.half,1,1) # p3 * n * p1 * p2
  E <- tensor(E,cov3.half,1,1) # n * p1 * p2 * p3
  aperm(E,c(2:4,1)) #p1 * p2 * n
}


#' generate X for matrix
#'
#' For simulation only, generate \eqn{X = bF \times_1 Q_1 \times_2 Q_2}
#'@name gen.X
#'@rdname gen.X
#'@aliases gen.X
#'@export
#'@param Q1 left loading matrix, p1*r1 orthonormal
#'@param Q2 right loading matrix, p2*r2 orthonormal
#'@param bF factor series, r1*r2*n
#'@return A tensor of size p1*p2*n
gen.X <- function(Q1,Q2,bF){
  X <- tensor(bF, Q1, 1, 2) # r2*n*p1
  X <- tensor(X, Q2, 1, 2) # n*p1*p2
  X <- aperm(X,c(2,3,1)) # p1*p2*n
  X
}


#' generate X for tensor
#'
#' For simulation only, generate \eqn{X = bF \times_1 Q_1 \times_2 Q_2 \times_3 Q_3}
#'@name gen.X.tensor3
#'@rdname gen.X.tensor3
#'@aliases gen.X.tensor3
#'@export
#'@param Q1 loading matrix, p1*r1 orthonormal
#'@param Q2 loading matrix, p2*r2 orthonormal
#'@param Q3 loading matrix, p3*r3 orthonormal
#'@param bF factor series, r1*r2*r3*n
#'@return A tensor of size p1*p2*p3*n
gen.X.tensor3 <- function(Q1,Q2,Q3,bF){
  X <- tensor(bF, Q1, 1, 2) # r2*r3*n*p1
  X <- tensor(X, Q2, 1, 2) # r3*n*p1*p2
  X <- tensor(X, Q3, 1, 2) # n*p1*p2*p3
  X <- aperm(X,c(2:4,1)) # p1*p2*p3*n
  X
}

