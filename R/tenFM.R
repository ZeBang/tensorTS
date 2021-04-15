#' Tensor factor estimation method
#'
#' Estimation function for tensor-valued time series factor model.
#'@name tenFM.est
#'@rdname tenFM.est
#'@aliases tenFM.est
#'@usage tenFM.est(x,r,method='TIPUP')
#'@export
#'@param x \eqn{T \times d_1 \times \cdots \times d_K} dimensional tensor-valued time series
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
#'@return return a list containing the following:\describe{
#'\item{\code{Ft}}{estimated factor matrices of dimension \eqn{T*r_1*r_2 *\cdots *r_k}}
#'\item{\code{Ft.all}}{Summation of factor matrices over time, dimension \eqn{r_1,r_2,\cdots,r_k} }
#'\item{\code{Q}}{loading matrices}
#'\item{\code{x.hat}}{fitted signal part}
#'\item{\code{niter}}{number of iterations}
#'\item{\code{fnorm.resid}}{Frobenius norm of residual, divide the Frobenius norm of the original tensor}
#'}
#'@examples
#'dims <- c(4,5,6) # dimensions of tensor time series
#'t <- 50
#'r <- c(2,2,2)  # dimensions of factor series
#'x = generate.FM(dims,array(c(runif(prod(r)),min=0.5,max=1),r),t)  # generate t*dims tensor
#'Et <- array(rnorm(t*prod(dims)),c(t,dims)) # add some iid noise to the tensor time series
#'x = x + Et
#'result <- tenFM.est(x,r,method='TOPUP')  # Estimation
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
  #x.hat <- get.hat(x.tnsr,ans.init$Q,d.seq)

  x.hat <- ttl(x.tnsr,lapply(ans.init$Q,t),d.seq)
  x.hat <- ttl(x.hat,ans.init$Q,d.seq)

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

        #x.hat <- get.hat(x.tnsr,ans.Q,d.seq)
        x.hat <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
        x.hat <- ttl(x.hat,ans.Q,d.seq)

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


#' Tensor factor model rank determination
#'
#' Function for rank determination of the tensor factor model
#'@name tenFM.rank
#'@rdname tenFM.rank
#'@aliases tenFM.rank
#'@usage tenFM.rank(x,r,rank='BIC',method='TIPUP')
#'@export
#'@param x \eqn{T \times d_1 \times \cdots \times d_K} dimensional tensor-valued time series
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
#'@param penalty takes value in {1,2,3,4,5}, decide which penalty function to use. \describe{
#'  \item{}{When \code{rank}= 'BIC':}
#'  \item{}{if \code{penalty}=1, \eqn{g_1= \frac{h_0 d^{2-2\nu}}{T}\log(\frac{dT}{d+T})}}
#'  \item{}{if \code{penalty}=2, \eqn{g_2= h_0 d^{2-2\nu}(\frac{1}{T}+\frac{1}{d})\log(\frac{dT}{d+T})}}
#'  \item{}{if \code{penalty}=3, \eqn{g_3= \frac{h_0 d^{2-2\nu}}{T} \log(\min{d,T})}}
#'  \item{}{if \code{penalty}=4, \eqn{g_4= h_0 d^{2-2\nu}(\frac{1}{T}+\frac{1}{d})\log(\min{d,T})}}
#'  \item{}{if \code{penalty}=5, \eqn{g_5= h_0 d^{2-2\nu}(\frac{1}{T}+\frac{1}{d})\log(\min{d_k,T})}}
#'  \item{}{When \code{rank}= 'ER':}
#'  \item{}{if \code{penalty}=1, \eqn{h_1= c_0 h_0}}
#'  \item{}{if \code{penalty}=2, \eqn{h_2= \frac{h_0 d^2}{T^2}}}
#'  \item{}{if \code{penalty}=3, \eqn{h_3= \frac{h_0 d^2}{T^2 d_k^2}}}
#'  \item{}{if \code{penalty}=4, \eqn{h_4= \frac{h_0 d^2}{T^2 d_k^2} + \frac{h_0 d_k^2}{T^2}}}
#'  \item{}{if \code{penalty}=5, \eqn{h_5= \frac{h_0 d^2}{T^2 d_k^2} + \frac{h_0 dd_k^2}{T^2}}}
#'}
#'@param delta1 weakest factor strength, a tuning parameter used for BIC method only
#'@param tol tolerance for each iteration
#'@param maxiter integer, the max number of iterations allowed
#'@return estimated number of factors \eqn{n_1,n_2,\cdots,n_k}
#'@examples
#'dims <- c(4,5,6) # dimensions of tensor time series
#'t <- 50
#'r <- c(2,2,2)  # dimensions of factor series
#'x = generate.FM(dims,array(c(runif(prod(r)),min=0.5,max=1),r),t)  # generate t*dims tensor
#'Et <- array(rnorm(t*prod(dims)),c(t,dims))  # add some iid noise to the tensor time series
#'x = x + Et
#'rank <- tenFM.rank(x,r,rank='BIC',method='TIPUP')  # Estimate the rank
tenFM.rank = function(x,r,h0=1,rank='BIC',method='TIPUP',inputr=FALSE,iter=TRUE,penalty=1,delta1=0,tol=1e-4,maxiter=100){
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
      factor.num[i,,1]=tensor.bic(ans.init$lambda[[i]]/sigmas.hat,h0,ddd[i],ddd[-i],n,delta1)
    }else if(rank=='ER'){
      factor.num[i,,1]=tensor.ratio(ans.init$lambda[[i]],ddd[i],ddd[-i],n)
    }else{
      stop('Wrong rank method !')  # should catch exception
    }
    eigen.gap[i,1:dd[i],1]=ans.init$lambda[[i]]
  }
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter)
  #x.hat <- get.hat(x.tnsr,ans.init$Q,d.seq)
  x.hat <- ttl(x.tnsr,lapply(ans.init$Q,t),d.seq)
  x.hat <- ttl(x.hat,ans.init$Q,d.seq)

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
          factor.num[i,,1+iiter]=tensor.bic(ans.iter$lambda[[1]]/sigmas.hat,h0,ddd[i],ddd[-i],n,delta1)
        }else if(rank=='ER'){
          factor.num[i,,1+iiter]=tensor.ratio(ans.iter$lambda[[1]],ddd[i],ddd[-i],n)
        }else{
          stop('Wrong rank method !')  # should be error
        }
        if(inputr==FALSE){
          r[i]=factor.num[i,penalty,1+iiter]
        }

        eigen.gap[i,1:dd[i],1+iiter]=ans.iter$lambda[[1]]
      }
      #x.hat <- get.hat(x.tnsr,ans.Q,d.seq)
      x.hat <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
      x.hat <- ttl(x.hat,ans.Q,d.seq)

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


#' Generate tensor time series from given factor structure
#'
#' For simulation only, generate tensor time series \eqn{X_t} from given factor structure, each factor entry follows an univariate AR(1) process
#'@name generate.FM
#'@rdname generate.FM
#'@aliases generate.FM
#'@export
#'@param dims dimensions of the output tensor time series,  \eqn{d_1\times d_2\cdots\times d_K}
#'@param ar1.coef ar coefficients of the factor process, same dimension with \eqn{F_t : r_1\times r_2\cdots\times r_K}
#'@param t length of time
#'@return A tensor of size  \eqn{T\times d_1\times d_2\cdots\times d_K}
#'@seealso \code{\link{generate.AR}}
#'@examples
#'dims <- c(4,5,6) # dimensions of tensor time series
#'t <- 50
#'r <- c(2,2,2)  # dimensions of factor series
#'ar1.coef = array(c(runif(prod(r)),min=0.5,max=1)),r)
#'x = generate.FM(dims,ar1.coef,t)  # generate t*dims tensor
generate.FM <- function(dims,ar1.coef,t){
  F.dims = dim(ar1.coef)
  Ft = array(NA,dim=c(F.dims[1],F.dims[2],F.dims[3],t))
  for(ir1 in 1:F.dims[1]){
    for(ir2 in 1:F.dims[2]){
      for(ir3 in 1:F.dims[3]){
        Ft[ir1,ir2,ir3,] = arima.sim(n=t, model=list(ar=ar1.coef[ir1,ir2,ir3]))
      }
    }
  }

  dd = length(dims)
  X <- Ft
  for(i in 1:dd){
    A = matrix(rnorm(dims[i]*r[i]),dims[i],r[i])
    Q = qr.Q(qr(A))
    X <- tensor(X,Q,1,2)
  }
  return(as.tensor(X))
}

tensor.bic<-function(reigen,h0=1,p1,p2,n,delta1=0){
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

tensor.ratio<-function(reigen,p1,p2,n){
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
    #x.hat <- get.hat(x.tnsr,ans.Q,1:k)
    x.hat <- ttl(x.tnsr,lapply(ans.Q,t),1:k)
    x.hat <- ttl(x.hat,ans.Q,1:k)
    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("M"=ans.M,"Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}

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
    #x.hat <- get.hat(x.tnsr,ans.Q,1:(d-1))
    x.hat <- ttl(x.tnsr,lapply(ans.Q,t),1:(d-1))
    x.hat <- ttl(x.hat,ans.Q,1:(d-1))

    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("M"=ans.M,"Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}
