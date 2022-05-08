#' Estimation for Tucker structure Factor Models of Tensor-Valued Time Series
#'
#' Estimation function for Tucker structure factor models of tensor-valued time series.
#' Two unfolding methods of the auto-covariance tensor, Time series Outer-Product Unfolding Procedure (TOPUP), Time series Inner-Product Unfolding Procedure (TIPUP),
#' are included, as determined by the value of \code{method}.
#'@details
#' Tensor factor model with Tucker structure has the following form,
#' \deqn{X_t = F_t \times_{1} A_1 \times_{2} \cdots \times_{K} A_k + E_t,}
#' where \eqn{A_k} is the deterministic loading matrix of size \eqn{d_k \times r_k} and \eqn{r_k \ll d_k},
#' the core tensor \eqn{F_t} itself is a latent tensor factor process of dimension \eqn{r_1 \times \cdots \times r_K},
#' and the idiosyncratic noise tensor \eqn{E_t} is uncorrelated (white) across time. Two estimation approaches, named TOPUP and TIPUP, are studied.
#' Time series Outer-Product Unfolding Procedure (TOPUP) are based on
#' \deqn{ {\rm{TOPUP}}_{k}(X_{1:T}) = \left(\sum_{t=h+1}^T \frac{{\rm{mat}}_{k}( X_{t-h}) \otimes {\rm{mat}}_k(X_t)} {T-h}, \ h=1,...,h_0 \right),}
#' where \eqn{h_0} is a predetermined positive integer, \eqn{\otimes} is tensor product. Note that
#' \eqn{ {\rm{TOPUP}}_k(\cdot)} is a function mapping a tensor time series to an order-5 tensor.
#' Time series Inner-Product Unfolding Procedure (TIPUP) replaces the tensor product in TOPUP with the inner product:
#' \deqn{ {\rm{TIPUP}}_k(X_{1:T})={\rm{mat}}_1\left(\sum_{t=h+1}^T \frac{{\rm{mat}}_k(X_{t-h}) {\rm{mat}}_k^\top(X_t)} {T-h}, \ h=1,...,h_0 \right).}
#'@name tenFM.est
#'@rdname tenFM.est
#'@aliases tenFM.est
#'@usage tenFM.est(x,r,h0=1,method='TIPUP',iter=TRUE,tol=1e-4,maxiter=100)
#'@export
#'@importFrom stats varimax
#'@param x \eqn{T \times d_1 \times \cdots \times d_K} tensor-valued time series.
#'@param r input rank of factor tensor.
#'@param h0 the number of lags used in auto-covariance tensor.
#'@param method character string, specifying the type of the estimation method to be used. \describe{
#'  \item{\code{"TIPUP",}}{TIPUP method.}
#'  \item{\code{"TOPUP",}}{TOPUP method.}
#'}
#'@param iter boolean, specifying using an iterative approach or an non-iterative approach.
#'@param tol tolerance in terms of the Frobenius norm.
#'@param maxiter maximum number of iterations if error stays above \code{tol}.
#'@return returns a list containing the following:\describe{
#'\item{\code{Ft}}{estimated factor processes of dimension \eqn{T \times r_1 \times r_2 \times \cdots \times r_k}.}
#'\item{\code{Ft.all}}{Summation of factor processes over time, of dimension \eqn{r_1,r_2,\cdots,r_k}.}
#'\item{\code{Q}}{a list of estimated factor loading matrices \eqn{Q_1,Q_2,\cdots,Q_K}. }
#'\item{\code{x.hat}}{fitted signal tensor, of dimension \eqn{T \times d_1 \times d_2 \times \cdots \times d_k}.}
#'\item{\code{niter}}{number of iterations.}
#'\item{\code{fnorm.resid}}{Frobenius norm of residuals, divide the Frobenius norm of the original tensor.}
#'}
#'@references
#'Chen, Rong, Dan Yang, and Cun-Hui Zhang. "Factor models for high-dimensional tensor time series." Journal of the American Statistical Association (2021): 1-59.
#'
#'Han, Yuefeng, Rong Chen, Dan Yang, and Cun-Hui Zhang. "Tensor factor model estimation by iterative projection." arXiv preprint arXiv:2006.02611 (2020).
#'
#'@examples
#'set.seed(333)
#'dims <- c(16,18,20) # dimensions of tensor time series
#'r <- c(3,3,3)  # dimensions of factor series
#'Ft <- tenAR.sim(t=100, dim=r, R=1, P=1, rho=0.9, cov='iid')
#'lambda <- sqrt(prod(dims))
#'x <- tenFM.sim(Ft,dims=dims,lambda=lambda,A=NULL,cov='iid') # generate t*dims tensor time series
#'result <- tenFM.est(x,r,h0=1,iter=TRUE,method='TIPUP')  # Estimation
#'Ft <- result$Ft
tenFM.est=function(x,r,h0=1,method='TIPUP',iter=TRUE,tol=1e-4,maxiter=100){
  x <- aperm(x,c(2:length(dim(x)),1))
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]

  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  if(method=="TIPUP"){
    ans.init <- tipup.init.tensor(x,r,h0,norm.true=TRUE)
  }else if(method=="TOPUP"){
    ans.init <- topup.init.tensor(x,r,h0,norm.true=TRUE)
  }else{
    stop('Wrong method !')
  }
  ddd=dd[-d]
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter+2)

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
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]
  fnorm.resid <- fnorm.resid^2
  x0 <- matrix(x,prod(dd[-d]))
  x0 <- t(scale(t(x0),scale=FALSE) )
  x0 <- array(x0,dd)
  return(list("Ft"=aperm(Ft@data,c(d,1:(d-1))),"Ft.all"=Ft.all,"Q"=ans.Q,"x.hat"=aperm(x.hat@data,c(d,1:(d-1))),"niter"=iiter,"fnorm.resid"=fnorm.resid[iiter]))
}

#' Rank Determination for Tensor Factor Models with Tucker Structure
#'
#' Function for rank determination of tensor factor models with Tucker Structure.
#' Two unfolding methods of the auto-covariance tensor, Time series Outer-Product Unfolding Procedure (TOPUP), Time series Inner-Product Unfolding Procedure (TIPUP),
#' are included, as determined by the value of \code{method}.
#' Different penalty functions for the information criterion (IC) and the eigen ratio criterion (ER) can be used,
#' which should be specified by the value of \code{rank} and \code{penalty}. The information criterion resembles BIC in the vector factor model literature.
#' And the eigen ratio criterion is similar to the eigenvalue ratio based methods in the vector factor model literature.
#'@details
#' Let \eqn{W} be a \eqn{p\times p} symmetric and non-negative definite matrix and \eqn{\widehat{W}} be its sample version, \eqn{{\hat\lambda}_j} be the eigenvalues of \eqn{\widehat{W}}
#' such that \eqn{{\hat\lambda}_1\geq {\hat\lambda}_2 \geq \cdots \hat{\lambda}_p}.
#' The rank determination methods using the information criterion ("IC") and the eigen ratio criterion ("ER") are defined as follows:
#' \deqn{IC(\widehat{W}) = \mathrm{argmin}_{0\leq m \leq m^{*}} \left\{ \sum_{j=m+1}^{p} {\hat\lambda}_j + mg(\widehat{W}) \right\},}
#' \deqn{ER(\widehat{W}) = \mathrm{argmin}_{0\leq m \leq m^{*}} \left\{ \frac{{\hat\lambda}_{m+1}+h(\widehat{W})}{ {\hat\lambda}_m +h(\widehat{W})} \right\},}
#' where \eqn{m^{*}} is a predefined upper bound, \eqn{g} and \eqn{h} are some appropriate positive penalty functions. We have provided 5 choices for \eqn{g} and \eqn{h};
#' see more details in the argument "\code{penalty}".
#' For non-iterative TOPUP and TIPUP methods, \eqn{\widehat{W}} is
#' \eqn{ {\rm mat}_1({\rm{TOPUP}}_{k}(X_{1:T})) {\rm mat}_1({\rm{TOPUP}}_{k}(X_{1:T}))^\top } or
#' \eqn{ ({\rm{TIPUP}}_{k}(X_{1:T})) ({\rm{TIPUP}}_{k}(X_{1:T}))^\top }, for each tensor mode \eqn{k}, \eqn{1\leq k \leq K},
#' where \eqn{{\rm{TOPUP}}_{k}(X_{1:T})} and \eqn{{\rm{TIPUP}}_{k}(X_{1:T})} are defined in the Details section of the function \code{\link{tenFM.est}}.
#' For iterative TOPUP and TIPUP methods, we refer to the literature in the References section for more information.
#'@name tenFM.rank
#'@rdname tenFM.rank
#'@aliases tenFM.rank
#'@usage tenFM.rank(x,r,h0=1,rank='IC',method='TIPUP',inputr=FALSE,iter=TRUE,penalty=1,
#'delta1=0,tol=1e-4,maxiter=100)
#'@export
#'@param x \eqn{T \times d_1 \times \cdots \times d_K} tensor-valued time series.
#'@param r initial guess of the rank of factor tensor.
#'@param h0 the number of lags used in auto-covariance tensor.
#'@param rank character string, specifying the type of the rank determination method to be used. \describe{
#'  \item{\code{"IC",}}{information criterion.}
#'  \item{\code{"ER",}}{eigen ratio criterion.}
#'}
#'@param method character string, specifying the type of the factor estimation method to be used. \describe{
#'  \item{\code{"TIPUP",}}{TIPUP method.}
#'  \item{\code{"TOPUP",}}{TOPUP method.}
#'}
#'@param inputr boolean, if TRUE, always use initial guess rank r in each iteration; if FLASE, the rank will be updated in each iteration.
#'@param iter boolean, specifying using an iterative approach or a non-iterative approach.
#'@param penalty takes value in {1,2,3,4,5}, decides which penalty function to use for each tesnor mode \eqn{k}. Here \eqn{\nu} is a tuning parameter defined in the argument "\code{delta1}", and \eqn{d=\prod_{i=1}^{K} d_k }.
#'  \describe{
#'  \item{}{When \code{rank}= '\code{IC}':}
#'  \item{}{if \code{penalty}=1, \eqn{g_1= \frac{h_0 d^{2-2\nu}}{T}\log(\frac{dT}{d+T})};}
#'  \item{}{if \code{penalty}=2, \eqn{g_2= h_0 d^{2-2\nu}(\frac{1}{T}+\frac{1}{d})\log(\frac{dT}{d+T})};}
#'  \item{}{if \code{penalty}=3, \eqn{g_3= \frac{h_0 d^{2-2\nu}}{T} \log(\min{(d,T)})};}
#'  \item{}{if \code{penalty}=4, \eqn{g_4= h_0 d^{2-2\nu}(\frac{1}{T}+\frac{1}{d})\log(\min{(d,T)})};}
#'  \item{}{if \code{penalty}=5, \eqn{g_5= h_0 d^{2-2\nu}(\frac{1}{T}+\frac{1}{d})\log(\min{(d_k,T)})}.}
#'  \item{}{When \code{rank}= '\code{ER}':}
#'  \item{}{if \code{penalty}=1, \eqn{h_1= c_0 h_0};}
#'  \item{}{if \code{penalty}=2, \eqn{h_2= \frac{h_0 d^2}{T^2}};}
#'  \item{}{if \code{penalty}=3, \eqn{h_3= \frac{h_0 d^2}{T^2 d_k^2}};}
#'  \item{}{if \code{penalty}=4, \eqn{h_4= \frac{h_0 d^2}{T^2 d_k^2} + \frac{h_0 d_k^2}{T^2}};}
#'  \item{}{if \code{penalty}=5, \eqn{h_5= \frac{h_0 d^2}{T^2 d_k^2} + \frac{h_0 dd_k^2}{T^2}}.}
#'}
#'@param delta1 weakest factor strength, a tuning parameter used for IC method only
#'@param tol tolerance in terms of the Frobenius norm.
#'@param maxiter maximum number of iterations if error stays above \code{tol}.
#'@return return a list containing the following:\describe{
#'\item{\code{path}}{a \eqn{K \times (\rm{niter}+1)} matrix of the estimated Tucker rank of the factor process as a path of the maximum number of iteration (\eqn{\rm{niter}}) used. The first row is the estimated rank under non-iterative approach, the \eqn{i+1}-th row is the estimated rank \eqn{\hat r_1, \hat r_2, \cdots, \hat r_K} at \eqn{(i)}-th iteration.}
#'\item{\code{factor.num}}{final solution of the estimated Tucker rank of the factor process \eqn{\hat r_1, \hat r_2, \cdots, \hat r_K}.}
#'}
#'@references
#'Han, Yuefeng, Cun-Hui Zhang, and Rong Chen. "Rank Determination in Tensor Factor Model." Available at SSRN 3730305 (2020).
#'@examples
#'set.seed(333)
#'dims <- c(16,18,20) # dimensions of tensor time series
#'r <- c(3,3,3)  # dimensions of factor series
#'Ft <- tenAR.sim(t=100, dim=r, R=1, P=1, rho=0.9, cov='iid')
#'lambda <- sqrt(prod(dims))
#'x <- tenFM.sim(Ft,dims=dims,lambda=lambda,A=NULL,cov='iid') # generate t*dims tensor time series
#'rank <- tenFM.rank(x,r=c(4,4,4),h0=1,rank='IC',iter=TRUE,method='TIPUP')  # Estimate the rank
tenFM.rank = function(x,r=NULL,h0=1,rank='IC',method='TIPUP',inputr=FALSE,iter=TRUE,penalty=1,delta1=0,tol=1e-4,maxiter=100){
  x <- aperm(x,c(2:length(dim(x)),1))
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]
  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  factor.num <- array(NA, c(d-1,5,maxiter+1))

  if(is.null(r)){
    r = rep(1,d-1)
  }
  
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
    if(rank=='BIC'|rank=='IC'){
      sigmas.hat=1
      factor.num[i,,1]=tensor.bic(ans.init$lambda[[i]]/sigmas.hat,h0,ddd[i],ddd[-i],n,delta1)
    }else if(rank=='ER'){
      factor.num[i,,1]=tensor.ratio(ans.init$lambda[[i]],ddd[i],ddd[-i],n)
    }else{
      stop('Wrong rank method !')
    }
  }
  r <- factor.num[,penalty,1]
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter+1)
  x.hat <- ttl(x.tnsr,lapply(ans.init$Q,t),d.seq)
  x.hat <- ttl(x.hat,ans.init$Q,d.seq)

  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q

  if(iter==TRUE){
    for(i in 1:d){r[i] = min(dd[i],r[i]+1)}
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
          stop('Wrong estimation method input !')
        }
        ddd=dd[-d]

        if(rank=='BIC'|rank=='IC'){
          factor.num[i,,1+iiter]=tensor.bic(ans.iter$lambda[[1]]/sigmas.hat,h0,ddd[i],ddd[-i],n,delta1)
        }else if(rank=='ER'){
          factor.num[i,,1+iiter]=tensor.ratio(ans.iter$lambda[[1]],ddd[i],ddd[-i],n)
        }else{
          stop('Wrong rank method !')
        }
        if(inputr==FALSE){
          r[i]=factor.num[i,penalty,1+iiter]
        }

      }
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

  # factor.num[,,maxiter]=factor.num[,,iiter]
  factor.num[,,maxiter]=factor.num[,,iiter-1]
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]

  # label the factor number path
  path = t(factor.num[,penalty,1:(iiter)])
  path.rowname = c()
  for(ii in 1:iiter){path.rowname <- c(path.rowname,paste('iteration ',ii-1,sep=''))}
  path.colname = c()
  for(ii in 1:(d-1)){path.colname <- c(path.colname,paste('mode ',ii,sep=''))}
  rownames(path)=path.rowname
  colnames(path)=path.colname
  
  # return(list("path"=t(factor.num[,penalty,1:(iiter)]),"factor.num"=factor.num[,penalty,maxiter]))
  return(list("path"=path,"factor.num"=factor.num[,penalty,maxiter]))
}

#' Generate Tensor Time series using given Factor Process and Factor Loading Matrices
#'
#' Simulate tensor time series \eqn{X_t} using a given factor process \eqn{F_t}. The factor process \eqn{F_t} can be generated by the function \code{\link{tenAR.sim}}.
#'@details
#' Simulate from the model :
#' \deqn{X_t = \lambda F_t \times_{1} A_1 \times_{2} \cdots \times_{K} A_k + E_t,}
#' where \eqn{A_k} is the deterministic loading matrix of size \eqn{d_k \times r_k} and \eqn{r_k \ll d_k},
#' the core tensor \eqn{F_t} itself is a latent tensor factor process of dimension \eqn{r_1 \times \cdots \times r_K},
#' \eqn{\lambda} is an additional signal strength parameter,
#' and the idiosyncratic noise tensor \eqn{E_t} is uncorrelated (white) across time. In this function, by default \eqn{A_k} are orthogonal matrices.
#'@name tenFM.sim
#'@rdname tenFM.sim
#'@aliases tenFM.sim
#'@usage tenFM.sim(Ft,dims=NULL,lambda=1,A=NULL,cov='iid',rho=0.2)
#'@export
#'@importFrom stats rnorm
#'@param Ft input of the factor process, of dimension \eqn{T \times r_1 \times r_2 \times \cdots \times r_k}. It can be TenAR(p) tensor time series generated by the function \link{tenAR.sim}.
#'@param dims dimensions of the output tensor at each time,  \eqn{d_1\times d_2\cdots\times d_K}.
#'@param A a list of the factor loading matrices \eqn{A_1, A_2, \cdots, A_K}. The default is random orthogonal matrices \eqn{A_k} of dimension \eqn{d_k \times r_k}.
#'@param lambda signal strength parameter of the tensor factor models, see Details section for more information.
#'@param cov covariance matrix of the error tensor: identity ("iid"), separable Kronecker structure ("separable"), random ("random").
#'@param rho a parameter only for "separable" covariance matrix of the error tensor. It is the off-diagonal element of the error matrices, with the diagonal being 1.
#'@return A tensor-valued time series of dimension \eqn{T\times d_1\times d_2\cdots\times d_K}.
#'@seealso \code{\link{tenAR.sim}}
#'@examples
#'set.seed(333)
#'dims <- c(16,18,20) # dimensions of tensor time series
#'r <- c(3,3,3)  # dimensions of factor series
#'Ft <- tenAR.sim(t=100, dim=r, R=1, P=1, rho=0.9, cov='iid')
#'lambda <- sqrt(prod(dims))
#'# generate t*dims tensor time series with iid error covaraince structure
#'x <- tenFM.sim(Ft,dims=dims,lambda=lambda,A=NULL,cov='iid')
#'# generate t*dims tensor time series with separable error covaraince structure
#'x <- tenFM.sim(Ft,dims=dims,lambda=lambda,A=NULL,cov='separable',rho=0.2)
tenFM.sim <- function(Ft,dims=NULL,lambda=1,A=NULL,cov='iid',rho=0.2){
  r <- dim(Ft)[-1] #dimensions of the factor process
  t <- dim(Ft)[1] #length of output series
  dd <- length(dims)
  if(length(r)!=dd & is.null(A)){
    stop("Incorrect length K of input dims or A, for order K tensor time series.")
  }
  if(is.null(A)){
    X <- aperm(Ft,c(2:(dd+1),1))
    for(i in 1:dd){
      Ai = matrix(rnorm(dims[i]*r[i]),dims[i],r[i])
      Q = qr.Q(qr(Ai))
      X <- tensor(X,Q,1,2)
    }
  } else{
    X <- aperm(Ft,c(2:(dd+1),1))
    for(i in 1:dd){
      Ai = A[[i]]
      Q = qr.Q(qr(Ai))
      X <- tensor(X,Q,1,2)
    }
  }
  if (cov == "iid"){
    E <- array(rnorm(prod(dims)*t),c(t,dims))
  } else if (cov == "separable"){
    E <- array(rnorm(prod(dims)*t),c(dims,t)) # d_1 * ... * d_K * t
    for(i in 1:dd){
      cov.matrix <- matrix(rho,dims[i],dims[i]) + diag(rep(1-rho,dims[i]))
      x.eig <- eigen(cov.matrix)
      values <- x.eig$values
      values.sqrt <- rep(0,nrow(cov.matrix))
      values.sqrt[values > 0] <- sqrt(values[values > 0])
      cov.half <- x.eig$vectors %*% diag(values.sqrt) %*% t(x.eig$vectors)
      E <- tensor(E,cov.half,1,1) # d_{i+1} * ... * d_K * t * d_1 *... *d_i
    }
  } else if (cov == "random"){
    E <- matrix(rnorm(prod(dims)*t),t,prod(dims)) %*% qr.Q(qr(matrix(rnorm(prod(dims)^2),prod(dims))))
    E <- array(E,c(t,dims))
  } else {
    stop("Please specify cov")
  }
  X <- lambda * X + E
  #return(as.tensor(X))
  return(X)
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
