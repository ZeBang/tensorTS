###Functions of Autoregressive Models

#' Estimation for Autoregressive Model of Tensor-Valued Time Series
#'
#' Estimation function for tensor-valued time series.
#' Projection method (proj), the Iterated least squares method (LSE), MLE under a Kronecker structured covariance matrix (MLE) and stacked vector AR(1) model (VAR), as determined by the value of \code{type}.
#'@name tenAR
#'@rdname tenAR
#'@aliases tenAR
#'@usage tenAR(xx, method)
#'@export
#'@param xx \eqn{T \times d_1 \times \cdots \times d_K} tensor-valued time series
#'@param R number of terms
#'@param P number of lag-one orders
#'@param method character string, specifying the type of the estimation method to be used. \describe{
#'  \item{\code{"PROJ",}}{Projection method.}
#'  \item{\code{"LSE",}}{Iterated least squares.}
#'  \item{\code{"MLE",}}{MLE under a Kronecker structured covariance matrix.}
#'  \item{\code{"VAR",}}{Stacked vector AR(1) Model.}
#'}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices \eqn{A_1,A_2,\cdots,A_K}}
#'\item{\code{res}}{residual of the model}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
#'@examples
#' dim <- c(3,3,3)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' tenAR(xx, R=1, P=1, method="LSE")
tenAR <- function(xx, R=1, P=1, method="LSE", k1=NULL, k2=NULL){
  if (identical("PROJ", method)) {
    tenAR.PROJ(xx, R, P)
  } else if (identical("LSE", method)) {
    tenAR.LS(xx, R, P)
  } else if (identical("MLE", method)) {
    tenAR.MLE(xx, R, P)
  } else if (identical("VAR", method)) {
    tenAR.VAR(xx, P)
  } else if (identical("RRLSE", method)) {
    MAR1.RR(xx, k1, k2)
  } else if (identical("RRMLE", method)) {
    MAR1.CC(xx, k1, k2)
  } else {
    stop("Please specify the type you want to use. See manuals or run ?tenAR for details.")
  }
}


#' Estimation for Autoregressive Model of Matrix-Valued Time Series
#'
#' Although MAR model is a special case in tenAR model. We still include this original function as reference.
#' Estimation function for matrix-valued time series, including the projection method (PROJ),
#' the Iterated least squares method (LSE) and MLE under a Kronecker structured covariance (MLE) and vectorized AR(1) model (VAR) as determined by the value of \code{method}.
#'@name MAR
#'@rdname MAR
#'@aliases MAR
#'@usage MAR(xx, method)
#'@export
#'@param xx T * p * q matrix-valued time series
#'@param method character string, specifying the method of the estimation to be used. \describe{
#'  \item{\code{"PROJ",}}{Projection method.}
#'  \item{\code{"LSE",}}{Iterated least squares.}
#'  \item{\code{"MLE",}}{MLE under a Kronecker structured covariance matrix.}
#'  \item{\code{"VAR",}}{Stacked vector AR(1) Model.}
#'}
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of LL, a p by p matrix}
#'\item{\code{RR}}{estimator of RR, a q by q matrix}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'}
#'@examples
#' dim <- c(3,3,3)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' MAR(xx, method="LSE")
MAR <- function(xx, method){
  if (identical("PROJ", method)) {
    MAR1.PROJ(xx)
  }
  if (identical("LSE", method)) {
    MAR1.LS(xx)
  }
  if (identical("MLE", method)) {
    MAR1.MLE(xx)
  }
  if (identical("VAR", method)) {
    tenAR.VAR(xx, P=1)
  }
  else {
    return("Please specify the type you want to use. See manuals for details.")
  }
}


#' Projection Method for MAR(1)
#'
#' MAR(1) one step projection estimation in the model \eqn{X_t = LL \times X_{t-1} \times {RR}^{\prime} + E_t}.
#'@name MAR1.PROJ
#'@rdname MAR1.PROJ
#'@aliases MAR1.PROJ
#'@export
#'@param xx T * p * q matrix-valued time series
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of LL, a p by p matrix}
#'\item{\code{RR}}{estimator of RR, a q by q matrix}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'}
MAR1.PROJ <- function(xx){
  # xx: T * p * q
  # X_t = LL X_{t-1} RR + E_t
  # Sig = cov(vec(E_t))
  # one-step projection estimation
  # Return LL, RR, and estimate of Sig
  dd <- dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  xx.mat <- matrix(xx,T,p*q)
  kroneck <- t(xx.mat[2:T,]) %*% xx.mat[1:(T-1),] %*% solve(t(xx.mat[1:(T-1),]) %*% xx.mat[1:(T-1),])
  ans.projection <- projection(kroneck,r=1,q,p,q,p)
  a <- svd(ans.projection$C,nu=0,nv=0)$d[1]
  LL <- ans.projection$C / a
  RR <- t(ans.projection$B) * a
  res=xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,1),LL,2,2),c(1,3,2))
  Sig <- matrix(tensor(res,res,1,1),p*q)/(T-1)
  return(list(LL=LL,RR=RR,res=res,Sig=Sig))
}





#' Least Squares Iterative Estimation for Matrix Time Series
#'
#' Iterated least squares estimation in the model \eqn{X_t = LL \times X_{t-1} \times {RR}^{\prime} + E_t}.
#'@name MAR1.LS
#'@rdname MAR1.LS
#'@aliases MAR1.LS
#'@export
#'@param xx T * p * q matrix-valued time series
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true printe LL and RR
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of LL, a p by p matrix}
#'\item{\code{RR}}{estimator of RR, a q by q matrix}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{dis}}{Frobenius norm difference of last update}
#'\item{\code{niter}}{number of iterations}
#'}
MAR1.LS <- function(xx,niter=50,tol=1e-6,print.true = FALSE){
  # xx: T * p * q
  # X_t = LL X_{t-1} RR + E_t
  # Sig = cov(vec(E_t))
  # LS criterion
  # iterative algorithm between LL <--> RR
  # Return LL, RR, and estimate of Sig
  dd=dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  LL.old <- diag(p)
  RR.old <- diag(q)
  dis <- 1
  iiter <- 1

  while(iiter <= niter & dis >= tol){
    # estimate RR0
    temp <- tensor(xx[1:(T-1),,,drop=FALSE],LL.old,2,2)  # (T-1) * q * p
    AA <- tensor(temp,temp,c(1,3),c(1,3))
    BB <- tensor(temp,xx[2:T,,,drop=FALSE],c(1,3),c(1,2))
    RR <- solve(AA,BB)
    # estimate LL0
    temp <- tensor(xx[1:(T-1),,,drop=FALSE],RR,3,1)  # (T-1) * p * q
    AA <- tensor(temp,temp,c(1,3),c(1,3))
    BB <- t(tensor(temp,xx[2:T,,,drop=FALSE],c(1,3),c(1,3)))
    LL <- t(solve(t(AA),t(BB)))
    a <- svd(LL,nu=0,nv=0)$d[1]
    LL <- LL / a
    RR <- RR * a
    # update for the next iteration
    dis <- sqrt(sum((kronecker(t(RR),LL)-kronecker(t(RR.old),LL.old))^2))
    LL.old <- LL
    RR.old <- RR
    iiter <- iiter + 1
    if(print.true==TRUE){
      print(LL)
      print(RR)
    }
  }
  res=xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,1),LL,2,2),c(1,3,2))
  Sig <- matrix(tensor(res,res,1,1),p*q)/(T-1)
  return(list(LL=LL,RR=RR,res=res,Sig=Sig,dis=dis,niter=iiter))
}


#' MLE under a structured covariance tensor
#'
#' MAR(1) iterative estimation with Kronecker covariance structure: \eqn{X_t = LL \times X_{t-1} \times {RR}^{\prime} + E_t} such that \eqn{\Sigma = cov(vec(E_t)) = \Sigma_r \otimes \Sigma_l}.
#'@name MAR1.MLE
#'@rdname MAR1.MLE
#'@aliases MAR1.MLE
#'@usage MAR1.MLE(xx,LL.init=NULL,Sigl.init=NULL,Sigr.init=NULL,niter=50,
#'tol=1e-6,print.true = FALSE)
#'@export
#'@param xx T * p * q matrix-valued time series
#'@param LL.init initial value of LL
#'@param Sigl.init initial value of Sigl
#'@param Sigr.init initial value of Sigr
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true print LL and RR
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of LL, a p by p matrix}
#'\item{\code{RR}}{estimator of RR, a q by q matrix}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sigl}}{one part of structured covariance matrix \eqn{\Sigma=\Sigma_r \otimes \Sigma_l}}
#'\item{\code{Sigr}}{one part of structured covariance matrix \eqn{\Sigma=\Sigma_r \otimes \Sigma_l}}
#'\item{\code{dis}}{Frobenius norm difference of the final update step}
#'\item{\code{niter}}{number of iterations}
#'}
MAR1.MLE <- function(xx,LL.init=NULL,Sigl.init=NULL,Sigr.init=NULL,niter=50,tol=1e-6,print.true = FALSE){
  # xx: T * p * q
  # X_t = LL X_{t-1} RR + E_t
  # Sig = cov(vec(E_t)) = Sigr otimes Sigl
  # optimization criterion is likelihood
  # iterative algorithm between LL <--> RR <--> Sig_r <--> Sig_l
  # Return LL, RR, Sigl, Sigr
  dd=dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  if(is.null(LL.init)){
    LL.old <- diag(p)
  } else{
    LL.old <- LL.init
  }
  RR.old <- diag(q)
  if(is.null(Sigl.init)){
    Sigl.old <- diag(p)
  } else{
    Sigl.old <- Sigl.init
  }
  if(is.null(Sigr.init)){
    Sigr.old <- diag(q)
  } else{
    Sigr.old <- Sigr.init
  }
  dis <- 1
  iiter <- 1

  while(iiter <= niter & dis >= tol){
    Sigl.inv.old <- ginv(Sigl.old)
    Sigr.inv.old <- ginv(Sigr.old)
    # estimate RR0
    temp1 <- tensor(xx[1:(T-1),,,drop=FALSE],LL.old,2,2)  # (T-1) * q * p
    temp2 <- tensor(temp1,Sigl.inv.old,3,2)  # (T-1) * q * p
    AA <- tensor(temp1,temp2,c(1,3),c(1,3))
    BB <- tensor(temp2,xx[2:T,,,drop=FALSE],c(1,3),c(1,2))
    RR <- solve(AA,BB)
    # estimate LL0
    temp1 <- tensor(xx[1:(T-1),,,drop=FALSE],RR,3,1)  # (T-1) * p * q
    temp2 <- tensor(temp1,Sigr.inv.old,3,1)  # (T-1) * p * q
    AA <- tensor(temp1,temp2,c(1,3),c(1,3))
    BB <- t(tensor(temp2,xx[2:T,,,drop=FALSE],c(1,3),c(1,3)))
    LL <- t(solve(t(AA),t(BB)))
    res=xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,1),LL,2,2),c(1,3,2))
    temp <- tensor(res,Sigl.inv.old,2,1) # (T-1) * q * p
    Sigr <- tensor(temp,res,c(1,3),c(1,2))/T/p
    temp <- tensor(res,ginv(Sigr),3,1) # (T-1) * p * q
    Sigl <- tensor(temp,res,c(1,3),c(1,3))/T/q
    a <- svd(LL,nu=0,nv=0)$d[1]
    LL <- LL / a
    RR <- RR * a
    a <- eigen(Sigl)$values[1]
    Sigl <- Sigl / a
    Sigr <- Sigr * a

    # update for the next iteration
    dis1 <- sqrt(sum((kronecker(t(RR),LL)-kronecker(t(RR.old),LL.old))^2))
    dis2 <- sqrt(sum((kronecker(Sigr,Sigl)-kronecker(Sigr.old,Sigl.old))^2))
    dis <- max(dis1,dis2)
    LL.old <- LL
    RR.old <- RR
    Sigr.old <- Sigr
    Sigl.old <- Sigl
    iiter <- iiter + 1
    if(print.true==TRUE){
      print(LL)
      print(RR)
    }
  }
  Sig <- kronecker(Sigr,Sigl)
  return(list(LL=LL,RR=RR,res=res,Sigl=Sigl,Sigr=Sigr,Sig=Sig,dis=dis,niter=iiter))
}


#' Asymptotic Covariance Matrix of \code{MAR1.MLE}
#'
#' Asymptotic covariance Matrix of \code{MAR1.MLE} for given A, B and matrix-valued time series xx, see Theory 3 in paper.
#'@name MAR.SE
#'@rdname MAR.SE
#'@aliases MAR.SE
#'@export
#'@param xx T * p * q matrix-valued time series
#'@param A p by p matrix in MAR(1) model
#'@param B q by q matrix in MAR(1) model
#'@param Sig covariance matrix cov(vec(E_t)) in MAR(1) model
#'@return asymptotic covariance matrix
#'@examples
#'# given T * p * q time series xx
#'out = MAR1.LS(xx)
#'FnormLL = sqrt(sum(out$LL))
#'xdim=3;ydim=4
#'out.Xi=MAR.SE(xx.nm,out$RR*FnormLL,out$LL/FnormLL,out$Sig)
#'out.SE=sqrt(diag(out.Xi))
#'SE.A=matrix(out.SE[1:xdim^2],nrow=xdim)
#'SE.B=t(matrix(out.SE[-(1:xdim^2)],nrow=ydim))
MAR.SE <- function(xx, B, A, Sigma){
  dd <- dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  BX <- tensor(xx,B,3,2) # Tpq
  AX <- tensor(xx,A,2,2) # Tqp
  BXI <- apply(BX,1,function(x){kronecker(t(x),diag(p))})
  AXI <- apply(AX,1,function(x){kronecker(diag(q),t(x))})
  BXI.array <- array(BXI,c(p*q,p^2,T)) #pq*p^2*T
  AXI.array <- array(AXI,c(p*q,q^2,T)) #pq*q^2*T
  # W transpose
  Wt <- abind(BXI.array,AXI.array,along=2) #pq*(p^2+q^2)*T
  EWWt <- tensor(Wt,Wt,c(1,3),c(1,3))/T #(p^2+q^2)*(p^2+q^2)
  alpha <- as.vector(A)
  beta <- as.vector(t(B))
  gamma <- c(alpha,rep(0,q^2))
  H <- EWWt + gamma %*% t(gamma) #(p^2+q^2)*(p^2+q^2)
  WSigma <- tensor(Wt,Sigma,1,1) #(p^2+q^2)*T*pq
  EWSigmaWt <- tensor(WSigma,Wt,c(3,2),c(1,3))/T
  Hinv <- solve(H)
  Xi <- Hinv %*% EWSigmaWt %*% Hinv
  Xi <- Xi/T
}


#' Stacked vector AR(p) Model
#'
#' vector AR(p) Model for vectorized Tensor time series
#'@name tenAR.VAR
#'@rdname tenAR.VAR
#'@aliases tenAR.VAR
#'@export
#'@param xx \eqn{T * d_1 * \cdots * d_K} tensor-valued time series
#'@param P number of orders
#'@return a list containing the following:\describe{
#'\item{\code{coef}}{coeficient of the fitted VAR(1) model}
#'\item{\code{res}}{residual of the VAR(1) model}
#'}
#'@examples
#'out.var = tenAR.VAR(xx)
#'sum(out.var$res**2)
tenAR.VAR <- function(xx, P){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  dd=xx@modes
  t <- dd[1]
  n <- prod(dd[-1])
  if (prod(dd[-1]) > t){stop("sample size T too small")}
  yy=apply(xx@data,MARGIN=1,as.vector)

  x = 0
  for (l in c(1:P)){
    if (l == 1){x = yy[,(P):(t-1)]} else {x = rbind(x, yy[,(P+1-l):(t-l)])}
  }
  y = yy[,(P+1):t]
  coef = t(solve(x %*% t(x)) %*% x %*% t(y))
  res= y - coef %*% x
  phi = list()
  for (l in c(1:P)){
    phi[[l]] = coef[,(n*(l-1)+1):(n*l)]
  }
  return(list(coef=phi, res=res))
}


#' Generate coefficient matrices
#'
#' For test only, generate coefficient matrices in TenAR model
#'@name tenAR.A
#'@rdname tenAR.A
#'@aliases tenAR.A
#'@export
#'@param dim an array of dimensions of matrices or tensors
#'@param R number of terms
#'@param P number of orders
#'@param rho spectral radius
#'@return a list containing coefficient matrices
#'@examples
#' dim <- c(3,3,3)
#' A <- tenAR.A(dim,R=2,P=1,rho=0.5)
#' xx <- tenAR.xx(t=500, A, setting='iid')
tenAR.A <- function(dim,R,P,rho){
  K <- length(dim)
  A <- lapply(1:P, function(p) {lapply(1:R, function(j) {lapply(1:K, function(i) {diag(dim[i])})})})
  for (p in c(1:P)){
    for (j in c(1:R)){
      for (i in c(1:K)){
        A[[p]][[j]][[i]] <- matrix(rnorm(dim[i]^2), c(dim[i],dim[i]))
      }
    }
  }

  eigen = M.eigen(A)
  while (eigen >= 1){
    for (l in c(1:P)){
      for (j in c(1:R)){
        A[[l]][[j]][[K]] <- rho * matrix(rnorm(dim[i]^2), c(dim[i],dim[i]))/ eigen
      }
      A[[l]] <- fro.order(fro.rescale(A[[l]]))
    }
    eigen = M.eigen(A)
  }
  stopifnot(M.eigen(A) < 1)
  return(A)
}


#' Generate an AR(1) tensor time series with given coefficient matrices
#'
#' For test only, with given coefficient matrices and time length t, generate an AR(1) tensor time series.
#'@name tenAR.xx
#'@rdname tenAR.xx
#'@aliases tenAR.xx
#'@export
#'@param t length of time
#'@param A true coefficient matrices (the dimensions, number of terms and orders are implied by structure of A)
#'@param setting the structure of random error, can be specified as "iid", "mle" and "svd", default value is "iid".
#'@return Tensor-valued time series generated by TenAR model
#'@seealso \code{\link{run.test}}
#'@examples
#' dim <- c(3,3,3)
#' A <- tenAR.A(dim,R=2,P=1,rho=0.5)
#' xx <- tenAR.xx(t=500, A, setting='iid')
tenAR.xx <- function(t,A, setting){
  P <- length(A)
  R <- length(A[[1]])
  K <- length(A[[1]][[1]])
  dim <- c()
  for (i in c(1:K)){
    dim[i] <- nrow(A[[1]][[1]][[i]])
  }
  x <- array(0, c(t,prod(dim)))
  for (l in c(1:P)){
    x[l,] <- rnorm(prod(dim))
  }
  for (i in c((P+1):t)){
    if (setting == "iid"){
      e <- rnorm(prod(dim))
    } else if (setting == "mle"){
      e <- kronecker_list(rev(Sigma.true)) %*% rnorm(prod(dim))
    } else if (setting == "svd"){
      e <- E %*% rnorm(prod(dim))
    } else {
      return("Please specify setting")
    }
    temp = 0
    for (l in c(1:P)){
      phi <-  Reduce("+", lapply(1:R, function(j) {kronecker_list(rev(A[[l]][[j]]))}))
      temp = temp + phi %*% x[i-l, ]
    }
    x[i,] <-  temp + e
  }
  return(array(x, c(t, dim)))
}


#' Projection estimation for TenAR(p) model
#'
#' projection estimation
#'@name tenAR.PROJ
#'@rdname tenAR.PROJ
#'@aliases tenAR.PROJ
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param R number of terms
#'@param P number of orders
#'@return a list containing the estimation of matrices
tenAR.PROJ <- function(xx,R,P){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  dim <- xx@modes[-1]
  mm <- tenAR.VAR(xx, P)$coef
  A = list()
  for (p in c(1:P)){
    tt <- trearrange(mm[[p]], rev(dim))
    tt <- as.tensor(aperm(tt@data))
    A[[p]] <- fro.order(fro.rescale(ten.proj(tt, dim, R)))
  }
  return(list(A=A))
}



#' Least Squares Iterative Estimation for TenAR(p)
#'
#' Iterated least squares estimation
#'@name tenAR.LS
#'@rdname tenAR.LS
#'@aliases tenAR.LS
#'@import tensor rTensor
#'@export
#'@param xx  \eqn{T * d_1 * \cdots * d_K} tensor-valued time series
#'@param R number of terms
#'@param P number of orders
#'@param init initial value of A
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true print \eqn{A_i}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices}
#'\item{\code{res}}{residual of the model}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
tenAR.LS <- function(xx,R, P, init=NULL, niter=500,tol=1e-6,print.true = FALSE){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  dim <- xx@modes[-1]
  K <- length(dim)
  t <- xx@modes[[1]]

  if (is.null(init)) {A.old <- tenAR.PROJ(xx,R,P)$A} else {A.old <- init}

  A.new <- A.old

  dis <- 1
  iiter <- 1
  a <- c()
  while(iiter <= niter & dis >= tol & dis <= 1e3){
    for (p in c(1:P)){
      for (r in c(1:R)){
        for (k in c(1:K)){
          s0 <- rTensor::ttl(xx, A.new[[p]][[r]][-k], c(2:(K+1))[-k])
          temp <- s0@data[(1+P-p):(t-p),,,,drop=FALSE]

          L1 <- 0
          for (l in c(1:P)){
            if (l == p){
              if (R > 1){
                L1 <- L1 + Reduce("+",lapply(c(1:R)[-r], function(n) {(rTensor::ttl(xx[(1+P-l):(t-l),,,], A.new[[l]][[n]], (c(1:K) + 1)))}))
              }
            } else {
              L1 <- L1 + Reduce("+",lapply(c(1:R), function(n) {(rTensor::ttl(xx[(1+P-l):(t-l),,,], A.new[[l]][[n]], (c(1:K) + 1)))}))
            }
          }


          L2 <-  xx[(1+P):t,,,,drop=FALSE] - L1
          temp2 <- L2@data

          RR <- tensor(temp,temp,c(1:4)[-(k+1)],c(1:4)[-(k+1)])
          LL <- tensor(temp2,temp,c(1:4)[-(k+1)],c(1:4)[-(k+1)])
          A.new[[p]][[r]][[k]] <- LL %*% solve(RR)
        }
      }
    }
    for (p in c(1:P)){
      for (r in c(1:R)){
        a <- c()
        for (k in c(1:K)){
          m <- A.new[[p]][[r]][[k]]
          if (k != K){
            a[k] <- svd(m,nu=0,nv=0)$d[1]
            A.new[[p]][[r]][[k]] <- m/a[k]
          } else {
            A.new[[p]][[r]][[k]] <- m * prod(a)
          }
        }
      }
    }
    phi.new <- list()
    phi.old <- list()
    for (p in c(1:P)){
      phi.new[[p]] <- Reduce("+", lapply(1:R, function(j) {rTensor::kronecker_list(rev(A.new[[p]][[j]]))}))
      phi.old[[p]] <- Reduce("+", lapply(1:R, function(j) {rTensor::kronecker_list(rev(A.old[[p]][[j]]))}))
    }

    dis <- Reduce("+", lapply(1:P, function(j) {sqrt(sum((phi.new[[j]] - phi.old[[j]])^2))}))
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  for (p in c(1:P)){
    A.new[[p]] <- fro.order(fro.rescale(A.new[[p]]))
  }
  res <- ten.res(xx,A.new,P,R,K)
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new,niter=iiter,Sig=Sig,res=res))
}



#' MLEs for TenAR(p) Model
#'
#' MLEs under Kronecker structures covariance matrix
#'@name tenAR.MLE
#'@rdname tenAR.MLE
#'@aliases tenAR.MLE
#'@import tensor rTensor
#'@export
#'@param xx  \eqn{T * d_1 * \cdots * d_K} tensor-valued time series
#'@param R number of terms
#'@param P number of orders
#'@param init.A initial value of A
#'@param init.sig initial value of Sigma
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true print \eqn{A_i}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices}
#'\item{\code{sigma}}{estimator of Sigma}
#'\item{\code{res}}{residual}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
tenAR.MLE <- function(xx, R, P, init.A=NULL, init.sig=NULL, niter=500,tol=1e-5, print.true = FALSE){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  dim <- xx@modes[-1]
  K <- length(dim)
  t <- xx@modes[[1]]
  if (is.null(init.sig)) {Sig.old <- lapply(1:K, function(i) {diag(dim[i])})} else {Sig.old <- init.sig}
  Sig.new <- Sig.old
  Sig.new.inv <- lapply(1:K, function (k) {solve(Sig.new[[k]])})
  if (is.null(init.A)) {A.old <- tenAR.PROJ(xx, R, P)$A} else {A.old <- init.A}
  A.new <- A.old
  dis <- 1
  iiter <- 1
  while(iiter <= niter & dis >= tol){
    stopifnot(dis <= 1e3)
    for (p in c(1:P)){
      for (r in c(1:R)){
        for (k in c(1:K)){

          sphi <-  lapply(1:K, function (k) {Sig.new.inv[[k]] %*% (A.new[[p]][[r]][[k]])})
          s0 <- rTensor::ttl(xx, A.new[[p]][[r]][-k], c(2:(K+1))[-k])
          temp <- s0@data[(1+P-p):(t-p),,,,drop=FALSE] # X_{t-1,k} * t(Phi_k^r)

          s1 <- rTensor::ttl(xx, sphi[-k],c(2:(K+1))[-k])
          temp1 <- s1@data[(1+P-p):(t-p),,,,drop=FALSE] # X_{t-1,k} * t(S_k^{-1}*Phi_k^r)

          L1 <- 0
          for (l in c(1:P)){
            if (l == p){
              if (R > 1){
                L1 <- L1 + Reduce("+",lapply(c(1:R)[-r], function(n) {(rTensor::ttl(xx[(1+P-l):(t-l),,,], A.new[[l]][[n]], (c(1:K) + 1)))}))
              }
            } else {
              L1 <- L1 + Reduce("+",lapply(c(1:R), function(n) {(rTensor::ttl(xx[(1+P-l):(t-l),,,], A.new[[l]][[n]], (c(1:K) + 1)))}))
            }
          }
          temp2 <-  (xx[(1+P):t,,,,drop=FALSE] - L1)@data


          RR <- tensor(temp,temp1,c(1:4)[-(k+1)],c(1:4)[-(k+1)])
          LL <- tensor(temp2,temp1,c(1:4)[-(k+1)],c(1:4)[-(k+1)])
          A.new[[p]][[r]][[k]] <- LL %*% solve(RR)
        }
      }
    }

    for (k in c(1:K)){
      res.old <- rTensor::as.tensor(ten.res(xx,A.new,P,R,K))
      rs <- rTensor::ttl(res.old, Sig.new.inv[-k], c(2:(K+1))[-k])
      Sig.new[[k]] <- tensor(res.old@data, rs@data, c(1:4)[-(k+1)],c(1:4)[-(k+1)])/(t-1)/prod(dim[-k])
      Sig.new.inv <- lapply(1:K, function (i) {solve(Sig.new[[i]])})
    }

    for (p in c(1:P)){
      A.new[[P]] <- svd.rescale(A.new[[P]])
    }
    Sig.new <- eigen.rescale(list(Sig.new))[[1]]

    phi.new <- list()
    phi.old <- list()
    for (p in c(1:P)){
      phi.new[[p]] <- Reduce("+", lapply(1:R, function(j) {rTensor::kronecker_list(rev(A.new[[p]][[j]]))}))
      phi.old[[p]] <- Reduce("+", lapply(1:R, function(j) {rTensor::kronecker_list(rev(A.old[[p]][[j]]))}))
    }
    dis <- Reduce("+", lapply(1:P, function(j) {sqrt(sum((phi.new[[j]] - phi.old[[j]])^2))}))

    Sig.old <- Sig.new
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }

  for (p in c(1:P)){
    A.new[[p]] <- fro.order(fro.rescale(A.new[[p]]))
  }
  res <- ten.res(xx,A.new,P,R,K)
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new, SIGMA=Sig.new, niter=iiter, Sig=Sig, res=res))
}


#' Reduced Rank MAR(1) iterative estimation
#'
#' LSE for the matrix reduced rank model \eqn{X_t = LL \times X_{t-1} \times {RR}^{\prime} + E_t}.
#'@name MAR1.RR
#'@rdname MAR1.RR
#'@aliases MAR1.RR
#'@import MASS
#'@usage MAR1.RR(xx, k1, k2, niter=200, tol=1e-4, print.true=FALSE,
#' LL.init=NULL, RR.init=NULL)
#'@param xx matrix-valued time series
#'@param k1 rank of A1
#'@param k2 rank of A2
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true print \eqn{A_i}
#'@param LL.init initial values of A1 in iterations, by default is diagonal matrices
#'@param RR.init initial values of A2 in iterations, by default is diagonal matrices
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of coefficient matrices \eqn{A_1}}
#'\item{\code{RR}}{estimator of coefficient matrices \eqn{A_2}}
#'\item{\code{res}}{residual}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{BIC}}{value of information criteria}
#'\item{\code{dis}}{Frobenius norm difference of last update}
#'\item{\code{niter}}{number of iterations}
#'}
#'@export
#'@examples
#'dim <- c(4,5) # dimension of matrix at time t is 4 * 5
#'A <- tenAR.A(dim, R=1, P=1, rho=0.5)
#'xx <- tenAR.xx(t=500, A, setting="iid")
#'model <- MAR1.RR(xx, k1=2, k2=2)
MAR1.RR <- function(xx, k1, k2, niter=200, tol=1e-4, print.true=FALSE, LL.init=NULL, RR.init=NULL){
  # xx: T * p * q
  # X_t = LL X_{t-1} RR' + E_t ### NOTE! This function is written with RR'.
  # Sig = cov(vec(E_t)) = Sigr \otimes Sigl
  # optimization criterion is likelihood
  # iterative algorithm between LL <--> Sigma_l <--> RR <--> Sig_r
  # Return LL, RR, Sigl, Sigr
  dd=dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  if(is.null(LL.init)){
    LL.old <- diag(p)
  } else{
    LL.old <- LL.init
  }
  if(is.null(RR.init)){
    RR.old <- diag(q)
  } else{
    RR.old <- RR.init
  }
  Tol=tol*sqrt(p^2+q^2)
  dis <- 1
  iiter <- 1

  while(iiter <= niter & dis >= Tol){


    ## Save old
    LL.oold=LL.old
    RR.oold=RR.old

    # estimate LL0
    temp1 <- tensor(xx[1:(T-1),,,drop=FALSE],RR.old,3,2)  # (T-1) * p * q
    AA <- tensor(temp1,temp1,c(1,3),c(1,3))
    BB <- tensor(xx[2:T,,,drop=FALSE],temp1,c(1,3),c(1,3))
    LL <- BB%*%ginv(AA)
    U <- svd(LL%*%t(BB))$u[,1:k1]
    LL <- U%*%t(U)%*%LL
    # update for next iteration
    a=svd(LL,nu=0,nv=0)$d[1]
    LL=LL/a
    dis3=sum((LL-LL.old)^2)
    LL.old <- LL

    # estimate RR0
    temp1 <- tensor(xx[1:(T-1),,,drop=FALSE],LL.old,2,2)  # (T-1) * q * p
    AA <- tensor(temp1,temp1,c(1,3),c(1,3))
    BB <- tensor(xx[2:T,,,drop=FALSE],temp1,c(1,2),c(1,3))
    RR <- BB%*%ginv(AA)
    U <- svd(RR%*%t(BB))$u[,1:k2]
    RR <- U%*%t(U)%*%RR
    # update for next iteration
    dis3=dis3+sum((RR-RR.old)^2)
    RR.old <- RR

    # update for the next iteration
    dis1 <- sqrt(sum((kronecker(t(RR),LL)-kronecker(t(RR.oold),LL.oold))^2))
    dis3 = sqrt(dis3)
    dis <- dis3
    iiter <- iiter + 1
    if(print.true==TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  a <- sqrt(sum(LL^2))
  LL <- LL / a
  RR <- RR * a
  res=xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,2),LL,2,2),c(1,3,2))
  Sig <- matrix(tensor(res,res,1,1),p*q)/(T-1)
  bic <- T*p*q*log(sum(res^2/(T*p*q))) ##+log(T*p*q)*(bic.penalty(p,k1)+bic.penalty(q,k2))
  return(list(LL=LL,RR=RR,res=res,Sig=Sig,BIC=bic,dis=dis,niter=iiter-1))
}


#' Reduced Rank MAR(1) iterative estimation with Kronecker covariance structure
#'
#' MLE for the matrix reduced rank model \eqn{X_t = LL \times X_{t-1} \times {RR}^{\prime} + E_t}.
#'@name MAR1.CC
#'@rdname MAR1.CC
#'@aliases MAR1.CC
#'@usage MAR1.CC(xx, k1, k2, LL.init=NULL,RR.init=LL,Sigl.init=NULL,Sigr.init=NULL,
#' niter=200,tol=1e-4,print.true = FALSE)
#'@export
#'@param xx  matrix-valued time series
#'@param k1  rank of A1
#'@param k2  rank of A2
#'@param niter maximum number of iterations if error stays above \code{tol}, by default niter=200
#'@param tol relative Frobenius norm error tolerance, by default tol=1e-4
#'@param print.true print \eqn{A_i}
#'@param LL.init  initial values of A1 in iterations, by default is diagonal matrices
#'@param RR.init  initial values of A2 in iterations, by default is diagonal matrices
#'@param Sigl.init  initial values of Sigma_1 in iterations, by default is diagonal matrices
#'@param Sigr.init  initial values of Sigma_2 in iterations, by default is diagonal matrices
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimated \eqn{A_1}}
#'\item{\code{RR}}{estimated \eqn{A_2}}
#'\item{\code{Sigl}}{estimated Sigma_1}
#'\item{\code{Sigr}}{estimated Sigma_2}
#'\item{\code{res}}{residual}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{dis}}{Frobenius norm difference of last update}
#'\item{\code{niter}}{number of iterations}
#'}
#'@examples
#'dim <- c(4,5) # dimension of matrix at time t is 4 * 5
#'A <- tenAR.A(dim, R=1, P=1, rho=0.5)
#'xx <- tenAR.xx(t=500, A, setting="iid")
#'model <- MAR1.RR(xx, k1=2, k2=2)
MAR1.CC <- function(xx,k1,k2,LL.init=NULL,RR.init=LL,Sigl.init=NULL,Sigr.init=NULL,niter=200,tol=1e-4,print.true = FALSE){
  # xx: T * p * q
  # X_t = LL X_{t-1} RR' + E_t ### NOTE! This function is written with RR'.
  # Sig = cov(vec(E_t)) = Sigr \otimes Sigl
  # optimization criterion is likelihood
  # iterative algorithm between LL <--> Sigma_l <--> RR <--> Sig_r
  # Return LL, RR, Sigl, Sigr
  dd=dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  if(is.null(LL.init)){
    LL.old <- diag(p)
  } else{
    LL.old <- LL.init
  }
  if(is.null(RR.init)){
    RR.old <- diag(q)
  } else{
    RR.old <- RR.init
  }
  if(is.null(Sigl.init)){
    Sigl.old <- diag(p)
  } else{
    Sigl.old <- Sigl.init
  }
  if(is.null(Sigr.init)){
    Sigr.old <- diag(q)
  } else{
    Sigr.old <- Sigr.init
  }
  Tol=tol*sqrt(p^2+q^2)
  dis <- 1
  iiter <- 1

  while(iiter <= niter & dis >= Tol){

    ## Save old
    LL.oold=LL.old
    RR.oold=RR.old
    Sigl.oold=Sigl.old
    Sigr.oold=Sigr.old

    # estimate LL0 and Sigl
    Sigr.inv.old <- ginv(Sigr.old)
    temp1 <- tensor(xx[1:(T-1),,,drop=FALSE],RR.old,3,2)  # (T-1) * p * q
    temp2 <- tensor(temp1,Sigr.inv.old,3,1)  # (T-1) * p * q
    AA <- tensor(temp1,temp2,c(1,3),c(1,3))
    BB <- tensor(xx[2:T,,,drop=FALSE],temp2,c(1,3),c(1,3))
    LL <- BB%*%ginv(AA)
    res <- xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR.old,3,2),LL,2,2),c(1,3,2)) # (T-1) * p * q
    temp <- tensor(res,Sigr.inv.old,3,1) # (T-1) * p * q
    Sigl <- tensor(temp,res,c(1,3),c(1,3))/T/q
    Sigl.spec <- eigen(Sigl)
    Sigl.root <- Sigl.spec$vectors%*%diag(sqrt(Sigl.spec$values))%*%t(Sigl.spec$vectors)
    Sigl.root.inv <- Sigl.spec$vectors%*%diag(1/sqrt(Sigl.spec$values))%*%t(Sigl.spec$vectors)
    U <- svd(Sigl.root.inv%*%LL%*%t(BB)%*%Sigl.root.inv)$u[,1:k1]
    LL <- Sigl.root%*%U%*%t(U)%*%Sigl.root.inv%*%LL
    res <- xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR.old,3,2),LL,2,2),c(1,3,2)) # (T-1) * p * q
    temp <- tensor(res,Sigr.inv.old,3,1) # (T-1) * p * q
    Sigl <- tensor(temp,res,c(1,3),c(1,3))/T/q
    # update for next iteration
    a=svd(LL,nu=0,nv=0)$d[1]
    LL=LL/a
    dis3=sum((LL-LL.old)^2)
    LL.old <- LL
    Sigl.old <- Sigl


    # estimate RR0 and Sigr
    Sigl.inv.old <- ginv(Sigl.old)
    temp1 <- tensor(xx[1:(T-1),,,drop=FALSE],LL.old,2,2)  # (T-1) * q * p
    temp2 <- tensor(temp1,Sigl.inv.old,3,1)  # (T-1) * q * p
    AA <- tensor(temp1,temp2,c(1,3),c(1,3))
    BB <- tensor(xx[2:T,,,drop=FALSE],temp2,c(1,2),c(1,3))
    RR <- BB%*%ginv(AA)
    res <- xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,2),LL.old,2,2),c(1,3,2)) # (T-1) * p * q
    temp <- tensor(res,Sigl.inv.old,2,1) # (T-1) * q * p
    Sigr <- tensor(temp,res,c(1,3),c(1,2))/T/p
    Sigr.spec <- eigen(Sigr)
    Sigr.root <- Sigr.spec$vectors%*%diag(sqrt(Sigr.spec$values))%*%t(Sigr.spec$vectors)
    Sigr.root.inv <- Sigr.spec$vectors%*%diag(1/sqrt(Sigr.spec$values))%*%t(Sigr.spec$vectors)
    U <- svd(Sigr.root.inv%*%RR%*%t(BB)%*%Sigr.root.inv)$u[,1:k2]
    RR <- Sigr.root%*%U%*%t(U)%*%Sigr.root.inv%*%RR
    res <- xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,2),LL.old,2,2),c(1,3,2)) # (T-1) * p * q
    temp <- tensor(res,Sigl.inv.old,2,1) # (T-1) * q * p
    Sigr <- tensor(temp,res,c(1,3),c(1,2))/T/p
    # update for next iteration
    dis3=dis3+sum((RR-RR.old)^2)
    RR.old <- RR
    Sigr.old <- Sigr


    a <- eigen(Sigl)$values[1]
    Sigl <- Sigl / a
    Sigr <- Sigr * a
    ### cat(eigen(Sigl)$values[1]," ",eigen(Sigr)$values[1], " ")

    # update for the next iteration
    dis1 <- sqrt(sum((kronecker(t(RR),LL)-kronecker(t(RR.oold),LL.oold))^2))
    ### cat(dis1," ")
    dis2 <- sqrt(sum((kronecker(Sigr,Sigl)-kronecker(Sigr.oold,Sigl.oold))^2))
    ### cat(dis2," ",dis3,"\n")
    dis3 = sqrt(dis3)
    ### dis <- max(dis1,dis2)
    dis <- dis3
    Sigr.old <- Sigr
    Sigl.old <- Sigl
    iiter <- iiter + 1
    if(print.true==TRUE){
      print(LL)
      print(RR)
    }
  }
  a <- sqrt(sum(LL^2))
  LL <- LL / a
  RR <- RR * a
  Sig <- kronecker(Sigr,Sigl)
  return(list(LL=LL,RR=RR,res=res,Sigl=Sigl,Sigr=Sigr,Sig=Sig,dis=dis,niter=iiter-1))
}


#' Asymptotic Covariance Matrix of LSE in TenAR(1)
#'
#' Asymptotic covariance Matrix of LSE in TenAR(1) for given a tensor-valued time series xx, see related Theorems in our paper.
#'@name tenAR.SE.LSE
#'@rdname tenAR.SE.LSE
#'@aliases tenAR.SE.LSE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param A.true coefficient matrices in TAR(1) model
#'@param Sigma covariance matrix cov(vec(E_t)) in TAR(1) model
#'@return asmptotic covariance matrix
#'@examples
#' dim <- c(2,2,2)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' Sigma <- diag(prod(dim))
#' SIGMA <- tenAR.SE.LSE(xx, A[[1]], Sigma)
tenAR.SE.LSE <- function(xx, A.true, Sigma){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  r <- length(A.true)
  dim <- xx@modes[-1]
  m1 <- dim[1]
  m2 <- dim[2]
  m3 <- dim[3]
  k <- length(dim)
  T <- xx@modes[[1]]
  ndim <- m1^2+m2^2+m3^2
  Gamma <- matrix(0,r*ndim,r*ndim)
  r1 <- matrix(0, r*ndim, 1)
  for (x in c(1:r)){
    for (y in c(1:(k-1))) {
      if (y == 1) {
        a <- as.matrix(as.vector(A.true[[x]][[y]]))
        r1[((x-1)*ndim + 1):((x-1)*ndim + m1^2),] <- a
      } else if (y == 2) {
        a <- as.matrix(as.vector(A.true[[x]][[y]]))
        r1[((x-1)*ndim + m1^2 + 1):((x-1)*ndim + m1^2 + m2^2),] <- a
      } else  {
        print("Only Support k=3 but now k>3")
      }
      Gamma <- Gamma + r1 %*% t(r1)

    }
  }
  Hdim <- r*ndim
  HT <- array(1,c(T,Hdim,Hdim))
  WT <- array(1,c(T,Hdim,prod(dim))) #T*r(m1^2+m2^2+m^3)*(m1m2m3)
  for (i in c(1:r)) {
    C <- A.true[[i]][[3]]
    B <- A.true[[i]][[2]]
    A <- A.true[[i]][[1]]
    for (t in c(1:T)){
      w1 <- k_unfold(xx[t,,,], m = 1)@data %*% t(kronecker(C,B))
      w2 <- k_unfold(xx[t,,,], m = 2)@data %*% t(kronecker(C,A))
      w3 <- k_unfold(xx[t,,,], m = 3)@data %*% t(kronecker(B,A))
      w <- rbind(kronecker(w1,diag(m1)) ,kronecker(w2,diag(m2)) %*% kronecker(diag(m3),pm(m2,m1)), kronecker(w3,diag(m3)) %*% pm(m3,m2*m1))
      WT[t, ((i-1)*ndim + 1):(i*ndim),] <-  w
    }
  }
  WSigma <- tensor(WT,Sigma,3,1) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  EWSigmaWt <- tensor(WSigma,WT,c(3,1),c(3,1))/T
  H <- tensor(WT,WT,c(3,1),c(3,1))/T + Gamma #r(m1^2+m2^2+m^3)*r(m1^2+m2^2+m^3)
  Hinv <- solve(H)
  Xi <- Hinv %*% EWSigmaWt %*% Hinv
  return(Xi)
}


#' Asymptotic Covariance Matrix of MLEs in TenAR(1)
#'
#' Asymptotic covariance Matrix of MLEs in TenAR(1) for given a tensor-valued time series xx, see related Theorems in our paper.
#'@name tenAR.SE.MLE
#'@rdname tenAR.SE.MLE
#'@aliases tenAR.SE.MLE
#'@import tensor rTensor
#'@export
#'@param xx  \eqn{T * d_1 * \cdots * d_K} tensor-valued time series
#'@param A.true coefficient matrices in TAR(1) model
#'@param Sigma covariance matrix cov(vec(E_t)) in TAR(1) model
#'@return asmptotic covariance matrix
#'@examples
#' dim <- c(2,2,2)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' Sigma <- diag(prod(dim))
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' SIGMA <- tenAR.SE.MLE(xx, A[[1]], Sigma)
tenAR.SE.MLE <- function(xx, A.true, Sigma){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  r <- length(A.true)
  dim <- xx@modes[-1]
  m1 <- dim[1]
  m2 <- dim[2]
  m3 <- dim[3]
  k <- length(dim)
  T <- xx@modes[[1]]
  ndim <- m1^2+m2^2+m3^2
  Gamma <- matrix(0,r*ndim,r*ndim)
  r1 <- matrix(0, r*ndim, 1)
  for (x in c(1:r)){
    for (y in c(1:(k-1))) {
      if (y == 1) {
        a <- as.matrix(as.vector(A.true[[x]][[y]]))
        r1[((x-1)*ndim + 1):((x-1)*ndim + m1^2),] <- a
      } else if (y == 2) {
        a <- as.matrix(as.vector(A.true[[x]][[y]]))
        r1[((x-1)*ndim + m1^2 + 1):((x-1)*ndim + m1^2 + m2^2),] <- a
      } else  {
        print("Only Support k=3 but now k>3")
      }
      Gamma <- Gamma + r1 %*% t(r1)
    }
  }
  Hdim <- r*ndim
  HT <- array(1,c(T,Hdim,Hdim))
  WT <- array(1,c(T,Hdim,prod(dim))) #T*r(m1^2+m2^2+m^3)*(m1m2m3)
  for (i in c(1:r)) {
    C <- A.true[[i]][[3]]
    B <- A.true[[i]][[2]]
    A <- A.true[[i]][[1]]
    for (t in c(1:T)){
      w1 <- k_unfold(xx[t,,,], m = 1)@data %*% t(kronecker(C,B))
      w2 <- k_unfold(xx[t,,,], m = 2)@data %*% t(kronecker(C,A))
      w3 <- k_unfold(xx[t,,,], m = 3)@data %*% t(kronecker(B,A))
      w <- rbind(kronecker(w1,diag(m1)) ,kronecker(w2,diag(m2)) %*% kronecker(diag(m3),pm(m2,m1)), kronecker(w3,diag(m3)) %*% pm(m3,m2*m1))
      WT[t, ((i-1)*ndim + 1):(i*ndim),] <-  w
    }
  }
  WSigma <- tensor(WT,solve(Sigma),3,1) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  EWSigmaWt <- tensor(WSigma,WT,c(3,1),c(3,1))/T
  H <- EWSigmaWt + Gamma
  Hinv <- solve(H)
  Xi <- Hinv %*% EWSigmaWt %*% Hinv
  return(Xi)
}

#' Asymptotic Covariance Matrix of Reduced rank MAR(1)
#'
#' Asymptotic covariance matrix of reduced rank MAR(1) for given a matrix-valued time series xx, see related Theorems in our paper.
#'@name MAR1.RRLS.SE
#'@rdname MAR1.RRLS.SE
#'@aliases MAR1.RRLS.SE
#'@usage MAR1.RRLS.SE(A1,A2,k1,k2,Sigma.e,RU1=diag(k1),RV1=diag(k1),RU2=diag(k2),
#'RV2=diag(k2),mpower=100)
#'@importFrom tensor tensor
#'@export
#'@param A1 coefficient matrix
#'@param A2 coefficient matrix
#'@param k1 rank of A1
#'@param k2 rank of A2
#'@param Sigma.e Cov(vec(\eqn{E_t})) = Sigma.e : of dimension \eqn{d_1 d_2 \times d_1 d_2}
#'@param RU1,RV1,RU2,RV2: rotation of U1,V1,U2,V2, e.g., new_U1=U1 RU1
#'@param mpower truncate vec(X_t) to mpower term, i.e. VMA(mpower). By default is 100.
#'@return a list containing the following:\describe{
#'\item{\code{Sigma.LS}}{asymptotic covariance matrix of (vec(\eqn{\hat{A}1}),vec(\eqn{\hat{A}2^{\prime}}))}
#'\item{\code{Theta1.LS.u}}{asymptotic covariance matrix of vec(\eqn{\hat{U}1})}
#'\item{\code{Theta1.LS.v}}{asymptotic covariance matrix of vec(\eqn{\hat{V}1})}
#'\item{\code{Theta2.LS.u}}{asymptotic covariance matrix of vec(\eqn{\hat{U}2})}
#'\item{\code{Theta2.LS.v}}{asymptotic covariance matrix of vec(\eqn{\hat{V}2})}
#'}
MAR1.RRLS.SE <- function(A1,A2,k1,k2,Sigma.e,RU1=diag(k1),RV1=diag(k1),RU2=diag(k2),RV2=diag(k2),mpower=100){
  # iterative least square
  # X_t = A1 X_{t-1} A2^T + E_t
  # RU1,RV1,RU2,RV2: rotation of U1,V1,U2,V2, e.g., new_U1=U1 RU1
  # Cov(vec(E_t)) = Sigma.e : of dimension d1d2\times d1d2
  # truncate vec(X_t) to mpower term, i.e. VMA(mpower)
  # return Sigma.LS: asymptotic covariance matrix of (vec(\hat A1),vec(\hat A2^T))
  #  Theta1.LS.u, Theta1.LS.v: asymptotic covariance matrix of vec(\hat U1), vec(\hat V1)
  #  Theta2.LS.u, Theta2.LS.v: asymptotic covariance matrix of vec(\hat U2), vec(\hat V2)
  d1 <- dim(A1)[1]
  d2 <- dim(A2)[1]
  Sigma.ee <- array(Sigma.e,c(d1,d2,d1,d2))
  B <- kronecker(A2,A1)
  #Sigma.x.vec, Sigma.xx: covariance matrix of vec(X_t)
  Sigma.x.vec <- Sigma.e
  for(i in 1:mpower){
    Sigma.x.vec <- Sigma.x.vec+(B%^%i)%*%Sigma.e%*%t(B%^%i)
  }
  Sigma.xx <- array(Sigma.x.vec,c(d1,d2,d1,d2))
  Gamma1 <- tensor(Sigma.xx,t(A1)%*%A1,c(1,3),c(1,2))
  Gamma2 <- tensor(Sigma.xx,t(A2)%*%A2,c(2,4),c(1,2))
  Gamma3 <- kronecker(diag(1,d1),A1)%*%matrix(aperm(Sigma.xx,c(3,1,4,2)),d1^2,d2^2)%*% kronecker(t(A2),diag(1,d2))
  H <- matrix(NA,d1^2+d2^2,d1^2+d2^2)
  H[1:d1^2,1:d1^2] <- as.vector(A1)%*%t(as.vector(A1)) + kronecker(Gamma2,diag(1,d1))
  H[1:d1^2,(d1^2+1):(d1^2+d2^2)] <- Gamma3
  H[(d1^2+1):(d1^2+d2^2),1:d1^2] <- t(Gamma3)
  H[(d1^2+1):(d1^2+d2^2),(d1^2+1):(d1^2+d2^2)] <- kronecker(diag(1,d2),Gamma1)
  H.inv <- solve(H)
  D1 <- Gamma2%*%t(A1)%*% ginv(A1%*%Gamma2%*%t(A1)) %*%A1
  D2 <- Gamma1%*%t(A2)%*% ginv(A2%*%Gamma1%*%t(A2)) %*%A2
  P1 <- A1%*%ginv(t(A1)%*%A1)%*%t(A1)
  P2 <- A2%*%ginv(t(A2)%*%A2)%*%t(A2)
  Sigma.Q <- matrix(NA,d1^2+d2^2,d1^2+d2^2)
  M1 <- kronecker(t(A2),P1)%*%Sigma.e%*%kronecker(A2,P1)
  M2 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e%*%kronecker(A2,P1)
  M3 <- kronecker(t(A2),P1)%*%Sigma.e%*%kronecker(A2,diag(d1)-P1)
  M4 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e%*%kronecker(A2,diag(d1)-P1)
  Sigma.Q1 <- Sigma.Q2 <- Sigma.Q3 <- Sigma.Q4 <- matrix(0,d1^2,d1^2)
  for(i in 1:d1){
    for(j in 1:d1){
      for(k in 1:d1){
        for(l in 1:d1){
          Sigma.Q1[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M1[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
          Sigma.Q2[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M2[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
          Sigma.Q3[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M3[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
          Sigma.Q4[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M4[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
        }
      }
    }
  }
  Sigma.Q[1:d1^2,1:d1^2] <- Sigma.Q1+kronecker(D1,diag(d1))%*%Sigma.Q2+Sigma.Q3%*%kronecker(t(D1),diag(d1))+
    kronecker(D1,diag(d1))%*%Sigma.Q4%*%kronecker(t(D1),diag(d1))
  M1 <- kronecker(P2,t(A1))%*%Sigma.e%*%kronecker(P2,A1)
  M2 <- kronecker(P2,t(A1))%*%Sigma.e%*%kronecker(diag(d2)-P2,A1)
  M3 <- kronecker(diag(d2)-P2,t(A1))%*%Sigma.e%*%kronecker(P2,A1)
  M4 <- kronecker(diag(d2)-P2,t(A1))%*%Sigma.e%*%kronecker(diag(d2)-P2,A1)
  Sigma.Q1 <- Sigma.Q2 <- Sigma.Q3 <- Sigma.Q4 <- matrix(0,d2^2,d2^2)
  for(i in 1:d2){
    for(j in 1:d2){
      for(k in 1:d2){
        for(l in 1:d2){
          Sigma.Q1[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M1[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
          Sigma.Q2[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M2[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
          Sigma.Q3[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M3[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
          Sigma.Q4[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M4[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
        }
      }
    }
  }
  Sigma.Q[(d1^2+1):(d1^2+d2^2),(d1^2+1):(d1^2+d2^2)] <- Sigma.Q1+Sigma.Q2%*%kronecker(diag(d2),t(D2))+
    kronecker(diag(d2),D2)%*%Sigma.Q3+kronecker(diag(d2),D2)%*%Sigma.Q4%*%kronecker(diag(d2),t(D2))
  M1 <- kronecker(t(A2),P1)%*%Sigma.e%*%kronecker(P2,A1)
  M2 <- kronecker(t(A2),P1)%*%Sigma.e%*%kronecker(diag(d2)-P2,A1)
  M3 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e%*%kronecker(P2,A1)
  M4 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e%*%kronecker(diag(d2)-P2,A1)
  Sigma.Q1 <- Sigma.Q2 <- Sigma.Q3 <- Sigma.Q4 <- matrix(0,d1^2,d2^2)
  for(i in 1:d1){
    for(j in 1:d1){
      for(k in 1:d2){
        for(l in 1:d2){
          Sigma.Q1[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M1[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
          Sigma.Q2[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M2[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
          Sigma.Q3[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M3[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
          Sigma.Q4[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M4[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
        }
      }
    }
  }
  Sigma.Q[1:d1^2,(d1^2+1):(d1^2+d2^2)] <- Sigma.Q1+Sigma.Q2%*%kronecker(diag(d2),t(D2))+
    kronecker(D1,diag(d1))%*%Sigma.Q3+kronecker(D1,diag(d1))%*%Sigma.Q4%*%kronecker(diag(d2),t(D2))
  Sigma.Q[(d1^2+1):(d1^2+d2^2),1:d1^2] <- t(Sigma.Q[1:d1^2,(d1^2+1):(d1^2+d2^2)])
  Sigma.LS <- H.inv%*%Sigma.Q%*%H.inv
  Sigma.LS11 <- Sigma.LS[1:d1^2,1:d1^2]
  Sigma.LS11.t <- array(Sigma.LS11,c(d1,d1,d1,d1))
  Sigma.LS11.t <- aperm(Sigma.LS11.t,c(2,1,4,3))
  Sigma.LS11.t <- matrix(Sigma.LS11.t, d1^2, d1^2)
  Sigma.LS22.t <- Sigma.LS[(d1^2+1):(d1^2+d2^2),(d1^2+1):(d1^2+d2^2)]
  Sigma.LS22 <- array(Sigma.LS22.t,c(d2,d2,d2,d2))
  Sigma.LS22 <- aperm(Sigma.LS22,c(2,1,4,3))
  Sigma.LS22 <- matrix(Sigma.LS22, d2^2, d2^2)
  svd.A1 <- svd(A1)
  #k1 <- length(svd.A1$d[svd.A1$d>1e-10])
  D1 <- diag(c(svd.A1$d[1:k1],1))[1:k1,1:k1]
  U1 <- svd.A1$u[,1:k1]
  U1c <- svd.A1$u[,(k1+1):d1]
  V1 <- svd.A1$v[,1:k1]
  V1c <- svd.A1$v[,(k1+1):d1]
  e1 <- diag(d1)
  J1 <- matrix(0,d1^2,d1^2)
  for(i in 1:d1){
    J1[,((i-1)*d1+1):(i*d1)] <- kronecker(diag(d1),e1[,i])
  }
  C1 <- kronecker(A1,diag(d1)) + kronecker(diag(d1),A1)%*%J1
  Sigma.LS.1u <- C1%*%Sigma.LS11%*%t(C1)
  C1.t <- kronecker(t(A1),diag(d1)) + kronecker(diag(d1),t(A1))%*%J1
  Sigma.LS.1v <- C1.t%*%Sigma.LS11.t%*%t(C1.t)
  e1 <- diag(k1)
  L1 <- matrix(0,k1^2,k1)
  for(i in 1:k1){
    L1[,i] <- kronecker(e1[,i],e1[,i])
  }
  R1u <- cbind(kronecker(diag(k1),U1), kronecker(diag(k1),U1c)) %*% rbind(solve(kronecker(D1^2,diag(k1)) -
                                                                                  kronecker(diag(k1),D1^2)+L1%*%t(L1)) %*% (diag(k1^2)-L1%*%t(L1)) %*% kronecker(t(U1),t(U1)),
                                                                          kronecker(diag(1/c(svd.A1$d[1:k1],1))[1:k1,1:k1]^2%*%t(U1), t(U1c))  )
  R1v <- cbind(kronecker(diag(k1),V1), kronecker(diag(k1),V1c)) %*% rbind(solve(kronecker(D1^2,diag(k1)) -
                                                                                  kronecker(diag(k1),D1^2)+L1%*%t(L1)) %*% (diag(k1^2)-L1%*%t(L1)) %*% kronecker(t(V1),t(V1)),
                                                                          kronecker(diag(1/c(svd.A1$d[1:k1],1))[1:k1,1:k1]^2%*%t(V1), t(V1c))  )
  Theta1.LS.u <- R1u%*% Sigma.LS.1u %*%t(R1u)
  Theta1.LS.v <- R1v%*% Sigma.LS.1v %*%t(R1v)
  Theta1.LS.u <- kronecker(t(RU1),diag(d1))%*%Theta1.LS.u
  Theta1.LS.v <- kronecker(t(RV1),diag(d1))%*%Theta1.LS.v
  svd.A2 <- svd(A2)
  #k2 <- length(svd.A2$d[svd.A2$d>1e-10])
  D2 <- diag(c(svd.A2$d[1:k2],1))[1:k2,1:k2]
  U2 <- svd.A2$u[,1:k2]
  U2c <- svd.A2$u[,(k2+1):d2]
  V2 <- svd.A2$v[,1:k2]
  V2c <- svd.A2$v[,(k2+1):d2]
  e2 <- diag(d2)
  J2 <- matrix(0,d2^2,d2^2)
  for(i in 1:d2){
    J2[,((i-1)*d2+1):(i*d2)] <- kronecker(diag(d2),e2[,i])
  }
  C2 <- kronecker(A2,diag(d2)) + kronecker(diag(d2),A2)%*%J2
  Sigma.LS.2u <- C2%*%Sigma.LS22%*%t(C2)
  C2.t <- kronecker(t(A2),diag(d2)) + kronecker(diag(d2),t(A2))%*%J2
  Sigma.LS.2v <- C2.t%*%Sigma.LS22.t%*%t(C2.t)
  e2 <- diag(k2)
  L2 <- matrix(0,k2^2,k2)
  for(i in 1:k2){
    L2[,i] <- kronecker(e2[,i],e2[,i])
  }
  R2u <- cbind(kronecker(diag(k2),U2), kronecker(diag(k2),U2c)) %*% rbind(solve(kronecker(D2^2,diag(k2)) -
                                                                                  kronecker(diag(k2),D2^2)+L2%*%t(L2)) %*% (diag(k2^2)-L2%*%t(L2)) %*% kronecker(t(U2),t(U2)),
                                                                          kronecker(diag(1/c(svd.A2$d[1:k2],1))[1:k2,1:k2]^2%*%t(U2), t(U2c))  )
  R2v <- cbind(kronecker(diag(k2),V2), kronecker(diag(k2),V2c)) %*% rbind(solve(kronecker(D2^2,diag(k2)) -
                                                                                  kronecker(diag(k2),D2^2)+L2%*%t(L2)) %*% (diag(k2^2)-L2%*%t(L2)) %*% kronecker(t(V2),t(V2)),
                                                                          kronecker(diag(1/c(svd.A2$d[1:k2],1))[1:k2,1:k2]^2%*%t(V2), t(V2c))  )
  Theta2.LS.u <- R2u%*% Sigma.LS.2u %*%t(R2u)
  Theta2.LS.v <- R2v%*% Sigma.LS.2v %*%t(R2v)
  Theta2.LS.u <- kronecker(t(RU2),diag(d2))%*%Theta2.LS.u
  Theta2.LS.v <- kronecker(t(RV2),diag(d2))%*%Theta2.LS.v
  return(list("Sigma.LS"=Sigma.LS,"Theta1.LS.u"=Theta1.LS.u,"Theta1.LS.v"=Theta1.LS.v,
              "Theta2.LS.u"=Theta2.LS.u,"Theta2.LS.v"=Theta2.LS.v))
}

#' Asymptotic Covariance Matrix of Reduced rank MAR(1) with Kronecker covariance structure
#'
#' Asymptotic covariance matrix of Reduced rank MAR(1) Kronecker covariance structure for given a matrix-valued time series xx, see related Theorems in our paper.
#'@name MAR1.RRCC.SE
#'@rdname MAR1.RRCC.SE
#'@aliases MAR1.RRCC.SE
#'@usage MAR1.RRCC.SE((A1,A2,k1,k2,Sigma1,Sigma2,RU1=diag(k1),RV1=diag(k1),RU2=diag(k2),
#'RV2=diag(k2),mpower=100))
#'@importFrom tensor tensor
#'@export
#'@param A1 coefficient matrix
#'@param A2 coefficient matrix
#'@param k1 rank of A1
#'@param k2 rank of A2
#'@param Sigma1,Sigma2 Cov(vec(\eqn{E_t})) = Sigma1 \eqn{\otimes} Sigma2
#'@param RU1,RV1,RU2,RV2: rotation of U1,V1,U2,V2, e.g., new_U1=U1 RU1
#'@param mpower truncate vec(X_t) to mpower term, i.e. VMA(mpower). By default is 100.
#'@return a list containing the following:\describe{
#'\item{\code{Sigma.CC}}{asymptotic covariance matrix of (vec(\eqn{\hat{A}1}),vec(\eqn{\hat{A}2^{\prime}}))}
#'\item{\code{Theta1.CC.u}}{asymptotic covariance matrix of vec(\eqn{\hat{U}1})}
#'\item{\code{Theta1.CC.v}}{asymptotic covariance matrix of vec(\eqn{\hat{V}1})}
#'\item{\code{Theta2.CC.u}}{asymptotic covariance matrix of vec(\eqn{\hat{U}2})}
#'\item{\code{Theta2.CC.v}}{asymptotic covariance matrix of vec(\eqn{\hat{V}2})}
#'}
MAR1.RRCC.SE <- function(A1,A2,k1,k2,Sigma1,Sigma2,RU1=diag(k1),RV1=diag(k1),RU2=diag(k2),RV2=diag(k2),mpower=100){
  # canonical correlation analysis
  # X_t = A1 X_{t-1} A2^T + E_t
  # RU1,RV1,RU2,RV2: rotation of U1,V1,U2,V2, e.g., new_U1=U1 RU1
  # Cov(vec(E_t)) = Sigma.e=Sigma2 \otimes \Sigma1 : of dimension d1d2\times d1d2
  # truncate vec(X_t) to mpower term, i.e. VMA(mpower)
  # return Sigma.CC: asymptotic covariance matrix of (vec(\hat A1),vec(\hat A2^T))
  #  Theta1.CC.u, Theta1.CC.v: asymptotic covariance matrix of vec(\hat U1), vec(\hat V1)
  #  Theta2.CC.u, Theta2.CC.v: asymptotic covariance matrix of vec(\hat U2), vec(\hat V2)
  d1 <- dim(A1)[1]
  d2 <- dim(A2)[1]
  Sigma1.inv=solve(Sigma1)
  Sigma2.inv=solve(Sigma2)
  Sigma.e=kronecker(Sigma2,Sigma1)
  Sigma.e.inv=kronecker(Sigma2.inv,Sigma1.inv)
  Sigma.ee <- array(Sigma.e,c(d1,d2,d1,d2))
  B <- kronecker(A2,A1)
  #Sigma.x.vec, Sigma.xx: covariance matrix of vec(X_t)
  Sigma.x.vec <- Sigma.e
  for(i in 1:mpower){
    Sigma.x.vec <- Sigma.x.vec+(B%^%i)%*%Sigma.e%*%t(B%^%i)
  }
  Sigma.xx <- array(Sigma.x.vec,c(d1,d2,d1,d2))
  Gamma1 <- tensor(Sigma.xx,t(A1)%*%Sigma1.inv%*%A1,c(1,3),c(1,2))
  Gamma2 <- tensor(Sigma.xx,t(A2)%*%Sigma2.inv%*%A2,c(2,4),c(1,2))
  Gamma3 <- kronecker(diag(1,d1),Sigma1.inv%*%A1)%*%matrix(aperm(Sigma.xx,c(3,1,4,2)),d1^2,d2^2)%*%
    kronecker(t(A2)%*%Sigma2.inv,diag(1,d2))
  H <- matrix(NA,d1^2+d2^2,d1^2+d2^2)
  H[1:d1^2,1:d1^2] <- as.vector(A1)%*%t(as.vector(A1)) + kronecker(Gamma2,Sigma1.inv)
  H[1:d1^2,(d1^2+1):(d1^2+d2^2)] <- Gamma3
  H[(d1^2+1):(d1^2+d2^2),1:d1^2] <- t(Gamma3)
  H[(d1^2+1):(d1^2+d2^2),(d1^2+1):(d1^2+d2^2)] <- kronecker(Sigma2.inv,Gamma1)
  H.inv <- solve(H)
  D1 <- Gamma2%*%t(A1)%*% ginv(A1%*%Gamma2%*%t(A1)) %*%A1
  D2 <- Gamma1%*%t(A2)%*% ginv(A2%*%Gamma1%*%t(A2)) %*%A2
  P1 <- Sigma1.inv%*%A1%*%ginv(t(A1)%*%Sigma1.inv%*%A1)%*%t(A1)
  P2 <- Sigma2.inv%*%A2%*%ginv(t(A2)%*%Sigma2.inv%*%A2)%*%t(A2)
  Sigma.Q <- matrix(NA,d1^2+d2^2,d1^2+d2^2)
  M1 <- kronecker(t(A2),P1)%*%Sigma.e.inv%*%kronecker(A2,P1)
  M2 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e.inv%*%kronecker(A2,P1)
  M3 <- kronecker(t(A2),P1)%*%Sigma.e.inv%*%kronecker(A2,diag(d1)-P1)
  M4 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e.inv%*%kronecker(A2,diag(d1)-P1)
  Sigma.Q1 <- Sigma.Q2 <- Sigma.Q3 <- Sigma.Q4 <- matrix(0,d1^2,d1^2)
  for(i in 1:d1){
    for(j in 1:d1){
      for(k in 1:d1){
        for(l in 1:d1){
          Sigma.Q1[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M1[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
          Sigma.Q2[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M2[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
          Sigma.Q3[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M3[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
          Sigma.Q4[(i-1)*d1+j,(k-1)*d1+l] <- sum(Sigma.xx[i,,k,] * M4[(1:d2-1)*d1+j,(1:d2-1)*d1+l])
        }
      }
    }
  }
  Sigma.Q[1:d1^2,1:d1^2] <- Sigma.Q1+kronecker(D1,diag(d1))%*%Sigma.Q2+Sigma.Q3%*%kronecker(t(D1),diag(d1))+
    kronecker(D1,diag(d1))%*%Sigma.Q4%*%kronecker(t(D1),diag(d1))
  M1 <- kronecker(P2,t(A1))%*%Sigma.e.inv%*%kronecker(P2,A1)
  M2 <- kronecker(P2,t(A1))%*%Sigma.e.inv%*%kronecker(diag(d2)-P2,A1)
  M3 <- kronecker(diag(d2)-P2,t(A1))%*%Sigma.e.inv%*%kronecker(P2,A1)
  M4 <- kronecker(diag(d2)-P2,t(A1))%*%Sigma.e.inv%*%kronecker(diag(d2)-P2,A1)
  Sigma.Q1 <- Sigma.Q2 <- Sigma.Q3 <- Sigma.Q4 <- matrix(0,d2^2,d2^2)
  for(i in 1:d2){
    for(j in 1:d2){
      for(k in 1:d2){
        for(l in 1:d2){
          Sigma.Q1[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M1[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
          Sigma.Q2[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M2[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
          Sigma.Q3[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M3[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
          Sigma.Q4[(i-1)*d2+j,(k-1)*d2+l] <- sum(Sigma.xx[,j,,l] * M4[(i-1)*d1+(1:d1),(k-1)*d1+(1:d1)])
        }
      }
    }
  }
  Sigma.Q[(d1^2+1):(d1^2+d2^2),(d1^2+1):(d1^2+d2^2)] <- Sigma.Q1+Sigma.Q2%*%kronecker(diag(d2),t(D2))+
    kronecker(diag(d2),D2)%*%Sigma.Q3+kronecker(diag(d2),D2)%*%Sigma.Q4%*%kronecker(diag(d2),t(D2))
  M1 <- kronecker(t(A2),P1)%*%Sigma.e.inv%*%kronecker(P2,A1)
  M2 <- kronecker(t(A2),P1)%*%Sigma.e.inv%*%kronecker(diag(d2)-P2,A1)
  M3 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e.inv%*%kronecker(P2,A1)
  M4 <- kronecker(t(A2),diag(d1)-P1)%*%Sigma.e.inv%*%kronecker(diag(d2)-P2,A1)
  Sigma.Q1 <- Sigma.Q2 <- Sigma.Q3 <- Sigma.Q4 <- matrix(0,d1^2,d2^2)
  for(i in 1:d1){
    for(j in 1:d1){
      for(k in 1:d2){
        for(l in 1:d2){
          Sigma.Q1[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M1[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
          Sigma.Q2[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M2[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
          Sigma.Q3[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M3[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
          Sigma.Q4[(i-1)*d1+j,(k-1)*d2+l] <- sum(Sigma.xx[i,,,l] * M4[(1:d2-1)*d1+j,(k-1)*d1+(1:d1)])
        }
      }
    }
  }
  Sigma.Q[1:d1^2,(d1^2+1):(d1^2+d2^2)] <- Sigma.Q1+Sigma.Q2%*%kronecker(diag(d2),t(D2))+
    kronecker(D1,diag(d1))%*%Sigma.Q3+kronecker(D1,diag(d1))%*%Sigma.Q4%*%kronecker(diag(d2),t(D2))
  Sigma.Q[(d1^2+1):(d1^2+d2^2),1:d1^2] <- t(Sigma.Q[1:d1^2,(d1^2+1):(d1^2+d2^2)])
  Sigma.CC <- H.inv%*%Sigma.Q%*%H.inv
  Sigma.CC11 <- Sigma.CC[1:d1^2,1:d1^2]
  Sigma.CC11.t <- array(Sigma.CC11,c(d1,d1,d1,d1))
  Sigma.CC11.t <- aperm(Sigma.CC11.t,c(2,1,4,3))
  Sigma.CC11.t <- matrix(Sigma.CC11.t, d1^2, d1^2)
  Sigma.CC22.t <- Sigma.CC[(d1^2+1):(d1^2+d2^2),(d1^2+1):(d1^2+d2^2)]
  Sigma.CC22 <- array(Sigma.CC22.t,c(d2,d2,d2,d2))
  Sigma.CC22 <- aperm(Sigma.CC22,c(2,1,4,3))
  Sigma.CC22 <- matrix(Sigma.CC22, d2^2, d2^2)
  svd.A1 <- svd(A1)
  #k1 <- length(svd.A1$d[svd.A1$d>1e-10])
  D1 <- diag(c(svd.A1$d[1:k1],1))[1:k1,1:k1]
  U1 <- svd.A1$u[,1:k1]
  U1c <- svd.A1$u[,(k1+1):d1]
  V1 <- svd.A1$v[,1:k1]
  V1c <- svd.A1$v[,(k1+1):d1]
  e1 <- diag(d1)
  J1 <- matrix(0,d1^2,d1^2)
  for(i in 1:d1){
    J1[,((i-1)*d1+1):(i*d1)] <- kronecker(diag(d1),e1[,i])
  }
  C1 <- kronecker(A1,diag(d1)) + kronecker(diag(d1),A1)%*%J1
  Sigma.CC.1u <- C1%*%Sigma.CC11%*%t(C1)
  C1.t <- kronecker(t(A1),diag(d1)) + kronecker(diag(d1),t(A1))%*%J1
  Sigma.CC.1v <- C1.t%*%Sigma.CC11.t%*%t(C1.t)
  e1 <- diag(k1)
  L1 <- matrix(0,k1^2,k1)
  for(i in 1:k1){
    L1[,i] <- kronecker(e1[,i],e1[,i])
  }
  R1u <- cbind(kronecker(diag(k1),U1), kronecker(diag(k1),U1c)) %*% rbind(solve(kronecker(D1^2,diag(k1)) -
                                                                                  kronecker(diag(k1),D1^2)+L1%*%t(L1)) %*% (diag(k1^2)-L1%*%t(L1)) %*% kronecker(t(U1),t(U1)),
                                                                          kronecker(diag(1/c(svd.A1$d[1:k1],1))[1:k1,1:k1]^2%*%t(U1), t(U1c))  )
  R1v <- cbind(kronecker(diag(k1),V1), kronecker(diag(k1),V1c)) %*% rbind(solve(kronecker(D1^2,diag(k1)) -
                                                                                  kronecker(diag(k1),D1^2)+L1%*%t(L1)) %*% (diag(k1^2)-L1%*%t(L1)) %*% kronecker(t(V1),t(V1)),
                                                                          kronecker(diag(1/c(svd.A1$d[1:k1],1))[1:k1,1:k1]^2%*%t(V1), t(V1c))  )
  Theta1.CC.u <- R1u%*% Sigma.CC.1u %*%t(R1u)
  Theta1.CC.v <- R1v%*% Sigma.CC.1v %*%t(R1v)
  Theta1.CC.u <- kronecker(t(RU1),diag(d1))%*%Theta1.CC.u
  Theta1.CC.v <- kronecker(t(RV1),diag(d1))%*%Theta1.CC.v
  svd.A2 <- svd(A2)
  #k2 <- length(svd.A2$d[svd.A2$d>1e-10])
  D2 <- diag(c(svd.A2$d[1:k2],1))[1:k2,1:k2]
  U2 <- svd.A2$u[,1:k2]
  U2c <- svd.A2$u[,(k2+1):d2]
  V2 <- svd.A2$v[,1:k2]
  V2c <- svd.A2$v[,(k2+1):d2]
  e2 <- diag(d2)
  J2 <- matrix(0,d2^2,d2^2)
  for(i in 1:d2){
    J2[,((i-1)*d2+1):(i*d2)] <- kronecker(diag(d2),e2[,i])
  }
  C2 <- kronecker(A2,diag(d2)) + kronecker(diag(d2),A2)%*%J2
  Sigma.CC.2u <- C2%*%Sigma.CC22%*%t(C2)
  C2.t <- kronecker(t(A2),diag(d2)) + kronecker(diag(d2),t(A2))%*%J2
  Sigma.CC.2v <- C2.t%*%Sigma.CC22.t%*%t(C2.t)
  e2 <- diag(k2)
  L2 <- matrix(0,k2^2,k2)
  for(i in 1:k2){
    L2[,i] <- kronecker(e2[,i],e2[,i])
  }
  R2u <- cbind(kronecker(diag(k2),U2), kronecker(diag(k2),U2c)) %*% rbind(solve(kronecker(D2^2,diag(k2)) -
                                                                                  kronecker(diag(k2),D2^2)+L2%*%t(L2)) %*% (diag(k2^2)-L2%*%t(L2)) %*% kronecker(t(U2),t(U2)),
                                                                          kronecker(diag(1/c(svd.A2$d[1:k2],1))[1:k2,1:k2]^2%*%t(U2), t(U2c))  )
  R2v <- cbind(kronecker(diag(k2),V2), kronecker(diag(k2),V2c)) %*% rbind(solve(kronecker(D2^2,diag(k2)) -
                                                                                  kronecker(diag(k2),D2^2)+L2%*%t(L2)) %*% (diag(k2^2)-L2%*%t(L2)) %*% kronecker(t(V2),t(V2)),
                                                                          kronecker(diag(1/c(svd.A2$d[1:k2],1))[1:k2,1:k2]^2%*%t(V2), t(V2c))  )
  Theta2.CC.u <- R2u%*% Sigma.CC.2u %*%t(R2u)
  Theta2.CC.v <- R2v%*% Sigma.CC.2v %*%t(R2v)
  Theta2.CC.u <- kronecker(t(RU2),diag(d2))%*%Theta2.CC.u
  Theta2.CC.v <- kronecker(t(RV2),diag(d2))%*%Theta2.CC.v
  return(list("Sigma.CC"=Sigma.CC,"Theta1.CC.u"=Theta1.CC.u,"Theta1.CC.v"=Theta1.CC.v,
              "Theta2.CC.u"=Theta2.CC.u,"Theta2.CC.v"=Theta2.CC.v))
}


#' select the number of terms by BIC
#'
#' select the number of terms by BIC
#'@name tenAR.bic
#'@rdname tenAR.bic
#'@aliases tenAR.bic
#'@import rTensor
#'@export
#'@param xx  \eqn{T \times d_1 \times \cdots \times d_K} tensor-valued time series
#'@examples
#' dim <- c(3,3,3)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' tenAR.bic(xx, rmax=5)
tenAR.bic <- function(xx, rmax=5){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  dim <- xx@modes[-1]
  t <- xx@modes[[1]]
  ans <- c()
  for (r in c(1:rmax)){
    est <- tenAR.LS(xx, R=r, P=1)
    res <- est$res
    ans[r] <- IC(xx,res,r,t, dim)
  }
  which.min(ans)
}

#' Plot Matrix-Valued Time Series
#'
#' Plot matrix-valued time series or Tensor-Valued time series by given mode.
#'@name mplot
#'@rdname mplot
#'@aliases mplot
#'@export
#'@param xx  \eqn{T \times d_1 \times d_2} matrix-valued time series. Note that the number of mode is 3, where the first mode is time.
#'@examples
#' dim <- c(3,3,3)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' mplot(xx[,,1])
mplot <- function(x){
  if (mode(x) == "S4"){x = x@data}
  dim = dim(x)
  time = array(c(1:dim[1]),dim[1])
  par(mfrow=c(dim[2],dim[3]),mai=0.05*c(1,1,1,1),oma=c(2,2,0,0))
  for(i in 1:dim[2]){
    for(j in 1:dim[3]){
      if(i!=dim[2] & j!=1){
        plot(time,x[,i,j],type='l',xaxt='n',yaxt='n',ylim=range(x[,i,]))
      }
      if(i!=dim[2] & j==1){
        plot(time,x[,i,j],type='l',xaxt='n',ylim=range(x[,i,]))
      }
      if(i==dim[2] & j!=1){
        plot(time,x[,i,j],type='l',yaxt='n',ylim=range(x[,i,]))
      }
      if(i==dim[2] & j==1){
        plot(time,x[,i,j],type='l',ylim=range(x[,i,]))
      }
    }
  }

}


#' Model Predictions
#'
#' predict.tenar is a function for predictions from the results of model fitting functions.
#' The function invokes particular methods which depend on the class of the first argument.
#' Our function is similar to the usage of classical `predict.ar` in package "stats".
#'@name predict.tenar
#'@rdname predict.tenar
#'@aliases predict.tenar
#'@import rTensor abind
#'@usage predict.tenar(object, xx, n.head, method="LSE")
#'@export
#'@param object a fit from TenAR(P) model
#'@param data data to which to apply the prediction.
#'@param n.head number of steps ahead at which to predict.
#'@param method method used by rolling forecast
#'@return predicted value
#'@seealso `predict.ar` or `predict.arima`
#'@examples
#' dim <- c(2,2,2)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' pred.xx <- predict.tenar(model, xx, n.head = 5)
predict.tenar <- function(object, xx, n.head, method="LSE"){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  A <- object$A
  P <- length(A)
  R <- length(A[[1]])
  K <- xx@num_modes - 1
  dim <- xx@modes
  ttt <- (dim[1]+1):(dim[1]+n.head)
  xx.pred <- array(0, c(n.head, dim[-1]))
  for(tt in ttt){
    tti <- tt - ttt[1] + 1
    L1 = 0
    for (l in c(1:P)){
      L1 <- L1 + Reduce("+",lapply(c(1:R), function(n) {(rTensor::ttl(xx[(tt-l),,,,drop=FALSE], A[[l]][[n]], (c(1:K) + 1)))}))
    }
    xx.pred[tti,,,] <- L1[1,,,]@data
    xx <- as.tensor(abind(xx@data, L1@data, along=1))
  }
  return(xx.pred)
}


#' Model Predictions by rolling forecast
#'
#' predict.rolling is a function for predictions from the results of model fitting functions.
#' The function invokes particular methods which depend on the class of the first argument.
#' Our function is similar to the usage of classical `predict.ar` in package "stats".
#'@name predict.rolling
#'@rdname predict.rolling
#'@aliases predict.rolling
#'@import rTensor abind
#'@usage predict.rolling(object, data, n.head, method="LSE")
#'@export
#'@param object a fit from TenAR(P) model
#'@param data data to which to apply the prediction.
#'@param n.head number of steps ahead at which to predict.
#'@param method method used by rolling forecast
#'@return predicted value
#'@seealso `predict.ar` or `predict.arima`
#'@examples
#' dim <- c(2,2,2)
#' A <- tenAR.A(dim,R=2,P=1,rho)
#' xx <- tenAR.xx(t=500, A, setting='iid')
#' pred.xx <- predict.rolling(model, xx, n.head = 5)
predict.rolling <- function(object, data, n.head, method="LSE"){
  if (mode(xx) != "S4") {xx <- rTensor::as.tensor(xx)}
  A <- object$A
  P <- length(A)
  R <- length(A[[1]])
  K <- xx@num_modes - 1
  dim <- xx@modes
  ttt <- (dim[1]+1):(dim[1]+n.head)
  xx.pred <- array(0, c(n.head, dim[-1]))
  for(tt in ttt){
    tti <- tt - ttt[1] + 1
    if (tti > 1){
      print(paste("==================complete",tti))
      xx.new <- array(x.mat[1:(tt-1)], c(dim[1] + tti - 1, dim[-1]))
      xx.nm <- .remove.mean(xx.new)
      model = tenAR(xx.nm, R, P, method)
      A <- model$A
    }
    L1 = 0
    for (l in c(1:P)){
      L1 <- L1 + Reduce("+",lapply(c(1:R), function(n) {(rTensor::ttl(xx[(tt-l),,,,drop=FALSE], A[[l]][[n]], (c(1:K) + 1)))}))
    }
    xx.pred[tti,,,] <- L1[1,,,]@data
    xx <- as.tensor(abind(xx@data, L1@data, along=1))
  }
  return(xx.pred)
}



.getpos <- function(mode, rank){
  pos = 0
  for (k in c(1:length(mode))){
    if (k > 1){mode[k] = mode[k] - 1}
    pos = pos + rank[k]*mode[k]
  }
  return(pos)
}

.getrank <- function(dim){
  rank = array(1, length(dim))
  for (k in c(1:length(dim))){
    if (k > 1){ for (q in c(1:(k-1))){rank[k] = rank[k]*(rev(dim)[q])}}
  }
  return(rank)
}

.remove.mean <- function(xx){
  dim <- xx@modes
  m <- apply(xx@data, c(2:xx@num_modes), mean)
  mm <- aperm(array(m, c(dim[-1],dim[1])), c(xx@num_modes,c(1:(xx@num_modes-1))))
  return(xx - mm)
}
