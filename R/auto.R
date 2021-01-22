###Functions of Autoregressive Models

#' Tensor Method
#'
#' Estimation function for tensor-valued time series.
#' Projection method, the Iterated least squares method, MLE under a structured covariance tensor and stacked vector AR(1) model, as determined by the value of \code{type}.
#'@name TAR
#'@rdname TAR
#'@aliases TAR
#'@usage TAR(xx, type = c("projection", "LS", "MLE", "ar"))
#'@export
#'@param xx \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param dim dimension of coefficient matrices
#'@param method character string, specifying the type of the estimation method to be used. \describe{
#'  \item{\code{"projection",}}{Projection method.}
#'  \item{\code{"lse",}}{Iterated least squares.}
#'  \item{\code{"mle",}}{MLE under a structured covariance tensor.}
#'  \item{\code{"ar",}}{Stacked vector AR(1) Model.}
#'}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices \eqn{A_1,A_2,\cdots,A_K}}
#'\item{\code{res}}{residual of the model}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
TAR <- function(xx, r=1, method="lse"){
  dim = xx@num_modes - 1
  if (dim == 1){
    var1(xx)
  } else if (dim == 2){
    if (r == 1){
      if (identical("projection", method)) {
        MAR1.projection(xx)
      } else if (identical("lse", method)) {
        MAR1.LS(xx)
      } else if (identical("mle", method)) {
        MAR1.otimes(xx)
      } else if (identical("ar", method)) {
        var1(xx)
      } else {
        stop("Please specify the type you want to use. See manuals or run ?TAR for details.")
      }
    } else if (r > 1){
      if (identical("projection", method)) {
        MAR1.projection(xx)
      } else if (identical("lse", method)) {
        MAR1.LS(xx)
      } else if (identical("mle", method)) {
        MAR1.otimes(xx)
      } else if (identical("ar", method)) {
        var1(xx)
      } else {
        stop("Please specify the type you want to use. See manuals or run ?TAR for details.")
      }
    }

  } else if (dim == 3){
    if (identical("projection", method)) {
      TAR1.projection(xx)
    } else if (identical("lse", method)) {
      TAR2.LS(xx,r)
    } else if (identical("mle", method)) {
      TAR2.MLE(xx,r)
    } else if (identical("ar", method)) {
      TAR1.VAR(xx)
    } else {
      stop("Please specify the type you want to use. See manuals or run ?TAR for details.")
    }
  } else {
    stop("lack of or wrong dimension (temporarily we only support dim = 2 or dim = 3)")
  }
}


#' Matrix Method
#'
#' Estimation function for matrix-valued time series, including the projection method,
#' the Iterated least squares method and MLE under a structured covariance tensor, as determined by the value of \code{type}.
#'@name MAR
#'@rdname MAR
#'@aliases MAR
#'@usage MAR(xx, type = c("projection", "LS", "MLE", "ar"))
#'@export
#'@param xx T * p * q matrix-valued time series
#'@param type character string, specifying the type of the estimation method to be used. \describe{
#'  \item{\code{"projection",}}{Projection method.}
#'  \item{\code{"LS",}}{Iterated least squares.}
#'  \item{\code{"MLE",}}{MLE under a structured covariance tensor.}
#'  \item{\code{"ar",}}{Stacked vector AR(1) Model.}
#'}
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of LL, a p by p matrix}
#'\item{\code{RR}}{estimator of RR, a q by q matrix}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'}
MAR <- function(xx, type){
  if (identical("projection", type)) {
    MAR1.projection(xx)
  }
  if (identical("LS", type)) {
    MAR1.LS(xx)
  }
  if (identical("MLE", type)) {
    MAR1.otimes(xx)
  }
  if (identical("ar", type)) {
    var1(xx)
  }
  else {
    return("Please specify the type you want to use. See manuals for details.")
  }
}


#'Rearrangement Operator
#'
#'Rearrangement Operator used for projection method.
#'@name rearrange
#'@rdname rearrange
#'@aliases rearrange
#'@param A m by n matrix such that \eqn{m = m1*n1} and \eqn{n = m2*n2}
#'@param m1 \code{ncol} of A
#'@param m2 \code{ncol} of B
#'@param n1 \code{nrow} of A
#'@param n2 \code{nrow} of B
#'@return rearengement matrix
#'@seealso \code{\link{MAR1.projection}}
#'@examples
#' A <- matrix(runif(6),ncol=2),
#' B <- matrix(runif(6),ncol=2),
#' M <- kronecker(B,A)
#' rearrange(M,3,3,2,2) == t(as.vector(A)) %*% as.vector(B)
#' 'TRUE'
rearrange <- function(A,m1,m2,n1,n2){
  # the inner function of "projection"
  # A: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # A \approx B otimes C
  # return RA

  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2 | m!=m1*m2){
    print("error")
    return
  }
  ans <- matrix(NA, m1*n1, m2*n2)
  for(i in 1:m1){
    for(j in 1:n1){
      ans[(j-1)*m1+i,] <- t(as.vector(A[(i-1)*m2+1:m2,(j-1)*n2+1:n2]))
    }
  }
  ans
}

#' Kronecker Product Approximation
#'
#' Kronecker product approximation used in Projection Method of matrix-value time series.
#'@name projection
#'@rdname projection
#'@aliases projection
#'@param A m by n matrix such that \eqn{m = m1*n1} and \eqn{n = m2*n2}
#'@param r number of terms
#'@param m1 \code{ncol} of A
#'@param m2 \code{ncol} of B
#'@param n1 \code{nrow} of A
#'@param n2 \code{nrow} of B
#'@return a list contaning two estimator (matrix)
#'@seealso \code{\link{MAR1.projection}}
#'@examples
#'A <- matrix(runif(6),ncol=2),
#'projection(A,3,3,2,2)
projection <- function(A,r,m1,m2,n1,n2){
  # the inner function of MAR1.projection
  # A: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # A \approx B \otimes C
  # return B and C
  dim = c(m1, n1)
  RA <- rearrange(A,m1,m2,n1,n2)
  A.hat <- lapply(1:r, function(j) {lapply(1:2, function(i) {diag(dim[i])})})
  if (r == 1){
    RA.svd <- svd(RA,nu=1,nv=1)
    B <- matrix(RA.svd$u * RA.svd$d[1], m1, n1)
    C <- matrix(RA.svd$v , m2, n2)
    list(B=B,C=C)
  } else {
    RA.svd <- svd(RA,nu=r,nv=r)
    for (i in c(1:r)){
      A.hat[[i]][[2]] <- matrix(RA.svd$u[,i] * RA.svd$d[i], m1, n1)
      A.hat[[i]][[1]] <- matrix(RA.svd$v[,i] , m2, n2)
    }
    return(A.hat)
  }
}


#' Projection Method
#'
#' MAR(1) one step projection estimation in the model \eqn{X_t = LL * X_{t-1} * RR + E_t}.
#'@name MAR1.projection
#'@rdname MAR1.projection
#'@aliases MAR1.projection
#'@export
#'@param xx T * p * q matrix-valued time series
#'@return a list containing the following:\describe{
#'\item{\code{LL}}{estimator of LL, a p by p matrix}
#'\item{\code{RR}}{estimator of RR, a q by q matrix}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'}
MAR1.projection <- function(xx){
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


#' Projection Method for matrix time series with multiple terms
#'
#' MAR(1) one step projection estimation in the model \eqn{X_t = \sum_{r=1}^{R} LL_r * X_{t-1} * RR_r + E_t}.
#'@name MAR2.projection
#'@rdname MAR2.projection
#'@aliases MAR2.projection
#'@export
#'@param xx T * p * q matrix-valued time series
#'@return a list containing the estimated coefficient matrices
MAR2.projection <- function(xx){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  xx.mat <- matrix(xx@data,t,dim[1]*dim[2])
  kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% solve(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),])
  A.old <- projection(kroneck, r, dim[1],dim[2],dim[1],dim[2])
  fnorm <- array(0,c(r,k)) # Rescale Result of PROJ
  for (j in c(1:r)){
    for (i in c(1:k)){
      if (i < k ){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]]/fnorm[j,i]
      } else if (i == k){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]] * prod(fnorm[j,1:(k-1)])
      } else {
        print("WRONG dimension")
      }
    }
  }
  A.new <- A.old
}


#' Least Squares Iterative Estimation for Matrix Time Series
#'
#' Iterated least squares estimation in the model \eqn{X_t = LL * X_{t-1} * RR + E_t}.
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
  #-------------
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


#' Least Squares Iterative Estimation for Matrix Time Series with Multiple terms
#'
#' Iterated least squares estimation in the model \eqn{X_t = \sum_{r=1}^{R} LL_r * X_{t-1} * RR_r + E_t}.
#'@name MAR2.LS
#'@rdname MAR2.LS
#'@aliases MAR2.LS
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
MAR2.LS <- function(xx,r,niter=80,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  xx.mat <- matrix(xx@data,t,dim[1]*dim[2])
  kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% solve(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),])
  A.old <- projection(kroneck, r, dim[1],dim[2],dim[1],dim[2])
  fnorm <- array(0,c(r,k)) # Rescale Result of PROJ
  for (j in c(1:r)){
    for (i in c(1:k)){
      if (i < k ){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]]/fnorm[j,i]
      } else if (i == k){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]] * prod(fnorm[j,1:(k-1)])
      } else {
        print("WRONG dimension")
      }
    }
  }
  A.new <- A.old
  dis <- 1
  iiter <- 1
  a <- c()
  while(iiter <= niter & dis >= tol & dis <= 1e3){ # stop when dis > 1e3
    for (j in c(1:r)){
      for (i in c(1:k)){
        s0 <- ttl(xx, A.new[[j]][-i], c(2:(k+1))[-i])
        temp <- s0@data[1:(t-1),,,drop=FALSE]

        L1 <- Reduce("+",lapply(c(1:r)[-j], function(n) {(ttl(xx[1:(t-1),,], A.new[[n]], (c(1:k) + 1)))}))
        if (r == 1){ L1 <- 0}
        L2 <-  xx[2:t,,,drop=FALSE] - L1
        temp2 <- L2@data[1:(t-1),,,drop=FALSE]

        RR <- tensor(temp,temp,c(1:3)[-(i+1)],c(1:3)[-(i+1)])
        LL <- tensor(temp2,temp,c(1:3)[-(i+1)],c(1:3)[-(i+1)])
        A.new[[j]][[i]] <- LL %*% solve(RR)
      }
    }

    for (j in c(1:r)){
      a <- c()
      for (i in c(1:k)){
        m <- A.new[[j]][[i]]
        if (i != k){
          a[i] <- svd(m,nu=0,nv=0)$d[1]
          A.new[[j]][[i]] <- m/a[i]
        } else {
          A.new[[j]][[i]] <- m * prod(a)
        }
      }
    }

    phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
    phi.old <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.old[[j]]))}))
    dis <- sqrt(sum((phi.new - phi.old)^2))
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
  disf <- sum((phi.new -phi)^2)
  res <- (xx[2:t,,,drop=FALSE] - Reduce("+",lapply(1:r, function(j) {(ttl(xx[1:(t-1),,], A.new[[j]], (c(1:k) + 1)))})))@data
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new,niter=iiter,Sig=Sig,res=res,disf=disf))
}


#' MLE under a structured covariance tensor for Matrix Time Series with Multiple Terms
#'
#' MAR(1) iterative estimation with Kronecker covariance structure: \eqn{X_t = \sum_{r=1}^{R} LL_r * X_{t-1} * RR_r + E_t} such that \eqn{\Sigma = cov(vec(E_t)) = \Sigma_r \otimes \Sigma_l}.
#'@name MAR2.otimes
#'@rdname MAR2.otimes
#'@aliases MAR1.otimes
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
#'\item{\code{SIGMA}}{structured covariance matrix \eqn{\Sigma=\Sigma_1, \cdots, \Sigma_R}}
#'\item{\code{dis}}{Frobenius norm difference of the final update step}
#'\item{\code{niter}}{number of iterations}
#'}
MAR2.otimes <- function(xx,r,niter=200,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  xx.mat <- matrix(xx@data,t,dim[1]*dim[2])
  kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% solve(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),])
  A.old <- projection(kroneck, r, dim[1],dim[2],dim[1],dim[2])
  Sig.old <- lapply(1:k, function(i) {diag(dim[i])})
  Sig.new <- Sig.old
  fnorm <- array(0,c(r,k))
  for (j in c(1:r)){
    for (i in c(1:k)){
      if (i < k ){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]]/fnorm[j,i]
      } else if (i == k){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]] * prod(fnorm[j,1:(k-1)])
      } else {
        print("WRONG dimension")
      }
    }
  }
  A.new <- A.old
  dis <- 1
  iiter <- 1
  while(iiter <= niter & dis >= tol & dis <= 1e3){ # stop when dis > 1e3
    for (j in c(1:r)){
      for (i in c(1:k)){
        Sig.new.inv <- lapply(1:k, function (i) {solve(Sig.new[[i]])})
        sphi <-  lapply(1:k, function (i) {Sig.new.inv[[i]] %*% (A.new[[j]][[i]])})

        s0 <- ttl(xx, A.new[[j]][-i], c(2:(k+1))[-i])
        temp <- s0@data[1:(t-1),,,drop=FALSE] # X_{t-1,k} * t(Phi_k^r)

        s1 <- ttl(xx, sphi[-i],c(2:(k+1))[-i])
        temp1 <- s1@data[1:(t-1),,,drop=FALSE] # X_{t-1,k} * t(S_k^{-1}*Phi_k^r)

        L1 <- Reduce("+",lapply(c(1:r)[-j], function (n) ttl(xx[1:(t-1),,], A.new[[n]], c(2:(k+1))))) # additional term
        if (r==1){L1 <- 0}
        L2 <- xx[2:t,,,drop=FALSE] - L1
        temp2 <- L2@data[1:(t-1),,,drop=FALSE]

        RR <- tensor(temp,temp1,c(1:3)[-(i+1)],c(1:3)[-(i+1)])
        LL <- tensor(temp2,temp1,c(1:3)[-(i+1)],c(1:3)[-(i+1)])
        A.new[[j]][[i]] <- LL %*% solve(RR)

        res.old <- xx[2:t,,,drop=FALSE] - Reduce("+",lapply(1:r, function(j) {(ttl(xx[1:(t-1),,], A.new[[j]], (c(1:k) + 1)))}))
        rs <- ttl(res.old, Sig.new.inv[-i], c(2:(k+1))[-i])
        Sig.new[[i]] <- tensor(res.old@data, rs@data, c(1:3)[-(i+1)],c(1:3)[-(i+1)])/(t-1)/prod(dim[-i])
      }
    }
    for (j in c(1:r)){
      a <- c()
      for (i in c(1:k)){
        m <- A.new[[j]][[i]]
        if (i != k){
          a[i] <- svd(m,nu=0,nv=0)$d[1]
          A.new[[j]][[i]] <- m/a[i]
        } else {
          A.new[[j]][[i]] <- m * prod(a)
        }
      }
    }
    b <- c()
    for (i in c(1:k)){
      s <- Sig.new[[i]]
      if (i != k){
        b[i] <- eigen(s)$values[1]
        Sig.new[[i]] <- s/b[i]
      } else {
        Sig.new[[i]] <- s * prod(b)
      }
    }
    phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
    phi.old <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.old[[j]]))}))
    dis <- sqrt(sum((phi.new-phi.old)^2))
    Sig.old <- Sig.new
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  disf <- sum((phi.new-phi)^2)
  res <- (xx[2:t,,,drop=FALSE] - Reduce("+",lapply(1:r, function(j) {(ttl(xx[1:(t-1),,], A.new[[j]], (c(1:k) + 1)))})))@data
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new, SIGMA=Sig.new, niter=iiter, Sig=Sig, res=res, disf=disf))
}

#' MLE under a structured covariance tensor
#'
#' MAR(1) iterative estimation with Kronecker covariance structure: \eqn{X_t = LL * X_{t-1} * RR + E_t} such that \eqn{\Sigma = cov(vec(E_t)) = \Sigma_r \otimes \Sigma_l}.
#'@name MAR1.otimes
#'@rdname MAR1.otimes
#'@aliases MAR1.otimes
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
MAR1.otimes <- function(xx,LL.init=NULL,Sigl.init=NULL,Sigr.init=NULL,niter=50,tol=1e-6,print.true = FALSE){
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
  #-------------
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


#' Stacked vector AR(1) Model
#'
#' vector AR(1) Model.
#'@name var1
#'@rdname var1
#'@aliases var1
#'@export
#'@param xx T * p * q matrix-valued time series
#'@return a list containing the following:\describe{
#'\item{\code{coef}}{coeficient of the fitted VAR(1) model}
#'\item{\code{res}}{residual of the VAR(1) model}
#'}
#'@examples
#'out.var1=var1(xx)
#'sum(out.var1$res**2)
var1 <- function(xx){
  dd=dim(xx)
  T <- dd[1]
  p <- dd[2]
  q <- dd[3]
  yy=apply(xx,MARGIN=1,as.vector)
  yy2=yy[,2:T];yy1=yy[,1:(T-1)]
  out=lm(t(yy2)~0+t(yy1))
  res=array(t(out$res),dim=c(dd[1]-1,dd[2],dd[3]))
  return(list(coef=out$coef, res=res))
}


#' Asymptotic Covariance Matrix of \code{MAR1.otimes}
#'
#' Asymptotic covariance Matrix of \code{MAR1.otimes} for given A, B and matrix-valued time series xx, see Theory 3 in paper.
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
#'out2=MAR1.LS(xx)
#'FnormLL=sqrt(sum(out2$LL))
#'xdim=p;ydim=q
#'out2Xi=MAR.SE(xx.nm,out2$RR*FnormLL,out2$LL/FnormLL,out2$Sig)
#'out2SE=sqrt(diag(out2Xi))
#'SE.A=matrix(out2SE[1:xdim^2],nrow=xdim)
#'SE.B=t(matrix(out2SE[-(1:xdim^2)],nrow=ydim))
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

#' Generate coefficient matrices
#'
#' For test only, to generate coefficient matrices
#'@name generateA
#'@rdname generateA
#'@aliases generateA
#'@export
#'@param dim an array of dimensions of matrices or tensors
#'@param R number of terms
#'@return a list containing coefficient matrices
#'@examples
#' A <- generateA(dim=c(2,2,2), R=1)
#' xx <- generate(c(m1,m2,m3), T=100)
generateA <- function(dim,R){
  K <- length(dim)
  A.true <- lapply(1:R, function(j) {lapply(1:K, function(i) {diag(dim[i])})})
  fnorm <- array(0, c(R,K))
  for (j in c(1:R)){
    for (i in c(1:(K))){
      if (i < K){
        A.true[[j]][[i]] <- matrix(rnorm(dim[i]^2), c(dim[i],dim[i]))
        fnorm[j,i] <- norm(A.true[[j]][[i]],'f')
        A.true[[j]][[i]] <- A.true[[j]][[i]]/fnorm[j,i]
      } else if (i == K){
        A.true[[j]][[i]] <- matrix(rnorm(dim[i]^2), c(dim[i],dim[i]))
        fnorm[j,i] <- norm(A.true[[j]][[i]],'f')
        A.true[[j]][[i]] <- A.true[[j]][[i]]/fnorm[j,i]
      } else {
        stop("sOMETHING WRONG in generate.A")
      }
    }
  }
  A.norm <- c()
  for (j in c(1:R)){
    A.norm[j] <- norm(A.true[[j]][[1]], 'f') * norm(A.true[[j]][[2]], 'f') * norm(A.true[[j]][[3]], 'f')
  }
  order.norm <- order(A.norm, decreasing=TRUE)
  A.temp <- A.true
  for (j in c(1:R)){
    A.true[[j]] <- A.temp[[order.norm[j]]]
  }
  return(A.true)
}


#' Generate an AR(1) tensor time series with given coefficient matrices
#'
#' For test only, with given coefficient matrices and time length t, generate an AR(1) tensor time series.
#'@name generate
#'@rdname generate
#'@aliases generate
#'@export
#'@param dim an array of dimensions of matrices or tensors
#'@param t length of time
#'@param setting the structure of random error, can be specified as "iid", "mle" and "svd", default value is "iid".
#'@return a list containing the following:\describe{
#'\item{\code{A}}{matrices \eqn{A_1,A_2,\cdots,A_k}}
#'\item{\code{X}}{AR(1) tensor time series}
#'}
#'@seealso \code{\link{run.test}}
#'@examples
#' A <- generateA(dim=c(2,2,2), R=1)
#' xx <- generate(c(m1,m2,m3), T=100)
generate <- function(A, t, setting="iid"){
  # to generate time series with given coefficient matrices and time length t
  # return time series xx
  r <- length(A)
  k <- length(A[[1]])
  if (k == 2){
    x <- list(rand_tensor(dim))  # initialize X1
    for (i in c(2:t)){
      e <- new("Tensor", as.integer(k), as.integer(dim), array(rnorm(prod(dim)), dim))
      x[[i]] <-  ttl(x[[i-1]], A, c(1:k)) + e
    }
    return(as.tensor(x))
  } else if (k == 3){
    x <- array(0, c(t,dim))
    x[1,,,] <- array(rnorm(prod(dim)), c(1,dim))
    for (i in c(2:t)){
      if (setting == "iid"){
        e <- array(rnorm(prod(dim)), dim)
      } else if (setting == "mle"){
        e <- as.vector(ttl(as.tensor(array(rnorm(prod(dim)), dim)), Sig.true, c(1:k)))
      } else if (setting == "svd"){
        e <- as.vector(array(mvrnorm(prod(dim), rep(0,prod(dim)), Sigma2), dim))
      } else {
        return("Please specify setting")
      }
      x[i,,,] <-  (Reduce("+", lapply(1:r, function (j) ttl(as.tensor(x[i-1,,,]), A[[j]], c(1:k)))) + e)@data  # sum terms 1:r and add error E, use reduce since ttl returns a list
    }
    return(as.tensor(x))
  }
}


#' Permutation matrix em
#'
#' Permutation matrix em.
#'@name em
#'@rdname em
#'@aliases em
#'@param m,n,i,j set \eqn{m \times n} zero matrix with \eqn{A_{ij} = 1}
#'@return Permutation matrix em such that \eqn{A_{ij} = 1} and other entries equals 0.
#'@seealso \code{\link{trearrange}}
#'@examples
#' A = em(3,3,1,1)
em <- function(m,n,i,j){
  mat <- matrix(0,m,n)
  mat[i,j] <- 1
  return(mat)
}


#' Permutation matrix pm
#'
#' Permutation matrix pm.
#'@name pm
#'@rdname pm
#'@aliases pm
#'@param m an array of dimensions of matrices \eqn{A_1,A_2,\cdots,A_k}
#'@param n length of time
#'@return Permutation matrix pm
#'@seealso \code{\link{trearrange}}
#'@examples
#' P = pm(3,3)
pm <- function(m,n){
  mat <- matrix(0,m*n,m*n)
  for (i in c(1:n)){
    for (j in c(1:m)){
      mat <- mat + kronecker(em(n,m,i,j),t(em(n,m,i,j)))
    }
  }
  return(mat)
}


#' (alpha version) rearrangement operator for tensor.
trearrange <- function(A,m1,m2,m3,n1,n2,n3){
  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2*n3 | m!=m1*m2*m3){
    stop("wrong dimention with your input Phi for rearrangement")
  }
  ans <- divide(A,m1,n1)
  dim <- c(m1*n1,m2*n2,m3*n3)
  t <- new("Tensor", as.integer(3), as.integer(dim), array(0, dim))
  for (i in c(1:m1)){
    for (j in c(1:n1)){
      t@data[(i-1)*n1+j,,] <- mrearrange(ans[[i]][[j]],m2,m3,n2,n3)
    }
  }
  return(t)
}


divide <- function(A,m,n){
  c <- dim(A)[1]/m
  l <- dim(A)[2]/n
  tmp <- lapply(1:m, function(i){
    lapply(1:n, function(j){
      A[((i-1)*c+1):(i*c),((j-1)*l+1):(j*l)]
    })
  })
  return(tmp)
}


mrearrange <- function(A,m1,m2,n1,n2){
  # the inner function of "projection"
  # A: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # A \approx B \otimes C
  # return RA
  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2 | m!=m1*m2){
    print("error m")
    return
  }
  ans <- matrix(NA, m1*n1, m2*n2)
  for(i in 1:m1){
    for(j in 1:n1){
      ans[(j-1)*m1+i,] <- t(as.vector(A[(i-1)*m2+1:m2,(j-1)*n2+1:n2]))
    }
  }
  ans
}


#' Projection Method for Tensor-Valued Time Series
#'
#' TAR(1) one step projection estimation in the model \eqn{X_t = X_{t-1} \times A_1 \times \cdots \times A_K + E_t}.
#'@name TAR1.projection
#'@rdname TAR1.projection
#'@aliases TAR1.projection
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param m1 dim(A1)
#'@param m2 dim(A2)
#'@param m3 dim(A3)
#'@return a list containing the estimation of matrices \eqn{A_1,A_2,\cdots,A_K}
TAR1.projection <- function(xx){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  m1 <- dim[1]; m2 <- dim[2]; m3 <- dim[3]
  n1 <- m1; n2 <- m2; n3 <- m3
  mm <- TAR1.VAR(xx)$coef
  tt <- trearrange(mm,m3,m2,m1,n3,n2,n1)
  cpd <- cp(tt,num_components = 1)
  u1 <- as.numeric(cpd$U[[1]])
  u2 <- as.numeric(cpd$U[[2]])
  u3 <- as.numeric(cpd$U[[3]])
  lam <- cpd$lambdas
  a <- u3/sqrt(sum(u3^2))
  b <- u2/sqrt(sum(u2^2))
  c <- u1*sqrt(sum(u3^2))*sqrt(sum(u2^2))*lam
  A <- list()
  A[[1]] <- list(matrix(a,m1,m1),matrix(b,m2,m2),matrix(c,m3,m3))
  return(list(A=A, disf=disf, res=res))
}


#' Projection Method for Tensor-Valued Time Series with Multiple Terms
#'
#' TAR(1) one step projection estimation in the model \eqn{X_t = \sum_{r=1}^{R} X_{t-1} \times A_1^{(r)} \times \cdots \times A_K^{(r)} + E_t}.
#'@name TAR2.projection
#'@rdname TAR2.projection
#'@aliases TAR2.projection
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param r number of terms
#'@return a list containing the estimation of matrices \eqn{A_1^{(1)},A_2^{(1)},\cdots,A_K^{(R)}}
TAR2.projection <- function(xx,r){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  m1 <- dim[1]; m2 <- dim[2]; m3 <- dim[3]
  n1 <- m1; n2 <- m2; n3 <- m3
  mm <- TAR1.VAR(xx)$coef
  tt <- trearrange(mm,m3,m2,m1,n3,n2,n1)
  cpd <- cp(tt,num_components = r)
  lam <- cpd$lambdas
  A <- list()
  for (j in c(1:r)){  # j is number of terms

    u1 <- cpd$U[[1]][,j]
    u2 <- cpd$U[[2]][,j]
    u3 <- cpd$U[[3]][,j]

    a <- u3/sqrt(sum(u3^2))
    b <- u2/sqrt(sum(u2^2))
    c <- u1*sqrt(sum(u3^2))*sqrt(sum(u2^2))*lam[j]

    A[[j]] <- list(matrix(a,dim[1],dim[1]),
                   matrix(b,dim[2],dim[2]),
                   matrix(c,dim[3],dim[3]))
  }
  return(list(A=A))
}


#' Least Squares Iterative Estimation for Tensor-Valued Time Series
#'
#' Iterated least squares estimation in the model \eqn{X_t = X_{t-1} \times A_1 \times \cdots \times A_K + E_t}.
#'@name TAR1.LS
#'@rdname TAR1.LS
#'@aliases TAR1.LS
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true printe \eqn{A_i}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices \eqn{A_1,A_2,\cdots,A_K}}
#'\item{\code{res}}{residual of the model}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
TAR1.LS <- function(xx,r=1,niter=80,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  A.old <- TAR1.projection(xx)$A
  A.new <- A.old
  dis <- 1
  iiter <- 1
  a <- c()
  while(iiter <= niter & dis >= tol){
    for (i in c(1:k)){
      s1 <- ttl(xx, A.new[[1]][-i], c(2:(k+1))[-i])
      temp <- s1@data[1:(t-1),,,,drop=FALSE]
      RR <- tensor(temp,temp,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
      LL <- tensor(xx@data[2:t,,,,drop=FALSE],temp,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
      A.new[[1]][[i]] <- LL %*% solve(RR)
    }
    for (i in c(1:k)){
      m <- A.new[[1]][[i]]
      if (i != k){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        A.new[[1]][[i]] <- m/a[i]
      } else {
        A.new[[1]][[i]] <- m * prod(a)
      }
    }
    phi.new <- kronecker_list(rev(A.new[[1]]))
    phi.old <- kronecker_list(rev(A.old[[1]]))
    dis <- sqrt(sum((phi.new - phi.old)^2))
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
  disf <- sum((phi.new -phi)^2)
  res <- (xx[2:t,,,,drop=FALSE] - ttl(xx[1:(t-1),,,], A.new[[1]], (c(1:k) + 1)))@data
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new, niter=iiter, Sig=Sig, res=res, disf=disf))
}


#' Least Squares Iterative Estimation for Tensor-Valued Time Series with Multiple Terms
#'
#' Iterated least squares estimation in the model \eqn{X_t = \sum_{r=1}^{R} X_{t-1} \times A_1^{(r)} \times \cdots \times A_K^{(r)} + E_t}.
#'@name TAR2.LS
#'@rdname TAR2.LS
#'@aliases TAR2.LS
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true print \eqn{A_i}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices \eqn{A_1^{(1)},A_2^{(1)},\cdots,A_K^{(R)}}}
#'\item{\code{res}}{residual of the model}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
TAR2.LS <- function(xx,r,niter=80,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  A.old <- TAR2.projection(xx,r)$A
  fnorm <- array(0,c(r,k)) # Rescale Result of PROJ
  for (j in c(1:r)){
    for (i in c(1:k)){
      if (i < k ){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]]/fnorm[j,i]
      } else if (i == k){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]] * prod(fnorm[j,1:(k-1)])
      } else {
        print("WRONG dimension")
      }
    }
  }
  A.new <- A.old
  phi <-  Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
  dis <- 1
  iiter <- 1
  a <- c()
  while(iiter <= niter & dis >= tol & dis <= 1e3){ # stop when dis > 1e3
    for (j in c(1:r)){
      for (i in c(1:k)){
        s0 <- ttl(xx, A.new[[j]][-i], c(2:(k+1))[-i])
        temp <- s0@data[1:(t-1),,,,drop=FALSE]

        L1 <- Reduce("+",lapply(c(1:r)[-j], function(n) {(ttl(xx[1:(t-1),,,], A.new[[n]], (c(1:k) + 1)))}))
        if (r == 1){ L1 <- 0}
        L2 <-  xx[2:t,,,,drop=FALSE] - L1
        temp2 <- L2@data[1:(t-1),,,,drop=FALSE]

        RR <- tensor(temp,temp,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
        LL <- tensor(temp2,temp,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
        A.new[[j]][[i]] <- LL %*% solve(RR)
      }
    }
    for (j in c(1:r)){
      a <- c()
      for (i in c(1:k)){
        m <- A.new[[j]][[i]]
        if (i != k){
          a[i] <- svd(m,nu=0,nv=0)$d[1]
          A.new[[j]][[i]] <- m/a[i]
        } else {
          A.new[[j]][[i]] <- m * prod(a)
        }
      }
    }
    phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
    phi.old <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.old[[j]]))}))
    dis <- sqrt(sum((phi.new - phi.old)^2))
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
  disf <- sum((phi.new -phi)^2)
  res <- (xx[2:t,,,,drop=FALSE] - Reduce("+",lapply(1:r, function(j) {(ttl(xx[1:(t-1),,,], A.new[[j]], (c(1:k) + 1)))})))@data
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new,niter=iiter,Sig=Sig,res=res,disf=disf))
}


#' MLE for Tensor-Valued Time Series with One Term Model Under a Structured Covariance Tensor
#'
#' MLE for the model \eqn{X_t = X_{t-1} \times A_1 \times \cdots \times A_K + E_t}.
#'@name TAR1.MLE
#'@rdname TAR1.MLE
#'@aliases TAR1.MLE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true printe \eqn{A_i}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coeficient matrices \eqn{A_1,A_2,\cdots,A_K}}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
TAR1.MLE <- function(xx, r=1,niter=80,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  Sig.old <- lapply(1:k, function(i) {diag(dim[i])})
  Sig.new <- Sig.old
  A.old <- TAR1.projection(xx)$A
  A.new <- A.old
  dis <- 1
  iiter <- 1
  a <- c()
  b <- c()
  while(iiter <= niter & dis >= tol & dis <= 1e3){ # stop when dis > 1e3
    for (i in c(1:k)){
      Sig.new.inv <- lapply(1:k, function (i) {solve(Sig.new[[i]])})
      sphi <-  lapply(1:k, function (i) {Sig.new.inv[[i]] %*% (A.new[[1]][[i]])})

      s0 <- ttl(xx, A.new[[1]][-i], c(2:(k+1))[-i])
      temp <- s0@data[1:(t-1),,,,drop=FALSE] # X_{t-1,k} * t(Phi_k^r)

      s1 <- ttl(xx, sphi[-i], c(2:(k+1))[-i])
      temp1 <- s1@data[1:(t-1),,,,drop=FALSE] # X_{t-1,k} * t(S_k^{-1}*Phi_k^r)

      L2 <- xx[2:t,,,,drop=FALSE]
      temp2 <- L2@data[1:(t-1),,,,drop=FALSE]

      RR <- tensor(temp,temp1,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
      LL <- tensor(temp2, temp1, c(1:4)[-(i+1)],c(1:4)[-(i+1)])
      A.new[[1]][[i]] <- LL %*% solve(RR)

      res.old <- xx[2:t,,,,drop=FALSE] - ttl(xx[1:(t-1),,,], A.new[[1]], (c(1:k) + 1))
      rs <- ttl(res.old, Sig.new.inv[-i], c(2:(k+1))[-i])
      Sig.new[[i]] <- tensor(res.old@data, rs@data, c(1:4)[-(i+1)],c(1:4)[-(i+1)])/(t-1)/prod(dim[-i])
    }
    for (i in c(1:k)){
      m <- A.new[[1]][[i]]
      s <- Sig.new[[i]]
      if (i != k){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        b[i] <- eigen(s)$values[1]
        A.new[[1]][[i]] <- m/a[i]
        Sig.new[[i]] <- s/b[i]
      } else {
        A.new[[1]][[i]] <- m * prod(a)
        Sig.new[[i]] <- s * prod(b)
      }
    }
    phi.new <- kronecker_list(rev(A.new[[1]]))
    phi.old <- kronecker_list(rev(A.old[[1]]))
    dis <- sqrt(sum((phi.new -phi.old)^2))
    Sig.old <- Sig.new
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
  disf <- sum((phi.new-phi)^2)
  res <- (xx[2:t,,,,drop=FALSE] - ttl(xx[1:(t-1),,,], A.new[[1]], (c(1:k) + 1)))@data
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new, SIGMA=Sig.new, niter=iiter, Sig=Sig, res=res, disf=disf))
}


#' MLE for Tensor-Valued Time Series with Multiple Terms Model Under a Structured Covariance Tensor
#'
#' MLE for the model \eqn{X_t = \sum_{r=1}^{R} X_{t-1} \times A_1^{(r)} \times \cdots \times A_K^{(r)} + E_t}.
#'@name TAR2.MLE
#'@rdname TAR2.MLE
#'@aliases TAR1.MLE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param niter maximum number of iterations if error stays above \code{tol}
#'@param tol relative Frobenius norm error tolerance
#'@param print.true print \eqn{A_i}
#'@return a list containing the following:\describe{
#'\item{\code{A}}{estimator of coefficient matrices \eqn{A_1,A_2,\cdots,A_K}}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
TAR2.MLE <- function(xx,r,niter=200,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  Sig.old <- lapply(1:k, function(i) {diag(dim[i])})
  Sig.new <- Sig.old
  A.old <- TAR2.projection(xx,r)$A
  fnorm <- array(0,c(r,k))
  for (j in c(1:r)){
    for (i in c(1:k)){
      if (i < k ){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]]/fnorm[j,i]
      } else if (i == k){
        fnorm[j,i] <- norm(A.old[[j]][[i]],"f")
        A.old[[j]][[i]] <- A.old[[j]][[i]] * prod(fnorm[j,1:(k-1)])
      } else {
        print("WRONG dimension")
      }
    }
  }
  A.new <- A.old
  dis <- 1
  iiter <- 1
  while(iiter <= niter & dis >= tol & dis <= 1e3){ # stop when dis > 1e3

    for (j in c(1:r)){
      for (i in c(1:k)){
        Sig.new.inv <- lapply(1:k, function (i) {solve(Sig.new[[i]])})
        sphi <-  lapply(1:k, function (i) {Sig.new.inv[[i]] %*% (A.new[[j]][[i]])})

        s0 <- ttl(xx, A.new[[j]][-i], c(2:(k+1))[-i])
        temp <- s0@data[1:(t-1),,,,drop=FALSE] # X_{t-1,k} * t(Phi_k^r)

        s1 <- ttl(xx, sphi[-i],c(2:(k+1))[-i])
        temp1 <- s1@data[1:(t-1),,,,drop=FALSE] # X_{t-1,k} * t(S_k^{-1}*Phi_k^r)

        L1 <- Reduce("+",lapply(c(1:r)[-j], function (n) ttl(xx[1:(t-1),,,], A.new[[n]], c(2:(k+1))))) # additional term
        if (r==1){L1 <- 0}
        L2 <- xx[2:t,,,,drop=FALSE] - L1
        temp2 <- L2@data[1:(t-1),,,,drop=FALSE]

        RR <- tensor(temp,temp1,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
        LL <- tensor(temp2,temp1,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
        A.new[[j]][[i]] <- LL %*% solve(RR)

        res.old <- xx[2:t,,,,drop=FALSE] - Reduce("+",lapply(1:r, function(j) {(ttl(xx[1:(t-1),,,], A.new[[j]], (c(1:k) + 1)))}))
        rs <- ttl(res.old, Sig.new.inv[-i], c(2:(k+1))[-i])
        Sig.new[[i]] <- tensor(res.old@data, rs@data, c(1:4)[-(i+1)],c(1:4)[-(i+1)])/(t-1)/prod(dim[-i])
      }
    }

    for (j in c(1:r)){
      a <- c()
      for (i in c(1:k)){
        m <- A.new[[j]][[i]]
        if (i != k){
          a[i] <- svd(m,nu=0,nv=0)$d[1]
          A.new[[j]][[i]] <- m/a[i]
        } else {
          A.new[[j]][[i]] <- m * prod(a)
        }
      }
    }

    b <- c()
    for (i in c(1:k)){
      s <- Sig.new[[i]]
      if (i != k){
        b[i] <- eigen(s)$values[1]
        Sig.new[[i]] <- s/b[i]
      } else {
        Sig.new[[i]] <- s * prod(b)
      }
    }

    phi.new <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.new[[j]]))}))
    phi.old <- Reduce("+", lapply(1:r, function(j) {kronecker_list(rev(A.old[[j]]))}))
    dis <- sqrt(sum((phi.new-phi.old)^2))
    Sig.old <- Sig.new
    A.old <- A.new
    iiter <- iiter + 1
    if (print.true == TRUE){
      print(dis)
      print(paste('iiter num=',iiter))
    }
  }
  disf <- sum((phi.new-phi)^2)
  res <- (xx[2:t,,,,drop=FALSE] - Reduce("+",lapply(1:r, function(j) {(ttl(xx[1:(t-1),,,], A.new[[j]], (c(1:k) + 1)))})))@data
  Sig <- matrix(tensor(res,res,1,1),prod(dim))/(t-1)
  return(list(A=A.new, SIGMA=Sig.new, niter=iiter, Sig=Sig, res=res, disf=disf))
}


#' Stacked vector AR(1) Model for Tensor-Valued Time Series
#'
#' vector AR(1) Model.
#'@name TAR1.VAR
#'@rdname TAR1.VAR
#'@aliases TAR1.VAR
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@return a list containing the following:\describe{
#'\item{\code{coef}}{coeficient of the fitted VAR(1) model}
#'\item{\code{res}}{residual of the VAR(1) model}
#'}
#'@examples
#' A <- generateA(dim=c(2,2,2), R=1)
#' xx <- generate(c(m1,m2,m3), T=100)
#' SIGMA <- TAR1.VAR(xx)
TAR1.VAR <- function(xx){
  dd=xx@modes
  T <- dd[1]
  d1 <- dd[2]
  d2 <- dd[3]
  d3 <- dd[4]
  yy=apply(xx@data,MARGIN=1,as.vector)
  yy2=yy[,2:T];yy1=yy[,1:(T-1)]
  out=lm(t(yy2)~0+t(yy1))
  # res=array(t(out$res),dim=c(dd[1]-1,dd[2],dd[3])) ## should NOT transpose!
  res=array((out$res),dim=c(dd[1]-1,dd[2],dd[3],dd[4]))
  return(list(coef=out$coef, res=res))
}


#' Asymptotic Covariance Matrix of \code{TAR1.LS}
#'
#' Asymptotic covariance Matrix of \code{TAR1.LS} for given A tensor-valued time series xx, see related Theorems in our paper.
#'@name TAR1.SE.LSE
#'@rdname TAR1.SE.LSE
#'@aliases TAR1.SE.LSE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param A.true coefficient matrices in TAR(1) model
#'@param Sigma covariance matrix \eqn{cov(vec(E_t))} in TAR(1) model
#'@return asymptotic covariance matrix
#'@examples
#' dim <- c(2,2,2)
#' A <- generateA(dim, R=1)
#' xx <- generate(dim, T=100)
#' SIGMA <- TAR1.SE.LSE(xx, A, Sigma=diag(prod(dim)))
TAR1.SE.LSE <- function(xx, A.true, Sigma){
  dim <- xx@modes[-1]
  m1 <- dim[1]
  m2 <- dim[2]
  m3 <- dim[3]
  k <- length(dim)
  T <- xx@modes[[1]]
  A <- A.true[[1]][[1]]
  B <- A.true[[1]][[2]]
  C <- A.true[[1]][[3]]
  a <- as.matrix(as.vector(A))
  b <- as.matrix(as.vector(B))
  c <- as.matrix(as.vector(C))
  r1 <- rbind(a,matrix(0,m2^2+m3^2,1))
  r2 <- rbind(matrix(0,m1^2,1),b,matrix(0,m3^2,1))
  Gamma <-  r1 %*% t(r1) +r2 %*% t(r2)

  Hdim <- m1^2+m2^2+m3^2
  HT <- array(1,c(T,Hdim,Hdim))
  WT <- array(1,c(T,Hdim,prod(dim))) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  for (t in c(1:T)){
    w1 <- k_unfold(xx[t,,,], m = 1)@data %*% t(kronecker(C,B))
    w2 <- k_unfold(xx[t,,,], m = 2)@data %*% t(kronecker(C,A))
    w3 <- k_unfold(xx[t,,,], m = 3)@data %*% t(kronecker(B,A))
    w <- rbind(kronecker(w1,diag(m1)) ,kronecker(w2,diag(m2)) %*% kronecker(diag(m3),pm(m2,m1)), kronecker(w3,diag(m3)) %*% pm(m3,m2*m1))
    WT[t,,] <-  w
  }
  EWWt <- tensor(WT,WT,c(3,1),c(3,1))/T #(m1^2+m2^2+m^3)*(m1^2+m2^2+m^3)
  WSigma <- tensor(WT,Sigma,3,1) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  EWSigmaWt <- tensor(WSigma,WT,c(3,1),c(3,1))/T
  H <- EWWt + Gamma
  Hinv <- solve(H)
  Xi <- Hinv %*% EWSigmaWt %*% Hinv
  return(Xi)
}


#' Asymptotic Covariance Matrix of \code{TAR1.MLE}
#'
#' Asymptotic covariance Matrix of \code{TAR1.MLE} for given A tensor-valued time series xx, see related Theorems in our paper.
#'@name TAR1.SE.MLE
#'@rdname TAR1.SE.MLE
#'@aliases TAR1.SE.MLE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param A.true coefficient matrices in TAR(1) model
#'@param Sigma covariance matrix cov(vec(E_t)) in TAR(1) model
#'@return asmptotic covariance matrix
#'@examples
#' dim <- c(2,2,2)
#' A <- generateA(dim, R=1)
#' xx <- generate(dim, T=100)
#' SIGMA <- TAR2.SE.MLE(xx, A, Sigma=diag(prod(dim)))
TAR1.SE.MLE <- function(xx, A.true, Sigma){
  dim <- xx@modes[-1]
  m1 <- dim[1]
  m2 <- dim[2]
  m3 <- dim[3]
  k <- length(dim)
  T <- xx@modes[[1]]
  A <- A.true[[1]][[1]]
  B <- A.true[[1]][[2]]
  C <- A.true[[1]][[3]]
  a <- as.matrix(as.vector(A))
  b <- as.matrix(as.vector(B))
  c <- as.matrix(as.vector(C))
  r1 <- rbind(a,matrix(0,m2^2+m3^2,1))
  r2 <- rbind(matrix(0,m1^2,1),b,matrix(0,m3^2,1))
  Gamma <-  r1 %*% t(r1) +r2 %*% t(r2)
  Hdim <- m1^2+m2^2+m3^2
  HT <- array(1,c(T,Hdim,Hdim))
  WT <- array(1,c(T,Hdim,prod(dim))) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  for (t in c(1:T)){
    w1 <- k_unfold(xx[t,,,], m = 1)@data %*% t(kronecker(C,B))
    w2 <- k_unfold(xx[t,,,], m = 2)@data %*% t(kronecker(C,A))
    w3 <- k_unfold(xx[t,,,], m = 3)@data %*% t(kronecker(B,A))
    w <- rbind(kronecker(w1,diag(m1)) ,kronecker(w2,diag(m2)) %*% kronecker(diag(m3),pm(m2,m1)), kronecker(w3,diag(m3)) %*% pm(m3,m2*m1))
    WT[t,,] <-  w
  }
  WSigma <- tensor(WT,solve(Sigma),3,1) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  EWSigmaWt <- tensor(WSigma,WT,c(3,1),c(3,1))/T
  H <- EWSigmaWt + Gamma
  Hinv <- solve(H)
  Xi <- Hinv %*% EWSigmaWt %*% Hinv
  return(Xi)
}


#' Asymptotic Covariance Matrix of \code{TAR2.LSE}
#'
#' Asymptotic covariance Matrix of \code{TAR2.LSE} for given A tensor-valued time series xx, see related Theorems in our paper.
#'@name TAR2.SE.LSE
#'@rdname TAR2.SE.LSE
#'@aliases TAR2.SE.LSE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param A.true coefficient matrices in TAR(1) model
#'@param Sigma covariance matrix cov(vec(E_t)) in TAR(1) model
#'@return asmptotic covariance matrix
#'@examples
#' dim <- c(2,2,2)
#' A <- generateA(dim, R=1)
#' xx <- generate(dim, T=100)
#' SIGMA <- TAR2.SE.LSE(xx, A, Sigma=diag(prod(dim)))
TAR2.SE.LSE <- function(xx, A.true, Sigma){
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


#' Asymptotic Covariance Matrix of \code{TAR2.MLE}
#'
#' Asymptotic covariance Matrix of \code{TAR2.MLE} for given A tensor-valued time series xx, see related Theorems in our paper.
#'@name TAR2.SE.MLE
#'@rdname TAR2.SE.MLE
#'@aliases TAR2.SE.MLE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param A.true coefficient matrices in TAR(1) model
#'@param Sigma covariance matrix cov(vec(E_t)) in TAR(1) model
#'@return asmptotic covariance matrix
#'@examples
#' dim <- c(2,2,2)
#' A <- generateA(dim, R=1)
#' xx <- generate(dim, T=100)
#' SIGMA <- TAR2.SE.MLE(xx, A, Sigma=diag(prod(dim)))
TAR2.SE.MLE <- function(xx, A.true, Sigma){
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


#' Test Function
#'
#' For test only
#'@name run.test
#'@rdname run.test
#'@aliases run.test
run.test <- function(m1,m2,m3,n=100,T){
  n1 <- m1
  n2 <- m2
  n3 <- m3
  err1 <- c()
  err12 <- c()
  err2 <- c()
  err22 <- c()
  err3 <- c()
  for (i in c(1:n)){
    tic("-----------------------------------complete one simulation")
    xx <- generate(c(m1,m2,m3),T)
    lse <- TAR1.LS(xx)
    pro <- TAR1.projection(xx)
    var <- TAR1.VAR(xx)
    err1 <- c(err1, log(lse[[4]]))
    err12 <- c(err12, log(lse[[5]]))
    err2 <- c(err2, log(pro[[1]]))
    err22 <- c(err22,log(pro[[2]]))
    err3 <- c(err3, log(var))
    print(paste("-----------------------------------complete",i))
    print(paste("-----------------------------------dim",m1))
    toc()
  }
  return(list(err1,err2,err3,err12,err22))
}

# devtools::load_all()
# devtools::document()
# devtools::build_manual()
