###Functions of autoregressive Models

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
#'#'@seealso \code{\link{MAR1.projection}}
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
#'@param m1 \code{ncol} of A
#'@param m2 \code{ncol} of B
#'@param n1 \code{nrow} of A
#'@param n2 \code{nrow} of B
#'@return a list contaning two estimator (matrix)
#'@seealso \code{\link{MAR1.projection}}
#'@examples
#'A <- matrix(runif(6),ncol=2),
#'projection(A,3,3,2,2)
projection <- function(A,m1,m2,n1,n2){
  # the inner function of MAR1.projection
  # A: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # A \approx B \otimes C
  # return B and C
  RA <- rearrange(A,m1,m2,n1,n2)
  RA.svd <- svd(RA,nu=1,nv=1)
  B <- matrix(RA.svd$u * RA.svd$d[1], m1, n1)
  C <- matrix(RA.svd$v , m2, n2)
  list(B=B,C=C)
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
  ans.projection <- projection(kroneck, q,p,q,p)
  a <- svd(ans.projection$C,nu=0,nv=0)$d[1]
  LL <- ans.projection$C / a
  RR <- t(ans.projection$B) * a
  res=xx[2:T,,,drop=FALSE] - aperm(tensor(tensor(xx[1:(T-1),,,drop=FALSE],RR,3,1),LL,2,2),c(1,3,2))
  Sig <- matrix(tensor(res,res,1,1),p*q)/(T-1)
  return(list(LL=LL,RR=RR,res=res,Sig=Sig))
}

#' Least Squares Iterative Estimation
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
#'@param print.true printe LL and RR
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
#'@return asmptotic covariance matrix
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


