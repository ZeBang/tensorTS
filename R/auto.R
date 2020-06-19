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

#' Generate an AR(1) tensor time series
#'
#' For test only, randomly generate matrices A1,A2,..Ak by iid standard gaussians with given dimensions, and then generate an AR(1) tensor time series.
#'@name generate
#'@rdname generate
#'@aliases generate
#'@param dim an array of dimentions of matrices \eqn{A_1,A_2,\cdots,A_k}
#'@param t length of time
#'@return a list containing the following:\describe{
#'\item{\code{A}}{matrices \eqn{A_1,A_2,\cdots,A_k}}
#'\item{\code{X}}{AR(1) tensor time series}
#'}
#'#'@seealso \code{\link{run.test}}
#'@examples
#' dim <- c(1,2,3)
#' T <- 100
#' xx <- generate(c(m1,m2,m3),T)
generate <- function(dim,t){
  # to generate matrices A1,A2,..Ak with given dimensions
  # return A1,A2,...Ak
  k <- length(dim)
  A <- list()
  for (i in c(1:k)){
    ma <- matrix(rnorm(dim[i]^2), c(dim[i],dim[i]))
    ma <- ma/norm(ma, type = "F")
    A <- append(A,list(ma))
  }
  X <- list(rand_tensor(dim))  # initialize X1
  for (i in c(2:t)){
    E <- new("Tensor", as.integer(k), as.integer(dim), array(rnorm(prod(dim)), dim))
    X[[i]] <-  ttl(X[[i-1]], A, c(1:k)) + E
  }
  return(list(A,X))
}

#' Permutation matrix em
#'
#' Permutation matrix em.
#'@name em
#'@rdname em
#'@aliases em
#'@param m
#'@param n
#'@param i
#'@param j
#'@return Permutation matrix em
#'#'@seealso \code{\link{trearrange}}
#'@examples
#' em(m,n,i,j)
em <- function(m,n,i,j){
  mat <- matrix(0,m,n)
  mat[i,j] <- 1
  return(mat)
}

#' Permutation matrix PM
#'
#' Permutation matrix PM.
#'@name PM
#'@rdname PM
#'@aliases PM
#'@param m an array of dimentions of matrices \eqn{A_1,A_2,\cdots,A_k}
#'@param n length of time
#'@return Permutation matrix PM
#'#'@seealso \code{\link{trearrange}}
#'@examples
#' PM(m,n)
PM <- function(m,n){
  mat <- matrix(0,m*n,m*n)
  for (i in c(1:n)){
    for (j in c(1:m)){
      mat <- mat + kronecker(em(n,m,i,j),t(em(n,m,i,j)))
    }
  }
  return(mat)
}


#' trearrange
#'
#' (alpha version) rearrangement operator for tensor.
#'@name trearrange
#'@rdname trearrange
#'@aliases trearrange
#'@export
trearrange <- function(A,m1,m2,m3,n1,n2,n3){
  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2*n3 | m!=m1*m2*m3){
    print("error t")
    return
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
#'@name projection
#'@rdname projection
#'@aliases projection
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param m1 dim(A1)
#'@param m2 dim(A2)
#'@param m3 dim(A3)
#'@return a list containing the estimation of matrices \eqn{A_1,A_2,\cdots,A_K}
projection <- function(xx,m1,m2,m3,n1,n2,n3){
  mm <- varlse(xx)
  tt <- trearrange(mm,m1,m2,m3,n1,n2,n3)
  cpd <- cp(tt,num_components = 1)
  a1 <- as.numeric(cpd$U[[1]])
  a2 <- as.numeric(cpd$U[[2]])
  a3 <- as.numeric(cpd$U[[3]])
  lam <- (cpd$lambdas)^(1/3)
  A.new <- list(lam*matrix(a1,m1,m1),lam*matrix(a2,m2,m2),lam*matrix(a3,m3,m3))
  return(A.new)
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
#'\item{\code{A}}{estimator of coeficient matrices \eqn{A_1,A_2,\cdots,A_K}}
#'\item{\code{res}}{residual of the MAR(1)}
#'\item{\code{Sig}}{covariance matrix cov(vec(E_t))}
#'\item{\code{niter}}{number of iterations}
#'}
TAR1.LS <- function(xx,niter=1000,tol=1e-6,print.true = FALSE){
  dim <- xx@modes[-1]
  k <- length(dim)
  t <- xx@modes[[1]]
  A.old <- projection(xx,dim[1],dim[2],dim[3],dim[1],dim[2],dim[3])
  # A.old <- lapply(1:k, function(i){diag(dim[i])})
  dis <- 1
  iiter <- 1
  a <- c()
  #-------------
  while(iiter <= niter & dis >= tol){
    A.new <- lapply(1:k, function(i) {
      s1 <- ttl(xx, A.old[-i], c(2:(k+1))[-i])
      temp <- s1@data[1:(t-1),,,,drop=FALSE]
      r <- tensor(temp,temp,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
      l <- tensor(xx@data[2:t,,,,drop=FALSE],temp,c(1:4)[-(i+1)],c(1:4)[-(i+1)])
      m <- l %*% solve(r)
      if (i != k){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        m/a[i]
      } else {
        m * prod(a)
      }
    })
    dis <- sqrt(sum((kronecker_list(A.new)-kronecker_list(A.old))^2))
    # dis2 <- sqrt(sum((kronecker_list(A.new)-kronecker_list(A))^2))
    A.old <- A.new
    iiter <- iiter + 1
    print(dis)
    print(paste('iiter num=',iiter))
  }
  # disf <- sum((kronecker_list(A.new)-kronecker_list(A))^2)
  res <- xx[2:t,,,,drop=FALSE] - ttl(xx[1:(t-1),,,], A.new, (c(1:k) + 1))
  Sig <- matrix(tensor(res@data,res@data,1,1),prod(dim))/(t-1)
  # disf2 <- sum(min((A.new[[1]] - A[[1]])^2,(-A.new[[1]] - A[[1]])^2)) +sum(min((A.new[[2]]-A[[2]])^2,(-A.new[[2]] - A[[2]])^2)) + sum(min((A.new[[3]] - A[[3]])^2,(-A.new[[3]] - A[[3]])^2))
  return(list(A=A.new,niter=iiter,Sig=Sig,res=res))
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
#'
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
#' Asymptotic covariance Matrix of \code{TAR1.LS} for given A tensor-valued time series xx, see Theory 2 in paper.
#'@name TAR.SE
#'@rdname TAR.SE
#'@aliases TAR.SE
#'@export
#'@param xx  \eqn{T * m_1 * \cdots * m_K} tensor-valued time series
#'@param A m1 by m1 matrix in MAR(1) model
#'@param B m2 by m2 matrix in MAR(1) model
#'@param C m3 by m3 matrix in MAR(1) model
#'@param Sig covariance matrix cov(vec(E_t)) in TAR(1) model
#'@return asmptotic covariance matrix
#'@examples
TAR.SE <- function(xx, C, B, A, Sigma){
  dim <- xx@modes[-1]
  m1 <- dim[1]
  m2 <- dim[2]
  m3 <- dim[3]
  k <- length(dim)
  T <- xx@modes[[1]]
  a <- as.matrix(as.vector(A))
  b <- as.matrix(as.vector(B))
  c <- as.matrix(as.vector(C))
  r1 <- rbind(a,matrix(0,m2^2+m3^2,1))
  r2 <- rbind(matrix(0,m1^2,1),b,matrix(0,m3^2,1))
  r <-  r1 %*% t(r1) +r2 %*% t(r2)
  Hdim <- m1^2+m2^2+m3^2
  HT <- array(1,c(T,Hdim,Hdim))
  WT <- array(1,c(T,Hdim,prod(dim))) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  for (t in c(1:T)){
    w1 <- k_unfold(xx[t,,,], m = 1)@data %*% t(kronecker(C,B))
    w2 <- k_unfold(xx[t,,,], m = 2)@data %*% t(kronecker(C,A))
    w3 <- k_unfold(xx[t,,,], m = 3)@data %*% t(kronecker(B,A))
    w <- rbind(kronecker(w1,diag(m1)) ,kronecker(w2,diag(m2)) %*% kronecker(diag(m3),PM(m2,m1)), kronecker(w3,diag(m3)) %*% PM(m3,m2*m1))
    WT[t,,] <-  w
  }
  EWWt <- tensor(WT,WT,c(3,1),c(3,1))/T #(m1^2+m2^2+m^3)*(m1^2+m2^2+m^3)
  WSigma <- tensor(WT,Sigma,3,1) #T*(m1^2+m2^2+m^3)*(m1m2m3)
  EWSigmaWt <- tensor(WSigma,WT,c(3,1),c(3,1))/T
  H <- EWWt + r
  Hinv <- solve(H)
  Xi <- Hinv %*% EWSigmaWt %*% Hinv
  Xi <- Xi/T
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
    pro <- projection(xx,m1,m2,m3,n1,n2,n3)
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
