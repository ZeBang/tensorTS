### Helper functions

# myslice <- function(xx, K, start, end){
#   if (K==2){
#     return(xx[start:end,,,drop=FALSE])
#   } else if (K==3){
#     return(xx[start:end,,,,drop=FALSE])
#   } else {
#     stop("not support tensor mode K > 3")
#   }
# }

# mat projection
matAR.PROJ <- function(xx, dim, r, t){
  xx.mat <- matrix(xx,t,dim[1]*dim[2])
  # kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% solve(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),])
  kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% chol2inv(chol(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),]))
  return(projection(kroneck, r, dim[1],dim[2],dim[1],dim[2]))
}

# Tensor Times List
tl <- function(x, list_mat, k = NULL){
  if (is.null(k)){
    tensor(tensor(tensor(x, list_mat[[1]], 2, 2), list_mat[[2]], 2, 2), list_mat[[3]], 2, 2)
  } else if (k == 1){
    tensor(tensor(x, list_mat[[1]], 3, 2), list_mat[[2]], 3, 2)
  } else if (k == 2){
    aperm(tensor(tensor(x, list_mat[[1]], 2, 2), list_mat[[2]], 3, 2),c(1,3,2,4))
  } else if (k == 3){
    aperm(tensor(tensor(x, list_mat[[1]], 2, 2), list_mat[[2]], 2, 2),c(1,3,4,2))
  } else {
    stop("not support tensor mode K > 3")
  }

}

# standard error extraction
covtosd <- function(cov, dim, R){
  K <- length(dim)
  P <- length(R)
  sd = list()
  for (p in c(1:P)){
    if (is.na(R[p])) stop("p != length(R)")
    if (R[p] == 0) next
    sd[[p]] <- lapply(1:R[p], function(j) {lapply(1:K, function(i) {list()})})
  }
  for (i in c(1:P)){
    for (j in c(1:R[i])){
      for (k in c(1:K)){
        left <- sum(dim^2)*sum(R[0:(i-1)]) + sum(dim^2)*(j-1) + sum((dim^2)[0:(k-1)]) + 1
        right <- sum(dim^2)*sum(R[0:(i-1)]) + sum(dim^2)*(j-1) + sum((dim^2)[0:k])
        sd[[i]][[j]][[k]] <- array(sqrt(diag(cov)[left:right]), c(dim[k], dim[k]))
      }
    }
  }
  return(sd)
}

# Permutation matrix em
em <- function(m,n,i,j){
  ## m,n,i,j set \eqn{m \times n} zero matrix with \eqn{A_{ij} = 1}
  ## return: Permutation matrix em such that \eqn{A_{ij} = 1} and other entries equals 0.
  mat <- matrix(0,m,n)
  mat[i,j] <- 1
  return(mat)
}

# Permutation matrix pm
pm <- function(m,n){
  ## m: an array of dimensions of matrices \eqn{A_1,A_2,\cdots,A_k}
  ## n: length of time
  ## return: Permutation matrix pm
  mat <- matrix(0,m*n,m*n)
  for (i in c(1:n)){
    for (j in c(1:m)){
      mat <- mat + kronecker(em(n,m,i,j),t(em(n,m,i,j)))
    }
  }
  return(mat)
}

# rearrangement operator for tensor
trearrange <- function(A,dim){
  m1 = dim[1]; m2 = dim[2]; m3 = dim[3]
  n1 = m1; n2 = m2; n3 = m3
  m <- nrow(A)
  n <- ncol(A)
  if(n!=n1*n2*n3 | m!=m1*m2*m3){
    stop("wrong dimension with your input Phi for rearrangement")
  }
  ans <- divide(A,m1,n1)
  dim <- c(m1*n1,m2*n2,m3*n3)
  t <- array(0, dim)
  for (i in c(1:m1)){
    for (j in c(1:n1)){
      t[(j-1)*m1+i,,] <- mrearrange(ans[[i]][[j]],m2,m3,n2,n3)
    }
  }
  return(t)
}


divide <- function(A,m,n){
  # the inner function of "trearrange"
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
    stop("error m")
  }
  ans <- matrix(NA_real_, m1*n1, m2*n2)
  for(i in 1:m1){
    for(j in 1:n1){
      ans[(j-1)*m1+i,] <- t(as.vector(A[(i-1)*m2+1:m2,(j-1)*n2+1:n2]))
    }
  }
  return(ans)
}

projection <- function(M,r,m1,m2,n1,n2){
  # the inner function of MAR1.projection
  # M: m1m2*n1n2
  # B: m1*n1
  # C: m2*n2
  # M \approx B \otimes C
  # return B and C
  RA <- mrearrange(M,m1,m2,n1,n2)
  RA.svd <- svd(RA,nu=r,nv=r)
  A <- list()
  for (i in c(1:r)){
    A[[i]] <- list(matrix(RA.svd$v[,i] * RA.svd$d[i], m2, n2), matrix(RA.svd$u[,i], m1, n1))
  }
  for (j in c(1:r)){
    A[[j]] <- rev(A[[j]])
    a <- c()
    for (i in c(1:2)){
      m <- A[[j]][[i]]
      if (i != 2){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        A[[j]][[i]] <- m/a[i]
      } else {
        A[[j]][[i]] <- m * prod(a)
      }
    }
  }
  return(A)
}

ten.proj <- function(tt, dim, R){
  ## inner func of "TenAR.proj"
  cpd <- rTensor::cp(rTensor::as.tensor(tt), num_components = R, max_iter = 100, tol = 1e-06)
  lam <- cpd$lambdas
  A.proj <- list()
  for (j in c(1:R)){
    u1 <- cpd$U[[1]][,j]
    u2 <- cpd$U[[2]][,j]
    u3 <- cpd$U[[3]][,j]
    f1 <- sqrt(sum(cpd$U[[1]][,j]^2))
    f2 <- sqrt(sum(cpd$U[[2]][,j]^2))
    f3 <- sqrt(sum(cpd$U[[3]][,j]^2))
    a1 <- u1/f1
    a2 <- u2/f2
    a3 <- u3*f1*f2*lam[j]
    A.proj[[j]] <- list(matrix(a1,dim[1],dim[1]),
                        matrix(a2,dim[2],dim[2]),
                        matrix(a3,dim[3],dim[3]))
  }
  return(fro.order(A.proj))
}


fro.rescale <- function(A){
  r <- length(A)
  k <- length(A[[1]])
  for (j in c(1:r)){
    a <- c()
    for (i in c(1:k)){
      m <- A[[j]][[i]]
      if (i < k ){
        a[i] <- norm(m,"f")
        A[[j]][[i]] <- m/a[i]
      } else if (i == k){
        A[[j]][[i]] <- m * prod(a)
      } else {
        stop("WRONG dimension")
      }
    }
  }
  return(A)
}

svd.rescale <- function(A){
  r <- length(A)
  k <- length(A[[1]])
  for (j in c(1:r)){
    a <- c()
    for (i in c(1:k)){
      m <- A[[j]][[i]]
      if (i < k ){
        a[i] <- svd(m,nu=0,nv=0)$d[1]
        A[[j]][[i]] <- m/a[i]
      } else if (i == k){
        A[[j]][[i]] <- m * prod(a)
      } else {
        stop("WRONG dimension")
      }
    }
  }
  return(A)
}

eigen.rescale <- function(A){
  r <- length(A)
  k <- length(A[[1]])
  for (j in c(1:r)){
    a <- c()
    for (i in c(1:k)){
      m <- A[[j]][[i]]
      if (i < k ){
        a[i] <- eigen(m)$values[1]
        A[[j]][[i]] <- m/a[i]
      } else if (i == k){
        A[[j]][[i]] <- m * prod(a)
      } else {
        stop("WRONG dimension")
      }
    }
  }
  return(A)
}

fro.order <- function(A){
  R <- length(A)
  K <- length(A[[1]])
  if (R == 1){return(A)}
  A.norm <- c()
  for (j in c(1:R)){
    A.norm[j] <- Reduce("*",lapply(c(1:K), function(k) { norm(A[[j]][[k]], 'f')}))
  }
  order.norm <- order(A.norm, decreasing=TRUE)
  A.temp <- A
  for (j in c(1:R)){
    A[[j]] <- A.temp[[order.norm[j]]]
  }
  return(A)
}

ten.dis.A <- function(A, B, R, K){
  P = length(R)
  dis <- 0
  for (p in c(1:P)){
    if (R[p] == 0) next
    for (r in c(1:R[p])){
      for (k in c(1:K)){
        dis <- dis + min(sum((A[[p]][[r]][[k]] - B[[p]][[r]][[k]])^2), sum((A[[p]][[r]][[k]] + B[[p]][[r]][[k]])^2))
      }
    }
  }
  return(sqrt(dis))
}

ten.dis.phi <- function(phi.A, phi.B){
  P <- length(phi.A)
  dis <- 0
  for (i in c(1:P)){
    dis <- dis + sqrt(sum((phi.A[[i]] - phi.B[[i]])^2))
  }
  return(dis)
}


ten.res <- function(subxx,A,P,R,K,ttlMode){
  L1 = 0
  for (l in c(1:P)){
    if (R[l] == 0) next
    L1 <- L1 + Reduce("+",lapply(c(1:R[l]), function(n) {rTensor::ttl(subxx[[l+1]], A[[l]][[n]], ttlMode)}))
  }
  res <- subxx[[1]] - L1
  return(res)
}


M.eigen <- function(A, R, P, dim){
  phi <- list()
  PP = P
  for (i in c(1:P)){
    if (sum(R[i:length(R)]) == 0){
      PP = i-1
      break
    }
    if (R[i] == 0){
      phi[[i]] = pracma::zeros(prod(dim))
    } else {
      phi[[i]] <- Reduce("+", lapply(1:R[i], function(j) {rTensor::kronecker_list(rev(A[[i]][[j]]))}))
    }
    if (i == 1){M <- phi[[1]]} else {M <- cbind(M, phi[[i]])}
  }
  K <- dim(phi[[1]])[[1]]

  M <- rbind(M, cbind(diag(K*(PP-1)), array(0,c(K*(PP-1),K))))
  return(max(Mod(eigen(M, only.values = TRUE)$values)))
}

specRadius <- function(M){
  return(max(Mod(eigen(M, only.values = TRUE)$values)))
}


initializer <- function(xx, k1=1, k2=1){
  PROJ = MAR1.PROJ(xx)
  if (specRadius(PROJ$A1)*specRadius(PROJ$A2) < 1){
    return(list(A1=PROJ$A1,A2=PROJ$A2))
  }
  MAR = MAR1.LS(xx)
  if (specRadius(MAR$A1)*specRadius(MAR$A2) < 1){
    return(list(A1=MAR$A1,A2=MAR$A2))
  }
  RRMAR = MAR1.RR(xx, k1, k2)
  if (specRadius(MAR1.RR$A1)*specRadius(MAR1.RR$A2) < 1){
    return(list(A1=MAR1.RR$A1,A2=MAR1.RR$A2))
  }
  stop('causality condition of initializer fails.')
}



# initializer <- function(xx, k1=1, k2=1){
#   dim = dim(xx)[-1]
#   p = dim(xx)[2]
#   q = dim(xx)[3]
# 
#   PROJ = MAR1.PROJ(xx)
#   if (specRadius(PROJ$A1)*specRadius(PROJ$A2) < 1){
#     A1 = PROJ$A1; A2 = PROJ$A2
#     eps1 = matrix(rnorm(p^2, sd=sqrt(sum(A1^2))/(p^2)), ncol=p)
#     eps2 = matrix(rnorm(q^2, sd=sqrt(sum(A2^2))/(q^2)), ncol=q)
#     return(list(A1=PROJ$A1+eps1,A2=PROJ$A2+eps2))
#   }
#   MAR = MAR1.LS(xx)
#   if (specRadius(MAR$A1)*specRadius(MAR$A2) < 1){
#     A1 = MAR$A1; A2 = MAR$A2
#     eps1 = matrix(rnorm(p^2, sd=sqrt(sum(A1^2))/(p^2)), ncol=p)
#     eps2 = matrix(rnorm(q^2, sd=sqrt(sum(A2^2))/(q^2)), ncol=q)
#     return(list(A1=MAR$A1+eps1,A2=MAR$A2+eps2))
#   }
#   RRMAR = MAR1.RR(xx, k1, k2)
#   if (specRadius(MAR1.RR$A1)*specRadius(MAR1.RR$A2) < 1){
#     A1 = RRMAR$A1; A2 = RRMAR$A2
#     eps1 = matrix(rnorm(p^2, sd=sqrt(sum(A1^2))/(p^2)), ncol=p)
#     eps2 = matrix(rnorm(q^2, sd=sqrt(sum(A2^2))/(q^2)), ncol=q)
#     return(list(A1=MAR1.RR$A1+eps1,A2=MAR1.RR$A2+eps2))
#   }
#   stop('causality condition of initializer fails.')
# }

initializer.sig <- function(xx){
  dim = dim(xx)[-1]
  t = dim(xx)[1]
  res = tenAR.VAR(xx, P=1)$res
  SIGMA = res %*% t(res) / (t-1)
  sig = projection(SIGMA, 1, dim[1],dim[2],dim[1],dim[2])[[1]]
  if (sig[[1]][1,1] < 0){sig[[1]] = - sig[[1]]}
  if (sig[[2]][1,1] < 0){sig[[2]] = - sig[[2]]}
  return(list(Sigl.init=sig[[1]], Sigr.init=sig[[2]]))
}


likelihood.lse <- function(fres, s, d, t){
  l1 <- fres/2/s^2
  l2 <- -(t - 1)*d*log(2*pi*s^2)/2
  return(l2 - l1)
}


IC <- function(xx,res,r,t,dim){
  N <- prod(dim)
  ic <- log(sum((res)^2)/(N*t))/2 + sum(r)*log(t)/t
  return(ic)
}

