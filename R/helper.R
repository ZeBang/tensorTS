### Helper functions

myslice <- function(xx, K, start, end){
  if (K==2){
    return(xx[start:end,,,drop=FALSE])
  } else if (K==3){
    return(xx[start:end,,,,drop=FALSE])
  } else {
    stop("not support tensor mode K > 3")
  }
}

# mat projection
matAR.PROJ <- function(xx, dim, r, t){
  xx.mat <- matrix(xx,t,dim[1]*dim[2])
  kroneck <- t(xx.mat[2:t,]) %*% xx.mat[1:(t-1),] %*% solve(t(xx.mat[1:(t-1),]) %*% xx.mat[1:(t-1),])
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
  sd <- lapply(1:R, function(n) {lapply(1:K, function(m) {list()})})
  for (j in c(1:R)){
    for (i in c(1:K)){
      left <- sum(dim^2)*R + sum((dim^2)[1:(i-1)])+1
      right <- sum(dim^2)*R + sum((dim^2)[1:i])
      sd[[j]][[i]] <- array(diag(cov)[left:right], c(dim[i], dim[i]))
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
    stop("wrong dimention with your input Phi for rearrangement")
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
  ans <- matrix(NA, m1*n1, m2*n2)
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
    A[[i]] <- list(matrix(RA.svd$v[,i], m2, n2),matrix(RA.svd$u[,i] * RA.svd$d[i], m1, n1))
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
  cpd <- rTensor::cp(rTensor::as.tensor(tt), num_components = R)
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
        print("WRONG dimension")
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
        print("WRONG dimension")
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
        print("WRONG dimension")
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

ten.dis.A <- function(A, B){
  P = length(A)
  R <- length(A[[1]])
  K <- length(A[[1]][[1]])
  dis <- 0
  for (p in c(1:P)){
    for (r in c(1:R)){
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


ten.res <- function(xx,A,P,R,K,t){
  L1 = 0
  for (l in c(1:P)){
    if (R[l] == 0) next
    L1 <- L1 + Reduce("+",lapply(c(1:R[l]), function(n) {rTensor::ttl(as.tensor(myslice(xx, K, 1+P-l, t-l)), A[[l]][[n]], (c(1:K) + 1))@data}))
  }
  res <- myslice(xx, K, 1+P, t) - L1
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
      phi[[i]] = zeroes(prod(dim))
    } else {
      phi[[i]] <- Reduce("+", lapply(1:R[i], function(j) {rTensor::kronecker_list(rev(A[[i]][[j]]))}))
    }
    if (i == 1){M <- phi[[1]]} else {M <- cbind(M, phi[[i]])}
  }
  K <- dim(phi[[1]])[[1]]

  M <- rbind(M, cbind(diag(K*(PP-1)), array(0,c(K*(PP-1),K))))
  return(max(Mod(eigen(M, only.values = TRUE)$values)))
}


likelihood <- function(xx, A, Sigma){
  r <- length(A[[1]])
  dd <- dim(xx)
  t <- dd[1]
  dim <- dd[-1]
  k <- length(dd[-1])
  i = 1
  res <- ten.res(xx,A,P=1,R=r,K=k,t=t)
  Sigma.inv <- lapply(1:k, function (i) {solve(Sigma[[i]])})
  ll <- tl(res, Sigma.inv)
  l1 <- sum(diag(tensor(ll, res, c(1:4)[-(i+1)],c(1:4)[-(i+1)])))
  l2 <- 0
  for (i in c(1:k)){
    l2 = l2 - prod(dim[-i]) * (t-1) * (log(det(Sigma[[i]])))
  }
  return((l2 - l1)/2)
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
