###Functions of autoregressive Models

#'Rearrangement Operator
#'
#'Rearrangement Operator used for projection method.
#'@name rearrange
#'@rdname rearrange
#'@aliases rearrange
#'@export
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
  # A \approx B \otimes C
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
#'@export
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
#'
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

