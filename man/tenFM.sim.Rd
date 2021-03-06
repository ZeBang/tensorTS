% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tenFM.R
\name{tenFM.sim}
\alias{tenFM.sim}
\title{Generate Tensor Time series using given Factor Process and Factor Loading Matrices}
\usage{
tenFM.sim(Ft,dims=NULL,lambda=1,A=NULL,cov='iid',rho=0.2)
}
\arguments{
\item{Ft}{input of the factor process, of dimension \eqn{T \times r_1 \times r_2 \times \cdots \times r_k}. It can be TenAR(p) tensor time series generated by the function \link{tenAR.sim}.}

\item{dims}{dimensions of the output tensor at each time,  \eqn{d_1\times d_2\cdots\times d_K}.}

\item{lambda}{signal strength parameter of the tensor factor models, see Details section for more information.}

\item{A}{a list of the factor loading matrices \eqn{A_1, A_2, \cdots, A_K}. The default is random orthogonal matrices \eqn{A_k} of dimension \eqn{d_k \times r_k}.}

\item{cov}{covariance matrix of the error tensor: identity ("iid"), separable Kronecker structure ("separable"), random ("random").}

\item{rho}{a parameter only for "separable" covariance matrix of the error tensor. It is the off-diagonal element of the error matrices, with the diagonal being 1.}
}
\value{
A tensor-valued time series of dimension \eqn{T\times d_1\times d_2\cdots\times d_K}.
}
\description{
Simulate tensor time series \eqn{X_t} using a given factor process \eqn{F_t}. The factor process \eqn{F_t} can be generated by the function \code{\link{tenAR.sim}}.
}
\details{
Simulate from the model :
\deqn{X_t = \lambda F_t \times_{1} A_1 \times_{2} \cdots \times_{K} A_k + E_t,}
where \eqn{A_k} is the deterministic loading matrix of size \eqn{d_k \times r_k} and \eqn{r_k \ll d_k},
the core tensor \eqn{F_t} itself is a latent tensor factor process of dimension \eqn{r_1 \times \cdots \times r_K},
\eqn{\lambda} is an additional signal strength parameter,
and the idiosyncratic noise tensor \eqn{E_t} is uncorrelated (white) across time. In this function, by default \eqn{A_k} are orthogonal matrices.
}
\examples{
set.seed(333)
dims <- c(16,18,20) # dimensions of tensor time series
r <- c(3,3,3)  # dimensions of factor series
Ft <- tenAR.sim(t=100, dim=r, R=1, P=1, rho=0.9, cov='iid')
lambda <- sqrt(prod(dims))
# generate t*dims tensor time series with iid error covaraince structure
x <- tenFM.sim(Ft,dims=dims,lambda=lambda,A=NULL,cov='iid')
# generate t*dims tensor time series with separable error covaraince structure
x <- tenFM.sim(Ft,dims=dims,lambda=lambda,A=NULL,cov='separable',rho=0.2)
}
\seealso{
\code{\link{tenAR.sim}}
}
