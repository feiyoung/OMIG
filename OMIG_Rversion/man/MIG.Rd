\name{OrMIG}
\alias{OrMIG}
\alias{OrMIG}
\title{Missing Value Imputation via Generalized Factor Models}
\description{
  Missing value imputation for high-dimensional mixed-types data via generalized factor models under missing at random.
}
\usage{
  \method{OrMIG}{default}(Xmis, group, types, q, algorithm= c("VEM", "AM"),
                             epsLogLike=1e-5, maxIter=30, offset = FALSE,
                             verbose=FALSE, parallel=FALSE, lambda=1e-5, seed=1)
  \method{OrMIG}{matrix}(Xmis, group, types, q, algorithm= c("VEM", "AM"),
                             epsLogLike=1e-5, maxIter=30, offset = FALSE,
                             verbose=FALSE, parallel=FALSE, lambda=1e-5, seed=1)
  \method{OrMIG}{list}(XmisList, types, q, algorithm= c("VEM", "AM"),
                       epsLogLike=1e-5, maxIter=30, offset = FALSE,
                       verbose=FALSE, parallel=FALSE, lambda=1e-5, seed=1)
  \method{OrMIG}{data.frame}(DFmis,types, q, ...)
}
\arguments{
 \item{Xmis}{a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-types data matrix, d is the types of variables, pj is the dimension of j-th types of variables.}
 \item{XmisList}{a matrix list, each component matrix is dimension of n*pk.}
 \item{DFmis}{a data.frame with missing values that requires to be imputed}
 \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
 \item{types}{a d-dimensional character vector, specify the types of variables in each group. For example, types=c('poisson', 'binomial'), and it is referred to the help file of \link{glm.fit} function for more details.}
 \item{q}{a positive integer, specify the number of factors.}
 \item{algorithm}{a string, specify the algorithm to be used for fitting model. Now it supports two algorithms: variational EM (VEM) and alternate maximization (AM) algorithm, default as VEM. Empirically, we observed that VEM is more robust than Alternate_Maximization to the high noise data.}
 \item{epsLogLike}{a positive real number, specify the tolerance of varing quantity of log-likelihood/ELBO function in the algorithm. Optional parameter with default as 1e-5.}
 \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 30.}
 \item{offset}{a logical value, whether add an offset term (the total counts for each row in the count matrix) when there are Poisson variables.}
 \item{verbose}{a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process, (optional) default as FALSE.}
 \item{parallel}{a logical value with TRUE or FALSE, indicates wheter to use prallel computating and is an optional parameter with default as FALSE. Moreover, only if n>1000 and p>1000, parallel schema will be used even though \code{parallel=TRUE}.}
 \item{seed}{an integer, the random seed.}
}



\details{
  If S3 method for matrix is used, user requires to specify the types of each variable by "group"" and "types". For example, when \code{p=100} and \code{types=c("gaussian", "binomial")}, then \code{group=c(rep(1,50), rep(2, 50))}, which means the first 50 variables are of continuous types and the last 50 variables are of categorical types.

  If S3 method for list is used, user needs to  specify the missing marix list and the varialbe type of each component matrix of list.

  If S3 method for data.frame is used, the MIG will automatically find the variable types for each variable from argument "types".
}
\value{
 return a list object including following components:
 \item{hX}{an imputed matrix if one uses S3 method for matrix}
 \item{DFimp}{an imputed data.frame if one uses S3 method for data.frame}
 \item{types}{a string vector, the used variable types in MIG}
 \item{gfmList}{a list including the estimated factor matrix, loading matrix and other quantities in generalized factor model}
}
\references{
   Liu, W., Lin, H., Zheng, S., & Liu, J. (2021). Generalized factor model for ultra-high dimensional correlated variables with mixed types. Journal of the American Statistical Association, (just-accepted), 1-42.
}
\author{
Wei Liu
}
\note{
nothing
}
\seealso{
 None
}
\examples{
q <- 3
dat <- gendata(q = q, n=100, p=100, rho=6, mis_vec= c(0.1, 0.2, 0.3))
dat$Xmis[1:5,1:5]
Xmis <- dat$Xmis
group <- dat$group
types <- 'gaussian'
misList <- OrMIG(Xmis, group, types, q, verbose=TRUE)
NAE(misList$hX, dat$X, dat$Xmis, group)
## Based on data.frame
misList2 <- OrMIG(as.data.frame(Xmis), types, q, verbose=TRUE)
NAE(as.matrix(misList2$DFimp), dat$X, dat$Xmis, group)

## based on list
misList3 <- OrMIG(dat$XmisList, types, q, verbose=TRUE)
NAE(misList3$hX, dat$X, dat$Xmis, group)

misList4 <- OrMIG(Xmis, group,types, q, algorithm='AM', verbose=TRUE)
NAE(misList4$hX, dat$X, dat$Xmis, group)



}
