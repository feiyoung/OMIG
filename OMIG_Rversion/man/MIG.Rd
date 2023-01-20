\name{OrMIG}
\alias{OrMIG}
\alias{OrMIG}
\title{Missing Value Imputation via Generalized Factor Models}
\description{
  Missing value imputation for high-dimensional mixed-type data via generalized factor models under missing at random.
}
\usage{
  \method{OrMIG}{matrix}(Xmis, group, type, q,epsLogLike=1e-5, maxIter=10,
                             verbose=FALSE, parallel=FALSE, lambda=1e-5, ...)
  \method{OrMIG}{data.frame}(DFmis,type, q, ...)
}
\arguments{
 \item{Xmis}{a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-type data matrix, d is the types of variables, pj is the dimension of j-th type of variables.}
 \item{DFmis}{a data.frame with missing values that requires to be imputed}
 \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
 \item{type}{a d-dimensional character vector, specify the type of variables in each group. For example, type=c('poisson', 'binomial'), and it is referred to the help file of \link{glm.fit} function for more details.}
 \item{q}{a positive integer, specify the number of factors.}
 \item{epsLogLike}{a positive real number, specify the tolerance of varing quantity of log-likelihood function in the algorithm. Optional parameter with default as 1e-4.}
 \item{maxIter}{a positive integer, specify the times of iteration. Optional parameter with default as 10.}
 \item{verbose}{a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process, (optional) default as FALSE.}
 \item{parallel}{a logical value with TRUE or FALSE, indicates wheter to use prallel computating and is an optional parameter with default as FALSE. Moreover, only if n>1000 and p>1000, parallel schema will be used even though \code{parallel=TRUE}.}
 \item{...}{Other arguments passed from S3 method for 'matrix' and 'default'.}
}



\details{
  If S3 method for matrix is used, user requires to specify the type of each variable by "group"" and "type". For example, when \code{p=100} and \code{type=c("gaussian", "binomial")}, then \code{group=c(rep(1,50), rep(2, 50))}, which means the first 50 variables are of continuous type and the last 50 variables are of categorical type.

  If S3 method for data.frame is used, the MIG will automatically find the variable type for each variable from argument "type".
}
\value{
 return a list object including following components:
 \item{hX}{an imputed matrix if one uses S3 method for matrix}
 \item{DFimp}{an imputed data.frame if one uses S3 method for data.frame}
 \item{type}{a string vector, the used variable types in MIG}
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
type <- 'gaussian'
misList <- OrMIG(Xmis, group, type, q,maxIter=3, verbose=TRUE)
## Based on data.frame
misList2 <- OrMIG(as.data.frame(Xmis), type, q, maxIter=3, verbose=TRUE)
}
