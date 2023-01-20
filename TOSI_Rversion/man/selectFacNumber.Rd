\name{selectFacNumber}
\alias{selectFacNumber}
\title{Select the number of fators}
\description{
  IC criteria is used to select the number of factors in the generalized factor model.
}
\usage{
  selectFacNumber(Xmis,group, type, q_set= 2:10, num_core=5,parallel=TRUE,...)
}
\arguments{
  \item{Xmis}{a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-type data matrix, d is the types of variables, pj is the dimension of j-th type of variables.}
 \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
 \item{type}{a d-dimensional character vector, specify the type of variables in each group. For example, type=c('poisson', 'binomial'), and it is referred to the help file of \link{glm.fit} function for more details.}
 \item{q}{a positive integer, specify the number of factors.}

}
\details{
  Nothing
}
\value{
  return a matrix with two columns with column names "ic" and "pc" and row names "q=xx", where "ic" means the IC criteria value and "pc" is the PC criteria value. We recommend using IC criteria.
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
  \link{MIG} and \link{OMIG}.
}
\examples{
# Mix of normal, poisson and binomial
  dat <- gendata(seed=3, n=200, p=200, type='npb',
               q=6, rho=7, mis_vec = 0.4)
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('gaussian','poisson', 'binomial')
  q_set <- 3:8

  # Select number of factors
  icMat <- selectFacNumber(dat$Xmis,group, type, q_set=q_set,
                           parallel = T, verbose=T)
  rownames(icMat)[which.min(icMat[,'ic'])]
  rownames(icMat)[which.min(icMat[,'pc'])]
}

