\name{selectFacNumber}
\alias{selectFacNumber}
\title{Select the number of fators}
\description{
  SVR or IC criterion is used to select the number of factors in the generalized factor model.
}
\usage{
  selectFacNumber(Xmis,group, types, select_method=c('SVR', 'IC'), q_set= 2:10,
                            num_core=5, parallel=TRUE, ...)
}
\arguments{
  \item{Xmis}{a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-type data matrix, d is the types of variables, pj is the dimension of j-th type of variables.}
 \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
 \item{types}{a d-dimensional character vector, specify the type of variables in each group. For example, type=c('poisson', 'binomial'), and it is referred to the help file of \link{glm.fit} function for more details.}
 \item{select_method}{a character, specify the method for selecting the number of factors. Here, two methods are provided: singular value ratio based methods (SVR) and the information criterion (IC), default as SVR due to more computational effeciency.}
 \item{q_set}{a positive integer, specify the number of factors.}
 \item{num_core}{a positive integer, the number of cores used when using the `IC` criterion.}
 \item{parallel}{a logical value, whether use the parallel computing when using the `IC` criterion.}
 \item{...}{Other arguments passed to \code{\link{OrMIG}}} function.
}
\details{
  Nothing
}
\value{
  return a matrix with two columns with column names "ic" and "pc" and row names "q=xx", where "ic" means the IC criteria value and "pc" is the PC criteria value. We recommend using IC criteria.
}
\references{
  Wei Liu, Lan Luo, Ling Zhou, Online missing value imputation for high-dimensional mixed-type data via generalized factor models, Computational Statistics & Data Analysis, 2023.
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
  hq <- selectFacNumber(dat$Xmis, group, type, q_set=q_set,select_method='SVR',
                      parallel = T, verbose=T)

}

