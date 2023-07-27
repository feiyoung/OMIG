\name{OMIG}
\alias{OMIG}
\alias{OMIG}
\title{Online Missing Value Imputation via Generalized Factor Models}
\description{
  Online missing value imputation for high-dimensional mixed-type data via generalized factor models under missing at random.
}
\usage{
  OMIG(Xmis_new, q,  group, type, hist.summary= NULL, lambda=1e-7,
                 verbose=TRUE)
}
\arguments{
  \item{Xmis_new}{a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-type data matrix, d is the types of variables, pj is the dimension of j-th type of variables.}
  \item{q}{a positive integer, specify the number of factors.}
  \item{group}{a vector with length equal to p, specify each column of X belonging to which group.}
 \item{type}{a d-dimensional character vector, specify the type of variables in each group. For example, type=c('poisson', 'binomial'), and it is referred to the help file of \link{glm.fit} function for more details.}
 \item{hist.summary}{a list, including the historical summary statistics to improve the estimation efficiency and imputation accuracy of current data batch. If \code{Xmis_new} is the first data batch, then \code{hist.summary= NULL}; see details for more.}
}
\details{
   User requires to specify the type of each variable by "group"" and "type". For example, when \code{p=100} and \code{type=c("gaussian", "binomial")}, then \code{group=c(rep(1,50), rep(2, 50))}, which means the first 50 variables are of continuous type and the last 50 variables are of categorical type.

  OMIG returns a list indcluding a component named hist.summary which includes the summary statistics such as loading matrix, hessian matrix and that for updating the observed probability of each variable.
}
\value{
  OMIG returns a list including following components:
 \item{hX}{an imputed matrix if one uses S3 method for matrix}
 \item{Hm}{the estimated factor matrix for current data batch}
 \item{hist.summary}{a list including the summary statistics that helps to improve the estimation efficiency and imputation accuracy of the later data batch.}
}
\references{
  None
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

# Gaussian data -----------------------------------------------------------

type <- 'heternorm'
q <- 4
dat <- gendata(q = q, n=100, p=100, type= type, rho=2, mis_vec= c(0.1, 0.2, 0.3))
str(dat)
dat$Xmis[1:5,1:5]
Xmis_new <- dat$Xmis;
group <- dat$group;

hist.summary= NULL;
lambda=1;
verbose=TRUE
types <- 'gaussian'
res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=lambda,
            verbose=TRUE)

NAE(res$hX, dat$X, dat$Xmis, dat$group)
## compare with the offline method MIG
res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = T)
NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)

res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = T)
NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)

## Batch 2
dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=0, verbose=TRUE)

## Normalized absolute error
NAE(res2$hX, dat$X, dat$Xmis, dat$group)
## compare with the offline method MIG
res_mig <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "VEM", verbose = T)
NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "AM", verbose = T)
NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)

## Batch 3
dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary)

## Normalized absolute error
NAE(res3$hX, dat$X, dat$Xmis, dat$group)
## compare with the offline method MIG
res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM")
NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM")
NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)


}

