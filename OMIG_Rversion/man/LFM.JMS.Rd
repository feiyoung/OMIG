% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/otherMethods.R
\name{LFM.JMS}
\alias{LFM.JMS}
\title{Missing data imputation using the methods (linear factor models)  in Jin (2021).}
\usage{
LFM.JMS(Xmis, q, maxIter = NULL, EM = TRUE)
}
\arguments{
\item{Xmis}{A n-by-p matrix with missing values, the missing data.}

\item{q}{A positive integer, the number of factors.}

\item{maxIter}{A positive integer, the maximum iterations if EM=TRUE.}

\item{EM}{A logical value, if it is TRUE, the LFM-EM is used, otherwise, the LFM-Init is used.}
}
\description{
The function \code{LFM.JMS} imputes missing data using the methods (linear factor models)  in Jin (2021).
}
\examples{
datList <- gendata(seed=1)
hX <- LFM.JMS(datList$Xmis, q=4)
}
