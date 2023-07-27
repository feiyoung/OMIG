# Others methods ----------------------------------------------------------

#' Missing data imputation using the methods (linear factor models)  in Jin (2021).
#' @importFrom MASS mvrnorm
#' @description
#' The function \code{LFM.JMS} imputes missing data using the methods (linear factor models)  in Jin (2021).
#' @export
#' @param Xmis A n-by-p matrix with missing values, the missing data.
#' @param q A positive integer, the number of factors.
#' @param maxIter A positive integer, the maximum iterations if EM=TRUE.
#' @param EM  A logical value, if it is TRUE, the LFM-EM is used, otherwise, the LFM-Init is used.
#'
#' @references Jin, S., Miao, K., and Su, L. (2021). On factor models with random missing: Em estimation, inference, and cross validation. Journal of Econometrics, 222(1):745–777
#' @examples
#' datList <- gendata(seed=1)
#' hX <- LFM.JMS(datList$Xmis, q=4)
LFM.JMS <- function(Xmis, q, maxIter=NULL, EM=TRUE){

  ## methods based on factor models in Jin (2021)
  # Jin, S., Miao, K., and Su, L. (2021). On factor models with random missing: Em estimation,
  # inference, and cross validation. Journal of Econometrics, 222(1):745–777.
  hmu <- colMeans(Xmis, na.rm = TRUE)
  mind <- is.na(Xmis)
  hr <- sum(mind) / prod(dim(Xmis))
  if(is.null(maxIter)) maxIter <- floor(log(0.001)/log(1-hr))
  X0 <- Xmis - matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE)
  X0[mind] <- 0
  svdX0 <- svd(X0 / hr, nu=q, nv=q)
  hH <- svdX0$u/nrow(X0)
  hB <- svdX0$v %*% diag(svdX0$d[1:q]) * nrow(X0)
  hX_center <- hH %*% t(hB)
  X0[mind] <- hX_center[mind]
  if(EM){
    for(i in 1:maxIter){
      cat('iter = ', i, '\n')
      svdX0 <- svd(X0, nu=q, nv=q)
      hX_center <- svdX0$u %*% diag(svdX0$d[1:q]) %*% t(svdX0$v)
      X0[mind] <- hX_center[mind]
    }
  }

  return(X0 + matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE))
}
# hX_final <- LFM.JMS(Xmis, q, 5, F)
# assessfun(hX_final[mind], dat$X[mind])
# cbind(hX_final[mind], dat$X[mind])[1:5, ]
# hX_final[1:4,1:4]


#' Missing data imputation using the methods (linear factor models)  in Xiong (2023).
#' @importFrom MASS mvrnorm
#' @description
#' The function \code{LFM.XP} imputes missing data using the methods (linear factor models)  in Xiong (2023).
#' @export
#' @param Xmis A n-by-p matrix with missing values, the missing data.
#' @param q A positive integer, the number of factors.
#' @param prop.weighted  A logical value, if it is TRUE, the LFM-PR is used, otherwise, the LFM-BR is used.
#'
#' @references Xiong, R., & Pelger, M. (2023). Large dimensional latent factor modeling with missing observations and applications to causal inference. Journal of Econometrics, 233(1), 271-301.
#' @examples
#' datList <- gendata(seed=1)
#' hX <- LFM.XP(datList$Xmis, q=4)

LFM.XP <- function(Xmis, q, prop.weighted=TRUE){
  #Xiong, R. and Pelger, M. (2021). Large dimensional latent factor modeling with missing
  # observations and applications to causal inference. arXiv preprint arXiv:1910.08273.
  n <- nrow(Xmis)
  hmu <- colMeans(Xmis, na.rm = TRUE)
  Oind <- (!is.na(Xmis)) *1
  if(prop.weighted) Oind <- Oind / matrix(colMeans(Oind),nrow=n, ncol=ncol(Xmis), byrow=TRUE )
  X0 <- Xmis - matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE)
  covX0 <- cov(X0, use="pairwise.complete.obs")
  eig <- eigen(covX0)
  hLam <- eig$vectors[,1:q] * sqrt(n)
  hH <- matrix(0, n, q)
  for(i in 1:n){
    hH[i,] <- coef(lm(X0[i,] ~ hLam+0, weights = Oind[i,]))
  }
  hX <- hH %*% t(hLam) + matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE)
  hX[!is.na(Xmis)] <- Xmis[!is.na(Xmis)]
  return(hX)
}

LFM_impute <- function(Xmis, q, group, type, prop.weighted=TRUE){
  n <- nrow(Xmis); p <- ncol(Xmis)

  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }

  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }

  gcell <- list()
  for(j in 1:ng){
    g1 <- which(group==j)
    gcell[[j]] <- g1
  }

  hmu <- colMeans(Xmis, na.rm = TRUE)
  Oind <- (!is.na(Xmis)) *1
  if(prop.weighted) Oind <- Oind / matrix(colMeans(Oind),nrow=n, ncol=ncol(Xmis), byrow=TRUE )
  X0 <- Xmis - matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE)
  covX0 <- cov(X0, use="pairwise.complete.obs")
  eig <- eigen(covX0)
  hLam <- eig$vectors[,1:q] * sqrt(n)
  hH <- matrix(0, n, q)
  for(i in 1:n){
    hH[i,] <- runif(q)
    try(hH[i,] <- coef(lm(X0[i,] ~ hLam+0, weights = Oind[i,])), silent = TRUE)
  }
  hHm=cbind(1, hH); hBm=cbind(hmu, hLam)
  c1 <- objMisfun(hHm, hBm, Xmis, gcell, type)
  criValue <- c("IC"=sum(log(c1+1e-12), q*(n+p)/(n*p)*log(n*p/(n+p))),
                #"IC"=sum(log(c1+1e-12), q/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)),
                "PC"=sum(c1, q * (n+p)/(n*p)*log(n*p/(n+p)))
  )
  hX <- hH %*% t(hLam) + matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE)
  hX[!is.na(Xmis)] <- Xmis[!is.na(Xmis)]

  return(list(hX=hX, hD = hX, hHm=hHm, hBm=hBm,cvVals=criValue))
}

