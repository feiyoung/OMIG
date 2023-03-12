Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}
cor.mat <- function (p, rho, type = "toeplitz"){
  if(p == 1) return(matrix(1,1,1))
  mat <- diag(p)
  if (type == "toeplitz") {
    for (i in 2:p) {
      for (j in 1:i) {
        mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
      }
    }
  }
  if (type == "identity") {
    mat[mat == 0] <- rho
  }
  return(mat)
}

#' Generate simulated missing mixed-type data.
#' @importFrom MASS mvrnorm
#' @description
#' The function \code{gendata} is used to generate simulated missing mixed-type data
#' @export
#' @param seed A non-negative integer, the random seed.
#' @param n A positive integer, the sample size.
#' @param p A positive integer, the variable dimension.
#' @param type A stering, the variable types to be generated.
#' @param q A positive integer, the number of factors.
#' @param mis_vec A real or real vector with values between 0 and 1, the missing rate of each variable.
#' @param rho  A positive real or postive vector with two components, the magnitude of loading matrix.
#' @param mechanism A string, the missing mechanism.
#'
#'
#' @examples
#' datList <- gendata(seed=1)
#'
gendata <- function(seed=1, n=300, p=50,  type='homonorm', q=6,mis_vec=0.3, rho=1, mechanism="MCAR"){
  #type = {'homonorm', 'heternorm', 'pois', 'norm_pois', 'pois_bino', 'npb'}
  if(length(rho)==1) rho <- c(rho, 1.5)
  factor_term <- rho[2]
  set.seed(1)
  Z <- matrix(rnorm(p*q), p, q)
  B <- qr(Z)
  eigvs <- sqrt(sort(eigen(t(Z)%*%Z)$values, decreasing = T))
  B0 <- rho[1]* qr.Q(B) %*% Diag(sqrt(eigvs))

  mu0 <- 0.4 * rnorm(p)
  Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  H0 <- MASS::mvrnorm(n, mu=rep(0,q), cor.mat(q, 0.5))

  if(type == 'homonorm'){
    X <- H0 %*% t(B0) + matrix(mu0, n,p, byrow=T) + 0.1*MASS::mvrnorm(n, rep(0,p), diag(p))
    group1 <- c(rep(1, ncol(X)))
  }else if(type == 'heternorm'){
    sigmas = 0.1 + 4* runif(p)
    X <- H0 %*% t(B0) + matrix(mu0, n,p, byrow=T) + MASS::mvrnorm(n, rep(0,p), diag(sigmas))
    group1 <- c(rep(1, ncol(X)))
  }else if(type == 'pois'){
    g1 <- 1:p
    B0[g1,] <- B0[g1,] /max(B0[g1,])* factor_term
    mu <- exp(H0 %*% t(B0) + matrix(mu0, n,p, byrow=T))
    X <- matrix(rpois(n*p,lambda=mu),n,p)
    group1 <- c(rep(1, ncol(X)))
  }else if(type == 'norm_pois'){
    g1 <- 1:floor(p/2)
    g2 <- (floor(p/2)+1):p
    Bm0[g2,-1] <- Bm0[g2,-1] /max(Bm0[g2,-1]) * factor_term/1.5
    mu1 <-  cbind(1,H0) %*% t(Bm0[g1,])
    mu2 <- exp(cbind(1, H0) %*% t(Bm0[g2,]) )
    X <- cbind(matrix(rnorm(prod(dim(mu1)), mu1,1), n, floor(p/2)),
               matrix(rpois(prod(dim(mu2)), mu2), n, ncol(mu2)))
    group1 <- c(rep(1, length(g1)), rep(2, length(g2)))
  }else if(type == 'pois_bino'){
    g1 <- 1:floor(p/2)
    g2 <- (floor(p/2)+1):p
    Bm0[g1,-1] <- Bm0[g1,-1] /max(Bm0[g1,-1])* factor_term
    mu1 <-  exp(cbind(1, H0) %*% t(Bm0[g1,]) )
    mu2 <- 1/(1+exp(-cbind(1, H0) %*% t(Bm0[g2,]) ))
    X <- cbind(matrix(rpois(prod(dim(mu1)), mu1), n, ncol(mu1)),
               matrix(rbinom(prod(dim(mu2)), 1, mu2), n, ncol(mu2)))
    group1 <- c(rep(1, length(g1)), rep(2, length(g2)))
  }else if(type == 'npb'){
    g1 <- 1:floor(p/3);
    g2 <- (floor(p/3)+1): floor(2*p/3)
    g3 <- (floor(2*p/3)+1):p
    Bm0[g2,-1] <- Bm0[g2,-1] /max(Bm0[g2,-1]) * factor_term
    mu1 <- cbind(1, H0) %*% t(Bm0[g1,])
    mu2 <- exp(cbind(1, H0) %*% t(Bm0[g2,]) )
    mu3 <- 1/(1+exp(-cbind(1, H0) %*% t(Bm0[g3,]) ))
    X <- cbind(matrix(rnorm(prod(dim(mu1)), mu1,1), n, ncol(mu1)),
               matrix(rpois(prod(dim(mu2)), mu2), n, ncol(mu2)),
               matrix(rbinom(prod(dim(mu3)), 1, mu3), n, ncol(mu3)))
    group1 <- c(rep(1, length(g1)), rep(2, length(g2)), rep(3, length(g3)))
  }else if (type == 'bino'){
    mu <- 1/(1+ exp(-H0 %*% t(B0) - matrix(mu0, n,p, byrow=T)))
    X <- matrix(rbinom(n*p,5, prob=mu),n,p)
    group1 <- c(rep(1, ncol(X)))
  }
  Hm0 <- cbind(1, H0)
  Xmis <- X
  hasBino <- FALSE
  if(is.element(type, c("binomial",'pois_bino', 'npb'))){
    hasBino <- TRUE
    idx_bino <- switch (type,
                        binomial = 1:p,
                        pois_bino = (floor(p/2)+1):p,
                        npb = (floor(2*p/3)+1):p
    )
  }


  # {
  #   mis_rate <- mis_vec
  #   mis.ind <- sample(n*p, floor(n*p*mis_rate))
  #   Xmis[mis.ind] <- NA
  # }else{
  if(mechanism=='MCAR'){
    n_mis <- length(mis_vec)
    for(kk in 1:n_mis){
      for(j in 1:p){
        if(j %% n_mis==kk-1){
          mis_rate <- mis_vec[kk]
          misj_ind <- sample(n, floor(n*mis_rate))
          Xmis[misj_ind, j] <- NA
        }
      }
    }
  }else if(mechanism=='MNAR'){
    n_mis <- length(mis_vec)
    for(kk in 1:n_mis){
      for(j in 1:p){
        if(j %% n_mis==kk-1){
          mis_rate <- mis_vec[kk]
          if(hasBino && j %in% idx_bino){
            misj_ind <- sample(n, floor(n*mis_rate))
          }else{
            misj_ind <- order(X[,j], decreasing = T)[1:floor(n*mis_rate)]
          }

          Xmis[misj_ind, j] <- NA
        }
      }
    }
  }else if(mechanism == "MAR"){
    n_mis <- length(mis_vec)
    p1 <- floor(p/2)
    XX <- cbind(X[,(p1+1):p],X[,1:p1])
    for(kk in 1:n_mis){
      for(j in 1:p){
        if(j %% n_mis==kk-1){
          mis_rate <- mis_vec[kk]
          if( j <= floor(p/2)){
            misj_ind <- order(XX[,j], decreasing = T)[1:floor(n*mis_rate)]
            Xmis[misj_ind, j] <- NA
          }else{
            misj_ind <- sample(n, floor(n*mis_rate))
            Xmis[misj_ind, j] <- NA
          }

        }
      }
    }
  }
  types <- c("gaussian", "poisson", "binomial")
  types <- types[unique(group1)]
  XmisList <- transferMat2List(Xmis, types, group1)


  return(list(XmisList=XmisList, Xmis=Xmis, X=X,  Bm0=Bm0, Hm0=Hm0, group=group1))
}
