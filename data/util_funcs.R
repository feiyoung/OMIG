
rm(list=ls())
library(MASS)
# Define functions --------------------------------------------------------
AE <- function(hX, X,Xmis, group){
  hmu <- colMeans(Xmis, na.rm=T)
  hmu[is.nan(hmu)] <- runif(1)
  Mu <- matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T)
  d <- length(unique(group))
  XList <- list()
  mea <- numeric(d)
  for(i in 1:d){
    mind1 <- is.na(Xmis[,group==i] )
    mea[i] <- mean(abs(hX[,group==i][mind1] - X[, group==i][mind1])) 
  }
  return(mea)
}


NAE <- function(hX, X,Xmis, group){
  hmu <- colMeans(Xmis, na.rm=T)
  hmu[is.nan(hmu)] <- runif(1)
  Mu <- matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T)
  d <- length(unique(group))
  XList <- list()
  mea <- numeric(d)
  for(i in 1:d){
    mind1 <- is.na(Xmis[,group==i] )
    mea[i] <- sum(abs(hX[,group==i][mind1] - X[, group==i][mind1])) /
      sum(abs(Mu[,group==i][mind1] - X[, group==i][mind1]))
  }
  return(mea)
}

imputeFun <- function(hD, Xmis, type, group){
  ind_set <- unique(group)
  ng <- length(ind_set)
  for(j in 1:ng){
    g1 <- which(group==j)
    if(type[j] == 'binomial'){
      N <- max(Xmis[, g1], na.rm = T)
    }
  }
  hX <- Xmis
  for(k in 1:ng){
    if(type[k]=='gaussian'){
      hX[,group==k] <- hD[,group==k]
    }else if(type[k] == 'binomial'){
      hX[,group==k] <- round(N/(1+exp(-hD[,group==k])) )
    }else if(type[k] == 'poisson'){
      hX[,group==k] <- round(exp(hD[,group==k]))
    }
    
  }
  hX[!is.na(Xmis)] <- Xmis[!is.na(Xmis)]
  return(hX)
}
getHessian <- function(Xmis1,Hm, Bm, type, group, rb_vec=NULL){
  n <- nrow(Xmis1)
  p <- ncol(Xmis1)
  O <- (!is.na(Xmis1))
  Xmis1[!O] <- 0
  ng <- length(type)
  if(is.null(rb_vec)){
    OdQ <- O / matrix(colMeans(O), n,p, byrow = T)
  }else{
    OdQ <- O / matrix(rb_vec, n,p, byrow = T)
  }
  
  gcell <- list()
  for(j in 1:ng){
    g1 <- which(group==j)
    gcell[[j]] <- g1
  }
  
  muMat<- matrix(NA, n, p)
  for(j in 1:ng){
    muMat[,group==j] <- switch(type[j],
                               gaussian = Hm %*% t(Bm[gcell[[j]],]),
                               poisson = exp(Hm %*% t(Bm[gcell[[j]],])),
                               binomial = 1/(1+exp(-Hm %*% t(Bm[gcell[[j]],])))
    )
  }
  
  hesList <- list()
  for(j in 1:p){
    hesList[[j]] <- switch (type[group[j]],
                            gaussian = t(Hm) %*% diag(OdQ[,j]) %*% Hm,
                            poisson = t(Hm) %*% diag(muMat[,j]*OdQ[, j]) %*% Hm,
                            binomial = t(Hm) %*%
                              diag(muMat[,j]*(1-muMat[,j])* OdQ[,j]) %*% Hm
    )
  }
  return(hesList)
}

getScore <- function(Xmisj, Hm, Bmj, typej, rb_vecj=NULL){
  n <- length(Xmisj)
  o <- (!is.na(Xmisj))
  Xmisj[!o] <- 0
  if(is.null(rb_vecj)){
    OdQ <- o / mean(o)
  }else{
    OdQ <- o / rb_vecj
  }
  
  ng <- length(type)
  muj <- switch(typej,
                gaussian = (Hm) %*% (Bmj),
                poisson = exp(Hm %*% (Bmj)),
                binomial = 1/(1+exp(-Hm %*% Bmj))
  )
  Uj <- t(Hm) %*% ((Xmisj - muj) * OdQ)
  return(Uj)
}
batchGroup <- function(n, nbatch1=0, nbatch=30, nleast=30){
  idxList <- list()
  n1 <- nbatch
  if(nbatch1==0){
    J <- floor(n/n1)
    for(j in 1:J){
      idxList[[j]] <- ((j-1)*n1+1):(j*n1)
    }
    if(n -n1*J >= nleast){
      J <- J + 1
      idxList[[J]] <- ((j*n1)+1): n
    }else if(n -n1*J > 0){
      idxList[[J]] <- c(idxList[[J]], ((j*n1)+1): n)
    }
  }else if(nbatch1 >0){
    idxList[[1]] <- 1: nbatch1
    n2 <- (n-nbatch1)
    J <- floor(n2/n1)+1
    for(j in 2:J){
      idxList[[j]] <- nbatch1 + (((j-2)*n1+1):((j-1)*n1))
    }
    
    if(n2 -n1*J >= nleast){
      J <- J + 1
      idxList[[J]] <- ((j*n1)+1): n
    }else if(nbatch1+((j-1)*n1) <n ){
      idxList[[J]] <- c(idxList[[J]], (nbatch1+((j-1)*n1)+1): n)
    }
  }
  return(idxList)
}

updateHs <- function(Xmis, Bm, Hm, gcell, type, rb_vec=NULL, parallel,lambda){
  # Xmis <- Xmisb; Bm <- Bms ; Hm <- Hmb
  n <- nrow(Xmis); p <- ncol(Xmis)
  O <- (!is.na(Xmis))
  Xmis[!O] <- 0
  if(is.null(rb_vec)){
    OdQ <- O / matrix(colMeans(O), n,p, byrow = T)
  }else{
    OdQ <- O / matrix(rb_vec, n,p, byrow = T)
  }
 
  if(ncol(Bm)-1 ==1){
    B <- matrix(Bm[,-1], ncol=1)
  }else{
    B <- Bm[,-1]
  }
  
  q <- ncol(Hm) - 1
  ng <- length(type)
  mucell <- list() # mean matrix
  for(j in 1:ng){
    muMat <- switch(type[j],
                          gaussian = Hm %*% t(Bm[gcell[[j]],]),
                          poisson = exp(Hm %*% t(Bm[gcell[[j]],])),
                          binomial = 1/(1+exp(-Hm %*% t(Bm[gcell[[j]],])))
                          
    )
    
    mucell[[j]] <- muMat
  }
  # Score matrix
  df2 <- matrix(0,n, q)
  for(j in 1:ng){
    df2 <- df2 + ((Xmis[,gcell[[j]]] - mucell[[j]]) *OdQ[,gcell[[j]]] ) %*% B[gcell[[j]],]
  }
  
  df2 <- df2 - lambda* Hm[,-1] # derivative of  regularized term
  
  
  
  # Hessian matrix or information matrix
  # d2f <- list()
  if(parallel){
    H2update <- single_parallel(parafunH ,iterable = 1:n, varlist= NULL, df2=df2, B=B, OdQ=OdQ, gcell=gcell,
                                mucell=mucell, type=type, lambda=lambda)
    H2update <- t(H2update)
  }else{
    H2update <- matrix(0, n, q)
    for(i in 1:n){
      Bng <- matrix(0, q, q)
      for(j in 1:ng){
        index <- gcell[[j]]
        dBng <- switch (type[j],
                        gaussian = t(B[index,]) %*% diag(OdQ[i, index]) %*% B[index,],
                        poisson = t(B[index,]) %*% diag(mucell[[j]][i,]*OdQ[i, index]) %*% B[index,],
                        binomial = t(B[index,]) %*%
                          diag(mucell[[j]][i,]*(1-mucell[[j]][i,])* OdQ[i, index]) %*% B[index,]
        )
        Bng <- Bng + dBng
      }
      H2update[i,] <- ginv(Bng + lambda*diag(q)) %*% df2[i,]
    }
  }
  
  Hm[,-1] <- Hm[,-1] + H2update
  return(Hm)
}


