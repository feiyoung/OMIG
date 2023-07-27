

OMIG <- function(Xmis_new, q,  group, type, hist.summary= NULL, lambda=1e-7,
                 verbose=TRUE){
  res <- streamGFMImpute(Xmis_new, q,  group, type, hist.summary= hist.summary, lambda=lambda,
                         verbose=verbose)
  return(res)
}



## imputation based on the estimates of common components hD and variable types.
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


## Evaluate the Hessian matrix in terms of loading-intercept vector for each data batch
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

## Evaluate the score vector in terms of loading-intercept vector for each data batch
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

## Estimate the factor matrix for each data batch
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
      H2update[i,] <- MASS::ginv(Bng + lambda*diag(q)) %*% df2[i,]
    }
  }

  Hm[,-1] <- Hm[,-1] + H2update
  return(Hm)
}



streamGFMImpute <- function(Xmis_new, q,  group, type, hist.summary= NULL, lambda=1e-7,
                             verbose=TRUE){
  ## compares the stastic GFMImpute for all data and the dynamic GFMImpute
  # lambda=0.1; verbose=T; X=NULL; Bm=NULL
  #Xmis_new <- dat$Xmis; group <- dat$group; hist.summary <- res$hist.summary


  X_imp <- Xmis_new
  n <- nrow(X_imp)
  p <- ncol(X_imp)


  Xmisb <- Xmis_new
  nrb_vec <- nrow(Xmisb) * colMeans(!is.na(Xmisb))
  ns <- nrow(Xmisb)
  rb_vec <- nrb_vec / ns
  if(is.null(hist.summary)){ # If this is the first data batch
    res <- try({
      misList1 <- gfmImpute(Xmisb, group, type, q=q, epsLogLike=1e-4,
                            maxIter=20,verbose=verbose, parallel=FALSE, lambda=lambda)
      Bms <- misList1$hBm
      Hmb <- misList1$hHm
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      hXb <- imputeFun(Hmb%*% t(Bms), Xmisb, type, group)
      X_imp[is.na(X_imp)] <- hXb[is.na(X_imp)]
      ## Evaluate the hessian matrix of current data batch
      hesLists <- getHessian(Xmisb,Hmb, Bms, type, group, rb_vec)

      hist.summary$nrb_vec <- nrb_vec
      hist.summary$ns <- ns
      hist.summary$Bms <- Bms
      hist.summary$hesLists <- hesLists
      }, silent = T)
    if(class(res) == 'try-error'){
      hmu <- colMeans(Xmisb, na.rm=T)
      X0 <- Xmisb - matrix(hmu, nrow(Xmisb), ncol(Xmisb), byrow=TRUE)
      X0[is.na(Xmisb)] <- 0
      Fac <- Factorm(X0, q, centered = TRUE)
      Hmb <- cbind(1, Fac$hH); Bms <- cbind(hmu, Fac$hB)
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      hXb <- Hmb%*% t(Bms)
      X_imp[is.na(X_imp)] <- hXb[is.na(X_imp)]
      hesLists <- getHessian(Xmisb,Hmb, Bms, type, group, rb_vec)

      hist.summary$nrb_vec <- nrb_vec
      hist.summary$ns <- ns
      hist.summary$Bms <- Bms
      hist.summary$hesLists <- hesLists
    }

  }else{

    ### update missing probability vector.
    nrb_vec <- hist.summary$nrb_vec + nrow(Xmisb) * colMeans(!is.na(Xmisb))
    hist.summary$nrb_vec <- nrb_vec
    ns <- nrow(Xmisb) + hist.summary$ns
    rb_vec <- nrb_vec / ns
    hist.summary$ns <- ns
    # rb_vec <- NULL

    ind_set <- unique(group)
    ng <- length(ind_set)
    gcell <- list()
    for(j in 1:ng){
      g1 <- which(group==j)
      gcell[[j]] <- g1
      if(type[j] == 'binomial'){
        N <- max(Xmisb[, g1], na.rm = T)
        Xmisb[, g1] <- Xmisb[, g1] / N
      }
    }
    Bms <- hist.summary$Bms
    ## Get initial values for estimating factor matrix of current data batch.
    hmu <- colMeans(Xmisb, na.rm=T)
    X0 <- Xmisb - matrix(hmu, nrow(Xmisb), p, byrow=T)
    X0[is.na(Xmisb)] <- 0
    Fac <- Factorm(X0, q, centered = TRUE)
    Hmb <- cbind(1, Fac$hH);
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    message("Estimate facor matrix: ")
    for(iter in 2:maxIter){
      try({Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)
      }, silent=T)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      message('dc=', format(dc,digits=4, scientific = TRUE))
      if(dc<1e-9) break
    }

    ## Given Hm and historical summary data, update Upsilon
    hesList2 <- getHessian(Xmisb,Hmb, Bms, type, group, rb_vec)
    hesLists1 <- hist.summary$hesLists
    hesLists <- lapply(1:p, function(j) hesLists1[[j]]+hesList2[[j]])
    rm(hesList2)
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    Bm1 <- Bms
    message("\nOnline updating for loading-intercept matrix: ")
    for(iter in 2:maxIter){
      for(j in 1:p){
        Ubj <- getScore(Xmisb[,j], Hmb, Bms[j,], type[group[j]], rb_vec[j])
        Uj <- hesLists1[[j]] %*% (Bm1[j,] - Bms[j,]) + Ubj
        Bms[j, ] <- Bms[j, ] + MASS::ginv(hesLists[[j]]) %*% Uj
      }
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      message('dc=', format(dc,digits=4, scientific = TRUE))
      if(dc<1e-9) break
    }

    hist.summary$hesLists <- hesLists
    rm(hesLists)
    # re-update H
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    message("\nRefit facor matrix: ")
    for(iter in 2:maxIter){
      try({Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)}, silent=T)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      message('dc=', format(dc,digits=4, scientific = TRUE))
      if(dc<1e-9) break
    }
    hist.summary$Bms <- Bms
    hXb <- imputeFun(Hmb%*% t(Bms), Xmisb, type, group)
    X_imp[is.na(X_imp)] <- hXb[is.na(X_imp)]

  }



  return(list(hX=X_imp, Hm=Hmb, hist.summary=hist.summary))
}

