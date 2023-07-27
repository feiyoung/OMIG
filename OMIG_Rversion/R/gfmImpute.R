

Diag <- function(x){
  if(length(x)==1){
    matrix(x,1,1)
  }else{
    diag(x)
  }

}

objMisfun <- function(Hm, Bm, Xmis,  gcell, type){

  mind <- !is.na(Xmis)
  probMis <- colMeans(mind)
  n <- nrow(Xmis); p <- ncol(Xmis)
  eps1 <- 1e-20
  BHm <- Hm %*% t(Bm)
  ng <- length(type)
  Q <- matrix(0, n, p)
  for(j in 1:ng){
    if(type[j]== 'gaussian'){
      Q[,gcell[[j]]] <- (Xmis[,gcell[[j]]] - BHm[,gcell[[j]]])^2
    }else if(type[j] == 'poisson'){
      me <- exp(BHm[,gcell[[j]]])
      Q[,gcell[[j]]] <- -log(ppois(Xmis[, gcell[[j]]], me)+eps1)
    }else if(type[j] == 'binomial'){
      me3 <- 1 / (1 + exp(-BHm[,gcell[[j]]]))
      Q[,gcell[[j]]] <- -Xmis[, gcell[[j]]] * log(me3+eps1) +
        (1-Xmis[, gcell[[j]]]) * log(1-me3 + eps1)
    }
  }
  obj <- sum((Q / matrix(probMis, n, p, byrow = TRUE))[mind]) / sum(mind) #
  return(obj)
}


family2func <- function(type){
  switch(type,
         gaussian= gaussian(link = "identity"),
         binomial = binomial(link = "logit"),
         poisson = poisson(link = "log"))
}
paraglmfit <- function(j, Xx, g1, XX, type1, mind, lambda){
  if(type1 == 'binomial' && sum(!is.na(unique(as.vector(XX[,g1])))) > 2){
    fun1 <- family2func(type1)
    glm.fit(x=Xx[mind[,g1[j]]==1,], y=XX[mind[,g1[j]]==1, g1[j]], family = fun1, intercept = FALSE)$coefficients
  }else{
    coef(glmnet(x=Xx[mind[,g1[j]]==1,], y=XX[mind[,g1[j]]==1, g1[j]],alpha=0,lambda=lambda, family = type1 , intercept = FALSE))[-1]
  }

}
localupdateB2 <- function(Xmis, g1, hH, type1, parallel=FALSE, lambda){
  # g1 <- gcell[[1]]; type1 <- type[1]; hH <- hHm[,-1]

  ## update loading-intercept matrix (mu, B)
  mind <- !is.na(Xmis)
  if(is.matrix(hH)) hH <- hH
  if(is.vector(hH)) hH <- matrix(hH, ncol=1)
  n <- nrow(Xmis);q <- ncol(hH)
  p1 <- length(g1); B1 <- matrix(0, q+1, p1)

  jg <- 1:p1

  if(parallel && length(jg)>= 1000){ # For a large number of variables, the parallel will be used.

    B1 <- single_parallel(paraglmfit,iterable = jg, Xx=cbind(1,hH), g1=g1, XX=Xmis, type1=type1, mind=mind, lambda=lambda)
  }else{
    if(type1 == 'binomial' && sum(!is.na(unique(as.vector(Xmis[,g1])))) > 2){
      fun1 <- family2func(type1)
      for(j in jg){
        B1[,j] <- glm.fit(x=cbind(1,hH[mind[,g1[j]]==1,]), y=Xmis[mind[,g1[j]]==1, g1[j]], family = fun1, intercept = FALSE)$coefficients
      }
      return(B1)
    }else{
      for(j in jg){
          fun1 <- family2func(type1)
          b1j <- glm.fit(x=cbind(1,hH[mind[,g1[j]]==1,]), y=Xmis[mind[,g1[j]]==1, g1[j]], family = fun1, intercept = FALSE)$coefficients

          # if(any(is.na(b1j))){
          #   x=cbind(1,hH[mind[,g1[j]]==1,]); y=Xmis[mind[,g1[j]]==1, g1[j]]
          #   b1j <- MASS::ginv(t(x)%*% x + 1e-7* diag(q+1))%*% t(x) %*% y
          # }
          B1[,j] <- b1j

      }
    }
  }

  return(B1)
}

updateH <- function(Xmis, Bm, Hm, gcell, type, parallel,lambda){
  # Bm <- hBm_init ; Hm <- hHm_init
  n <- nrow(Xmis); p <- ncol(Xmis)
  O <- (!is.na(Xmis))
  Xmis[!O] <- 0
  OdQ <- O / matrix(colMeans(O), n,p, byrow = TRUE)
  if(ncol(Bm)-1 ==1){
    B <- matrix(Bm[,-1], ncol=1)
  }else{
    B <- Bm[,-1]
  }

  q <- ncol(Hm) - 1
  ng <- length(type)
  mucell <- list() # mean matrix
  for(j in 1:ng){
    mucell[[j]] <- switch(type[j],
                          gaussian = Hm %*% t(Bm[gcell[[j]],]),
                          poisson = exp(Hm %*% t(Bm[gcell[[j]],])),
                          binomial = 1/(1+exp(-Hm %*% t(Bm[gcell[[j]],])))
    )
  }
  # Score matrix
  df2 <- matrix(0,n, q)
  for(j in 1:ng){
    df2 <- df2 + ((Xmis[,gcell[[j]]] - mucell[[j]]) *OdQ[,gcell[[j]]] ) %*% B[gcell[[j]],]
  }

  df2 <- df2 - lambda* Hm[,-1] # derivative of  regularized term



  # Hessian matrix or information matrix
  # d2f <- list()
  if(parallel && n >= 1000){
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
      try({
        H2update[i,] <- qr.solve(Bng + lambda*diag(q)) %*% df2[i,]
      }, silent=TRUE)

    }
  }
  H2update[is.na(H2update)] <- 0
  Hm[,-1] <- Hm[,-1] + H2update
  return(Hm)
}
parafunH <- function(i,df2, B, OdQ, gcell, mucell,type,  lambda){


  q <- ncol(B)
  ng <- length(gcell)
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
  return(qr.solve(Bng + lambda*diag(q)) %*% df2[i,])
}


gfmImpute <- function(Xmis, group, type, q,
                      epsLogLike=1e-5, maxIter=10,
                      verbose=0, parallel=FALSE, lambda=1e-20){
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
    if(type[j] == 'binomial'){
      N <- max(Xmis[, g1], na.rm = TRUE)
      Xmis[, g1] <- Xmis[, g1] / N
    }
  }
  n <- nrow(Xmis); p <- ncol(Xmis)
  if(length(group) != p){
    stop("The length of group must match with column of Xmis!")
  }
  #initialize
  hmu <- colMeans(Xmis, na.rm=TRUE)
  X0 <- Xmis - matrix(hmu, n, p, byrow=TRUE)
  X0[is.na(Xmis)] <- 0
  Fac <- Factorm(X0, q, centered = TRUE)
  hHm <- cbind(1, Fac$hH); hBm <- cbind(hmu, Fac$hB)
  eps1 <- 1e-4
  dBm <- Inf; dH <- Inf; dc =Inf; dOmega <- max(dBm, dH)
  tmpBm <- matrix(0, p, q+1); tmpHm <- hHm; tmpc <- 1e7
  k <- 1; history <- list()
  tic <- proc.time()
  while(k <= maxIter && dOmega > eps1 && dc >epsLogLike){
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(Xmis, gcell[[j]], hHm[,-1], type[j],parallel, lambda)
      hhB <- cbind(hhB, B1)
    }
    hBm <- t(hhB)
    # message("number of NA in Bm is: ",sum(is.na(hBm)), '\n')
    dB <- norm(hBm-tmpBm, "F") / norm(hBm, 'F')
    tmpBm <- hBm

    hHm <- updateH(Xmis, hBm, hHm, gcell, type, parallel,lambda)
    # message("number of NA in Hm is: ",sum(is.na(hHm)), '\n')
    dH <- norm(hHm- tmpHm, 'F')/norm(hHm, 'F')
    tmpH <- hHm

    dOmega <- max(dB, dH)
    c1 <- objMisfun(hHm, hBm, Xmis, gcell, type)
    dc <- abs(c1 - tmpc)/abs(tmpc)
    tmpc <- c1
    if(verbose){
      message('Iter=', k, ', dB=',format(dB,digits=4, scientific = TRUE), ', dH=', format(dH,digits=4),
          ',dc=', format(dc,digits=4, scientific = TRUE), ', negLogLike=',
          format(c1,digits=4, scientific = TRUE))
    }
    history$dB[k] <- dB; history$dH[k] <- dH; history$dc[k] <- dc;
    history$negLogLik[k] <- c1
    k <- k + 1
  }
  # Reupdate hBm once
  hhB <- NULL
  for(j in 1:ng){
    B1 <- localupdateB2(Xmis, gcell[[j]], hHm[,-1], type[j], parallel, lambda)
    hhB <- cbind(hhB, B1)
  }
  hBm <- t(hhB)
  stoc <- proc.time() - tic
  history$realIter <- k -1
  history$maxIter <- maxIter
  history$elapsedTime <- stoc

  criValue <- c("IC"=sum(log(c1+1e-12), q*(n+p)/(n*p)*log(n*p/(n+p))),
                #"IC"=sum(log(c1+1e-12), q/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)),
                "PC"=sum(c1, q * (n+p)/(n*p)*log(n*p/(n+p)))
  )

  hD <- hHm %*% t(hBm)
  ## add identifiability conditions
  hmu <- colMeans(hD)
  tD <- hD - matrix(hmu, n, p, byrow=TRUE)
  svdtD <- svd(tD, nu=q, nv=q)
  hHm <- cbind(1, svdtD$u*sqrt(n))
  hBm <- cbind(hmu, svdtD$v %*% Diag(svdtD$d[1:q])/sqrt(n))
  ng <- length(type)
  hX <- Xmis;
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

  return(list(hX=hX, hD=hD, hHm=hHm, hBm=hBm, history=history, cvVals=criValue))
}
parafun1 <- function(q, Xmis, group, type, ...){
  out <- try({
    gfmImpute(Xmis, group, type, q,parallel = FALSE,...)
  }, silent = FALSE)
  if(class(out) == 'try-error'){
    res <- rep(Inf,2)
  }else{
    res <- out$cvVals
  }

  return(res)
}



selectFacNumber <- function(Xmis,group, types, select_method=c('SVR', 'IC'), q_set= 2:10,
                            num_core=5,parallel=TRUE,
                             ...){
  nq <- length(q_set)

  if(select_method=="IC" || select_method=="PC"){
    if(parallel){
      # require(parallel)
      if (num_core > 1) {
        if (num_core > parallel::detectCores()) {
          warning("selectClustNumber:: the number of cores you're setting is larger than detected cores!")
          num_core = detectCores()
        }
      }

      require(doSNOW)
      cl <- makeSOCKcluster(num_core)
      registerDoSNOW(cl)

      ## set Prgogress Bar
      pb <- txtProgressBar(min=1, max=nq, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      k <- 1
      icMat <- foreach::foreach(q = 1:nq,.packages= "OMIG" ,.options.snow=opts,
                                .combine='rbind') %dopar% {
                                  # out <- gfmImpute(Xmis, group, types, q_set[q],parallel = FALSE,...)
                                  # out$cvVals
                                  out <- try({
                                    gfmImpute(Xmis, group, types, q_set[q],parallel = FALSE)
                                  }, silent = TRUE)
                                  if(class(out) == 'try-error'){
                                    rep(Inf,2)
                                  }else{

                                    out$cvVals
                                  }
                                }

      close(pb)
      stopCluster(cl)

    }else{
      icMat <- matrix(NA, nq, 2)
      pb <- txtProgressBar()
      for(k in 1:nq){
        out <- try({
          gfmImpute(Xmis, group, types, q_set[k])
        }, silent = TRUE)
        if(class(out) == 'try-error'){
          res <- rep(Inf,2)
        }else{

          res <- out$cvVals
        }

        setTxtProgressBar(pb, k)
        icMat[k, ] <- res
      }
      close(pb)
    }
    colnames(icMat) <- c("IC", 'PC')
    row.names(icMat) <- paste0("q=", q_set)
    if(select_method=="IC"){
      q <- q_set[which.min(icMat[,'IC'])]
      message('IC criterion estimates the factor number q  as ', q, '. \n')
    }
    if(select_method=="PC"){
      q <- q_set[which.min(icMat[,'PC'])]
      message('IC criterion estimates the factor number q  as ', q, '. \n')
    }



  }

  if(select_method == 'SVR'){
    q_max <- max(q_set)
    resmisList <- OrMIG(Xmis, group, types, q=q_max, ...)


    svalues <- svd(resmisList$fitList$B)$d
    svalues_use <- svalues[svalues>1e-2]
    if(length(svalues_use)<2) svalues_use <- svalues[1:2]
    q_max <- length(svalues_use)

    q <- which.max(svalues[-q_max] / svalues[-1])

    message('SVR estimates the factor number q  as ', q, '. \n')
  }


  return(q)

}

selectFacNumber2 <- function(Xmis, qmax=15){
  mnlamjFun <- function(eigvals, j){
    p <- length(eigvals)
    lamj <- eigvals[j]
    Sum <- 0
    for(l in (j+1):p){
      Sum <- Sum + 1/(eigvals[l] - lamj)
    }
    res <- Sum + 1/ ((3*lamj + eigvals[j+1])/4 - lamj)
    return(res/(p-j))
  }
  mtnlamjFun <- function(n, eigvals, j){
    p <- length(eigvals)
    rhojn <-  (p-j)/(n-1)
    res <- -(1-rhojn)/ eigvals[j] + rhojn * mnlamjFun(eigvals, j)
    return(res)
  }
  ##Reference: Fan, J., Guo, J., & Zheng, S. (2020). Estimating number of factors by adjusted eigenvalues thresholding. Journal of the American Statistical Association, 1-10.
  n <- nrow(Xmis)
  p <- ncol(Xmis)
  corMat <- cor(Xmis, use="pairwise.complete.obs")
  evalues <- eigen(corMat)$values
  hq1 <- sum(evalues>1+sqrt(p/(n-1)))
  if(hq1 < 15){
    hq <- hq1
  }else{ # ajdust the eigvalues
    adj.eigvals <- sapply(1:(p-1), function(j) -1/mtnlamjFun(n, evalues, j))
    hq <- sum(adj.eigvals >1) # overselect
  }
  propvar <- sum(evalues[1:hq]) / sum(evalues)
  res <- list()
  res$q <- hq
  res$propvar <- sum(evalues[1:hq]) / sum(evalues)

  return(res)
}
#selectFacNumber2(dat$Xmis)
