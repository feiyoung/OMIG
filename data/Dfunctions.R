

rm(list = ls())
library(MASS)
######################################################################
#### Missing value imputatation via generalized factor models ######
#######################################################################
##

# Main function of our methods --------------------------------------------

##
## main function: gfmImpute
##
##
## ---- Input Variables ----
## Xmis: a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-types data matrix, d is the types of variables, pj is the dimension of j-th types of variables.
## group: a vector with length equal to p, specify each column of X belonging to which group.
## type: a d-dimensional character vector, specify the types of variables in each group. For example, types=c('poisson', 'binomial'), and it is referred to the help file of glm.fit function for more details.
## q: a positive integer, specify the number of factors.
## epsLogLike: a positive real number, specify the tolerance of varing quantity of log-likelihood function in the algorithm. 
## maxIter: 	a positive integer, specify the times of iteration.
## verbose: a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process
## parallel: whether use parallel computing
## lambda: a ridge penalized parameter to increase the computational stability.

## 
##
## main function: gfmImpute
##

gfmImpute <- function(Xmis, group, type, q,
                      epsLogLike=1e-5, maxIter=10,
                      verbose=0, parallel=F, lambda=1e-20){
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
      N <- max(Xmis[, g1], na.rm = T)
      Xmis[, g1] <- Xmis[, g1] / N
    }
  }
  n <- nrow(Xmis); p <- ncol(Xmis)
  if(length(group) != p){
    stop("The length of group must match with column of Xmis!")
  }
  #initialize
  hmu <- colMeans(Xmis, na.rm=T)
  X0 <- Xmis - matrix(hmu, n, p, byrow=T)
  X0[is.na(Xmis)] <- 0
  Fac <- Factorm(X0, q, centered = T)
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
    # cat("number of NA in Bm is: ",sum(is.na(hBm)), '\n')
    dB <- norm(hBm-tmpBm, "F") / norm(hBm, 'F')
    tmpBm <- hBm
    
    hHm <- updateH(Xmis, hBm, hHm, gcell, type, parallel,lambda)
    # cat("number of NA in Hm is: ",sum(is.na(hHm)), '\n')
    dH <- norm(hHm- tmpHm, 'F')/norm(hHm, 'F')
    tmpH <- hHm
    
    dOmega <- max(dB, dH)
    c1 <- objMisfun(hHm, hBm, Xmis, gcell, type)
    dc <- abs(c1 - tmpc)/abs(tmpc)
    tmpc <- c1
    if(verbose){
      cat('Iter=', k, ', dB=',round(dB,4), ', dH=', round(dH,4),
          ',dc=', round(dc,4), ', c1=', round(c1,4), '\n')
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
  tD <- hD - matrix(hmu, n, p, byrow=T)
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

Factorm <- function (X, q = NULL, centered=T)
{
  if(!centered){
    X <- scale(X, scale=F)
  }
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
    svdX <- eigen(X %*% t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2 <- n^(-1) * t(X) %*% hatF
    sB <- sign(B2[1, ])
    hB <- B2 * matrix(sB, nrow = p, ncol = q, byrow = T)
    hH <- sapply(1:q, function(k) hatF[, k] * sign(B2[1,
                                                      ])[k])
  }
  else {
    svdX <- eigen(t(X) %*% X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])
    hH1 <- n^(-1) * X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u * sqrt(n), hH1)
    if (q == 1) {
      hB1 <- hB1 %*% svdH$d[1:q] * sqrt(n)
    }
    else {
      hB1 <- hB1 %*% diag(svdH$d[1:q]) * sqrt(n)
    }
    sB <- sign(hB1[1, ])
    hB <- hB1 * matrix(sB, nrow = p, ncol = q, byrow = T)
    hH <- sapply(1:q, function(j) hH2[, j] * sB[j])
  }
  sigma2vec <- colMeans((X - hH %*% t(hB))^2)
  res <- list()
  res$hH <- hH
  res$hB <- hB
  res$q <- q
  res$sigma2vec <- sigma2vec
  res$propvar <- sum(evalues[1:q])/sum(evalues)
  res$egvalues <- evalues
  attr(res, "class") <- "fac"
  return(res)
}

signrevise <- function(A1, A2){
  nzid1 <- which(rowSums(A1^2)> 1e-5)[1]
  q <- ncol(A1)
  A <- sapply(1:q, function(k){
    if(sign(A1[nzid1,k]) != sign(A2[nzid1,k]))
      return(-A1[,k])
    return(A1[,k])
  })
  return(A)
}
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
  obj <- sum((Q / matrix(probMis, n, p, byrow = T))[mind]) / sum(mind) #
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
    glm.fit(x=Xx[mind[,g1[j]]==1,], y=XX[mind[,g1[j]]==1, g1[j]], family = fun1, intercept = F)$coefficients
  }else{
    coef(glmnet(x=Xx[mind[,g1[j]]==1,], y=XX[mind[,g1[j]]==1, g1[j]],alpha=0,lambda=lambda, family = type1 , intercept = F))[-1]
  }
  
}
localupdateB2 <- function(Xmis, g1, hH, type1, parallel=F, lambda){
  # g1 <- gcell[[1]]; type1 <- type[1]; hH <- hHm[,-1]
  require(glmnet)
  mind <- !is.na(Xmis)
  if(is.matrix(hH)) hH <- hH
  if(is.vector(hH)) hH <- matrix(hH, ncol=1)
  n <- nrow(Xmis);q <- ncol(hH)
  p1 <- length(g1); B1 <- matrix(0, q+1, p1)
  
  jg <- 1:p1
  
  if(parallel){
    
    B1 <- single_parallel(paraglmfit,iterable = jg, Xx=cbind(1,hH), g1=g1, XX=Xmis, type1=type1, mind=mind, lambda=lambda)
  }else{
    if(type1 == 'binomial' && sum(!is.na(unique(as.vector(Xmis[,g1])))) > 2){
      fun1 <- family2func(type1)
      for(j in jg){
        B1[,j] <- glm.fit(x=cbind(1,hH[mind[,g1[j]]==1,]), y=Xmis[mind[,g1[j]]==1, g1[j]], family = fun1, intercept = F)$coefficients
      }
      return(B1)
    }else{
      for(j in jg){
        if(lambda==0){
          fun1 <- family2func(type1)
          b1j <- glm.fit(x=cbind(1,hH[mind[,g1[j]]==1,]), y=Xmis[mind[,g1[j]]==1, g1[j]], family = fun1, intercept = F)$coefficients
          
          # if(any(is.na(b1j))){
          #   x=cbind(1,hH[mind[,g1[j]]==1,]); y=Xmis[mind[,g1[j]]==1, g1[j]]
          #   b1j <- MASS::ginv(t(x)%*% x + 1e-7* diag(q+1))%*% t(x) %*% y
          # }
          B1[,j] <- b1j
        }else{
          B1[,j] <- coef(glmnet(x=cbind(1,hH[mind[,g1[j]]==1,]), y=Xmis[mind[,g1[j]]==1, g1[j]],alpha=0,
                                lambda=lambda, family = type1 , intercept = F))[-1]
          
        }
        
        
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
  OdQ <- O / matrix(colMeans(O), n,p, byrow = T)
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
      try({
        H2update[i,] <- qr.solve(Bng + lambda*diag(q)) %*% df2[i,]
      }, silent=T)
      
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


parafun1 <- function(q, Xmis, group, type, ...){
  out <- try({
    gfmImpute(Xmis, group, type, q,parallel = F,...)
  }, silent = F)
  if(class(out) == 'try-error'){
    res <- rep(Inf,2)
  }else{
    res <- out$cvVals
  }
  
  return(res)
}

### select the number of factors
selectFacNumber <- function(Xmis,group, type, q_set= 2:10, num_core=5,parallel=T, par.type='parallel', ...){
  nq <- length(q_set)
  if(parallel){
    require(parallel)
    if (num_core > 1) {
      if (num_core > detectCores()) {
        warning("selectClustNumber:: the number of cores you're setting is larger than detected cores!")
        num_core = detectCores()
      }
    }
    if(par.type=='doSNOW'){
      require(doSNOW)
      cl <- makeSOCKcluster(num_core)
      registerDoSNOW(cl)
      
      ## set Prgogress Bar
      pb <- txtProgressBar(min=1, max=nq, style=3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
      k <- 1
      icMat <- foreach::foreach(q = 1:nq,.packages= "gfmImpute" ,.options.snow=opts,
                                .combine='rbind') %dopar% {
                                  # out <- gfmImpute(Xmis, group, type, q_set[q],parallel = F,...)
                                  # out$cvVals
                                  out <- try({
                                    gfmImpute(Xmis, group, type, q_set[q],parallel = F,...)
                                  }, silent = T)
                                  if(class(out) == 'try-error'){
                                    rep(Inf,2)
                                  }else{
                                    
                                    out$cvVals
                                  }
                                }
      
      close(pb)
      stopCluster(cl)
    }else if(par.type=='parallel'){
      
      
      cl <- makeCluster(num_core)
      # clusterExport(cl, list("gfmImpute", "updateH"))
      cat("Starting parallel computing...")
      clusterCall(cl, function(){
        library(gfmImpute);library(glmnet)
      } )
      # Run
      icMat <- parallel::parSapply(cl, X=q_set, parafun1, Xmis=Xmis, group=group, type=type, ...)
      stopCluster(cl)
      icMat <- t(icMat)
    }
    
  }else{
    icMat <- matrix(NA, nq, 2)
    pb <- txtProgressBar()
    for(k in 1:nq){
      out <- try({
        gfmImpute(Xmis, group, type, q_set[k], ...)
      }, silent = T)
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
  colnames(icMat) <- c("ic", 'pc')
  row.names(icMat) <- paste0("q=", q_set)
  return(icMat)
  
}







# Online imputation via genelized factor models (OMIG) ---------------------------

##
## main function: streamGFMImpute
##
##
## ---- Input Variables ----
## Xmis: a matrix with mssing values and dimension of n*p(p=(p1+p2+..+pd)),observational mixed-types data matrix, d is the types of variables, pj is the dimension of j-th types of variables.
## group: a vector with length equal to p, specify each column of X belonging to which group.
## type: a d-dimensional character vector, specify the types of variables in each group. For example, types=c('poisson', 'binomial'), and it is referred to the help file of glm.fit function for more details.
## q: a positive integer, specify the number of factors.
## idxList: a list, the index set of each data batch.
## verbose: a logical value with TRUE or FALSE, specify whether ouput the mediate information in iteration process
## lambda: a ridge penalized parameter to increase the computational stability.
## X: the complete observed matrix, used for computing the imputation error.
## Bm: the intercept-loading matrix, used for computing the estimation error.

## 
##
## main function: streamGFMImpute
##
streamGFMImpute <- function(Xmis, q,  group, type, idxList,  lambda=0,
                            verbose=T, X=NULL, Bm=NULL){
  # compare offline using current data and the streaming gfmImpute.
  # lambda=0.1; verbose=T; X=NULL; Bm=NULL
  
  J <- length(idxList)
  X_imp <- Xmis
  b <- 1
  n <- nrow(Xmis)
  p <- ncol(Xmis)
  if(is.null(X)){
    NAE_mat <- NULL
    AE_mat <- NULL
  }else{
    NAE_mat <- matrix(NA, J, length(type))
    AE_mat <- matrix(NA, J, length(type))
  }
  if(is.null(Bm)){
    corB_vec <- NULL
  }else{
    corB_vec <- rep(NA, J)
  }
  Xmisb <- Xmis[idxList[[b]], ]
  nrb_vec <- nrow(Xmisb) * colMeans(!is.na(Xmisb)) 
  ns <- nrow(Xmisb)
  rb_vec <- nrb_vec / ns
  
  # misList1 <- gfmImpute(Xmisb, group, type, q=q, epsLogLike=1e-4, 
  #                       maxIter=20,verbose=verbose, parallel=F, lambda=lambda)
  res <- try({
    misList1 <- gfmImpute(Xmisb, group, type, q=q, epsLogLike=1e-4, 
                          maxIter=20,verbose=verbose, parallel=F, lambda=lambda)
    Bms <- misList1$hBm
    Hm <- misList1$hHm
    Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
    hXb <- imputeFun(Hm%*% t(Bms), Xmisb, type, group)
    X_imp[idxList[[b]], ] <- hXb
    1}, silent = T
  )
  
  if(class(res) == 'try-error'){
    hmu <- colMeans(Xmisb, na.rm=T)
    X0 <- Xmisb - matrix(hmu, nrow(Xmisb), ncol(Xmisb), byrow=T)
    X0[is.na(Xmisb)] <- 0
    Fac <- Factorm(X0, q, centered = T)
    Hm <- cbind(1, Fac$hH); Bms <- cbind(hmu, Fac$hB)
    Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
    hXb <- Hm%*% t(Bms)
  }
  Hms <- Hm
  hesLists <- getHessian(Xmisb,Hms, Bms, type, group, rb_vec)
  hXb <- imputeFun(Hms%*% t(Bms), Xmisb, type, group)
  X_imp[idxList[[b]],] <- hXb
  if(!is.null(X)){
    NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
    AE_mat[b,] <- AE(hXb, X[idxList[[b]],], Xmisb, group)
  }
  if(!is.null(Bm)){
    corB_vec[b] <- cancor(Bms, Bm)$cor[q]
  }
  #lambda <- 0
  
  for(b in 2:J){
    # b <- 2
    idx <- idxList[[b]]
    # cat(idx, '\n')
    Xmisb <- Xmis[idx, ]
    
    ### update missing probability vector.
    nrb_vec <- nrb_vec + nrow(Xmisb) * colMeans(!is.na(Xmisb)) 
    ns <- nrow(Xmisb) + ns
    rb_vec <- nrb_vec / ns
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
    hmu <- colMeans(Xmisb, na.rm=T)
    X0 <- Xmisb - matrix(hmu, nrow(Xmisb), p, byrow=T)
    X0[is.na(Xmisb)] <- 0
    Fac <- Factorm(X0, q, centered = T)
    Hmb <- cbind(1, Fac$hH); 
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    for(iter in 2:maxIter){
      Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    Hms <- rbind(Hms, Hmb)
    
    ## Given Hm, update Upsilon
    hesList2 <- getHessian(Xmisb,Hmb, Bms, type, group, rb_vec)
    hesLists1 <- hesLists
    hesLists <- lapply(1:p, function(j) hesLists1[[j]]+hesList2[[j]])
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    Bm1 <- Bms
    for(iter in 2:maxIter){
      for(j in 1:p){
        Ubj <- getScore(Xmisb[,j], Hmb, Bms[j,], type[group[j]], rb_vec[j])
        Uj <- hesLists1[[j]] %*% (Bm1[j,] - Bms[j,]) + Ubj
        Bms[j, ] <- Bms[j, ] + qr.solve(hesLists[[j]]) %*% Uj
      }
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    
    if(!is.null(Bm)){
      corB_vec[b] <- cancor(Bms, Bm)$cor[q]
    }
    hXb <- imputeFun(Hmb%*% t(Bms), Xmisb, type, group)
    if(!is.null(X)){
      NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
      AE_mat[b,] <- AE(hXb, X[idxList[[b]],], Xmisb, group)
    }
    X_imp[idxList[[b]],] <- hXb
  }
  #hXall <- imputeFun(Hms%*% t(Bms), Xmis, type, group)
  return(list(hX=X_imp, Hm=Hms, Bm=Bms, NAE_mat=NAE_mat,AE_mat=AE_mat, corB_vec=corB_vec))
}


streamGFMImpute2 <- function(Xmis, q,  group, type,idxList,  lambda=0,
                             verbose=T, X=NULL, Bm0=NULL, H0=NULL){
  ## compares the stastic GFMImpute for all data and the dynamic GFMImpute 
  # lambda=0.1; verbose=T; X=NULL; Bm=NULL
  # lambda =0; verbose =T;X=dat$X; Bm0=dat$Bm0; H0=dat$Hm0[,-1];group=group
  
  J <- length(idxList)
  errD <- NULL
  if(is.null(X)){
    NAE_mat <- NULL
  }else{
    NAE_mat <- matrix(NA, J, length(type))
  }
  if(is.null(Bm0)){
    corB_vec <- NULL
  }else{
    corB_vec <- rep(NA, J)
    errD <- NULL
  }
  if(is.null(H0)){
    corH_vec <- NULL
  }else{
    corH_vec <- rep(NA, J)
    errD <- NULL
  }
  
  
  
  X_imp <- Xmis
  b <- 1
  n <- nrow(Xmis)
  p <- ncol(Xmis)
  if(!is.null(H0) && !is.null(Bm0)){
    errD <- numeric(J)
    D0 <- cbind(1, H0) %*% t(Bm0)
  }
  
  Xmisb <- Xmis[idxList[[b]], ]
  nrb_vec <- nrow(Xmisb) * colMeans(!is.na(Xmisb)) 
  ns <- nrow(Xmisb)
  rb_vec <- nrb_vec / ns
  
  # misList1 <- gfmImpute(Xmisb, group, type, q=q, epsLogLike=1e-4, 
  #                       maxIter=20,verbose=verbose, parallel=F, lambda=lambda)
  res <- try({
    misList1 <- gfmImpute(Xmisb, group, type, q=q, epsLogLike=1e-4, 
                          maxIter=20,verbose=verbose, parallel=F, lambda=lambda)
    Bms <- misList1$hBm
    Hm <- misList1$hHm
    Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
    hXb <- imputeFun(Hm%*% t(Bms), Xmisb, type, group)
    X_imp[idxList[[b]], ] <- hXb
  }, silent = T
  )
  
  if(class(res) == 'try-error'){
    hmu <- colMeans(Xmisb, na.rm=T)
    X0 <- Xmisb - matrix(hmu, nrow(Xmisb), ncol(Xmisb), byrow=T)
    X0[is.na(Xmisb)] <- 0
    Fac <- Factorm(X0, q, centered = T)
    Hm <- cbind(1, Fac$hH); Bms <- cbind(hmu, Fac$hB)
    Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
    hXb <- Hm%*% t(Bms)
  }
  if(!is.null(X)){
    NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
  }
  if(!is.null(Bm0)){
    corB_vec[b] <- cancor(Bms, Bm0)$cor[q+1]
  }
  if(!is.null(H0)){
    corH_vec[b] <- cancor(Hm[,-1], H0[idxList[[b]],])$cor[q]
  }
  
  
  Hms <- Hm
  hesLists <- getHessian(Xmisb,Hms, Bms, type, group, rb_vec)
  hXb <- imputeFun(Hms%*% t(Bms), Xmisb, type, group)
  X_imp[idxList[[b]],] <- hXb
  
  if(!is.null(H0) && !is.null(Bm0)){
    idx <- 1: tail(idxList[[b]],1)
    errD[b] <- sum(abs((Hms %*% t(Bms) - D0[idx,])))/ prod(dim( D0[idx,]))
  }
  
  
  for(b in 2:J){
    # b <- 2
    idx <- idxList[[b]]
    # cat(idx, '\n')
    Xmisb <- Xmis[idx, ]
    
    ### update missing probability vector.
    nrb_vec <- nrb_vec + nrow(Xmisb) * colMeans(!is.na(Xmisb)) 
    ns <- nrow(Xmisb) + ns
    rb_vec <- nrb_vec / ns
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
    hmu <- colMeans(Xmisb, na.rm=T)
    X0 <- Xmisb - matrix(hmu, nrow(Xmisb), p, byrow=T)
    X0[is.na(Xmisb)] <- 0
    Fac <- Factorm(X0, q, centered = T)
    Hmb <- cbind(1, Fac$hH); 
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    for(iter in 2:maxIter){
      # iter <- 2
      try({Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)
      }, silent=T)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    Hms_int <- Hms
    Hms_int <- rbind(Hms_int, Hmb)
    
    ## Given Hm, update Upsilon
    hesList2 <- getHessian(Xmisb,Hmb, Bms, type, group, rb_vec)
    hesLists1 <- hesLists
    hesLists <- lapply(1:p, function(j) hesLists1[[j]]+hesList2[[j]])
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    Bm1 <- Bms
    for(iter in 2:maxIter){
      # iter <- 2
      for(j in 1:p){
        # j <- 1
        Ubj <- getScore(Xmisb[,j], Hmb, Bms[j,], type[group[j]], rb_vec[j])
        Uj <- hesLists1[[j]] %*% (Bm1[j,] - Bms[j,]) + Ubj
        Bms[j, ] <- Bms[j, ] + MASS::ginv(hesLists[[j]]) %*% Uj
      }
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    # re-update H
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    for(iter in 2:maxIter){
      # iter <- 2
      try({Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)}, silent=T)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    Hms <- rbind(Hms, Hmb)
    
    hXb <- imputeFun(Hmb%*% t(Bms), Xmisb, type, group)
    X_imp[idxList[[b]],] <- hXb
    
    if(!is.null(X)){
      NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
    }
    if(!is.null(Bm0)){
      corB_vec[b] <- cancor(Bms, Bm0)$cor[q+1]
    }
    if(!is.null(H0)){
      corH_vec[b] <- cancor(Hmb[,-1], H0[idxList[[b]],])$cor[q]
    }
    
    if(!is.null(H0) && !is.null(Bm0)){
      idx <- 1: tail(idxList[[b]],1)
      errD[b] <- sum(abs((Hms %*% t(Bms) - D0[idx,])))/ prod(dim( D0[idx,]))
    }
  }
  hXall <- imputeFun(Hms%*% t(Bms), Xmis, type, group)
  return(list(hX=X_imp, hXall=hXall, Hm=Hms, Bm=Bms, NAE_mat=NAE_mat,
              corB_vec=corB_vec, corH_vec=corH_vec, errD_vec =errD))
}


#### Offline version of OMIG
offlineGFMImpute <- function(Xmis, q, group, type, idxList, lambda=0, verbose=T, X=NULL,
                             Bm=NULL){
  # compare offline using current data and the streaming gfmImpute.
  X_imp <- Xmis
  b <- 1
  n <- nrow(Xmis)
  J <- length(idxList)
  Hms <- NULL
  if(is.null(X)){
    NAE_mat <- NULL
    AE_mat <- NULL
  }else{
    NAE_mat <- matrix(NA, J, length(type))
    AE_mat <- matrix(NA, J, length(type))
  }
  if(is.null(Bm)){
    corB_vec <- NULL
  }else{
    corB_vec <- rep(NA, J)
  }
  for(b in 1:J){
    Xmisb <- Xmis[idxList[[b]], ]
    res <- try({
      misList1 <- gfmImpute(Xmisb, group, type, q=q, epsLogLike=1e-4, 
                            maxIter=20,verbose=verbose, parallel=F, lambda=lambda)
      Bms <- misList1$hBm
      Hm <- misList1$hHm
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      hXb <- imputeFun(Hm%*% t(Bms), Xmisb, type, group)
      X_imp[idxList[[b]], ] <- hXb
      1}, silent = T
    )
    if(class(res) == 'try-error'){
      hmu <- colMeans(Xmisb, na.rm=T)
      X0 <- Xmisb - matrix(hmu, nrow(Xmisb), ncol(Xmisb), byrow=T)
      X0[is.na(Xmisb)] <- 0
      Fac <- Factorm(X0, q, centered = T)
      Hm <- cbind(1, Fac$hH); Bms <- cbind(hmu, Fac$hB)
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      hXb <- Hm%*% t(Bms)
      X_imp[idxList[[b]], ] <- hXb
    }
    if(!is.null(X)){
      NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
      AE_mat[b,] <- AE(hXb, X[idxList[[b]],], Xmisb, group)
    }
    if(!is.null(Bm)){
      
      corB_vec[b] <- cancor(Bms, Bm)$cor[q]
    }
    Hms <- rbind(Hms, Hm)
    
  }
  return(list(hX=X_imp, hHm=Hms, hBm=Bms, NAE_mat=NAE_mat,AE_mat=AE_mat, corB_vec=corB_vec))
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



# Compared methods --------------------------------------------------------

## methods based on factor models in Jin (2021)
# Jin, S., Miao, K., and Su, L. (2021). On factor models with random missing: Em estimation,
# inference, and cross validation. Journal of Econometrics, 222(1):745鈥?77.
LFM.JMS <- function(Xmis, q, maxIter=NULL, EM=T){
  
  hmu <- colMeans(Xmis, na.rm = T)
  mind <- is.na(Xmis)
  hr <- sum(mind) / prod(dim(Xmis))
  if(is.null(maxIter)) maxIter <- floor(log(0.001)/log(1-hr))
  X0 <- Xmis - matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T)
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
  
  return(X0 + matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T))
}

#Xiong, R. and Pelger, M. (2019). Large dimensional latent factor modeling with missing
# observations and applications to causal inference. arXiv preprint arXiv:1910.08273.
LFM.XP <- function(Xmis, q, prop.weighted=T){
  n <- nrow(Xmis)
  hmu <- colMeans(Xmis, na.rm = T)
  Oind <- (!is.na(Xmis)) *1
  if(prop.weighted) Oind <- Oind / matrix(colMeans(Oind),nrow=n, ncol=ncol(Xmis), byrow=T )
  X0 <- Xmis - matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T)
  covX0 <- cov(X0, use="pairwise.complete.obs")
  eig <- eigen(covX0)
  hLam <- eig$vectors[,1:q] * sqrt(n)
  hH <- matrix(0, n, q)
  for(i in 1:n){
    hH[i,] <- coef(lm(X0[i,] ~ hLam+0, weights = Oind[i,]))
  }
  hX <- hH %*% t(hLam) + matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T)
  hX[!is.na(Xmis)] <- Xmis[!is.na(Xmis)]
  return(hX)
}

# Data generation functions -----------------------------------------------

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
  
  return(list(Xmis=Xmis, X=X, Bm0=Bm0, Hm0=Hm0, group=group1))
}



# Define functions --------------------------------------------------------

## absolute error to measure the imputation error
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

## normalized absolute error to measure the imputation error
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



