
add_identifiability <- function(H, B, mu){
  # Load the irlba library
  library(irlba)

  # Perform SVD decomposition with rank k = 10

  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  out <- try({
    svdHB <- irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
  }, silent = TRUE)
  if(class(out)=='try-error'){
    svdHB <- svd((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q, nu=q)
  }
  signB1 <- sign(svdHB$v[1,])
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)

  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)

  return(list(H=H, B=B, mu=mu))
}

transferMat2List <- function(Xmis, types, group){


  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types])
  if(length(unique(group)) != length(types)) stop("transferMat2List: The unique group must be the same length as that of types!")

  XmisList <- list()
  for(i in seq_along(typeID)){
    XmisList[[i]] <- Xmis[,group==i]
  }

  return(XmisList)

}

transferList2Mat <- function(XmisList, types){


  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types])

  group <- NULL
  Xmis <- NULL
  for(i in seq_along(typeID)){
    Xmis <- cbind(Xmis, XmisList[[i]])
    group <- c(group, rep(i, ncol(XmisList[[i]])))
  }

  return(list(Xmis=Xmis, group=group))

}


OrMIG_vb.fit <- function(XList, types, q, offset=FALSE, epsELBO=1e-5, maxIter=30, verbose=TRUE, seed = 1){

  ## epsELBO=1e-4; maxIter=30; verbose=TRUE


  n <- nrow(XList[[1]]);
  pvec <- sapply(XList, ncol); p <- sum(pvec)
  A <- matrix(0, n, p);
  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types])



  if(offset && any(typeID==2)){

    AList <- XList
    A <- NULL
    for(i in seq_along(typeID)){
      dim1 <- dim(AList[[i]])
      AList[[i]] <- matrix(0, dim1[1], dim1[2])
      if(typeID[i] == 2){
        tmpvec <- rowSums(AList[[i]], na.rm = TRUE)
        AList[[i]] <- matrix(log(tmpvec+(tmpvec==0)), dim1[1], dim1[2])
      }
      A <- cbind(A, AList[[i]]) # re-obtain A.
    }

    rm(AList)
  }

  library(Matrix)
  A <- as(A, "sparseMatrix")

  Mu_y_int <- NULL
  for(i in seq_along(typeID)){
    if(typeID[i]!=2){
      Mu_y_int <- cbind(Mu_y_int,XList[[i]])
    }else if(typeID[i] == 2){
      Mu_y_int <- cbind(Mu_y_int, log(1+ XList[[i]]))
    }

  }
  S_y_int = matrix(1, n, p);
  set.seed(seed)
  B_int = matrix(rnorm(p*q), p, q)*0.1;
  mu_int <- rep(0, p)
  Mu_h_int <-  matrix(rnorm(n*q), n, q)
  S_h_int  <- diag(rep(1, q))
  Sigma_h_int <- diag(rep(1, q))
  invLambda_int = rep(1, p);
  reslist <- VB_OrMIGcpp(XList, typeID, A,  Mu_y_int, S_y_int, invLambda_int,
                         B_int, mu_int, Mu_h_int, S_h_int, Sigma_h_int, epsELBO=epsELBO, maxIter=maxIter,
                         verbose=verbose)

  tmp_list <- transferList2Mat(XList, types)
  Xmis <- tmp_list$Xmis; group <- tmp_list$group
  rm(tmp_list)
  ng <- length(types)
  gcell <- list()
  for (j in 1:ng) {
    g1 <- which(group == j)
    gcell[[j]] <- g1
    if (types[j] == "binomial") {
      N <- max(Xmis[, g1], na.rm = TRUE)
      Xmis[, g1] <- Xmis[, g1]/N
    }
  }

  hHm <- cbind(1, reslist$H)
  hBm <- cbind(as.vector(reslist$mu), reslist$B)

  ## try
  try({
    hhB <- NULL
    for (j in 1:ng) {
      B1 <- localupdateB2(Xmis, gcell[[j]], hHm[, -1],
                          types[j], FALSE, 0)
      hhB <- cbind(hhB, B1)
    }
    hBm <- t(hhB)
    hHm <- updateH(Xmis, hBm, hHm, gcell, types, FALSE, 0)
    hhB <- NULL
    for (j in 1:ng) {
      B1 <- localupdateB2(Xmis, gcell[[j]], hHm[, -1],
                          types[j], FALSE, 0)
      hhB <- cbind(hhB, B1)
    }
    hBm <- t(hhB)
  }, silent=TRUE)

  hD <- hHm %*% t(hBm)


  res_idents <- add_identifiability(reslist$H, reslist$B, as.vector(reslist$mu))
  reslist$H <- res_idents$H
  reslist$B <- res_idents$B
  reslist$mu <- res_idents$mu




  hX <- Xmis
  for (k in 1:ng) {
    if (types[k] == "gaussian") {
      hX[, group == k] <- hD[, group == k]
    }
    else if (types[k] == "binomial") {
      N_vec <- apply(Xmis[,group==k], 2,  max, na.rm=TRUE)
      hX[, group == k] <- round(matrix(N_vec, nrow=n, ncol=length(N_vec), byrow = TRUE)/(1 + exp(-hD[, group == k])))
    }
    else if (types[k] == "poisson") {
      hX[, group == k] <- round(exp(hD[, group == k]))
    }
  }
  hX[!is.na(Xmis)] <- Xmis[!is.na(Xmis)]
  return(list(hX = hX, fitList=reslist))
}


