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

Factorm <- function (X, q = NULL, centered=TRUE)
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
    hB <- B2 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
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
    hB <- hB1 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
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

ortheB <- function(Z){
  B <- qr(Z)
  eigvs <- sqrt(sort(eigen(t(Z)%*%Z)$values, decreasing = TRUE))
  B1 <- qr.Q(B) %*% Diag(eigvs)
  B0 <- B1 %*% Diag(sign(B1[1,]))
  return(B0)
}
ortheH <- function(H){
  H1 <- qr.Q(qr(H)) * sqrt(nrow(H))
  hH <- H1 %*% Diag(sign(H[1,]) * sign(H1[1,]))
  return(hH)
}


normvec <- function(x) sqrt(sum(x^2))


##  The normalized absolute error for each types of variables
NAE <- function(hX, X,Xmis, group){
  hmu <- colMeans(Xmis, na.rm=TRUE)
  Mu <- matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=TRUE)
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



# parallel computing for the single dynamic parameter.
single_parallel <- function(func,iterable,varlist=NULL,...){
  "
  :param func: function to be parallel computing
  :param iteralbe:a dynamic parameter(vector銆乴ist) of func.
  :param ...: all static paramters of func.
  :return list whose length is equivalent to that of iterable.
  "
  #1.load package
  library(parallel)
  #2.count the number of cores
  cores <- 10#detectCores(logical = FALSE)
  #3.open parallel computing nodes
  cl <- makeCluster(cores)
  #4.pass objects for each node.
  funcname <- deparse(substitute(func))
  varlist <- c(funcname,varlist)
  parallel::clusterExport(cl, varlist = varlist, envir = environment())
  # Put the reqiured functions in GFM package into all nodes.
  parallel::clusterCall(cl, function() library(glmnet))
  #5.start to parallel computing
  result <- parallel::parSapply(cl=cl,X=iterable,FUN=func,...)
  #6.close parallel computing
  stopCluster(cl)
  return(result)
}


preprocess_data <- function(Xmis, group, types){

  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match types!")
  }

  if(ng != length(types)){
    stop("The number of groups must match with length of types!")
  }
  # check zero-variance variables
  stdVec <- apply(Xmis, 2, sd, na.rm=TRUE)
  stdVec[is.na(stdVec)] <- 0 # remove the column that all entries are zeros.
  var_id <- which(stdVec < 1e-5)
  if(length(var_id)){
    message("There are ", length(var_id), " zero-variance variables in Xmis that will be removed!")
    meassage("They are ", var_id, " -th variables");

    Xmis <- Xmis[, -var_id]
    group <- group[-var_id]
  }
  type_set <- 1:ng
  for(s in 1:ng){
    if(sum(group==s)==0){
      type_set <- setdiff(type_set, s)
    }
  }
  types <- types[type_set]
  Xmiss <- Xmis


  n <- nrow(Xmis)
  type_scale <- list()
  for(s in 1:ng){
    switch (types[s],
      poisson = {
        cutoff <- 50
        id_types <- which(group==s)
        maxVec <- apply(Xmis[,id_types], 2, max, na.rm=TRUE)
        id <- id_types[maxVec > cutoff]
        if(!is.null(id)){
          Xmis[,id] <- Xmis[,id] / matrix(maxVec[maxVec> cutoff], n, length(id), byrow=TRUE)
        }
        type_scale[[types[s]]] <- cbind(id, maxVec[maxVec> cutoff])
      },
      gaussian = {
        cutoff <- 10
        id_types <- which(group==s)
        stdVec <- apply(Xmis[,id_types], 2, sd, na.rm=TRUE)
        id <- id_types[stdVec > cutoff]
        if(!is.null(id)){
          Xmis[,id] <- Xmis[,id] / matrix(stdVec[stdVec> cutoff], n, length(id), byrow=TRUE)
        }
        type_scale[[types[s]]] <- cbind(id, stdVec[stdVec> cutoff])
      }
    )

  }

  return(list(Xmis=Xmis, Xmiss = Xmiss, group=group, types=types, type_scale=type_scale))
}


## Find the group based on types
# find_group <- function(Xdf, types){
#
#   p <- ncol(Xdf)
#   ng <- length(types)
#   if(ng==1){
#     group <- rep(1, p)
#   }
#   if(ng >1){
#     group_type <- rep(NA, p)
#     id_res <- 1:p
#     for(k in 1:ng){
#     switch (types[k],
#           binomial = {
#             bin.flag1 <- apply(Xdf[,id_res], 2, function(x) length(unique(x))< 10)
#             bin.flag2 <- apply(Xdf[,id_res], 2, function(x) class(x)%in% c("character", "factor"))
#             idk <- id_res[bin.flag1 | bin.flag2]
#             group_type[idk] <- 'binomial'
#             id_res <- setdiff(id_res, idk)
#             if(length(id_res)==0){
#                break;
#             }
#
#           },
#           gaussian = {
#             cont.flag1 <- apply(Xdf[,id_res], 2, function(x) length(unique(x)) >= 10)
#             cont.flag2 <- apply(Xdf[,id_res], 2, function(x) class(x) =='numeric')
#             # cont.flag3 <- apply(Xdf[,id_res], 2, function(x) sum(x < 0, na.rm = TRUE) )
#             idk <-
#             group[cont.flag1 & cont.flag2] <- 'gaussian'
#           }
#     )
#
#     }
#   }
#
# }
