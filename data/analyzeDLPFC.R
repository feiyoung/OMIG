


# Run each compared method --------------------------------------


### Evaluate all missing cases:
setwd("./")
source("./util_funcs.R")
library(R.matlab)
dat <- readMat("./genBrainMatDat/misRateAllData.mat")
str(dat)
n_mis <- length(dat$idna10.cell)

X <- dat$X
q <- 4
library(OMIG)
library(missMDA)
group <- dat$group
N <- ncol(dat$idna10.cell[[1]][[1]])
naeArray <- array(NA, dim = c(6,2, N, n_mis))
aeArray <- naeArray
timeArray <- array(NA, dim = c(6, N, n_mis))
for(imis in 1: n_mis){
  # imis <- 1
  cat('imis = ', imis, '\n')
  idnaMat <- dat$idna10.cell[[imis]][[1]]
  for(i in 1:N){
    # i <- 1
    cat('i = ', i, '\n')
    Xmis <- X
    Xmis[idnaMat[,i]] <- NA
    # FAMD and iterative FAMD
    Xdat <- as.data.frame(Xmis)
    tic <- proc.time()
    hX_famd <- imputeFAMD(Xdat, ncp=q, maxiter = 1)$completeObs
    hX_famd_mat <- as.matrix(hX_famd)
    toc <-  proc.time() 
    time_famd <- toc[3] - tic[3]
    nae_famd <- NAE(hX_famd_mat, X, Xmis, group)
    ae_famd <- AE(hX_famd_mat, X, Xmis, group)
    
    tic <- proc.time()
    hX_iterfamd <- imputeFAMD(Xdat, ncp=q)$completeObs
    hX_iterfamd_mat <- as.matrix(hX_iterfamd)
    toc <-  proc.time() 
    time_iterfamd <- toc[3] - tic[3]
    nae_iterfamd <- NAE(hX_iterfamd_mat,X, Xmis, group)
    ae_iterfamd <- AE(hX_iterfamd_mat,X, Xmis, group)
    
    
    ## JMS method
    tic <- proc.time()
    hX_JMS_int <- LFM.JMS(Xmis, q, EM=F)
    nae_JMS_int <- NAE(hX_JMS_int, X, Xmis, group)
    ae_JMS_int <- AE(hX_JMS_int, X, Xmis, group)
    toc <-  proc.time() 
    time_JMS_int <- toc[3] - tic[3]
    tic <- proc.time()
    hX_JMS_final <- LFM.JMS(Xmis, q,maxIter = 5, EM=T)
    nae_JMS_final <- NAE(hX_JMS_final, X, Xmis, group)
    ae_JMS_final <- AE(hX_JMS_final, X, Xmis, group)
    toc <-  proc.time() 
    time_JMS_final <- toc[3] - tic[3]
    
    # PX method
    tic <- proc.time()
    hX_XP_int <- LFM.XP(Xmis, q, prop.weighted = F)
    toc <-  proc.time() 
    time_XP_int <- toc[3] - tic[3]
    nae_XP_int <- NAE(hX_XP_int, X, Xmis, group)
    ae_XP_int <- AE(hX_XP_int, X, Xmis, group)
    
    tic <- proc.time()
    hX_XP_final <- LFM.XP(Xmis, q,prop.weighted = T)
    toc <-  proc.time() 
    time_XP_final <- toc[3] - tic[3]
    nae_XP_final <- NAE(hX_XP_final, X, Xmis, group)
    ae_XP_final <- AE(hX_XP_final, X, Xmis, group)
    
    naeArray[,,i, imis] <- rbind(nae_famd, nae_iterfamd, nae_JMS_int, nae_JMS_final, nae_XP_int, nae_XP_final)
    
    aeArray[,,i, imis] <- rbind(ae_famd, ae_iterfamd, ae_JMS_int, ae_JMS_final, ae_XP_int, ae_XP_final)
    
    timeArray[,i, imis] <- c(time_famd, time_iterfamd, time_JMS_int, time_JMS_final, time_XP_int, time_XP_final)
  }
}

apply(naeArray, c(1,2), mean)

save(naeArray, aeArray, file='./resBrainDat/otherMethod_brain_NAE_AE.Rdata')



# Compare OMIG with OfMIG-------------------------------------------------------
#setwd('/home/ligz/LiuWei_files/OMIG')
source("./util_funcs.R")
library(OMIG)
setwd("./genBrainMatDat")
library(R.matlab)
dat <- readMat("misRateAllData.mat")
str(dat)
X <- dat$X
q <- 4
library(gfmImpute)
library(missMDA)
group <- dat$group
X[,group==2] <- X[,group==2] / 10 # scale the raw data for stability
N <- 10
naeArray <- array(NA, dim = c(6,2, N))
aeArray <- array(NA, dim = c(6,2, N))
for(imis in 3){
  # imis <- 3
  cat('imis = ', imis, '\n')
  idnaMat <- dat$idna10.cell[[imis]][[1]]
  for(i in 1:N){
    # i <- 1
    cat('i = ', i, '\n')
    Xmis <- X
    Xmis[idnaMat[,i]] <- NA
  }
}
mis_rate <- 0.3
jj <- 1
q <- 4
n <- nrow(X); p <- ncol(X); type1 = 'norm_pois'
N <- 10
d <- 2 # two types

MethodNames <- c("OfMIG", "OMIG")

##
idxList <- batchGroup(n, nbatch1=1000, nbatch = 400)
J <- length(idxList)
errArray <- array(NA, dim=c(J,d, 2,  N) )
colnames(errArray) <- c("Normal", "Poisson")
AE_errArray <- errArray
lambda <- 0
## get the data
imis <- 3
cat('imis = ', imis, '\n')
idnaMat <- dat$idna10.cell[[imis]][[1]]

for(i in 1:N){
  # i <- 1
  cat('i = ', i, '\n')
  Xmis <- X
  Xmis[idnaMat[,i]] <- NA
  length(idnaMat[,i]) / prod(dim(Xmis))
  group <- dat$group
  type <- c('gaussian', 'poisson')
  
  NAE_mat <- matrix(NA, J, d)
  AE_mat <-  matrix(NA, J, d)
  for(j in 1:J){
    # j <- 2
    idx <- idxList[[j]]
    out <- try({
      misList <- OrMIG(Xmis[idx,], group, type,q,lambda=lambda, verbose=T,
                           maxIter=20,epsLogLike=1e-4)
     
      hX <- misList$hX
      NAE_mat[j, ] <- NAE(hX, X[idxList[[j]],], Xmis[idxList[[j]],], group)
      AE_mat[j, ] <- AE(hX, X[idxList[[j]],], Xmis[idxList[[j]],], group)
      
      
    })
    errArray[,,1,i] <- NAE_mat
    AE_errArray[,, 1,i] <- AE_mat
    out <- try({
      sList <- dynamicGFMImpute_real(Xmis, q,  group, type, idxList, lambda=lambda,
                                     verbose=T, X=X)
    }, silent = T)
    
    errArray[,,2,i] <- sList$NAE_mat
    AE_errArray[,, 2,i] <- sList$AE_mat
  }
}



errArray[,,,2]
apply(errArray, c(1,2,3), mean, na.rm=T)
apply(AE_errArray, c(1,2,3), mean, na.rm=T)
save(errArray, AE_errArray, file='./RealData/resBrainDat/AE_dynamic_static_brain_misRate03.Rdata')


dynamicGFMImpute_real <- function(Xmis, q,  group, type, idxList, lambda=lambda,
                                  verbose=T, X=X){
  
  
  require(OMIG)
  
  J <- length(idxList)
  NAE_mat <- matrix(NA, J, length(type))
  AE_mat <- NAE_mat
  X_imp <- Xmis
  b <- 1
  n <- nrow(Xmis)
  
  Xmisb <- Xmis[idxList[[b]], ]
  nrb_vec <- nrow(Xmisb) * colMeans(!is.na(Xmisb)) 
  ns <- nrow(Xmisb)
  rb_vec <- nrb_vec / ns
  # verbose <- T
  
  res <- try({
    misList1 <- OrMIG(Xmisb, group, type, q=q, epsLogLike=1e-4, 
                          maxIter=20,verbose=verbose, parallel=F, lambda=lambda)
    Bms <- misList1$hBm
    Hm <- misList1$hHm
    Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
    hXb <- imputeFun(Hm%*% t(Bms), Xmisb, type, group)
    X_imp[idxList[[b]], ] <- hXb
    1}, silent = T
  )
  
  
  Hms <- Hm
  hesLists <- getHessian(Xmisb,Hms, Bms, type, group, rb_vec)
  hXb <- imputeFun(Hms%*% t(Bms), Xmisb, type, group)
  X_imp[idxList[[b]],] <- hXb
  
  NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
  AE_mat[b,] <- AE(hXb, X[idxList[[b]],], Xmisb, group)
  
  for(b in 2:J){
    # b <- 3
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
    maxIter <- 30
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    for(iter in 2:maxIter){
      
      # try({Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)
      # }, silent=T)
      Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('iter=',iter,', dc=', dc, '\n')
      if(dc<1e-9) break
    }
    
    
    
    Hms_int <- Hms
    Hms_int <- rbind(Hms_int, Hmb)
    
    ## Given Hm, update Upsilon
    hesList2 <- getHessian(Xmisb,Hmb, Bms, type, group, rb_vec)
    hesLists1 <- hesLists
    hesLists <- lapply(1:p, function(j) hesLists1[[j]]+hesList2[[j]])
    maxIter <- 30
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    Bm1 <- Bms
    for(iter in 2:maxIter){
      for(j in 1:p){
        Ubj <- getScore(Xmisb[,j], Hmb, Bms[j,], type[group[j]], rb_vec[j])
        Uj <- hesLists1[[j]] %*% (Bm1[j,] - Bms[j,]) + Ubj
        #try({Bms[j, ] <- Bms[j, ] + qr.solve(hesLists[[j]]) %*% Uj}, silent=T)
        #Bms[j, ] <- Bms[j, ] + qr.solve(hesLists[[j]]) %*% Uj
        Bms[j, ] <- Bms[j, ] + MASS::ginv(hesLists[[j]]) %*% Uj
      }
      Bms[is.na(Bms)] <- runif(sum(is.na(Bms)))
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    # sum(abs(Bms[group==2,])>10)
    # Bms[group==2,][abs(Bms[group==2,])>10] <- 0.01
    # re-update H
    maxIter <- 15
    negloglike_seq <- numeric(maxIter)
    negloglike_seq[1] <- 1e10
    for(iter in 2:maxIter){
      #try({Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)}, silent=T)
      Hmb <- updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec,parallel=F,lambda=lambda)
      negloglike_seq[iter] <- objMisfun(Hmb, Bms, Xmisb, gcell, type)
      dc <- abs(negloglike_seq[iter]-negloglike_seq[iter-1])/abs(negloglike_seq[iter-1])
      cat('dc=', dc, '\n')
      if(dc<1e-9) break
    }
    # Hmss <- Hms
    Hms <- rbind(Hms, Hmb)
    
    hXb <- imputeFun(Hmb%*% t(Bms), Xmisb, type, group)
    X_imp[idxList[[b]],] <- hXb
    
    NAE_mat[b,] <- NAE(hXb, X[idxList[[b]],], Xmisb, group)
    AE_mat[b,] <- AE(hXb, X[idxList[[b]],], Xmisb, group)
  }
  hXall <- imputeFun(Hms%*% t(Bms), Xmis, type, group)
  return(list(hX=X_imp, hXall=hXall, Hm=Hms, Bm=Bms, AE_mat=AE_mat, NAE_mat=NAE_mat))
}



