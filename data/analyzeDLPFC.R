


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

