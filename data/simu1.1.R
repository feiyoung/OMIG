# # Ex1. mix of Poisson and Binomial -------------------------------------------
setwd("./")
source("D:/LearnFiles/Research paper/idea/2. OMIGimpute/Rcode/simu/submit_simucode/Dfunctions.R")


# Given q -----------------------------------------------------------------
## Screen ex1.1 run: jj=1, 2 and jme=1(MCAR)
## Screen ex1.2 run: jj=1,2 and jme=2(MAR)
library(missMDA)
misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1), seq(0.4-0.1, 0.4+0.1, by=0.1))
jj <- 2
q <- 3
n <- 200; p <- 200; type1 = "pois_bino"
rho <- 8
N <- 5 ## try 5 runs; more runs will produce more stable results
d <- 2 # two types
mecha_set <- c("MCAR", "MAR") ## 4 cases
jme <- 2
MethodNames <- c("gfmImpute", "FAMD", "iterFAMD", " LFM-Init", " LFM-EM", "LFM-BR", "LFM-PR")
errMat <- array(0, dim=c(N, length(MethodNames),d) )
colnames(errMat) <- MethodNames

for(i in 1:N){
  # i <- 1
  cat("i=",i, ", n=", n, ",p=",p, '\n')
  dat <- gendata(seed=i, n=n, p=p, type=type1, 
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('poisson', 'binomial')
  
  out <- try({
    misList <- gfmImpute(dat$Xmis, group, type, q=q, epsLogLike=1e-4, 
                         maxIter=20,verbose=1, parallel=F, lambda=0.1)
    hX <- misList$hX
    mind<- is.na(dat$Xmis)
    AE(hX, dat$X, dat$Xmis, group)
  }, silent = T)
  if(class(out) == 'try-error'){
    errMat[i, 1, ] <- NA
  }else{
    errMat[i, 1, ] <-  out
  }
  
  # FAMD and iterative FAMD
  Xdat <- as.data.frame(dat$Xmis)
  id_bino <- which(group == 2)
  for(j in id_bino) Xdat[,j] <- as.factor(Xdat[,j])
  out <- try({
    hX_famd <- imputeFAMD(Xdat, ncp=q, maxiter = 1)$completeObs
    for(j in id_bino) hX_famd[,j] <- as.numeric(hX_famd[, j]) -1
    hX_famd_mat <- as.matrix(hX_famd)
    AE(hX_famd_mat, dat$X, dat$Xmis, group)
  }, silent = T)
  if(class(out) == 'try-error'){
    errMat[i, 2, ] <- NA
  }else{
    errMat[i, 2, ] <-   out
  }
  
  out <- try({
    hX_iterfamd <- imputeFAMD(Xdat, ncp=q)$completeObs
    for(j in id_bino) hX_iterfamd[,j] <- as.numeric(hX_iterfamd[, j]) -1
    hX_iterfamd_mat <- as.matrix(hX_iterfamd)
    AE(hX_iterfamd_mat, dat$X, dat$Xmis, group)
  }, silent = T)
  if(class(out) == 'try-error'){
    errMat[i, 3, ] <- NA
  }else{
    errMat[i, 3, ] <-   out
  }
  
  ## JMS method
  hX_JMS_int <- LFM.JMS(dat$Xmis, q, EM=F)
  errMat[i, 4, ] <- AE(hX_JMS_int, dat$X, dat$Xmis, group)
  hX_JMS_final <- LFM.JMS(dat$Xmis, q,maxIter = 5, EM=T)
  errMat[i, 5, ] <- AE(hX_JMS_final, dat$X, dat$Xmis, group)
  
  # PX method
  hX_XP_int <- LFM.XP(dat$Xmis, q, prop.weighted = F)
  errMat[i, 6,] <- AE(hX_XP_int, dat$X, dat$Xmis, group)
  hX_XP_final <- LFM.XP(dat$Xmis, q,prop.weighted = T)
  errMat[i, 7,] <- AE(hX_XP_final, dat$X, dat$Xmis, group)
  
}
errMat[1:2,,]
colMeans(errMat)
stmp <- paste0('./AE_FixQsimu1_misRate',jj,"mech", mecha_set[jme])
save(errMat, file=paste0(stmp,"errMat.Rdata"))


# Select q ----------------------------------------------------------------

misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1), seq(0.4-0.1, 0.4+0.1, by=0.1))
jj <- 1
q <- 3
n <- 200; p <- 200; type1 = "pois_bino"
rho <- 8
N <- 5 ## N=500 is set when run formally
d <- 2 # two types
mecha_set <- c("MCAR", "MAR")
jme <- 1
MethodNames <- c("gfmImpute")
errMat <- array(0, dim=c(N, length(MethodNames),d) )
colnames(errMat) <- MethodNames
# XList <- list()
# XmisList <- list()
q_set <- 1:6
qVec <- numeric(N)
for(i in 1:N){
  # i <- 1
  tic <- proc.time()
  cat("i=",i, ", n=", n, ",p=",p, '\n')
  dat <- gendata(seed=i, n=n, p=p, type=type1, 
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  # XList[[i]] <- dat$X
  # XmisList[[i]] <- dat$Xmis;
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('poisson', 'binomial')
  
  out <- try({
    icMat <- selectFacNumber(dat$Xmis,group, type, q_set=q_set, par.type = 'doSNOW', 
                             parallel = T, verbose=T, lambda=0) # 0.01
    q_set[which.min(icMat[,'ic'])]
  }, silent = T)
  if(class(out) == 'try-error'){
    errMat[i, 1, ] <- NA
  }else{
    hq <-  out
    qVec[i] <- hq
    cat(" ---------------select q=", hq, '-----------\n')
    out <- try({
      misList <- gfmImpute(dat$Xmis, group, type, q=hq, epsLogLike=1e-4, maxIter=20,
                           verbose=1, parallel=F, lambda=0) # 1e-20
      hX <- misList$hX
      mind<- is.na(dat$Xmis)
      AE(hX, dat$X, dat$Xmis, group)
    }, silent = T)
    if(class(out) == 'try-error'){
      errMat[i, 1, ] <- NA
    }else{
      errMat[i, 1, ] <-  out
    }
  }
  toc <- proc.time()
}
errMat[1:2,,]
colMeans(errMat)
stmp <- paste0('./Rdata_R1/AE_chooseQsimu1_misRate',jj,"mech", mecha_set[jme])
save(qVec, errMat, file=paste0(stmp,"errMat.Rdata"))






