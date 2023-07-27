# # Ex2. mix of normal and Poisson------------------------------------------
setwd("./")
source("./Dfunctions.R")

## Screen ex2.1 run: jj=1, 2 and jme=1(MCAR)
## Screen ex2.2 run: jj=1,2 and jme=2(MAR)
# Given q -----------------------------------------------------------------
library(missMDA)
misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1), seq(0.4-0.1, 0.4+0.1, by=0.1))
jj <- 2
q <- 4
n <- 200; p <- 200; type1 = 'norm_pois'
rho <- 8
N <- 5   ## N=500 is set when run formally
d <- 2 # two types
mecha_set <- c("MCAR", "MAR")
jme <- 2
MethodNames <- c("gfmImpute", "FAMD", "iterFAMD", " LFM-Init", " LFM-EM", "LFM-BR", "LFM-PR")
errMat <- array(0, dim=c(N, length(MethodNames),d) )
colnames(errMat) <- MethodNames
XList <- list()
XmisList <- list()
for(i in 1:N){
  # i <- 10
  cat("i=",i, ", n=", n, ",p=",p, '\n')
  dat <- gendata(seed=i, n=n, p=p, type=type1, 
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  XList[[i]] <- dat$X
  XmisList[[i]] <- dat$Xmis;
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('gaussian','poisson')
  
  # jme = 3, jj=2: lambda=0.01
  # jme = 2, jj=2: lambda =0.01
  out <- try({
    misList <- gfmImpute(dat$Xmis, group, type, q=q, epsLogLike=1e-4, 
                         maxIter=20,verbose=1, parallel=F, lambda=0)
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
  out <- try({
    hX_famd <- imputeFAMD(Xdat, ncp=q, maxiter = 1)$completeObs
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
errMat[1:3,,]
colMeans(errMat)
apply(errMat, c(2,3), mean)
stmp <- paste0('./Rdata_R1/AE_FixQsimu2_misRate',jj,"mech", mecha_set[jme])
save(errMat, file=paste0(stmp,"errMat.Rdata"))

names(XmisList) <- paste0("sample", 1:N)
names(XList) <- paste0("sample", 1:N)
R.matlab::writeMat(paste0(stmp, 'Xlist.mat'), 
                   fixNames=F, XmisList=XmisList,XList=XList,group=group)





# Choose q ----------------------------------------------------------------

misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1), seq(0.4-0.1, 0.4+0.1, by=0.1))
jj <- 2
q <- 4
n <- 200; p <- 200; type1 = 'norm_pois'
rho <- 8
N <- 5  ## N=500 is set when run formally
d <- 2 # two types
mecha_set <- c("MCAR", "MAR")
jme <- 2
MethodNames <- c("gfmImpute")
errMat <- array(0, dim=c(N, length(MethodNames),d) )
colnames(errMat) <- MethodNames
XList <- list()
XmisList <- list()
q_set <- 1:10
qVec <- numeric(N)
for(i in 1:N){
  # i <- 1
  cat('mis=mis',jj,", i=",i, ", n=", n, ",p=",p, '\n')
  dat <- gendata(seed=i, n=n, p=p, type=type1, 
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  # XList[[i]] <- dat$X
  # XmisList[[i]] <- dat$Xmis;
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('gaussian','poisson')
  
  out <- try({
    icMat <- selectFacNumber(dat$Xmis,group, type, q_set=q_set, par.type = 'doSNOW', 
                             parallel = T, verbose=T, lambda=0)
    q_set[which.min(icMat[,'ic'])]
  }, silent = T)
  if(class(out) == 'try-error'){
    errMat[i, 1, ] <- NA
  }else{
    hq <-  out
    qVec[i] <- hq
    cat(" ---------------select q=", hq, '-----------\n')
    out <- try({
      misList <- gfmImpute(dat$Xmis, group, type, q=hq, epsLogLike=1e-4,
                           maxIter=20,verbose=1, lambda=0,parallel=F)
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
  stmp <- paste0('./Rdata_R1/AE_chooseQsimu2_misRate',jj,"mech", mecha_set[jme])
  # save(qVec, errMat, file=paste0(stmp,"errMat.Rdata"))
  
}
errMat[1:2,,]
colMeans(errMat)

stmp <- paste0('./Rdata_R1/AE_chooseQsimu2_misRate',jj,"mech", mecha_set[jme])
save(qVec, errMat, file=paste0(stmp,"errMat.Rdata"))





# Consider more complex case ----------------------------------------------
#  involving complex features such as multiple groups, skewness and heavy tails
gendata_complex <- function(seed=1, n=300, p=50, q=6,mis_vec=0.3, rho=1, mechanism="MCAR"){
    
  #  seed=1; n=300; p=50; q=6;mis_vec=0.3; rho=1; mechanism="MCAR"
  
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
    
    g1 <- 1:floor(p/2)
    g2 <- (floor(p/2)+1):p
    Bm0[g2,-1] <- Bm0[g2,-1] /max(Bm0[g2,-1]) * factor_term/1.5
    mu1 <-  cbind(1,H0) %*% t(Bm0[g1,])
    mu2 <- exp(cbind(1, H0) %*% t(Bm0[g2,]) )
    igroup_index1 <- 1:floor(n/2); igroup_index2 <- (floor(n/2)+1):n
    X_igroup1_p1 <- matrix(rnorm(prod(dim(mu1[igroup_index1,])), mu1[igroup_index1,],1), 
                           length(igroup_index1), length(g1))
    X_igroup2_p2 <- mu1[igroup_index2,] + matrix(rt(prod(dim(mu1[igroup_index2,])), 6), 
                                                 length(igroup_index2), length(g1))
    X_group <- rbind(X_igroup1_p1, X_igroup2_p2)
    X <- cbind(X_group,
               matrix(rpois(prod(dim(mu2)), mu2), n, ncol(mu2)))
    
    
    group1 <- c(rep(1, length(g1)), rep(2, length(g2)))
    Hm0 <- cbind(1, H0)
    Xmis <- X
    
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

# dat <- gendata_complex(n=200, p=200, q=4, rho=c(1, 0.3))
# str(dat)
# Xmis <- dat$Xmis
# group <- dat$group

library(missMDA)
misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1), seq(0.4-0.1, 0.4+0.1, by=0.1))
jj <- 1
q <- 4
n <- 200; p <- 200
rho <- c(8, 0.2)
N <- 5
d <- 2 # two types
mecha_set <- c("MCAR", "MAR")
jme <- 2
MethodNames <- c("gfmImpute", "FAMD", "iterFAMD", " LFM-Init", " LFM-EM", "LFM-BR", "LFM-PR")
errMat <- array(0, dim=c(N, length(MethodNames),d) )
colnames(errMat) <- MethodNames
# XList <- list()
# XmisList <- list()
for(i in 1:N){
  # i <- 10
  cat("i=",i, ", n=", n, ",p=",p, '\n')
  dat <-  gendata_complex(seed=i, n=n, p=p,
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  # XList[[i]] <- dat$X
  # XmisList[[i]] <- dat$Xmis;
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('gaussian','poisson')
  
  # jme = 3, jj=2: lambda=0.01
  # jme = 2, jj=2: lambda =0.01
  out <- try({
    misList <- gfmImpute(dat$Xmis, group, type, q=q, epsLogLike=1e-4, 
                         maxIter=20,verbose=1, parallel=F, lambda=0)
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
  out <- try({
    hX_famd <- imputeFAMD(Xdat, ncp=q, maxiter = 1)$completeObs
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
errMat[1:3,,]
colMeans(errMat)
apply(errMat, c(2,3), mean)
stmp <- paste0('./simu_R1/Rdata_R1/Complex_AE_FixQsimu2_misRate',jj,"mech", mecha_set[jme])
save(errMat, file=paste0(stmp,"errMat.Rdata"))



