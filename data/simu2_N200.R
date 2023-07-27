# Comparison of the imputation accuracy of OMIG and OfMIG for each data batch--------


setwd("./")
source("./Dfunctions.R")

misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1)) # , seq(0.4-0.1, 0.4+0.1, by=0.1)
jj <- 1
q <- 4
n <- 200; p <- 200; type1 = 'norm_pois'
rho <- c(8, 0.2) 
N <- 500
d <- 2 # two types
mecha_set <- c("MCAR", "MAR")
jme <- 2
MethodNames <- c("OfMIG", "OMIG")

idxList <- batchGroup(n, nbatch1=0, nbatch = 40)
J <- length(idxList)
errArray <- array(NA, dim=c(J,d, 2,  N) )
colnames(errArray) <- c("Normal", "Poisson")
AE_errArray <- errArray
corBArray <- array(NA, dim=c(J,2, N) )
colnames(corBArray) <-  c("Normal", "Poisson")
lambda <- 0
for(i in 1:N){
  # i <- 1
  cat("i=",i, ", n=", n, ",p=",p, '\n')
  dat <- gendata(seed=i, n=n, p=p, type=type1, 
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  Xmis <- dat$Xmis
  group <- dat$group
  type <- c('gaussian', 'poisson')
  
  
  # jme = 3, jj=2: lambda=0.01
  # jme = 2, jj=2: lambda =0.01
  out <- try({
    misList <- offlineGFMImpute(dat$Xmis, q,  group, type,idxList ,lambda=lambda, verbose=T,
                                X=dat$X, Bm=dat$Bm0)
    errArray[,,1,i] <- misList$NAE_mat
    AE_errArray[,,1,i] <- misList$AE_mat
    corBArray[,1,i] <- misList$corB_vec
    
    
  }, silent = T)
  
  out <- try({
    sList <- streamGFMImpute(dat$Xmis, q,  group, type, idxList, lambda=lambda,
                             verbose=T, X= dat$X, Bm=dat$Bm0)
    errArray[,,2,i] <- sList$NAE_mat
    AE_errArray[,,2,i] <- sList$AE_mat
    corBArray[,2,i] <- abs(sList$corB_vec)
    
  }, silent = T)
  
  AE_errArray[,,,i]
}

row.names(errArray) <- paste0('b',1:J)
row.names(AE_errArray) <- paste0('b',1:J)

stmp <- paste0('./simu_R1/Rdata_R1/AE_StreamFixQsimu2_misRate',jj,"mech", mecha_set[jme])
save(errArray, AE_errArray, corBArray, file=paste0(stmp, '.Rdata') )
AE_errArray[AE_errArray>1e5] <- NA
apply(AE_errArray, c(1:3), mean, na.rm=T)
apply(AE_errArray, c(1:3), mean)
