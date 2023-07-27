
####
##Comparison of the estimation accuracy of OMIG and OrMIG. 

setwd("./")
source("D:/LearnFiles/Research paper/idea/2. OMIGimpute/Rcode/simu/submit_simucode/Dfunctions.R")


SMAE <- function(hX, X,Xmis, group){
  hmu <- colMeans(Xmis, na.rm=T)
  Mu <- matrix(hmu, nrow(Xmis), ncol(Xmis), byrow=T)
  d <- length(unique(group))
  XList <- list()
  mea <- numeric(d)
  for(i in 1:d){
    mind1 <- is.na(Xmis[,group==i] )
    tmp <- 0
    for(j in 1:ncol(mind1)){
      tmp <- tmp + sum(abs((hX[,group==i][mind1[,j]==1,j] - X[, group==i][mind1[,j]==1,j])))/
        sum(abs((Mu[,group==i][mind1[,j]==1,j] - X[, group==i][mind1[,j]==1,j])))
    }
    
    mea[i] <- tmp/ ncol(mind1)
  }
  return(mea)
}
MAE <- function(hX, X, Xmis, group){
  
  d <- length(unique(group))
  mea <- numeric(d)
  for(i in 1:d){
    mind1 <- is.na(Xmis[,group==i] )
    tmp <- 0
    for(j in 1:ncol(mind1)){
      tmp <- tmp + sum(abs((hX[,group==i][mind1[,j]==1,j] - X[, group==i][mind1[,j]==1,j])))
    }
    
    mea[i] <- tmp/ sum(mind1)
  }
  return(mea)
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
          misj_ind <- order(XX[,j], decreasing = T)[1:floor(n*mis_rate)]
          Xmis[misj_ind, j] <- NA
        }
      }
    }
  }
  
  return(list(Xmis=Xmis, X=X, Bm0=Bm0, Hm0=Hm0, group=group1))
}
# Simulation --------------------------------------------------------------


misrateList<- list(seq(0.2-0.1, 0.2+0.1,by=0.1), seq(0.4-0.1, 0.4+0.1, by=0.1))
jj <- 1
q <- 4
n <- 500; p <- 200; type1 = 'norm_pois'
rho <- 8
N <- 5 # run 5 times for example
d <- 2 # two types
mecha_set <- c("MCAR", "MAR")
jme <- 1
MethodNames <- c("GFM-OL", "GFM-ST")

## consider nbatch=60 or 30.
idxList <- batchGroup(n, nbatch1=300, nbatch = 30)
J <- length(idxList)
errArray <- array(NA, dim=c(J,d, 2,  N) )
colnames(errArray) <- c("Normal", "Poisson")
corBArray <- array(NA, dim=c(J,2, N) )
colnames(corBArray) <-  c("OrMIG", "OMIG")
errDArray <- array(NA, dim=c(J,2, N) )
colnames(errDArray) <-  c("OrMIG", "OMIG")
lambda <- 0
for(i in 1:N){
  # i <- 1
  cat("i=",i, ", n=", n, ",p=",p, '\n')
  dat <- gendata(seed=i, n=n, p=p, type=type1, 
                 q=q, rho=rho, mis_vec = misrateList[[jj]], mechanism = mecha_set[jme])
  Xmis <- dat$Xmis
  apply(Xmis, 2, function(x) sum(is.na(x))/n)
  group <- dat$group
  type <- c('gaussian', 'poisson')
  D0 <- dat$Hm0%*% t(dat$Bm0)
  
  # jme = 3, jj=2: lambda=0.01
  # jme = 2, jj=2: lambda =0.01
  corB_vec <- numeric(J)
  
  errD_vec <- numeric(J)
  NAE_mat <- matrix(NA, J, d)
  MAE_mat <- matrix(NA, J, d)
  for(j in 1:J){
    # j <- 1
    idx <- 1: tail(idxList[[j]],1)
    out <- try({
      misList <- gfmImpute(dat$Xmis[idx,], group, type,q,lambda=lambda, verbose=T,
                           maxIter=20,epsLogLike=1e-4)
      
      
      corB_vec[j] <- cancor(misList$hBm, dat$Bm0)$cor[q+1]
      
      NAE_mat[j, ] <- NAE(misList$hX[idxList[[j]],], dat$X[idxList[[j]],],
                           dat$Xmis[idxList[[j]],], group)
      
      MAE_mat[j, ] <- MAE(misList$hX, dat$X[idxList[[j]],],
          dat$Xmis[idxList[[j]],], group)
      errD_vec[j] <- sum(abs((misList$hHm[idx,] %*% t(misList$hBm) - D0[idx,])))/ prod(dim( D0[idx,]))
    }, silent = T)
  }
  corBArray[,1,i] <- corB_vec
  errArray[,,1,i] <- NAE_mat
  errDArray[,1,i] <- errD_vec
  out <- try({ ## Add AE to this function
    sList <- streamGFMImpute2(dat$Xmis, q,  group, type, idxList, lambda=lambda,
                              verbose=T, X=dat$X, Bm0=dat$Bm0, H0=dat$Hm0[,-1])
  }, silent = T)
  
  corBArray[,2,i] <- sList$corB_vec
  errArray[,,2,i] <- sList$NAE_mat
  errDArray[,2,i] <- sList$errD_vec
}
apply(errArray, c(1,2,3), mean, na.rm=T)
apply(errDArray, c(1,2), mean, na.rm=T)
apply(corBArray, c(1,2), mean, na.rm=T)
row.names(errArray) <- paste0('b',1:J)
boxplot(t(errArray[1,1,,]))

row.names(errArray) <- paste0('b',1:J)
stmp <- paste0('./Rdata/StreamOffLineBHm_nbatch60_simu2_misRate',jj,"mech", mecha_set[jme])
save(errDArray, errArray, corBArray, file=paste0(stmp, '.Rdata') )
