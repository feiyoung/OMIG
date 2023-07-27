# library(OMIG)
# example(OrMIG)
#
# q=2;
# epsLogLike=1e-5; maxIter=10;
#
# verbose=FALSE; parallel=FALSE; lambda=1e-5
# misList2 <- OrMIG(Xmis, group,types, q, verbose=TRUE)
#
# misList2 <- OrMIG(Xmis, group,types, q, algorithm='AM', verbose=TRUE)
#
#
#
# # Mix of normal, poisson and binomial
# dat <- gendata(seed=3, n=200, p=200, type='npb',
#                q=6, rho=7, mis_vec = 0.4)
#
# Xmis <- dat$Xmis
# group <- dat$group
# type <- c('gaussian','poisson', 'binomial')
#
# q_set <- 3:8
#
# hq <- selectFacNumber(dat$Xmis, group, type, q_set=q_set,select_method='ratio_test',
#                       parallel = T, verbose=T)
#
# # Select number of factors
# hq <- selectFacNumber(dat$Xmis,group, type, q_set=q_set, select_method= 'IC',
#                          parallel = T, verbose=T)
#
# hq <- selectFacNumber(dat$Xmis,group, type, q_set=q_set, select_method= 'PC',
#                       parallel = T, verbose=T)
#
#
#
# library(OMIG)
# q <- 3
# dat <- gendata(q = q, n=40, p=100, rho=2, mis_vec= c(0.1, 0.2, 0.3))
# dat$Xmis[1:5,1:5]
# type <- 'gaussian'
#
# dat <- gendata(seed=3, n=20, p=20, type='npb',
#                q=6, rho=7, mis_vec = 0.4)
#
# Xmis <- dat$Xmis
# type <- c('gaussian','poisson', 'binomial')
#
# res_mig <- OrMIG(dat$Xmis, group=dat$group, types=type, q=q, algorithm = "VEM", verbose = T)
# NAE(res_mig$hX, dat$X, dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, group=dat$group, types=type, q=q, algorithm = "AM", verbose = T)
#
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type, hist.summary= NULL, lambda=0,
#             verbose=TRUE)
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
#
#
#
# # Gaussian data -----------------------------------------------------------
#
# type <- 'heternorm'
# dat <- gendata(q = q, n=100, p=100, type= type, rho=2, mis_vec= c(0.1, 0.2, 0.3))
# str(dat)
# dat$Xmis[1:5,1:5]
# Xmis_new <- dat$Xmis;
# q <- 4
# group <- dat$group;
#
# hist.summary= NULL;
# lambda=1;
# verbose=TRUE
# types <- 'gaussian'
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=lambda,
#             verbose=TRUE)
#
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 2
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
# res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=0,
#              verbose=TRUE)
#
# ## Normalized absolute error
# NAE(res2$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "VEM", verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "AM", verbose = T)
# NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 3
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
# res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary)
#
# ## Normalized absolute error
# NAE(res3$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
#
#
#
#
# #### Binomial data#######
# type <- 'bino'
#
# dat <- gendata(q = q, n=500, p=500, type= type, rho=2, mis_vec= c(0.1, 0.2, 0.3))
# str(dat)
# dat$Xmis[1:5,1:5]
# Xmis_new <- dat$Xmis;
# q <- 4
# group <- dat$group;
# types = types;
# hist.summary= NULL;
# lambda=1;
# verbose=TRUE
# types <- 'binomial'
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=1,
#             verbose=TRUE)
# GFM::measurefun(res$Hm, dat$Hm0)
# GFM::measurefun(res$hist.summary$Bms, dat$Bm0)
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# library(MASS)
# GFM::measurefun(cbind(1,res_mig$fitList$H), dat$Hm0)
# GFM::measurefun(cbind(res_mig$fitList$mu,res_mig$fitList$B), dat$Bm0)
#
#
# ## Batch 2
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
# res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=0,
#                         verbose=TRUE)
#
# ## Normalized absolute error
# NAE(res2$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 3
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
# res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary)
#
# ## Normalized absolute error
# NAE(res3$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM")
#
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
#
#
#
# #### Poisson data#######
# type <- 'pois'
# dat <- gendata(q = q, n=100, p=100, type= type, rho=2, mis_vec= c(0.1, 0.2, 0.3))
# str(dat)
# dat$Xmis[1:5,1:5]
# Xmis_new <- dat$Xmis;
# q <- 4
# group <- dat$group;
#
# #### Which indicates the advantage of online missing imputation!!!
# hist.summary= NULL;
# lambda=1;
# verbose=TRUE
# types <- 'poisson'
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=lambda,
#             verbose=TRUE)
#
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 2
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
# res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=0,
#              verbose=TRUE)
#
# ## Normalized absolute error
# NAE(res2$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "VEM", verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "AM", verbose = T)
# NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 3
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
# res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary)
#
# ## Normalized absolute error
# NAE(res3$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# #### Normal and Poisson mixed data#######
# q <- 4
# type <- "norm_pois"
# dat <- gendata(q = q, n=100, p=100, type= type, rho=2, mis_vec= c(0.1, 0.2, 0.3))
# str(dat)
# dat$Xmis[1:5,1:5]
# Xmis_new <- dat$Xmis;
# q <- 4
# group <- dat$group;
#
# #### Which indicates the advantage of online missing imputation!!!
# hist.summary= NULL;
# lambda=1;
# verbose=TRUE
# types <- c('gaussian', 'poisson')
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=lambda,
#             verbose=TRUE)
#
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = F)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = F)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 2
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
# res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=0,
#              verbose=TRUE)
#
# ## Normalized absolute error
# NAE(res2$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "VEM", verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# # res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "AM", verbose = T)
# # NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 3
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
# res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary)
#
# ## Normalized absolute error
# NAE(res3$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM")
# NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
#
# #### Poisson and binomial mixed data#######
# type <- "pois_bino"
# dat <- gendata(q = q, n=100, p=100, type= type, rho=2, mis_vec= c(0.1, 0.2, 0.3))
# str(dat)
# dat$Xmis[1:5,1:5]
# Xmis_new <- dat$Xmis;
# q <- 4
# group <- dat$group;
#
# #### Which indicates the advantage of online missing imputation!!!
# hist.summary= NULL;
# lambda=1;
# verbose=TRUE
# types <- c('poisson', 'binomial')
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=lambda,
#             verbose=TRUE)
#
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = F)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = F)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 2
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
# res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=0,
#              verbose=TRUE)
#
# ## Normalized absolute error
# NAE(res2$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "VEM", verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# # res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "AM", verbose = T)
# # NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 3
# dat <- gendata(q = q, n=40, p=100, rho=2,type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
# res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary)
#
# ## Normalized absolute error
# NAE(res3$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM")
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM")
# NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
#
#
# #### Normal, Poisson and binomial mixed data#######
# type <- "npb"
# dat <- gendata(q = q, n=100, p=100, type= type, rho=c(2, 0.1), mis_vec= c(0.1, 0.2, 0.3))
# str(dat)
# dat$Xmis[1:5,1:5]
# Xmis_new <- dat$Xmis;
# q <- 4
# group <- dat$group;
#
# #### Which indicates the advantage of online missing imputation!!!
# hist.summary= NULL;
# lambda=1;
# verbose=TRUE
# types <- c('gaussian','poisson', 'binomial')
# res <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, type= types, hist.summary= NULL, lambda=lambda,
#             verbose=TRUE)
#
# NAE(res$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'AM', verbose = F)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# res_mig <- OrMIG(dat$Xmis, dat$group, type= types, q, algorithm = 'VEM', verbose = F)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 2
# dat <- gendata(q = q, n=40, p=100, rho=c(2, 0.1),type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 2)
# res2 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res$hist.summary, lambda=1,
#              verbose=TRUE)
#
# ## Normalized absolute error
# NAE(res2$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "VEM", verbose = T)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, algorithm = "AM", verbose = T)
# NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
# ## Batch 3
# dat <- gendata(q = q, n=40, p=100, rho= c(2, 0.1),type= type, mis_vec= c(0.1, 0.2, 0.3), seed = 5)
# res3 <- OMIG(Xmis_new=dat$Xmis, q,  dat$group, types, hist.summary= res2$hist.summary, lambda=1)
#
# ## Normalized absolute error
# NAE(res3$hX, dat$X, dat$Xmis, dat$group)
# ## compare with the offline method MIG
# res_mig <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "VEM", lambda = 1)
# NAE(res_mig$hX, dat$X,  dat$Xmis, dat$group)
# res_mig2 <- OrMIG(dat$Xmis, dat$group, types, q, verbose = T, algorithm = "AM", lambda=1)
# NAE(res_mig2$hX, dat$X,  dat$Xmis, dat$group)
#
#
