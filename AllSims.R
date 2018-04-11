
runSims <- function(param, nrep, delta, n, B, MR, b){
  shape <- as.numeric(param[1])
  cenprob <- as.numeric(param[2])
  usehaz <- param[3]
  path <- "~/../Dropbox/EfficientLinYing/Code/"
  #library(Rcpp)
  library(matrixStats)
  Rcpp::sourceCpp(paste(path,"fillmatKunequal.cpp", sep=""))
  source(paste(path,"estBetaSmoothUnequal3.R", sep=""))
  
  prepdata <- function(T0, T1, ev0, ev1){
    n0 <- length(T0)
    n1 <- length(T1)
    n <- n0+n1
    X <- c(T0, T1)
    ev <- c(ev0, ev1)
    numcen <- sum(1-ev)
    X[duplicated(X)] <- X[duplicated(X)] + runif(sum(duplicated(X)), 0, min(X/1e6))
    Z <- c(rep(0, n0), rep(1, n1))
    
    obsdf <- data.frame(X = X, Z= Z, ev=ev)
    
    sortlist <- obsdf[order(obsdf$X),]
    return(sortlist )
  }
  K <- 10
  
  #delta <- 1
  alpha <- 0.05
  #source(paste(path,"speedEx.R", sep=""))
  
  betas <- rep(0, B)
  betahats <- rep(0, nrep)
  LYbeta <- vhat <- rep(0, nrep)
  CIs <- matrix(nrow=nrep, ncol=2)
  CIsAH <- CIs2 <- CIs
  Bshapes <- Bscales <- matrix(nrow=B, ncol=K)
  
  for(i in 1:nrep){
    
    if(usehaz=="weibull"){
      T0 <- rweibull(n, shape=shape, scale=1)
      T1 <- pmin(rweibull(n, shape=shape, scale=1), rexp(n,delta))
      C0 <- rexp(n,1/qweibull(cenprob, shape=shape, scale=1))
      C1 <- rexp(n,1/qweibull(cenprob, shape=shape, scale=1))
    }else if(usehaz=="lognormal"){
      T0 <- rlnorm(n, sdlog=shape, meanlog=0)
      T1 <- pmin(rlnorm(n, sdlog=shape, meanlog=0), rexp(n,delta))
      C0 <- rexp(n,1/qlnorm(cenprob, sdlog=shape, meanlog=0))
      C1 <- rexp(n,1/qlnorm(cenprob, sdlog=shape, meanlog=0))
    }else if(usehaz=="sine"){
      T0 <- rsin(n)
      T1 <- pmin(rsin(n), rexp(n, delta))
      meanCen <- optimize(function(t) (cenprob-Fsin(t))^2, lower=0, 5)$minimum
      C0 <- rexp(n, 1/meanCen)
      C1 <- rexp(n, 1/meanCen)
    }
    ev0 <- as.numeric(T0 < C0)
    ev1 <- as.numeric(T1 < C1)
    #ev0 <- ev1 <- rep(1, n)
    sortlist <- prepdata(pmin(T0,C0), pmin(T1,C1), ev0, ev1)
    groups <- rep(1:K, times=(2*n/K))
    sortlist$groups <- groups
    betas <- rep(0, B)
    
    ahfit <-  estBetaLY(sortlist)
    ahbeta <- ahfit$beta
    
    b <- 2*n
    for(j in 1:B){
      subsmp <- sort(sample(1:(2*n), b, replace=TRUE))
      
      sl1 <- sortlist[subsmp,]
      dupes <- duplicated(sl1$X)
      sl1$X[dupes] <- sl1$X[dupes]+sort(runif(dupes, 0, 1e-6))
      
      groups <- rep((1:K), ceiling(b/K))[1:b]
      
      betas[j] <- estBetaSmoothU(sl1, groups, rep(1,b), MR)
      
      
      if(is.nan(betas[j])) failSL <- sl1
      
    }
    betahat <- estBetaSmoothU(sortlist, sortlist$groups, rep(1, (2*n)), MR)
    
    betahats[i] <- betahat
    
    #crit.val <- quantile(betas-betahat, probs=c(1-alpha/2, alpha/2))#*(sqrt(2*n)/(sqrt(2*n) + sqrt(SN)))
    #CIs[i,] <- betahat - crit.val
    CIs[i,] <- quantile(betas, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE)
    
    #ahse <- ahaz:::summary.ahaz(ahfit)$coefficients[2]
    #CIsAH[i,1] <- ahbeta-1.96*ahse
    #CIsAH[i,2] <- ahbeta+1.96*ahse
    LYbeta[i] <- ahbeta
    vhat[i] <- ahfit$B/ahfit$A^2
  }
  
  return(cbind( betahats, CIs, LYbeta, vhat))
}

library(parallel)

#shape <-c(1/2, 2/3, 3/4, 1, 1.5, 2, 3)
#shape <-c(1/2,3)
cenprob <- c(0.2,  0.4,  0.6, 0.8)
#cenprob <- c(0.4, 0.6)
shape <- c(1/3, 1/2, 2/3, 1, 1.5, 2, 2.5, 3)

#cenprob <- c(0.4, 0.6, 0.8)
param <- cbind(rep(shape, length(cenprob)) ,rep(cenprob, each=length(shape)), rep("weibull", length(shape)*length(cenprob)))
param <- rbind(param, cbind(rep(shape, length(cenprob)) ,rep(cenprob, each=length(shape)), rep("lognormal", length(shape)*length(cenprob))))
param <- rbind(param, cbind(rep(1, length(cenprob)), cenprob, rep("sine", length(cenprob))))
paramlist <- vector("list", length(cenprob)*length(shape))
for(i in 1:dim(param)[1]){
  paramlist[[i]] <- param[i,]
}


n <-  100
#nrep <- 2000
#B <- 2000
nrep <- 10
delta <- 0.1
B <- 1000
MR <- 10
b <- floor((2*n)^.75)

nclus <- detectCores()
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 567)
clusterEvalQ(cl, library(MASS))
test1 <- parLapply(cl=cl, paramlist[1:12], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test1) <- as.character(paramlist)[1:12]
# 
save(test1, file="n100Del01part1.RData")


nclus <- detectCores()
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 678)
clusterEvalQ(cl, library(MASS))
test2 <- parLapply(cl=cl, paramlist[13:24], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test2) <- as.character(paramlist)[13:24]
# 
save(test2, file="n100Del01part2.RData")

nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 789)
clusterEvalQ(cl, library(MASS))
test3 <- parLapply(cl=cl, paramlist[25:35], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test3) <- as.character(paramlist)[25:35]
# 
save(test3, file="n100Del01part3.RData")


nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 891)
clusterEvalQ(cl, library(MASS))
test4 <- parLapply(cl=cl, paramlist[36:46], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test4) <- as.character(paramlist)[36:46]
# 
save(test4, file="n100Del01part4.RData")


nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 912)
clusterEvalQ(cl, library(MASS))
test5 <- parLapply(cl=cl, paramlist[47:57], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test5) <- as.character(paramlist)[47:57]
# 
save(test5, file="n100Del01part5.RData")

nclus <- 11
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 123)
clusterEvalQ(cl, library(MASS))
test6 <- parLapply(cl=cl, paramlist[58:68], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test6) <- as.character(paramlist)[58:68]
# 
save(test6, file="n100Del01part6.RData")

#########################################################################
#########################################################################
#########################################################################
#########################################################################

runSims <- function(param, nrep, delta, n, B, MR, b){
  shape <- as.numeric(param[1])
  cenprob <- as.numeric(param[2])
  usehaz <- param[3]
  path <- "~/../Dropbox/EfficientLinYing/Code/"
  #library(Rcpp)
  library(matrixStats)
  Rcpp::sourceCpp(paste(path,"fillmatKunequal.cpp", sep=""))
  source(paste(path,"estBetaSmoothUnequal3.R", sep=""))
  
  prepdata <- function(T0, T1, ev0, ev1){
    n0 <- length(T0)
    n1 <- length(T1)
    n <- n0+n1
    X <- c(T0, T1)
    ev <- c(ev0, ev1)
    numcen <- sum(1-ev)
    X[duplicated(X)] <- X[duplicated(X)] + runif(sum(duplicated(X)), 0, min(X/1e6))
    Z <- c(rep(0, n0), rep(1, n1))
    
    obsdf <- data.frame(X = X, Z= Z, ev=ev)
    
    sortlist <- obsdf[order(obsdf$X),]
    return(sortlist )
  }
  K <- 10
  
  #delta <- 1
  alpha <- 0.05
  #source(paste(path,"speedEx.R", sep=""))
  
  betas <- rep(0, B)
  betahats <- rep(0, nrep)
  LYbeta <- vhat <- rep(0, nrep)
  CIs <- matrix(nrow=nrep, ncol=2)
  CIsAH <- CIs2 <- CIs
  Bshapes <- Bscales <- matrix(nrow=B, ncol=K)
  
  for(i in 1:nrep){
    
    if(usehaz=="weibull"){
      T0 <- rweibull(n, shape=shape, scale=1)
      T1 <- pmin(rweibull(n, shape=shape, scale=1), rexp(n,delta))
      C0 <- rexp(n,1/qweibull(cenprob, shape=shape, scale=1))
      C1 <- rexp(n,1/qweibull(cenprob, shape=shape, scale=1))
    }else if(usehaz=="lognormal"){
      T0 <- rlnorm(n, sdlog=shape, meanlog=0)
      T1 <- pmin(rlnorm(n, sdlog=shape, meanlog=0), rexp(n,delta))
      C0 <- rexp(n,1/qlnorm(cenprob, sdlog=shape, meanlog=0))
      C1 <- rexp(n,1/qlnorm(cenprob, sdlog=shape, meanlog=0))
    }else if(usehaz=="sine"){
      T0 <- rsin(n)
      T1 <- pmin(rsin(n), rexp(n, delta))
      meanCen <- optimize(function(t) (cenprob-Fsin(t))^2, lower=0, 5)$minimum
      C0 <- rexp(n, 1/meanCen)
      C1 <- rexp(n, 1/meanCen)
    }
    ev0 <- as.numeric(T0 < C0)
    ev1 <- as.numeric(T1 < C1)
    #ev0 <- ev1 <- rep(1, n)
    sortlist <- prepdata(pmin(T0,C0), pmin(T1,C1), ev0, ev1)
    groups <- rep(1:K, times=(2*n/K))
    sortlist$groups <- groups
    betas <- rep(0, B)
    
    ahfit <-  estBetaLY(sortlist)
    ahbeta <- ahfit$beta
    
    b <- 2*n
    for(j in 1:B){
      subsmp <- sort(sample(1:(2*n), b, replace=TRUE))
      
      sl1 <- sortlist[subsmp,]
      dupes <- duplicated(sl1$X)
      sl1$X[dupes] <- sl1$X[dupes]+sort(runif(sum(dupes), 0, 1e-6))
      
      groups <- rep((1:K), ceiling(b/K))[1:b]
      
      betas[j] <- estBetaSmoothU(sl1, groups, rep(1,b), MR)
      
      
      if(is.nan(betas[j])) failSL <- sl1
      
    }
    betahat <- estBetaSmoothU(sortlist, sortlist$groups, rep(1, (2*n)), MR)
    
    betahats[i] <- betahat
    
    #crit.val <- quantile(betas-betahat, probs=c(1-alpha/2, alpha/2))#*(sqrt(2*n)/(sqrt(2*n) + sqrt(SN)))
    #CIs[i,] <- betahat - crit.val
    CIs[i,] <- quantile(betas, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE)
    
    #ahse <- ahaz:::summary.ahaz(ahfit)$coefficients[2]
    #CIsAH[i,1] <- ahbeta-1.96*ahse
    #CIsAH[i,2] <- ahbeta+1.96*ahse
    LYbeta[i] <- ahbeta
    vhat[i] <- ahfit$B/ahfit$A^2
  }
  
  return(cbind( betahats, CIs, LYbeta, vhat))
}

library(parallel)

#shape <-c(1/2, 2/3, 3/4, 1, 1.5, 2, 3)
#shape <-c(1/2,3)
cenprob <- c(0.2,  0.4,  0.6, 0.8)
#cenprob <- c(0.4, 0.6)
shape <- c(1/3, 1/2, 2/3, 1, 1.5, 2, 2.5, 3)

#cenprob <- c(0.4, 0.6, 0.8)
param <- cbind(rep(shape, length(cenprob)) ,rep(cenprob, each=length(shape)), rep("weibull", length(shape)*length(cenprob)))
param <- rbind(param, cbind(rep(shape, length(cenprob)) ,rep(cenprob, each=length(shape)), rep("lognormal", length(shape)*length(cenprob))))
param <- rbind(param, cbind(rep(1, length(cenprob)), cenprob, rep("sine", length(cenprob))))
paramlist <- vector("list", length(cenprob)*length(shape))
for(i in 1:dim(param)[1]){
  paramlist[[i]] <- param[i,]
}


n <-  100
#nrep <- 2000
#B <- 2000
nrep <- 2000
delta <- 0.3
B <- 1000
MR <- 10
b <- floor((2*n)^.75)

nclus <- detectCores()
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 567)
clusterEvalQ(cl, library(MASS))
test1 <- parLapply(cl=cl, paramlist[1:12], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test1) <- as.character(paramlist)[1:12]
# 
save(test1, file="n100Del03part1.RData")


nclus <- detectCores()
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 678)
clusterEvalQ(cl, library(MASS))
test2 <- parLapply(cl=cl, paramlist[13:24], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test2) <- as.character(paramlist)[13:24]
# 
save(test2, file="n100Del03part2.RData")

nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 789)
clusterEvalQ(cl, library(MASS))
test3 <- parLapply(cl=cl, paramlist[25:35], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test3) <- as.character(paramlist)[25:35]
# 
save(test3, file="n100Del03part3.RData")


nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 891)
clusterEvalQ(cl, library(MASS))
test4 <- parLapply(cl=cl, paramlist[36:46], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test4) <- as.character(paramlist)[36:46]
# 
save(test4, file="n100Del03part4.RData")


nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 912)
clusterEvalQ(cl, library(MASS))
test5 <- parLapply(cl=cl, paramlist[47:57], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test5) <- as.character(paramlist)[47:57]
# 
save(test5, file="n100Del03part5.RData")

nclus <- 11
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 123)
clusterEvalQ(cl, library(MASS))
test6 <- parLapply(cl=cl, paramlist[58:68], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test6) <- as.character(paramlist)[58:68]
# 
save(test6, file="n100Del03part6.RData")


#########################################################################
#########################################################################
#########################################################################
#########################################################################

runSims <- function(param, nrep, delta, n, B, MR, b){
  shape <- as.numeric(param[1])
  cenprob <- as.numeric(param[2])
  usehaz <- param[3]
  path <- "~/../Dropbox/EfficientLinYing/Code/"
  #library(Rcpp)
  library(matrixStats)
  Rcpp::sourceCpp(paste(path,"fillmatKunequal.cpp", sep=""))
  source(paste(path,"estBetaSmoothUnequal3.R", sep=""))
  
  prepdata <- function(T0, T1, ev0, ev1){
    n0 <- length(T0)
    n1 <- length(T1)
    n <- n0+n1
    X <- c(T0, T1)
    ev <- c(ev0, ev1)
    numcen <- sum(1-ev)
    sl1$X[dupes] <- sl1$X[dupes]+sort(runif(sum(dupes), 0, 1e-6))
    Z <- c(rep(0, n0), rep(1, n1))
    
    obsdf <- data.frame(X = X, Z= Z, ev=ev)
    
    sortlist <- obsdf[order(obsdf$X),]
    return(sortlist )
  }
  K <- 10
  dloglogis <- function(t, scale, shape){
    ((shape/scale)*(t/scale)^(shape-1))/(1+(t/scale)^shape)^2
  }
  
  ploglogis <- function(t, scale, shape){
    1/(1+(t/scale)^(-shape))
  }
  
  qloglogis <- function(p, scale, shape){
    scale*(p/(1-p))^(1/shape)
  }
  
  rloglogis <- function(n, scale, shape){
    u <- runif(n)
    scale*(u/(1-u))^(1/shape)
  }
  #delta <- 1
  alpha <- 0.05
  #source(paste(path,"speedEx.R", sep=""))
  
  betas <- rep(0, B)
  betahats <- rep(0, nrep)
  LYbeta <- vhat <- rep(0, nrep)
  CIs <- matrix(nrow=nrep, ncol=2)
  CIsAH <- CIs2 <- CIs
  Bshapes <- Bscales <- matrix(nrow=B, ncol=K)
  
  for(i in 1:nrep){
    
    if(usehaz=="weibull"){
      T0 <- rweibull(n, shape=shape, scale=1)
      T1 <- pmin(rweibull(n, shape=shape, scale=1), rexp(n,delta))
      C0 <- rexp(n,1/qweibull(cenprob, shape=shape, scale=1))
      C1 <- rexp(n,1/qweibull(cenprob, shape=shape, scale=1))
    }else if(usehaz=="lognormal"){
      T0 <- rlnorm(n, sdlog=shape, meanlog=0)
      T1 <- pmin(rlnorm(n, sdlog=shape, meanlog=0), rexp(n,delta))
      C0 <- rexp(n,1/qlnorm(cenprob, sdlog=shape, meanlog=0))
      C1 <- rexp(n,1/qlnorm(cenprob, sdlog=shape, meanlog=0))
    }else if(usehaz=="sine"){
      T0 <- rsin(n)
      T1 <- pmin(rsin(n), rexp(n, delta))
      meanCen <- optimize(function(t) (cenprob-Fsin(t))^2, lower=0, 5)$minimum
      C0 <- rexp(n, 1/meanCen)
      C1 <- rexp(n, 1/meanCen)
    }else if(usehaz=="loglogis"){
      T0 <- rloglogis(n,1,shape)
      T1 <- pmin(rloglogis(n,1,shape), rexp(n,delta))
      C0 <- rexp(n, 1/qloglogis(cenprob, 1, shape))
      C1 <- rexp(n, 1/qloglogis(cenprob, 1, shape))
    }
    ev0 <- as.numeric(T0 < C0)
    ev1 <- as.numeric(T1 < C1)
    #ev0 <- ev1 <- rep(1, n)
    sortlist <- prepdata(pmin(T0,C0), pmin(T1,C1), ev0, ev1)
    groups <- rep(1:K, times=(2*n/K))
    sortlist$groups <- groups
    betas <- rep(0, B)
    
    ahfit <-  estBetaLY(sortlist)
    ahbeta <- ahfit$beta
    
    b <- 2*n
    for(j in 1:B){
      subsmp <- sort(sample(1:(2*n), b, replace=TRUE))
      
      sl1 <- sortlist[subsmp,]
      dupes <- duplicated(sl1$X)
      sl1$X[dupes] <- sl1$X[dupes]+sort(runif(sum(dupes), 0, 1e-6))
      
      groups <- rep((1:K), ceiling(b/K))[1:b]
      
      betas[j] <- estBetaSmoothU(sl1, groups, rep(1,b), MR)
      
      
      if(is.nan(betas[j])) failSL <- sl1
      
    }
    betahat <- estBetaSmoothU(sortlist, sortlist$groups, rep(1, (2*n)), MR)
    
    betahats[i] <- betahat
    
    #crit.val <- quantile(betas-betahat, probs=c(1-alpha/2, alpha/2))#*(sqrt(2*n)/(sqrt(2*n) + sqrt(SN)))
    #CIs[i,] <- betahat - crit.val
    CIs[i,] <- quantile(betas, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE)
    
    #ahse <- ahaz:::summary.ahaz(ahfit)$coefficients[2]
    #CIsAH[i,1] <- ahbeta-1.96*ahse
    #CIsAH[i,2] <- ahbeta+1.96*ahse
    LYbeta[i] <- ahbeta
    vhat[i] <- ahfit$B/ahfit$A^2
  }
  
  return(cbind( betahats, CIs, LYbeta, vhat))
}

library(parallel)

#shape <-c(1/2, 2/3, 3/4, 1, 1.5, 2, 3)
#shape <-c(1/2,3)
cenprob <- c(0.2,  0.4,  0.6, 0.8)
#cenprob <- c(0.4, 0.6)
shape <- c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)

#cenprob <- c(0.4, 0.6, 0.8)
param <- cbind(rep(shape, length(cenprob)) ,rep(cenprob, each=length(shape)), rep("loglogis", length(shape)*length(cenprob)))
#param <- rbind(param, cbind(rep(shape, length(cenprob)) ,rep(cenprob, each=length(shape)), rep("lognormal", length(shape)*length(cenprob))))
#param <- rbind(param, cbind(rep(1, length(cenprob)), cenprob, rep("sine", length(cenprob))))
param <- rbind(param, param)
paramlist <- vector("list", length(cenprob)*length(shape))
for(i in 1:dim(param)[1]){
  paramlist[[i]] <- param[i,]
}


n <-  100
#nrep <- 2000
#B <- 2000
nrep <- 2000
delta <- 0.1
B <- 1000
MR <- 10
b <- floor((2*n)^.75)

nclus <- detectCores()
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 567)
clusterEvalQ(cl, library(MASS))
test1 <- parLapply(cl=cl, paramlist[1:12], runSims, nrep,   delta, n, B, MR=MR, b)

stopCluster(cl)

names(test1) <- as.character(paramlist)[1:12]
# 
save(test1, file="n100LLDel01part1.RData")


nclus <- detectCores()
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 678)
clusterEvalQ(cl, library(MASS))
test2 <- parLapply(cl=cl, paramlist[13:24], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test2) <- as.character(paramlist)[13:24]
# 
save(test2, file="n100LLDel01part2.RData")

nclus <- 8
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 789)
clusterEvalQ(cl, library(MASS))
test3 <- parLapply(cl=cl, paramlist[25:32], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test3) <- as.character(paramlist)[25:32]
# 
save(test3, file="n100LLDel01part3.RData")

delta <- 0.3

nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 891)
clusterEvalQ(cl, library(MASS))
test4 <- parLapply(cl=cl, paramlist[33:43], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test4) <- as.character(paramlist)[33:43]
# 
save(test4, file="n100LLDel03part4.RData")


nclus <- detectCores()-1
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 912)
clusterEvalQ(cl, library(MASS))
test5 <- parLapply(cl=cl, paramlist[44:54], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test5) <- as.character(paramlist)[44:54]
# 
save(test5, file="n100LLDel03part5.RData")

nclus <- 10
cl <- makeCluster(nclus)
clusterSetRNGStream(cl, iseed = 123)
clusterEvalQ(cl, library(MASS))
test6 <- parLapply(cl=cl, paramlist[55:64], runSims, nrep,   delta, n, B, MR=MR, b)
stopCluster(cl)

names(test6) <- as.character(paramlist)[55:64]
# 
save(test6, file="n100LLDel03part6.RData")