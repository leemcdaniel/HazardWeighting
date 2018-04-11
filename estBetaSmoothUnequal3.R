estBetaSmoothU <- function(sortlist, groups, usevec, maxRat){
  K <- length(unique(groups))
  N <- dim(sortlist)[1]
  useable <- sum(usevec)
  gsize <- ceiling(useable - useable/K)
  
  tm <- fillMat(sortlist$X, groups, gsize, K)
  zm <- fillMat(sortlist$Z, groups, gsize, K)
  em <- fillMat(sortlist$ev, groups, gsize, K)
  usemat <- fillMat(usevec, groups, gsize, K)
  usematBig <- matrix(c(rep(1, K*useable), rep(0, K*(N-useable))), ncol=K, byrow=TRUE)
  
  # Calculate derived values (risk set, covariates at risk,
  # time difference, average covariates)
  #rmat <- matrix(rep((N - N/K):1 , times=K), ncol=K, byrow=FALSE)
  totmat <- matrix(rep(apply(usemat, 2, sum), dim(usemat)[1]+1), ncol=K, byrow=TRUE)
  rmat <- (totmat - rbind(rep(0, K), colCumsums(usemat)))[-dim(totmat)[1],]
  cmat <- sumCovsK(sortlist$Z, colSums(zm), groups, gsize, K)
  
  
  
  
  
  tm1 <- (rbind(rep(0, K), tm[1:(dim(tm)[1]-1), ]))
  tdiff <- tm-tm1
  zbar <- cmat/rmat
  zbar[is.nan(zbar)] <- 0
  
  # Calculate initial estimates from LY model
  num <- colSums(usemat*(zm-zbar)*em)
  denom <- colSums(usemat*(cmat*(1-zbar)^2 + usemat*(rmat-cmat)*(zbar^2))*tdiff)
  results <- num/denom
  resmat <- matrix(rep(results, N), ncol=K, byrow=TRUE)
  resmatS <- resmat[1:gsize,]
  
  hazraw <- colCumsums(usemat*em/rmat) - usemat*resmatS*tm*cmat/rmat
  
  regList <- fixHazards(hazraw, em, K, FALSE)
  
  # Get Estimated Scale/Shape for Weibull
  
  timesmat <- haztimeMat(tm, K, groups)
  pieceExp <- timesmat
  bw.pilot <- (apply(tm, 2, max) )/(8*apply(em, 2, sum)^0.2)
  for(i in 1:K){
    pieceExp[,i] <- smoothHaz(timesmat[,i], tm[usemat[,i]*regList$useEv[,i]==1,i], diff(c(0, regList$hazfix[usemat[,i]*regList$useEv[,i]==1,i])), bw.pilot[i])
  }
  
  
  cvec <- c(cmat[1,1]+sortlist$Z[1],cmat[,1])
  cvec <- c(sum(sortlist$Z), sum(sortlist$Z)-cumsum(sortlist$Z)[1:(N-1)])
  rvec <- N:1
  
  hazGroup1 <- pieceExp
  hazGroup0 <- pieceExp
  for(i in 1:K){
    if(any((hazGroup1[,i] + results[i]) < 0)){
      hazGroup1[,i] <- hazGroup1[,i] - min(hazGroup1[,i])
    }else{
      hazGroup1[,i] <- hazGroup1[,i] + results[i]
    }
  }
  
  hazGroup0[hazGroup0 < 0] <- 0
  hazGroup1[hazGroup1 < 0] <- 0
  mashup <- rbind(hazGroup0, hazGroup1)
  minvec <- apply(hazGroup0, 2, min)
  maxvec <- apply(hazGroup0, 2, max)
  meanvec <- apply(hazGroup0, 2, mean)
  
  if(any(meanvec == 0)){
    ind <- which(meanvec == 0)
    hazGroup0[,ind] <- 1
    meanvec[ind] <- 1
  }
  
  
  if(any(meanvec/minvec > maxRat)){
    for(i in 1:K){
      hazGroup0[which(meanvec[i]/hazGroup0[,i] > maxRat), i] <- mean(hazGroup0[,i])/maxRat
    }
  }
  
  if(any(maxvec/meanvec > maxRat)){
    for(i in 1:K){
      hazGroup0[which(hazGroup0[,i]/meanvec[i] > maxRat), i] <- mean(hazGroup0[,i])*maxRat 
    }
  }
  
  hazGroup1 <- hazGroup0 + resmat
  for(i in 1:K){
    hazGroup1[,i] <- pmax(hazGroup1[,i], meanvec[i]/maxRat)
  }
  
  zt <- (cvec/hazGroup1)/(cvec/hazGroup1 + (rvec-cvec)/hazGroup0)
  
  numer <- addMatVals(sortlist$ev*(sortlist$Z*(1-zt)/(hazGroup1) + (1-sortlist$Z)*(-zt)/hazGroup0), groups, useable)
  denom <- addMatVals((((sortlist$X - c(0,sortlist$X[1:(N-1)]))*((1/(hazGroup1))*cvec*(1 - zt)^2 + (1/(hazGroup0))*(rvec-cvec)*zt^2))), groups, useable)
  
  return(numer/denom)  
}




estBetaLY <- function(sortlist){
  #K <- length(unique(groups))
  N <- dim(sortlist)[1]
  #tm <- fillMat(sortlist$X, rep(2,N), 1)
  #zm <- fillMat(sortlist$Z, groups, K)
  #em <- fillMat(sortlist$ev, groups, K)
  
  # Calculate derived values (risk set, covariates at risk,
  # time difference, average covariates)
  riskvec <- N:1
  covrisk <- c(sum(sortlist$Z), sum(sortlist$Z)-cumsum(sortlist$Z)[1:(N-1)])
  
  tdiff <- diff(c(0, sortlist$X))
  
  zbar <- covrisk/riskvec
  
  # Calculate initial estimates from LY model
  num <- sum((sortlist$Z-zbar)*sortlist$ev)
  denom <- sum((covrisk*(1-zbar)^2 + (riskvec-covrisk)*(zbar^2))*tdiff)
  results <- num/denom
  #resmat <- matrix(rep(results, N), ncol=K, byrow=TRUE)
  #resmatS <- resmat[1:(N-(N/K)),]
  
  hazraw <- cumsum(sortlist$ev/riskvec) - results*sortlist$X*covrisk/riskvec
  
  B <- sum(((sortlist$Z-zbar)^2)*sortlist$ev)
  
  return(list(beta=results, hazraw=hazraw, B=B, A=denom))  
}



vb <- function(t, b, kern, haz, cumhaz, chtimes, t.grid){
  if(any((t-b*t.grid) >= 0)){
    chaz <- cumhaz[findInterval(t - b*t.grid, c(-Inf, chtimes))]
  }else{
    chaz <- 1
  }
  integrand <- (kern(t.grid,b)^2)*haz(t - b*t.grid, b)/exp(-chaz)
  (1/(length(cumhaz-1)*b))*sum(diff(c(0,t.grid[length(t.grid)-1]))*(integrand[1:(length(integrand)-1)] + integrand[2:length(integrand)]))/2
}

biasnum <- function(t, b, kern, haz, t.grid){
  integrand <- kern(t.grid,b)*haz(t - b*t.grid, b)
  sum(diff(c(0,t.grid[length(t.grid)-1]))*(integrand[1:(length(integrand)-1)] + integrand[2:length(integrand)]))/2 - haz(t, b)
}


epkern <- function(t, b){
  (3/4) * (1 - (t/b)^2)*(abs(t/b) <= 1)
}

# n.bw <- 51
# n.time <- 51
# min.mat <- matrix(0, nrow=n.time, ncol=n.bw)
# for(i in 1:K){
#   bw.grid <- seq(from=0.2*bw.pilot[i], to=5*bw.pilot[i], length.out=n.bw)
#   time.grid <- seq(from=min(tm[,i]), to=max(tm[,i]), length.out=n.time)
#   hf <- function(t, b) smoothHaz(t, tm[,i], diff(c(0, regList$hazfix[regList$useEv[,i]==1,i])), b)
#   for(j in 1:n.bw){
#     for(k in 1:n.time){
#       min.mat[j,k] <- biasnum(time.grid[k], bw.grid[j], epkern, hf, time.grid)^2 + 
#         vb(time.grid[k], bw.grid[j], epkern, hf, diff(c(0, regList$hazfix[regList$useEv[,i]==1,i])), tm[,i], time.grid)
#     }
#   }
#   
# }

rsin <- function(n){
  u <- rexp(n)
  return(floor(u/4)*pi + floor(0.5 + u/4)*pi + acos(-(u-2*floor(u/2) - 1)))
}

Hsin <- function(t){
  2*floor(t/pi) + 1 + cos(t)*(-1)^(floor(t/pi)+1)
}
Fsin <- function(t){
  1-exp(-Hsin(t))
}