path <- "~/../Dropbox/EfficientLinYing/Code/"
#library(Rcpp)
library(matrixStats)
#library(ahaz)
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

########## ALL CAUSE ################
library(ISwR)
dat <- nickel

set.seed(1)
dat$event <- as.numeric(dat$icd > 0)



ahfit <- ahaz(Surv((dat$ageout+runif(dim(dat)[1], 0 ,1e-6)-dat$agein), dat$event), as.numeric(dat$exposure > 0))

T0 <- (dat$ageout-dat$agein)[dat$exposure == 0]
T1 <- (dat$ageout-dat$agein)[dat$exposure > 0]
ev0 <- dat$event[dat$exposure==0]
ev1 <- dat$event[dat$exposure > 0]
sortlist <- (prepdata(T0,T1, ev0, ev1))
B <- 1000
K <- 10
N <- dim(sortlist)[1]
MR <- 10
groups <- rep(1:K, times=ceiling(N/K))[1:N]
sortlist$groups <- groups
betas <- rep(0, B)
b <- N

t1 <- proc.time()
for(j in 1:B){
  
  subsmp <- sort(sample(1:N, N, replace=TRUE))
  
  sl1 <- sortlist[subsmp,]
  dupes <- duplicated(sl1$X)
  sl1$X[dupes] <- sl1$X[dupes]+sort(runif(sum(dupes), 0, 1e-6))  
  groups <- rep((1:K), ceiling(b/K))[1:b]
  
  
  
  
  
  betas[j] <- estBetaSmoothU(sl1, groups, rep(1,b), MR)
  
  
  if(is.nan(betas[j])) failSL <- sl1
  
  
}
t2 <- proc.time()
t2-t1
betahat <- estBetaSmoothU(sortlist, sortlist$groups, rep(1,dim(sortlist)[1]), MR)


alpha <- 0.05
#CI <- betahat - quantile(sqrt(b)*(betas-betahat), probs=c(1-alpha/2, alpha/2), na.rm=TRUE)/sqrt(N)
CI <- quantile(betas, probs=c(alpha/2, 1-alpha/2))

CIah <- c(summary(ahfit)$coefficients[1] - qnorm(0.975)*summary(ahfit)$coefficients[2], summary(ahfit)$coefficients[1] + qnorm(0.975)*summary(ahfit)$coefficients[2] )
diff(CI)/diff(CIah)

allcause <- sortlist

############ NASAL CANCER ##############

set.seed(2)

dat <- nickel

dat$event <- as.numeric(dat$icd == 160)



ahfit <- ahaz(Surv((dat$ageout+runif(dim(dat)[1], 0 ,1e-6)-dat$agein), dat$event), as.numeric(dat$exposure > 0))

T0 <- (dat$ageout-dat$agein)[dat$exposure == 0]
T1 <- (dat$ageout-dat$agein)[dat$exposure > 0]
ev0 <- dat$event[dat$exposure==0]
ev1 <- dat$event[dat$exposure > 0]
sortlist <- (prepdata(T0,T1, ev0, ev1))
B <- 1000
K <- 10
N <- dim(sortlist)[1]
MR <- 10
groups <- rep(1:K, times=ceiling(N/K))[1:N]
sortlist$groups <- groups
betas <- rep(0, B)

t1 <- proc.time()
for(j in 1:B){
  
  subsmp <- sort(sample(1:N, N, replace=TRUE))
  
  sl1 <- sortlist[subsmp,]
  dupes <- duplicated(sl1$X)
  sl1$X[dupes] <- sl1$X[dupes]+sort(runif(sum(dupes), 0, 1e-6))  
  groups <- rep((1:K), ceiling(b/K))[1:b]
  
  
  
  
  
  betas[j] <- estBetaSmoothU(sl1, groups, rep(1,b), MR)
  
  
  if(is.nan(betas[j])) failSL <- sl1
  
  
}
t2 <- proc.time()
t2-t1
betahat <- estBetaSmoothU(sortlist, sortlist$groups, rep(1,dim(sortlist)[1]), MR)


alpha <- 0.05
#CI <- betahat - quantile(sqrt(b)*(betas-betahat), probs=c(1-alpha/2, alpha/2), na.rm=TRUE)/sqrt(N)
CI <- quantile(betas, probs=c(alpha/2, 1-alpha/2))

CIah <- c(summary(ahfit)$coefficients[1] - qnorm(0.975)*summary(ahfit)$coefficients[2], summary(ahfit)$coefficients[1] + qnorm(0.975)*summary(ahfit)$coefficients[2] )
diff(CI)/diff(CIah)

nasal <- sortlist

############ LUNG CANCER ##############

set.seed(3)

dat <- nickel

dat$event <- as.numeric(dat$icd == 162 | dat$icd == 163)



ahfit <- ahaz(Surv((dat$ageout+runif(dim(dat)[1], 0 ,1e-6)-dat$agein), dat$event), as.numeric(dat$exposure > 0))

T0 <- (dat$ageout-dat$agein)[dat$exposure == 0]
T1 <- (dat$ageout-dat$agein)[dat$exposure > 0]
ev0 <- dat$event[dat$exposure==0]
ev1 <- dat$event[dat$exposure > 0]
sortlist <- (prepdata(T0,T1, ev0, ev1))
B <- 1000
K <- 10
N <- dim(sortlist)[1]
MR <- 10
groups <- rep(1:K, times=ceiling(N/K))[1:N]
sortlist$groups <- groups
betas <- rep(0, B)

t1 <- proc.time()
for(j in 1:B){
  
  subsmp <- sort(sample(1:N, N, replace=TRUE))
  
  sl1 <- sortlist[subsmp,]
  dupes <- duplicated(sl1$X)
  sl1$X[dupes] <- sl1$X[dupes]+sort(runif(sum(dupes), 0, 1e-6))  
  groups <- rep((1:K), ceiling(b/K))[1:b]
  
  
  
  
  
  betas[j] <- estBetaSmoothU(sl1, groups, rep(1,b), MR)
  
  
  if(is.nan(betas[j])) failSL <- sl1
  
  
}
t2 <- proc.time()
t2-t1
betahat <- estBetaSmoothU(sortlist, sortlist$groups, rep(1,dim(sortlist)[1]), MR)


alpha <- 0.05
#CI <- betahat - quantile(sqrt(b)*(betas-betahat), probs=c(1-alpha/2, alpha/2), na.rm=TRUE)/sqrt(N)
CI <- quantile(betas, probs=c(alpha/2, 1-alpha/2))

CIah <- c(summary(ahfit)$coefficients[1] - qnorm(0.975)*summary(ahfit)$coefficients[2], summary(ahfit)$coefficients[1] + qnorm(0.975)*summary(ahfit)$coefficients[2] )
diff(CI)/diff(CIah)

lung <- sortlist

###################################
###################################
###################################
###################################
library(muhaz)

ac1 <- muhaz(allcause$X[sortlist$Z==1], allcause$ev[sortlist$Z==1])
ac0 <- muhaz(allcause$X[sortlist$Z==0], allcause$ev[sortlist$Z==0])

n1 <- muhaz(nasal$X[sortlist$Z==1], nasal$ev[sortlist$Z==1])
n0 <- muhaz(nasal$X[sortlist$Z==0], nasal$ev[sortlist$Z==0])

l1 <- muhaz(lung$X[sortlist$Z==1], lung$ev[sortlist$Z==1])
l0 <- muhaz(lung$X[sortlist$Z==0], lung$ev[sortlist$Z==0])

library(ggplot2)

npl <- data.frame(Time=c(ac0$est.grid, ac1$est.grid, n0$est.grid, n1$est.grid, l0$est.grid, l1$est.grid),
                  Hazard = c(ac0$haz.est, ac1$haz.est, n0$haz.est, n1$haz.est, l0$haz.est, l1$haz.est),
                  Exposed=rep(c(rep("No", 101), rep("Yes", 101)), times=3),
                  Cause = rep(c("All Cause", "Nasal Cancer", "Lung Cancer"), each=202))

pl <- ggplot(npl) + geom_line(aes(Time, Hazard, color=Exposed), size=2) + xlab("Time (years)") + ylab("Hazard (events/year)")

#pl <- pl + geom_text(aes(x=25, y=0.035, label="Not Exposed")) + geom_text(aes(x=15, y=0.065, label="Exposed")) + theme(legend.position="none")
pl <- pl + facet_wrap("Cause", ncol=1, scales="free_y")
pl


ExportPlot <- function(gplot, filename, width=2, height=1.5) {
  # Export plot in PDF and EPS.
  # Notice that A4: width=11.69, height=8.27
  ggsave(paste(filename, '.pdf', sep=""), gplot, width = width, height = height)
  postscript(file = paste(filename, '.eps', sep=""), width = width, height = height)
  print(gplot)
  dev.off()
  png(file = paste(filename, '_.png', sep=""), width = width * 100, height = height * 100)
  print(gplot)
  dev.off()
}


ExportPlot(pl, "Nickel", width=6, height=8)

ggsave("nickel.png", path="~/../Dropbox/EfficientLinYing/Plots")
