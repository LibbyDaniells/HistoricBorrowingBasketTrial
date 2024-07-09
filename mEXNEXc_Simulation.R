library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)


#mEXNEXc with Hellinger Weights--------------------------------------------------
diff <- function(vec,n){  #Compute pairwise differences in response rates
  len <- 1:length(vec)
  diff_mat <- matrix(,nrow=length(vec),ncol=length(vec))
  for(i in 1:length(vec)){
    for(j in len[-i]){
      diff_mat[i,j] <- abs(vec[i]-vec[j])
    }
  }
  return(diff_mat)
}

pi_fun_H <- function(vec,n,c){
  data <- data.frame(Responses=vec,No.Patients=n)
  y <- data$Responses
  x <- seq(0,1,by=0.0001)
  p <- matrix(,nrow=length(vec),ncol=length(x)) #Get posterior densities
  for(i in 1:length(vec)){
    p[i,] <- dbeta(x,y[i]+1,n[i]-y[i]+1)
  }
  x <- p
  est <- vec/n
  diff_mat <- diff(est,n)
  rem <- c()
  keep <- c()
  pi <- rep('NA',length(vec))
  for(i in 1:length(vec)){
    min <- as.numeric(min(diff_mat[i,],na.rm=TRUE))
    if(min>(c+0.00001)){ #If min pairwise difference is greater than c_e then treat as independent
      rem[i] <- i
      pi[i] <- 0
    }else{
      keep[i] <- i
    }
  }
  keep <- keep[!is.na(keep)]
  rem <- rem[!is.na(rem)]
  if(length(rem)==length(vec)){ #If all 'extreme' treat all independent
    pi <- rep(0,length(vec))
  }else if(length(rem)==0){ #If all not 'extreme' compute all pairwise Hellinger distances and take average
    mat <- 1-CalcHellingerDist(x)
    for(i in 1:length(vec)){
      pi[i] <- mean(mat[i,-i])
    }
  }else{ #Compute pairwise Hellinger distances between those not deemed 'extreme' and treat the rest independent
    y <- x[-rem,]
    mat2 <- 1-CalcHellingerDist(y)
    for(i in 1:length(keep)){
      pi[keep[i]] <- mean(mat2[i,-i])
    }
  }
  pi <- as.numeric(pi) #Probability vector for the EXNEX model
  return(pi)
}



mEXNEXhell <- function(K,y,yh,n,nh,q0,pw,H,hist.scale,pi){
  prob <- pi
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  yh <- yh[-which(nh==0)]
  nh <- nh[-which(nh==0)]
  prob <- pi_fun_H(yh,nh,1)
  if(length(prob)==K){
    pi <- prob
  }else{
    pi <- c(prob,rep(mean(prob)*hist.scale,K-H))
  }
  pi <- cbind(pi,1-pi)
  exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
  jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=100000,n.chains=4)
  samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p','delta'),n.iter=100000,silent=T)
  mexnexhell <- as.data.frame(samples.exnex2[[1]])
  return(mexnexhell)
}


#Calibration---------------------------------------
mEXNEXc_Cal <- function(K,p,n,ph,nh,q0,pw,pi,hist.scale,run){
  Delta.mat <- matrix(NA,nrow=run,ncol=K)
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  H <- sum(ph!=0)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    yh <- rbinom(K,nh,ph)
    y.hist <- yh[-which(nh==0)]
    n.hist <- nh[-which(nh==0)]
    prob <- pi_fun_H(y.hist,n.hist,1)
    if(length(prob)==K){
      pi <- prob
    }else{
      pi <- c(prob,rep(mean(prob)*hist.scale,K-H))
    }
    pi <- cbind(pi,1-pi)
    exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
    jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=1000,n.chains=4,quiet=T)
    samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=100000,silent=T)
    mexnexhell <- as.data.frame(samples.exnex2[[1]])
    Delta.mat[i,] <- apply(mexnexhell,2,post.prob.fun)
    print(i)
  }
  Delta <- colQuantiles(Delta.mat,probs=0.9)
  return(Delta)
}


# K <- 5
# p <- rep(0.1,K)
# n <- rep(34,K)
# ph <- c(0.1,0.1,0.1,0,0)
# nh <- c(13,13,13,0,0)
# q0 <- 0.1
# run <- 10000
# pw <- 0.2
# hist.scale <- 0.8
# pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# 
# mexnexc <- mEXNEXc_Cal(K,p,n,ph,nh,q0,pw,pi,hist.scale,run)
# 
# save(mexnexc,file='mEXNEXc_Delta.RData')


# 
# #Simulation-----------------------------
# mEXNEXc_Sim <- function(pmat,n,nh,K,q0,pw,pi,run,mexnexc.delta,hist.scale){
#   H <- sum(nh!=0)
#   p <- pmat[1:K]
#   ph <- pmat[K+1:H]
#   hypo <- matrix(NA,nrow=run,ncol=K)
#   true <- rep(0,K)
#   pointests <- matrix(NA,nrow=run,ncol=K)
#   nexmu <- rep(log(pw/(1-pw)),K)
#   nexsigma <- rep((1/pw)+(1/(1-pw)),K)
#   for(t in 1:K){
#     if(p[t]<=q0){true[t]<-0}else{true[t]<-1}
#   }
#   post.prob.fun <- function(x){
#     fun<-sum(x>q0)/100000
#     return(fun)
#   }
#   for(i in 1:run){
#     y <- rbinom(K,n,p)
#     yh <- rbinom(K,nh,ph)
#     y.hist <- yh[-which(nh==0)]
#     n.hist <- nh[-which(nh==0)]
#     prob <- pi_fun_H(y.hist,n.hist,1)
#     if(length(prob)==K){
#       pi <- prob
#     }else{
#       pi <- c(prob,rep(mean(prob)*hist.scale,K-H))
#     }
#     pi <- cbind(pi,1-pi)
#     exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
#     jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=1000,n.chains=4,quiet=T)
#     samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=100000,silent=T)
#     mexnexhell <- as.data.frame(samples.exnex2[[1]])
#     pointests[i,] <- colMeans(mexnexhell)
#     hypo[i,] <- as.integer(apply(mexnexhell,2,post.prob.fun)>mexnexc.delta)
#     print(i)
#   }
#   perfect <- 0 
#   for(j in 1:run){
#     if(all(hypo[j,]==true)){perfect<-perfect+1}
#   }
#   fwer_true <- which(true==0)
#   if(sum(true)==K){fwer <- rep('NA',K)}else{
#     fwer <- 0
#     for(a in 1:run){
#       error <- rep(0,length(fwer_true))
#       for(b in 1:length(fwer_true)){
#         if(hypo[a,fwer_true[b]]==1){
#           error[b] <- 1
#         }else{error[b] <- 0}
#       }
#       if(sum(error)!=0){fwer <- fwer+1}
#     }
#     fwer <- fwer/run
#   }
#   Pointmeans <- colMeans(pointests)
#   Pointsds <- apply(pointests,2,sd)
#   my_list <- list('Hypo'=hypo,'Point Estimates'=pointests,'Point Estimate Means'=Pointmeans,'Point Estimate Sds'=Pointsds,'Error Rates'=colMeans(hypo)*100,'Perfect'=perfect/run,'FWER'=fwer)
#   return(my_list)
# }


# K <- 5
# n <- rep(34,K)
# nh <- c(13,13,13,0,0)
# q0 <- 0.1
# mexnexc.delta <- 0.8340992
# run <- 10000
# pw <- 0.2
# hist.scale <- 0.8
# pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# 
# p1 <- c(c(0.1,0.1,0.1,0.1,0.1),c(0.1,0.1,0.1))
# p2 <- c(c(0.1,0.1,0.1,0.1,0.1),c(0.25,0.1,0.1))
# p3 <- c(c(0.1,0.1,0.1,0.1,0.1),c(0.25,0.25,0.25))
# p4 <- c(c(0.25,0.1,0.1,0.1,0.1),c(0.25,0.1,0.1))
# p5 <- c(c(0.25,0.1,0.1,0.1,0.1),c(0.25,0.25,0.1))
# p6 <- c(c(0.25,0.1,0.1,0.1,0.1),c(0.25,0.25,0.25))
# p7 <- c(c(0.25,0.25,0.1,0.1,0.1),c(0.25,0.25,0.1))
# p8 <- c(c(0.25,0.25,0.1,0.1,0.1),c(0.25,0.1,0.1))
# p9 <- c(c(0.25,0.25,0.1,0.1,0.1),c(0.25,0.25,0.25))
# p10 <- c(c(0.25,0.25,0.25,0.1,0.1),c(0.25,0.25,0.25))
# p11 <- c(c(0.25,0.25,0.25,0.1,0.1),c(0.25,0.25,0.1))
# p12 <- c(c(0.25,0.25,0.25,0.1,0.1),c(0.1,0.1,0.1))
# p13 <- c(c(0.25,0.25,0.25,0.25,0.1),c(0.25,0.25,0.25))
# p14 <- c(c(0.25,0.25,0.25,0.25,0.1),c(0.25,0.25,0.1))
# p15 <- c(c(0.25,0.25,0.25,0.25,0.1),c(0.1,0.1,0.1))
# p16 <- c(c(0.25,0.25,0.25,0.25,0.25),c(0.25,0.25,0.25))
# p17 <- c(c(0.25,0.25,0.25,0.25,0.25),c(0.25,0.25,0.1))
# p18 <- c(c(0.25,0.25,0.25,0.25,0.25),c(0.1,0.1,0.1))
# 
# list.scenarios1 <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
# list.scenarios2 <- list(p10,p11,p12,p13,p14,p15,p16,p17,p18)
# 
# mexnexc.sim1_9 <- mclapply(list.scenarios1,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,mexnexc.delta,hist.scale,mc.cores=10)
# save(mexnexc.sim1_9,file='mEXNEXc_Sim1_9.RData')
# 
# mexnexc.sim10_18 <- mclapply(list.scenarios2,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,mexnexc.delta,hist.scale,mc.cores=10)
# save(mexnexc.sim10_18,file='mEXNEXc_Sim10_18.RData')




#Simulation Across Scenarios-----------------------------
mEXNEXc_Sim <- function(pmat,n,nh,K,q0,pw,pi,run,hist.scale){
  H <- sum(nh!=0)
  p <- pmat[1:K]
  yh <- pmat[K+1:H]
  post.prob <- matrix(NA,nrow=run,ncol=K)
  true <- rep(0,K)
  pointests <- matrix(NA,nrow=run,ncol=K)
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  for(t in 1:K){
    if(p[t]<=q0){true[t]<-0}else{true[t]<-1}
  }
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    n.hist <- nh[-which(nh==0)]
    prob <- pi_fun_H(yh,n.hist,1)
    if(length(prob)==K){
      pi <- prob
    }else{
      pi <- c(prob,rep(mean(prob)*hist.scale,K-H))
    }
    pi <- cbind(pi,1-pi)
    exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
    jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=1000,n.chains=4,quiet=T)
    samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=100000,silent=T)
    mexnexhell <- as.data.frame(samples.exnex2[[1]])
    pointests[i,] <- colMeans(mexnexhell)
    post.prob[i,] <- apply(mexnexhell,2,post.prob.fun)
    print(i)
  }
  Pointmeans <- colMeans(pointests)
  Pointsds <- apply(pointests,2,sd)
  my_list <- list('Point Estimates'=pointests,'Point Estimate Means'=Pointmeans,'Point Estimate Sds'=Pointsds,'Posterior Probabilities'=post.prob)
  return(my_list)
}




K <- 5
n <- rep(34,K)
nh <- c(13,13,13,0,0)
q0 <- 0.1
run <- 5000
pw <- 0.2
hist.scale <- 0.8
pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))

p1 <- c(0.1,0.1,0.1,0.1,0.1,1,1,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,1,1,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,1,1,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,1,1,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,1,1,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,1,1,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,1,1,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,1,1,1)

list.scenarios1 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
mexnexc.sim0 <- mclapply(list.scenarios1,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=9)
save(mexnexc.sim0,file='mEXNEXc_Sim0.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,1,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,1,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,1,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,1,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,1,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,1,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,1,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,1,1)

list.scenarios2 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
mexnexc.sim1 <- mclapply(list.scenarios2,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=9)
save(mexnexc.sim1,file='mEXNEXc_Sim1.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,1)

list.scenarios3 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
mexnexc.sim2 <- mclapply(list.scenarios3,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=9)
save(mexnexc.sim2,file='mEXNEXc_Sim2.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,3)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,3)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,3)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,3)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,3)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,3)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,3)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,3)

list.scenarios4 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
mexnexc.sim3 <- mclapply(list.scenarios4,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=9)
save(mexnexc.sim3,file='mEXNEXc_Sim3.RData')



#Extra Scenarios
p9 <- c(0.25,0.25,0.1,0.25,0.1,1,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,1,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,1,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,1,1,1)

extra.list.scenarios1 <- list(p9,p10,p11,p12)
extra.mexnexc.sim0 <- mclapply(extra.list.scenarios1,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=5)
save(extra.mexnexc.sim0,file='mEXNEXc_Sim0_Extra.RData')

p9 <- c(0.25,0.25,0.1,0.25,0.1,3,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,1,1)

extra.list.scenarios2 <- list(p9,p10,p11,p12)
extra.mexnexc.sim1 <- mclapply(extra.list.scenarios2,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=5)
save(extra.mexnexc.sim1,file='mEXNEXc_Sim1_Extra.RData')


p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,1)

extra.list.scenarios3 <- list(p9,p10,p11,p12)
extra.mexnexc.sim2 <- mclapply(extra.list.scenarios3,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=5)
save(extra.mexnexc.sim2,file='mEXNEXc_Sim2_Extra.RData')


p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,3)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,3)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,3)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,3)

extra.list.scenarios4 <- list(p9,p10,p11,p12)
extra.mexnexc.sim3 <- mclapply(extra.list.scenarios4,mEXNEXc_Sim,n,nh,K,q0,pw,pi,run,hist.scale,mc.cores=5)
save(extra.mexnexc.sim3,file='mEXNEXc_Sim3_Extra.RData')


#Results Analysis---------------------------------------------------
# calmat <- 0
# for(i in 1:5){
#   calmat <- c(calmat,mexnexc.sim1[[i]]$`Posterior Probabilities`[,5])
# }
# cut.off <- quantile(calmat,0.9)

#cut.off <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9))



#Sim 0------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mexnexc.sim0[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,mexnexc.sim0[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)

cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mexnexc.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios1[[j]][1:5]>q0)
  post.prob <- mexnexc.sim0[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios1[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- mexnexc.sim0[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=8,ncol=K)
point.sd <- matrix(NA,nrow=8,ncol=K)
for(i in 1:8){
  point.mean[i,] <- mexnexc.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mexnexc.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mexnexc.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mexnexc.results0,file='mEXNEXc_Results0.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mexnexc.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
  post.prob <- extra.mexnexc.sim0[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- extra.mexnexc.sim0[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=4,ncol=K)
point.sd <- matrix(NA,nrow=4,ncol=K)
for(i in 1:4){
  point.mean[i,] <- extra.mexnexc.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mexnexc.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mexnexc.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mexnexc.results0,file='mEXNEXc_Results0_Extra.RData')



#Sim 1------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mexnexc.sim1[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat23 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat23 <- c(calmat23,mexnexc.sim1[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.23 <- quantile(calmat23,0.9)

calmat1 <- 0
for(i in 1:8){
  if(i==1|i==7){
    calmat1 <- c(calmat1,mexnexc.sim1[[i]]$`Posterior Probabilities`[,1])}
}
cut.off.1 <- quantile(calmat1,0.9)

cut.off <- c(cut.off.1,rep(cut.off.23,2),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mexnexc.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios2[[j]][1:5]>q0)
  post.prob <- mexnexc.sim1[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios2[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- mexnexc.sim1[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=8,ncol=K)
point.sd <- matrix(NA,nrow=8,ncol=K)
for(i in 1:8){
  point.mean[i,] <- mexnexc.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mexnexc.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mexnexc.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mexnexc.results1,file='mEXNEXc_Results1.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mexnexc.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
  post.prob <- extra.mexnexc.sim1[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- extra.mexnexc.sim1[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=4,ncol=K)
point.sd <- matrix(NA,nrow=4,ncol=K)
for(i in 1:4){
  point.mean[i,] <- extra.mexnexc.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mexnexc.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mexnexc.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mexnexc.results1,file='mEXNEXc_Results1_Extra.RData')


#Sim 2------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mexnexc.sim2[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat3 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat3 <- c(calmat3,mexnexc.sim2[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.3 <- quantile(calmat3,0.9)

calmat12 <- 0
for(i in 1:8){
  if(i==1|i==7|i==2|i==8){
    calmat12 <- c(calmat12,mexnexc.sim2[[i]]$`Posterior Probabilities`[,2])}
}
cut.off.12 <- quantile(calmat12,0.9)


cut.off <- c(rep(cut.off.12,2),cut.off.3,rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mexnexc.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios3[[j]][1:5]>q0)
  post.prob <- mexnexc.sim2[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios3[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- mexnexc.sim2[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=8,ncol=K)
point.sd <- matrix(NA,nrow=8,ncol=K)
for(i in 1:8){
  point.mean[i,] <- mexnexc.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mexnexc.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mexnexc.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mexnexc.results2,file='mEXNEXc_Results2.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mexnexc.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
  post.prob <- extra.mexnexc.sim2[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- extra.mexnexc.sim2[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=4,ncol=K)
point.sd <- matrix(NA,nrow=4,ncol=K)
for(i in 1:4){
  point.mean[i,] <- extra.mexnexc.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mexnexc.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mexnexc.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mexnexc.results2,file='mEXNEXc_Results2_Extra.RData')


#Sim 3------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mexnexc.sim3[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,mexnexc.sim3[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)


cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mexnexc.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios4[[j]][1:5]>q0)
  post.prob <- mexnexc.sim3[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios4[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- mexnexc.sim3[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=8,ncol=K)
point.sd <- matrix(NA,nrow=8,ncol=K)
for(i in 1:8){
  point.mean[i,] <- mexnexc.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mexnexc.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mexnexc.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mexnexc.results3,file='mEXNEXc_Results3.RData')


#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mexnexc.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios4[[j]][1:5]>q0)
  post.prob <- extra.mexnexc.sim3[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)

fwer <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios4[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- extra.mexnexc.sim3[[j]]$`Posterior Probabilities`[,1:5]
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
    }
    if(sum(hypo[fwer.true])!=0){
      fwer[j] <- fwer[j]+1
    }
  }
}
FWER <- 100*(fwer/run)

point.mean <- matrix(NA,nrow=4,ncol=K)
point.sd <- matrix(NA,nrow=4,ncol=K)
for(i in 1:4){
  point.mean[i,] <- extra.mexnexc.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mexnexc.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mexnexc.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mexnexc.results3,file='mEXNEXc_Results3_Extra.RData')










