library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)

#EXNEX--------------------------------------------------------------------------
EXNEX <- function(K,y,n,q0,pw,pi){
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
  jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=100000,n.chains=4)
  samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p','delta'),n.iter=100000,silent=T)
  exnex <- as.data.frame(samples.exnex2[[1]])
  return(exnex)
}

#Calibration---------------------------------------
EXNEX_Cal <- function(K,p,n,q0,pw,pi,run){
  Delta.mat <- matrix(NA,nrow=run,ncol=K)
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
    jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=1000,n.chains=4,quiet=T)
    samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=100000,silent=T)
    exnex <- as.data.frame(samples.exnex2[[1]])
    Delta.mat[i,] <- apply(exnex,2,post.prob.fun)
    print(i)
  }
  Delta <- colQuantiles(Delta.mat,probs=0.9)
  return(Delta)
}


# K <- 5
# p <- rep(0.1,K)
# n <- rep(34,K)
# q0 <- 0.1
# run <- 10000
# pw <- 0.2
# pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# 
# exnex.current <- EXNEX_Cal(K,p,n,q0,pw,pi,run)
# 
# save(exnex.current,file='EXNEX_Current_Delta.RData')



#Simulation-----------------------------
# EXNEX_Sim <- function(p,n,K,q0,pw,pi,run,exnex.delta){
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
#     exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
#     jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=1000,n.chains=4,quiet=T)
#     samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=100000,silent=T)
#     exnex <- as.data.frame(samples.exnex2[[1]])
#     pointests[i,] <- colMeans(exnex)
#     hypo[i,] <- as.integer(apply(exnex,2,post.prob.fun)>exnex.delta)
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
# 
# 
# K <- 5
# n <- rep(34,K)
# q0 <- 0.1
# exnex.delta <- 0.8426812
# run <- 10000
# pw <- 0.2
# pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# 
# p1 <- c(0.1,0.1,0.1,0.1,0.1)
# p2 <- c(0.25,0.1,0.1,0.1,0.1)
# p3 <- c(0.25,0.25,0.1,0.1,0.1)
# p4 <- c(0.25,0.25,0.25,0.1,0.1)
# p5 <- c(0.25,0.25,0.25,0.25,0.1)
# p6 <- c(0.25,0.25,0.25,0.25,0.25)
# 
# list.scenarios <- list(p1,p2,p3,p4,p5,p6)
# exnex.sim <- mclapply(list.scenarios,EXNEX_Sim,n,K,q0,pw,pi,run,exnex.delta,mc.cores=7)
# 
# save(exnex.sim,file='EXNEX_Sim.RData')




#Calibrate Across Scenarios
EXNEX_Sim <- function(p,n,K,pw,q0,run,pi){
  post.prob <- matrix(NA,nrow=run,ncol=K)
  pointests <- matrix(NA,nrow=run,ncol=K)
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    exnex2.data <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
    jags.exnex2 <- jags.model(file='EXNEX.txt',data=exnex2.data,n.adapt=1000,n.chains=4,quiet=T)
    samples.exnex2 <- coda.samples(jags.exnex2,variable.names=c('p'),n.iter=100000,silent=T)
    exnex <- as.data.frame(samples.exnex2[[1]])
    pointests[i,] <- colMeans(exnex)
    post.prob[i,] <- apply(exnex,2,post.prob.fun)
    print(i)
  }
  Pointmeans <- colMeans(pointests)
  Pointsds <- apply(pointests,2,sd)
  my_list <- list('Point Estimates'=pointests,'Point Estimate Means'=Pointmeans,'Point Estimate Sds'=Pointsds,'Posterior Probabilities'=post.prob)
  return(my_list)
}





K <- 5
n <- rep(34,K)
q0 <- 0.1
run <- 5000
pw <- 0.2
pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))


p1 <- c(0.1,0.1,0.1,0.1,0.1)
p2 <- c(0.25,0.1,0.1,0.1,0.1)
p3 <- c(0.25,0.25,0.1,0.1,0.1)
p4 <- c(0.25,0.25,0.25,0.1,0.1)
p5 <- c(0.25,0.25,0.25,0.25,0.1)
p6 <- c(0.25,0.25,0.25,0.25,0.25)
p7 <- c(0.1,0.1,0.1,0.25,0.1)
p8 <- c(0.25,0.1,0.1,0.25,0.1)

list.scenarios <- list(p1,p2,p3,p4,p5,p6,p7,p8)
exnex.sim <- mclapply(list.scenarios,EXNEX_Sim,n,K,pw,q0,run,pi,mc.cores=9)
save(exnex.sim,file='EXNEX_Sim.RData')

p9 <- c(0.25,0.25,0.1,0.25,0.1)
p10 <- c(0.1,0.1,0.1,0.25,0.25)
p11 <- c(0.25,0.1,0.1,0.25,0.25)
p12 <- c(0.25,0.25,0.1,0.25,0.25)

extra.list.scenarios <- list(p9,p10,p11,p12)
extra.exnex.sim <- mclapply(extra.list.scenarios,EXNEX_Sim,n,K,pw,q0,run,pi,mc.cores=5)
save(extra.exnex.sim,file='EXNEX_Sim_Extra.RData')


#Results Analysis---------------------------------------------------


calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
calmat <- c(calmat,exnex.sim[[i]]$`Posterior Probabilities`[,5])}
}
cut.off <- quantile(calmat,0.9)


#cut.off <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9))
# cut.off.123 <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9)[1:3])
# cut.off.45 <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9)[4:5])

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
post.probs <- exnex.sim[[i]]$`Posterior Probabilities`
for(k in 1:5){
  ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off))*100
}
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios[[j]]>q0)
  post.prob <- exnex.sim[[j]]$`Posterior Probabilities`
for(i in 1:run){
  hypo <- rep(NA,K)
  for(k in 1:5){
    hypo[k] <- as.numeric(post.prob[i,k]>cut.off)
  }
  if(all(hypo==true)){
    perfect[j] <- perfect[j]+1
  }
}
}
Perfect <- 100*(perfect/run)


fwer <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios[[j]]>q0)
  fwer.true <- which(true==0)
  post.prob <- exnex.sim[[j]]$`Posterior Probabilities`
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off)
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
  point.mean[i,] <- exnex.sim[[i]]$`Point Estimate Means`
  point.sd[i,] <- exnex.sim[[i]]$`Point Estimate Sds`
}
point.mean
point.sd

exnex.results <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(exnex.results,file='EXNEX_Results.RData')


#Extra Scenarios
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.exnex.sim[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off))*100
  }
}


perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios[[j]]>q0)
  post.prob <- extra.exnex.sim[[j]]$`Posterior Probabilities`
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off)
    }
    if(all(hypo==true)){
      perfect[j] <- perfect[j]+1
    }
  }
}
Perfect <- 100*(perfect/run)


fwer <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios[[j]]>q0)
  fwer.true <- which(true==0)
  post.prob <- extra.exnex.sim[[j]]$`Posterior Probabilities`
  for(i in 1:run){
    hypo <- rep(NA,K)
    for(k in 1:5){
      hypo[k] <- as.numeric(post.prob[i,k]>cut.off)
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
  point.mean[i,] <- extra.exnex.sim[[i]]$`Point Estimate Means`
  point.sd[i,] <- extra.exnex.sim[[i]]$`Point Estimate Sds`
}
point.mean
point.sd

extra.exnex.results <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.exnex.results,file='EXNEX_Results_Extra.RData')


# #Results + Extra-----------------------------------
# calmat45 <- 0
# for(i in 1:8){
#   if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
#     calmat45 <- c(calmat45,exnex.sim[[i]]$`Posterior Probabilities`[,5])}
# }
# for(j in 1:4){
#   if(j==1){
#     calmat45 <- c(calmat45,extra.exnex.sim[[j]]$`Posterior Probabilities`[,5])
#   }
# }
# cut.off.45 <- quantile(calmat45,0.9)
# 
# calmat <- 0
# for(i in 1:8){
#   if(i==1|i==2|i==3|i==7|i==8){
#     calmat <- c(calmat,exnex.sim[[i]]$`Posterior Probabilities`[,3])}
# }
# for(j in 1:4){
#   if(j==1){
#     calmat <- c(calmat,extra.exnex.sim[[j]]$`Posterior Probabilities`[,3])
#   }
# }
# cut.off.123 <- quantile(calmat,0.9)
# 
# 
# 
# ErrorMat <- matrix(NA,ncol=K,nrow=8)
# for(i in 1:8){
#   post.probs <- exnex.sim[[i]]$`Posterior Probabilities`
#   for(k in 1:3){
#     ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off.123))*100
#   }
#   for(j in 4:5){
#     ErrorMat[i,j] <- mean(as.integer(post.probs[,j]>cut.off.45))*100
#   }
# }
# ErrorMatExtra <- matrix(NA,ncol=K,nrow=4)
# for(i in 1:4){
#   post.probs <- extra.exnex.sim[[i]]$`Posterior Probabilities`
#   for(k in 1:3){
#     ErrorMatExtra[i,k] <- mean(as.integer(post.probs[,k]>cut.off.123))*100
#   }
#   for(j in 4:5){
#     ErrorMatExtra[i,j] <- mean(as.integer(post.probs[,j]>cut.off.45))*100
#   }
# }
# ErrorCombined <- rbind(ErrorMat,ErrorMatExtra)
# 
# perfect <- rep(0,8)
# for(j in 1:8){
#   true <- as.numeric(list.scenarios[[j]]>q0)
#   post.prob <- exnex.sim[[j]]$`Posterior Probabilities`
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:3){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off.123)
#     }
#     for(t in 4:5){
#       hypo[t] <- as.numeric(post.prob[i,t]>cut.off.45)
#     }
#     if(all(hypo==true)){
#       perfect[j] <- perfect[j]+1
#     }
#   }
# }
# Perfect <- 100*(perfect/run)
# perfectExtra <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios[[j]]>q0)
#   post.prob <- extra.exnex.sim[[j]]$`Posterior Probabilities`
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:3){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off.123)
#     }
#     for(t in 4:5){
#       hypo[t] <- as.numeric(post.prob[i,t]>cut.off.45)
#     }
#     if(all(hypo==true)){
#       perfectExtra[j] <- perfectExtra[j]+1
#     }
#   }
# }
# PerfectExtra <- 100*(perfect/run)
# PerfectCombined <- c(Perfect,PerfectExtra)
# 
# fwer <- rep(0,8)
# for(j in 1:8){
#   true <- as.numeric(list.scenarios[[j]]>q0)
#   fwer.true <- which(true==0)
#   post.prob <- exnex.sim[[j]]$`Posterior Probabilities`
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:3){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off.123)
#     }
#     for(t in 4:5){
#       hypo[t] <- as.numeric(post.prob[i,t]>cut.off.45)
#     }
#     if(sum(hypo[fwer.true])!=0){
#       fwer[j] <- fwer[j]+1
#     }
#   }
# }
# FWER <- 100*(fwer/run)
# fwerExtra <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios[[j]]>q0)
#   fwer.true <- which(true==0)
#   post.prob <- extra.exnex.sim[[j]]$`Posterior Probabilities`
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:3){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off.123)
#     }
#     for(t in 4:5){
#       hypo[t] <- as.numeric(post.prob[i,t]>cut.off.45)
#     }
#     if(sum(hypo[fwer.true])!=0){
#       fwerExtra[j] <- fwer[j]+1
#     }
#   }
# }
# FWERExtra <- 100*(fwer/run)
# FWERCombined <- c(FWER,FWERExtra)
# 
# point.mean <- matrix(NA,nrow=12,ncol=K)
# point.sd <- matrix(NA,nrow=12,ncol=K)
# for(i in 1:8){
#   point.mean[i,] <- exnex.sim[[i]]$`Point Estimate Means`
#   point.sd[i,] <- exnex.sim[[i]]$`Point Estimate Sds`
# }
# for(j in 9:12){
#   point.mean[j,] <- extra.exnex.sim[[j-8]]$`Point Estimate Means`
#   point.sd[j,] <- extra.exnex.sim[[j-8]]$`Point Estimate Sds`
# }
# point.mean
# point.sd
# 
# extra.exnex.results <- list('Reject'=ErrorCombined,'FWER'=FWERCombined,'Perfect'=PerfectCombined,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)
# 
# save(extra.exnex.results,file='EXNEX_Results_Extra.RData')
