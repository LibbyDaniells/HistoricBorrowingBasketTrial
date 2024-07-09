library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)


#Calibrate Across Scenarios
EXNEX_Sim <- function(p,n,yh,nh,K,pw,q0,run,pi){
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
    y.new <- y+yh
    n.new <- n+nh
    exnex2.data <- list('K'=K,'n'=n.new,'y'=y.new,'q0'=q0,'nexmu'=nexmu,'nexsigma'=nexsigma,'pi'=pi)
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
nh <- c(13,13,13,0,0)


p1 <- c(0.1,0.1,0.1,0.1,0.1)
p2 <- c(0.25,0.1,0.1,0.1,0.1)
p3 <- c(0.25,0.25,0.1,0.1,0.1)
p4 <- c(0.25,0.25,0.25,0.1,0.1)
p5 <- c(0.25,0.25,0.25,0.25,0.1)
p6 <- c(0.25,0.25,0.25,0.25,0.25)
p7 <- c(0.1,0.1,0.1,0.25,0.1)
p8 <- c(0.25,0.1,0.1,0.25,0.1)
list.scenarios <- list(p1,p2,p3,p4,p5,p6,p7,p8)

yh <- c(1,1,1,0,0)
exnex.sim0 <- mclapply(list.scenarios,EXNEX_Sim,n,yh,nh,K,pw,q0,run,pi,mc.cores=2)
save(exnex.sim0,file='EXNEXAllV2_Sim0.RData')

yh <- c(3,1,1,0,0)
exnex.sim1 <- mclapply(list.scenarios,EXNEX_Sim,n,yh,nh,K,pw,q0,run,pi,mc.cores=2)
save(exnex.sim1,file='EXNEXAllV2_Sim1.RData')

yh <- c(3,3,1,0,0)
exnex.sim2 <- mclapply(list.scenarios,EXNEX_Sim,n,yh,nh,K,pw,q0,run,pi,mc.cores=2)
save(exnex.sim2,file='EXNEXAllV2_Sim2.RData')

yh <- c(3,3,3,0,0)
exnex.sim3 <- mclapply(list.scenarios,EXNEX_Sim,n,yh,nh,K,pw,q0,run,pi,mc.cores=2)
save(exnex.sim3,file='EXNEXAllV2_Sim3.RData')




#Sim 0------------------------------------------------
exnexall.sim0 <- exnex.sim0
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,exnexall.sim0[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,exnexall.sim0[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)

cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- exnexall.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  post.prob <- exnexall.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- exnexall.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- exnexall.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- exnexall.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


exnexall.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(exnexall.results0,file='EXNEXAllV2_Results0.RData')



#Sim 1------------------------------------------------
exnexall.sim1 <- exnex.sim1
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,exnexall.sim1[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat23 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat23 <- c(calmat23,exnexall.sim1[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.23 <- quantile(calmat23,0.9)

calmat1 <- 0
for(i in 1:8){
  if(i==1|i==7){
    calmat1 <- c(calmat1,exnexall.sim1[[i]]$`Posterior Probabilities`[,1])}
}
cut.off.1 <- quantile(calmat1,0.9)

cut.off <- c(cut.off.1,rep(cut.off.23,2),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- exnexall.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  post.prob <- exnexall.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- exnexall.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- exnexall.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- exnexall.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


exnexall.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(exnexall.results1,file='EXNEXAllV2_Results1.RData')



#Sim 2------------------------------------------------
exnexall.sim2 <- exnex.sim2
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,exnexall.sim2[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat3 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat3 <- c(calmat3,exnexall.sim2[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.3 <- quantile(calmat3,0.9)

calmat12 <- 0
for(i in 1:8){
  if(i==1|i==7|i==2|i==8){
    calmat12 <- c(calmat12,exnexall.sim2[[i]]$`Posterior Probabilities`[,2])}
}
cut.off.12 <- quantile(calmat12,0.9)


cut.off <- c(rep(cut.off.12,2),cut.off.3,rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- exnexall.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  post.prob <- exnexall.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- exnexall.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- exnexall.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- exnexall.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


exnexall.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(exnexall.results2,file='EXNEXAllV2_Results2.RData')


#Sim 3------------------------------------------------
exnexall.sim3 <- exnex.sim3
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,exnexall.sim3[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,exnexall.sim3[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)


cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- exnexall.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  post.prob <- exnexall.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(list.scenarios[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- exnexall.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- exnexall.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- exnexall.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


exnexall.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(exnexall.results3,file='EXNEXAllV2_Results3.RData')


