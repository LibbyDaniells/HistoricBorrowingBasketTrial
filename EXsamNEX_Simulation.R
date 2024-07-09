library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)

#Simulation-----------------------------
SAMPrior_Sim <- function(pmat,n,nh,K,q0,pw,pi,run,clin.diff,a,b){
  H <- sum(nh!=0)
  p <- pmat[1:K]
  yh <- c(pmat[K+1:H],rep(0,K-H))
  hist <- as.numeric(nh>0)
  post.prob <- matrix(NA,nrow=run,ncol=K)
  pointests <- matrix(NA,nrow=run,ncol=K)
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(j in 1:run){
    y <- rbinom(K,n,p)
    w <- rep(0,K)
    for(i in 1:K){
      theta_hat <- (a+yh[i])/(a+b+nh[i])
      m1 <- ((theta_hat+clin.diff)^y[i])*((1-theta_hat-clin.diff)^(n[i]-y[i]))
      m2 <- ((theta_hat-clin.diff)^y[i])*((1-theta_hat+clin.diff)^(n[i]-y[i]))
      R <- ((theta_hat^y[i])*(1-theta_hat)^(n[i]-y[i]))/max(m1,m2) 
      if(hist[i]==1){
        w[i] <- R/(1+R)}else{
          w[i] <- 0
        }
    }
    pi.sam <- cbind(w,1-w)
    jags.data.sam <- list('K'=K,'n'=n,'y'=y,'q0'=q0,'pi'=pi,'pi.sam'=pi.sam,'hist'=hist,'a'=a,'b'=b,'yh'=yh,'nh'=nh)
    jags.fit.sam <- jags.model(file='SAMPriorinNEX.txt',data=jags.data.sam,n.adapt=1000,n.chains=4,quiet=T)
    samples.sam <- coda.samples(jags.fit.sam,variable.names=c('p'),n.iter=100000,silent=T)
    sam <- as.data.frame(samples.sam[[1]])
    pointests[j,] <- colMeans(sam)
    post.prob[j,] <- apply(sam,2,post.prob.fun)
    print(j)
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
clin.diff <- 0.15
a <- 1
b <- 1
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
samprior.sim0 <- mclapply(list.scenarios1,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=9)
save(samprior.sim0,file='SAMPrior_Sim0.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,1,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,1,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,1,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,1,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,1,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,1,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,1,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,1,1)

list.scenarios2 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
samprior.sim1 <- mclapply(list.scenarios2,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=9)
save(samprior.sim1,file='SAMPrior_Sim1.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,1)

list.scenarios3 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
samprior.sim2 <- mclapply(list.scenarios3,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=9)
save(samprior.sim2,file='SAMPrior_Sim2.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,3)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,3)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,3)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,3)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,3)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,3)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,3)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,3)

list.scenarios4 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
samprior.sim3 <- mclapply(list.scenarios4,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=9)
save(samprior.sim3,file='SAMPrior_Sim3.RData')


#Extra Scenarios
p9 <- c(0.25,0.25,0.1,0.25,0.1,1,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,1,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,1,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,1,1,1)

extra.list.scenarios1 <- list(p9,p10,p11,p12)
extra.samprior.sim0 <- mclapply(extra.list.scenarios1,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=5)
save(extra.samprior.sim0,file='SAMPrior_Sim0_Extra.RData')

p9 <- c(0.25,0.25,0.1,0.25,0.1,3,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,1,1)

extra.list.scenarios2 <- list(p9,p10,p11,p12)
extra.samprior.sim1 <- mclapply(extra.list.scenarios2,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=5)
save(extra.samprior.sim1,file='SAMPrior_Sim1_Extra.RData')

p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,1)

extra.list.scenarios3 <- list(p9,p10,p11,p12)
extra.samprior.sim2 <- mclapply(extra.list.scenarios3,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=5)
save(extra.samprior.sim2,file='SAMPrior_Sim2_Extra.RData')


p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,3)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,3)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,3)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,3)

extra.list.scenarios4 <- list(p9,p10,p11,p12)
extra.samprior.sim3 <- mclapply(extra.list.scenarios4,SAMPrior_Sim,n,nh,K,q0,pw,pi,run,clin.diff,a,b,mc.cores=5)
save(extra.samprior.sim3,file='SAMPrior_Sim3_Extra.RData')




#Results Analysis---------------------------------------------------
# calmat <- 0
# for(i in 1:5){
#   calmat <- c(calmat,samprior.sim1[[i]]$`Posterior Probabilities`[,5])
# }
# cut.off <- quantile(calmat,0.9)

#cut.off <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9))

#Sim 0------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,samprior.sim0[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,samprior.sim0[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)

cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- samprior.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios1[[j]][1:5]>q0)
  post.prob <- samprior.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- samprior.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- samprior.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- samprior.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


samprior.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(samprior.results0,file='SAMprior_Results0.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.samprior.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
  post.prob <- extra.samprior.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.samprior.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.samprior.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.samprior.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.samprior.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.samprior.results0,file='SAMprior_Results0_Extra.RData')



#Sim 1------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,samprior.sim1[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat23 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat23 <- c(calmat23,samprior.sim1[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.23 <- quantile(calmat23,0.9)

calmat1 <- 0
for(i in 1:8){
  if(i==1|i==7){
    calmat1 <- c(calmat1,samprior.sim1[[i]]$`Posterior Probabilities`[,1])}
}
cut.off.1 <- quantile(calmat1,0.9)

cut.off <- c(cut.off.1,rep(cut.off.23,2),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- samprior.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios2[[j]][1:5]>q0)
  post.prob <- samprior.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- samprior.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- samprior.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- samprior.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


samprior.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(samprior.results1,file='SAMprior_Results1.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.samprior.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
  post.prob <- extra.samprior.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.samprior.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.samprior.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.samprior.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.samprior.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.samprior.results1,file='SAMprior_Results1_Extra.RData')




#Sim 2------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,samprior.sim2[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat3 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat3 <- c(calmat3,samprior.sim2[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.3 <- quantile(calmat3,0.9)

calmat12 <- 0
for(i in 1:8){
  if(i==1|i==7|i==2|i==8){
    calmat12 <- c(calmat12,samprior.sim2[[i]]$`Posterior Probabilities`[,2])}
}
cut.off.12 <- quantile(calmat12,0.9)


cut.off <- c(rep(cut.off.12,2),cut.off.3,rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- samprior.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios3[[j]][1:5]>q0)
  post.prob <- samprior.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- samprior.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- samprior.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- samprior.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


samprior.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(samprior.results2,file='SAMprior_Results2.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.samprior.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
  post.prob <- extra.samprior.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.samprior.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.samprior.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.samprior.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.samprior.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.samprior.results2,file='SAMprior_Results2_Extra.RData')



#Sim 3------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,samprior.sim3[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,samprior.sim3[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)


cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- samprior.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios4[[j]][1:5]>q0)
  post.prob <- samprior.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- samprior.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- samprior.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- samprior.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


samprior.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(samprior.results3,file='SAMprior_Results3.RData')


#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.samprior.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios4[[j]][1:5]>q0)
  post.prob <- extra.samprior.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.samprior.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.samprior.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.samprior.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.samprior.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.samprior.results3,file='SAMprior_Results3_Extra.RData')





