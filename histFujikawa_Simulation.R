library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)

#Adapted Power Prior/Fujikawa's Design---------------------------------------------------
adapt.fujikawa <- function(y,n,yh,nh,tau,hist.scalar,epsilon){
  post <- list()
  hist <- list()
  post.list <- c()
  x <- seq(0,1,0.00001)
  wh <- matrix(NA,nrow=K,ncol=K)
  for(i in 1:K){
    post[[i]] <- dbeta(x,1+y[i],1+n[i]-y[i])/sum(dbeta(x,1+y[i],1+n[i]-y[i]))
    hist[[i]] <- dbeta(x,1+yh[i],1+nh[i]-yh[i])/sum(dbeta(x,1+yh[i],1+nh[i]-yh[i]))
    post.list <- rbind(post.list,post[[i]])
  }
  wc <- (1-JSD(post.list,unit='log2'))^epsilon
  for(k in 1:K){
    for(i in 1:K){
      if(wc[k,i]<tau){
        wc[k,i] <- 0
      }
      if(nh[i]!=0){
        wh[k,i] <- (1-JSD(rbind(post[[k]],hist[[i]]),unit='log2'))^epsilon}else{
          wh[k,i] <- 0
        }
      if(wh[k,i]<tau){wh[k,i] <- 0}
    }
  }
  wh
  wc
  
  a.par <- c()
  b.par <- c()
  for(k in 1:K){
    a.par[k] <- 1+hist.scalar*sum(wh[k,]*yh)+sum(wc[k,]*y)
    b.par[k] <- 1+hist.scalar*sum(wh[k,]*(nh-yh))+sum(wc[k,]*(n-y))
  }
  my_list <- list('Parameters'=cbind(a.par,b.par),'Historic weights'=wh,'Current weights'=wc)
  return(my_list)
}


#Simulation-----------------------------
Fujikawa_Sim <- function(pmat,K,n,nh,tau,hist.scalar,epsilon,q0,run){
  H <- sum(nh!=0)
  p <- pmat[1:K]
  yh <- c(pmat[K+1:H],rep(0,K-H))
  post.prob <- matrix(NA,nrow=run,ncol=K)
  pointests <- matrix(NA,nrow=run,ncol=K)

  for(i in 1:run){
    y <- rbinom(K,n,p)
    Par <- adapt.fujikawa(y,n,yh,nh,tau,hist.scalar,epsilon)$Parameters
    for(k in 1:K){
      post.prob[i,k] <- pbeta(q0,Par[k,1],Par[k,2],lower.tail=F)
      pointests[i,k] <- Par[k,1]/(Par[k,1]+Par[k,2])
    }
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
tau <- 0.2
hist.scalar <- 0.8
epsilon <- 2
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
fujikawa.sim0 <- mclapply(list.scenarios1,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=9)
save(fujikawa.sim0,file='Fujikawa_Sim0.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,1,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,1,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,1,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,1,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,1,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,1,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,1,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,1,1)

list.scenarios2 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
fujikawa.sim1 <- mclapply(list.scenarios2,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=9)
save(fujikawa.sim1,file='Fujikawa_Sim1.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,1)

list.scenarios3 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
fujikawa.sim2 <- mclapply(list.scenarios3,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=9)
save(fujikawa.sim2,file='Fujikawa_Sim2.RData')

p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,3)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,3)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,3)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,3)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,3)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,3)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,3)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,3)

list.scenarios4 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
fujikawa.sim3 <- mclapply(list.scenarios4,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=9)
save(fujikawa.sim3,file='Fujikawa_Sim3.RData')


#Extra Scenarios
p9 <- c(0.25,0.25,0.1,0.25,0.1,1,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,1,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,1,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,1,1,1)

extra.list.scenarios1 <- list(p9,p10,p11,p12)
extra.fujikawa.sim0 <- mclapply(extra.list.scenarios1,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=5)
save(extra.fujikawa.sim0,file='Fujikawa_Sim0_Extra.RData')


p9 <- c(0.25,0.25,0.1,0.25,0.1,3,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,1,1)

extra.list.scenarios2 <- list(p9,p10,p11,p12)
extra.fujikawa.sim1 <- mclapply(extra.list.scenarios2,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=5)
save(extra.fujikawa.sim1,file='Fujikawa_Sim1_Extra.RData')

p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,1)

extra.list.scenarios3 <- list(p9,p10,p11,p12)
extra.fujikawa.sim2 <- mclapply(extra.list.scenarios3,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=5)
save(extra.fujikawa.sim2,file='Fujikawa_Sim2_Extra.RData')


p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,3)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,3)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,3)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,3)

extra.list.scenarios4 <- list(p9,p10,p11,p12)
extra.fujikawa.sim3 <- mclapply(extra.list.scenarios4,Fujikawa_Sim,K,n,nh,tau,hist.scalar,epsilon,q0,run,mc.cores=5)
save(extra.fujikawa.sim3,file='Fujikawa_Sim3_Extra.RData')



#Results Analysis---------------------------------------------------
# calmat <- 0
# for(i in 1:5){
#   calmat <- c(calmat,fujikawa.sim3[[i]]$`Posterior Probabilities`[,5])
# }
# cut.off <- quantile(calmat,0.9)

#cut.off <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9))

#Sim 0------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,fujikawa.sim0[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,fujikawa.sim0[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)

cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- fujikawa.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios1[[j]][1:5]>q0)
  post.prob <- fujikawa.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- fujikawa.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- fujikawa.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- fujikawa.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


fujikawa.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(fujikawa.results0,file='Fujikawa_Results0.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.fujikawa.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
  post.prob <- extra.fujikawa.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.fujikawa.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.fujikawa.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.fujikawa.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.fujikawa.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.fujikawa.results0,file='Fujikawa_Results0_Extra.RData')


#Sim 1------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,fujikawa.sim1[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat23 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat23 <- c(calmat23,fujikawa.sim1[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.23 <- quantile(calmat23,0.9)

calmat1 <- 0
for(i in 1:8){
  if(i==1|i==7){
    calmat1 <- c(calmat1,fujikawa.sim1[[i]]$`Posterior Probabilities`[,1])}
}
cut.off.1 <- quantile(calmat1,0.9)

cut.off <- c(cut.off.1,rep(cut.off.23,2),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- fujikawa.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios2[[j]][1:5]>q0)
  post.prob <- fujikawa.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- fujikawa.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- fujikawa.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- fujikawa.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


fujikawa.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(fujikawa.results1,file='Fujikawa_Results1.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.fujikawa.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
  post.prob <- extra.fujikawa.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.fujikawa.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.fujikawa.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.fujikawa.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.fujikawa.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.fujikawa.results1,file='Fujikawa_Results1_Extra.RData')


#Sim 2------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,fujikawa.sim2[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat3 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat3 <- c(calmat3,fujikawa.sim2[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.3 <- quantile(calmat3,0.9)

calmat12 <- 0
for(i in 1:8){
  if(i==1|i==7|i==2|i==8){
    calmat12 <- c(calmat12,fujikawa.sim2[[i]]$`Posterior Probabilities`[,2])}
}
cut.off.12 <- quantile(calmat12,0.9)


cut.off <- c(rep(cut.off.12,2),cut.off.3,rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- fujikawa.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios3[[j]][1:5]>q0)
  post.prob <- fujikawa.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- fujikawa.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- fujikawa.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- fujikawa.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


fujikawa.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(fujikawa.results2,file='Fujikawa_Results2.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.fujikawa.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
  post.prob <- extra.fujikawa.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.fujikawa.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.fujikawa.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.fujikawa.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.fujikawa.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.fujikawa.results2,file='Fujikawa_Results2_Extra.RData')


#Sim 3------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,fujikawa.sim3[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,fujikawa.sim3[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)


cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- fujikawa.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios4[[j]][1:5]>q0)
  post.prob <- fujikawa.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- fujikawa.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- fujikawa.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- fujikawa.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


fujikawa.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(fujikawa.results3,file='Fujikawa_Results3.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.fujikawa.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios4[[j]][1:5]>q0)
  post.prob <- extra.fujikawa.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.fujikawa.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.fujikawa.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.fujikawa.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.fujikawa.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.fujikawa.results3,file='Fujikawa_Results3_Extra.RData')


