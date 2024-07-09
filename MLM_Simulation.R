library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)



#Multi-Level Mixture Model-------------------------------------------------------
MLMixture <- function(K,y,yh,n,nh,q0,H,pi.delta,pi.epsilon1,pi.epsilon2,a,b){
  m <- matrix(1,nrow=K,ncol=K)
  diag(m) <- 4
  mod1 <- cbind(m,matrix(1,nrow=K,ncol=H),rep(3,K))
  mod2 <- cbind(m,rep(3,K))
  yh.new <- yh[which(nh!=0)]
  nh.new <- nh[which(nh!=0)]
  ycomb <- cbind(matrix(rep(c(y,yh.new),K),nrow=K,ncol=K+H,byrow=T),y,matrix(rep(y,K),nrow=K,ncol=K,byrow=T),y)
  ncomb <- cbind(matrix(rep(c(n,nh.new),K),nrow=K,ncol=K+H,byrow=T),n,matrix(rep(n,K),nrow=K,ncol=K,byrow=T),n)
  yhcomb <- cbind(matrix(rep(c(yh,rep(0,H)),K),nrow=K,ncol=K+H,byrow=T),yh)
  nhcomb <- cbind(matrix(rep(c(nh,rep(0,H)),K),nrow=K,ncol=K+H,byrow=T),nh)
  exnex.data <- list('pi.delta'=pi.delta,'K'=K,'H'=H,'q0'=q0,'y'=ycomb,'n'=ncomb,'pi.epsilon1'=pi.epsilon1,'pi.epsilon2'=pi.epsilon2,'mod1'=mod1,'mod2'=mod2,'a'=a,'b'=b,'yh'=yhcomb,'nh'=nhcomb)
  jags.exnex <- jags.model(file='MultiLevelMixture.txt',data=exnex.data,n.adapt=100000,n.chains=4)
  samples.exnex <- coda.samples(jags.exnex,variable.names=c('p.extract','Delta','epsilon.extract'),n.iter=100000,silent=T)
  mix.ex <- as.data.frame(samples.exnex[[1]])
  return(mix.ex)
}


#Calibration---------------------------------------
MLMixture_Cal <- function(K,p,n,ph,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run){
  Delta.mat <- matrix(NA,nrow=run,ncol=K)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    yh <- rbinom(K,nh,ph)
    H <- sum(ph!=0)
    m <- matrix(1,nrow=K,ncol=K)
    diag(m) <- 4
    mod1 <- cbind(m,matrix(1,nrow=K,ncol=H),rep(3,K))
    mod2 <- cbind(m,rep(3,K))
    yh.new <- yh[which(nh!=0)]
    nh.new <- nh[which(nh!=0)]
    ycomb <- cbind(matrix(rep(c(y,yh.new),K),nrow=K,ncol=K+H,byrow=T),y,matrix(rep(y,K),nrow=K,ncol=K,byrow=T),y)
    ncomb <- cbind(matrix(rep(c(n,nh.new),K),nrow=K,ncol=K+H,byrow=T),n,matrix(rep(n,K),nrow=K,ncol=K,byrow=T),n)
    yhcomb <- cbind(matrix(rep(c(yh,rep(0,H)),K),nrow=K,ncol=K+H,byrow=T),yh)
    nhcomb <- cbind(matrix(rep(c(nh,rep(0,H)),K),nrow=K,ncol=K+H,byrow=T),nh)
    exnex.data <- list('pi.delta'=pi.delta,'K'=K,'H'=H,'q0'=q0,'y'=ycomb,'n'=ncomb,'pi.epsilon1'=pi.epsilon1,'pi.epsilon2'=pi.epsilon2,'mod1'=mod1,'mod2'=mod2,'a'=a,'b'=b,'yh'=yhcomb,'nh'=nhcomb)
    jags.exnex <- jags.model(file='MultiLevelMixture.txt',data=exnex.data,n.adapt=1000,n.chains=4,quiet=T)
    samples.exnex <- coda.samples(jags.exnex,variable.names=c('p.extract'),n.iter=100000,silent=T)
    mix.ex <- as.data.frame(samples.exnex[[1]])
    Delta.mat[i,] <- apply(mix.ex,2,post.prob.fun)
    print(i)
  }
  Delta <- colQuantiles(Delta.mat,probs=0.9)
  return(Delta)
}


# K <- 5
# p <- rep(0.1,K)
# ph <- c(0.1,0.1,0.1,0,0)
# nh <- c(13,13,13,0,0)
# n <- rep(34,K)
# q0 <- 0.1
# run <- 10000
# pw <- 0.2
# a <- 1
# b <- 1
# 
# pi.delta <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# pi.epsilon1 <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# pi.epsilon2 <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# 
# 
# 
# mlmixture <- MLMixture_Cal(K,p,n,ph,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run)
# print(mlmixture)
# 
# save(mlmixture.current,file='MLMixture_Delta.RData')


#MLM Simulation--------------------------------------------------------------------
MLMixture_Sim <- function(pmat,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run){
  H <- sum(nh!=0)
  p <- pmat[1:K]
  yh <- c(pmat[K+1:H],rep(0,K-H))
  post.prob <- matrix(NA,nrow=run,ncol=K)
  pointests <- matrix(NA,nrow=run,ncol=K)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    m <- matrix(1,nrow=K,ncol=K)
    diag(m) <- 4
    mod1 <- cbind(m,matrix(1,nrow=K,ncol=H),rep(3,K))
    mod2 <- cbind(m,rep(3,K))
    yh.new <- yh[which(nh!=0)]
    nh.new <- nh[which(nh!=0)]
    ycomb <- cbind(matrix(rep(c(y,yh.new),K),nrow=K,ncol=K+H,byrow=T),y,matrix(rep(y,K),nrow=K,ncol=K,byrow=T),y)
    ncomb <- cbind(matrix(rep(c(n,nh.new),K),nrow=K,ncol=K+H,byrow=T),n,matrix(rep(n,K),nrow=K,ncol=K,byrow=T),n)
    yhcomb <- cbind(matrix(rep(c(yh,rep(0,H)),K),nrow=K,ncol=K+H,byrow=T),yh)
    nhcomb <- cbind(matrix(rep(c(nh,rep(0,H)),K),nrow=K,ncol=K+H,byrow=T),nh)
    exnex.data <- list('pi.delta'=pi.delta,'K'=K,'H'=H,'q0'=q0,'y'=ycomb,'n'=ncomb,'pi.epsilon1'=pi.epsilon1,'pi.epsilon2'=pi.epsilon2,'mod1'=mod1,'mod2'=mod2,'a'=a,'b'=b,'yh'=yhcomb,'nh'=nhcomb)
    jags.exnex <- jags.model(file='MultiLevelMixture.txt',data=exnex.data,n.adapt=1000,n.chains=4,quiet=T)
    samples.exnex <- coda.samples(jags.exnex,variable.names=c('p.extract'),n.iter=100000,silent=T)
    mix.ex <- as.data.frame(samples.exnex[[1]])
    post.prob[i,] <- apply(mix.ex,2,post.prob.fun)
    pointests[i,] <- colMeans(mix.ex)
    print(i)
  }
  Pointmeans <- colMeans(pointests)
  Pointsds <- apply(pointests,2,sd)
  my_list <- list('Point Estimates'=pointests,'Point Estimate Means'=Pointmeans,'Point Estimate Sds'=Pointsds,'Posterior Probabilities'=post.prob)
  return(my_list)
}


K <- 5
p <- rep(0.1,K)
ph <- c(0.1,0.1,0.1,0,0)
nh <- c(13,13,13,0,0)
n <- rep(34,K)
q0 <- 0.1
run <- 5000
pw <- 0.2
a <- 1
b <- 1

pi.delta <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
pi.epsilon1 <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
pi.epsilon2 <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))


p1 <- c(0.1,0.1,0.1,0.1,0.1,1,1,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,1,1,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,1,1,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,1,1,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,1,1,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,1,1,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,1,1,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,1,1,1)

list.scenarios1 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
mlmixture.sim0 <- mclapply(list.scenarios1,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=9)
save(mlmixture.sim0,file='MLMixture_Sim0.RData')
# 
p1 <- c(0.1,0.1,0.1,0.1,0.1,3,1,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,1,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,1,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,1,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,1,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,1,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,1,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,1,1)

list.scenarios2 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# mlmixture.sim1 <- mclapply(list.scenarios2,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=9)
# save(mlmixture.sim1,file='MLMixture_Sim1.RData')
# 
p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,1)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,1)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,1)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,1)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,1)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,1)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,1)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,1)

list.scenarios3 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# mlmixture.sim2 <- mclapply(list.scenarios3,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=9)
# save(mlmixture.sim2,file='MLMixture_Sim2.RData')
# 
p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,3)
p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,3)
p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,3)
p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,3)
p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,3)
p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,3)
p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,3)
p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,3)

list.scenarios4 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# mlmixture.sim3 <-mclapply(list.scenarios4,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=9)
# save(mlmixture.sim3,file='MLMixture_Sim3.RData')



#Extra Scenarios
p9 <- c(0.25,0.25,0.1,0.25,0.1,1,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,1,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,1,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,1,1,1)

extra.list.scenarios1 <- list(p9,p10,p11,p12)
extra.mlmixture.sim0 <-mclapply(extra.list.scenarios1,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=5)
save(extra.mlmixture.sim0,file='MLMixture_Sim0_Extra.RData')
# 
p9 <- c(0.25,0.25,0.1,0.25,0.1,3,1,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,1,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,1,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,1,1)

extra.list.scenarios2 <- list(p9,p10,p11,p12)
# extra.mlmixture.sim1 <-mclapply(extra.list.scenarios2,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=5)
# save(extra.mlmixture.sim1,file='MLMixture_Sim1_Extra.RData')
# 
# 
p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,1)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,1)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,1)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,1)

extra.list.scenarios3 <- list(p9,p10,p11,p12)
# extra.mlmixture.sim2 <-mclapply(extra.list.scenarios3,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=5)
# save(extra.mlmixture.sim2,file='MLMixture_Sim2_Extra.RData')
# 
# 
# 
p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,3)
p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,3)
p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,3)
p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,3)

extra.list.scenarios4 <- list(p9,p10,p11,p12)
# extra.mlmixture.sim3 <-mclapply(extra.list.scenarios4,MLMixture_Sim,K,n,nh,q0,pi.delta,pi.epsilon1,pi.epsilon2,a,b,run,mc.cores=5)
# save(extra.mlmixture.sim3,file='MLMixture_Sim3_Extra.RData')






#Sim 0------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mlmixture.sim0[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,mlmixture.sim0[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)

cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mlmixture.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios1[[j]][1:5]>q0)
  post.prob <- mlmixture.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- mlmixture.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- mlmixture.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mlmixture.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mlmixture.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mlmixture.results0,file='MLMixture_Results0.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mlmixture.sim0[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
  post.prob <- extra.mlmixture.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.mlmixture.sim0[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.mlmixture.sim0[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mlmixture.sim0[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mlmixture.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mlmixture.results0,file='MLMixture_Results0_Extra.RData')


#Sim 1------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mlmixture.sim1[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat23 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat23 <- c(calmat23,mlmixture.sim1[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.23 <- quantile(calmat23,0.9)

calmat1 <- 0
for(i in 1:8){
  if(i==1|i==7){
    calmat1 <- c(calmat1,mlmixture.sim1[[i]]$`Posterior Probabilities`[,1])}
}
cut.off.1 <- quantile(calmat1,0.9)

cut.off <- c(cut.off.1,rep(cut.off.23,2),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mlmixture.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios2[[j]][1:5]>q0)
  post.prob <- mlmixture.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- mlmixture.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- mlmixture.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mlmixture.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mlmixture.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mlmixture.results1,file='MLMixture_Results1.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mlmixture.sim1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
  post.prob <- extra.mlmixture.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.mlmixture.sim1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.mlmixture.sim1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mlmixture.sim1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mlmixture.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mlmixture.results1,file='MLMixture_Results1_Extra.RData')



#Sim 2------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mlmixture.sim2[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat3 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat3 <- c(calmat3,mlmixture.sim2[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.3 <- quantile(calmat3,0.9)

calmat12 <- 0
for(i in 1:8){
  if(i==1|i==7|i==2|i==8){
    calmat12 <- c(calmat12,mlmixture.sim2[[i]]$`Posterior Probabilities`[,2])}
}
cut.off.12 <- quantile(calmat12,0.9)


cut.off <- c(rep(cut.off.12,2),cut.off.3,rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mlmixture.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios3[[j]][1:5]>q0)
  post.prob <- mlmixture.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- mlmixture.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- mlmixture.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mlmixture.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mlmixture.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mlmixture.results2,file='MLMixture_Results2.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mlmixture.sim2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
  post.prob <- extra.mlmixture.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.mlmixture.sim2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.mlmixture.sim2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mlmixture.sim2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mlmixture.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mlmixture.results2,file='MLMixture_Results2_Extra.RData')


#Sim 3------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,mlmixture.sim3[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,mlmixture.sim3[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)


cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- mlmixture.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios4[[j]][1:5]>q0)
  post.prob <- mlmixture.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- mlmixture.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- mlmixture.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- mlmixture.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


mlmixture.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(mlmixture.results3,file='MLMixture_Results3.RData')

#Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.mlmixture.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios4[[j]][1:5]>q0)
  post.prob <- extra.mlmixture.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.mlmixture.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.mlmixture.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.mlmixture.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.mlmixture.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.mlmixture.results3,file='MLMixture_Results3_Extra.RData')

