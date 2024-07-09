library(rjags)
library(textmineR)
library(matrixStats)
library(philentropy)
library(LaplacesDemon)
library(parallel)
library(MASS)

#Power prior in the NEX component--------------------------------------------------
PowerPrior <- function(K,y,n,q0,pi,pw,yh,nh,alpha0,a,b){
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  hist <- as.numeric(nh>0)
  jags.data <- list('n'=n,'K'=K,'y'=y,'pi'=pi,'nexmu'=nexmu,'nexsigma'=nexsigma,
                    'yh'=yh,'nh'=nh,'alpha0'=alpha0,'a'=a,'b'=b,'hist'=hist,'q0'=q0)
  jags.fit <- jags.model(file='PowerPriorinNEX.txt',data=jags.data,n.adapt=100000,n.chains=4)
  samples <- coda.samples(jags.fit,variable.names = c('p','weight'),n.iter=100000,silent=TRUE) #Fit the model
  powerp <- as.data.frame(samples[[1]])
  return(powerp)
}

#Calibration---------------------------------------
PowerPrior_Cal <- function(K,p,n,q0,pi,pw,ph,nh,alpha0,a,b,run){
  Delta.mat <- matrix(NA,nrow=run,ncol=K)
  nexmu <- rep(log(pw/(1-pw)),K)
  nexsigma <- rep((1/pw)+(1/(1-pw)),K)
  hist <- as.numeric(nh>0)
  post.prob.fun <- function(x){
    fun<-sum(x>q0)/100000
    return(fun)
  }
  for(i in 1:run){
    y <- rbinom(K,n,p)
    yh <- rbinom(K,nh,ph)
    jags.data <- list('n'=n,'K'=K,'y'=y,'pi'=pi,'nexmu'=nexmu,'nexsigma'=nexsigma,
                      'yh'=yh,'nh'=nh,'alpha0'=alpha0,'a'=a,'b'=b,'hist'=hist,'q0'=q0)
    jags.fit <- jags.model(file='PowerPriorinNEX.txt',data=jags.data,n.adapt=1000,n.chains=4,quiet=T)
    samples <- coda.samples(jags.fit,variable.names = c('p'),n.iter=100000,silent=TRUE) #Fit the model
    powerp <- as.data.frame(samples[[1]])
    Delta.mat[i,] <- apply(powerp,2,post.prob.fun)
    print(i)
  }
  Delta <- colQuantiles(Delta.mat,probs=0.9)
  return(Delta)
}

# 
# K <- 5
# p <- rep(0.1,K)
# ph <- c(0.1,0.1,0.1,0,0)
# nh <- c(13,13,13,0,0)
# n <- rep(34,K)
# q0 <- 0.1
# run <- 10000
# pw <- 0.2
# pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
# alpha0 <- 0.5
# a <- 1
# b <- 1
# 
# power.prior <- PowerPrior_Cal(K,p,n,q0,pi,pw,ph,nh,alpha0,a,b,run)
# print(power.prior)
# 
# save(power.prior,file='Power_Prior_Delta.RData')


# 
# 
# #Simulation-----------------------------
# PowerPrior_Sim <- function(pmat,n,nh,K,q0,pw,pi,run,alpha0,a,b,powerprior.delta){
#   H <- sum(nh!=0)
#   p <- pmat[1:K]
#   ph <- pmat[K+1:H]
#   hist <- as.numeric(nh>0)
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
#     jags.data <- list('n'=n,'K'=K,'y'=y,'pi'=pi,'nexmu'=nexmu,'nexsigma'=nexsigma,
#                       'yh'=yh,'nh'=nh,'alpha0'=alpha0,'a'=a,'b'=b,'hist'=hist,'q0'=q0)
#     jags.fit <- jags.model(file='PowerPriorinNEX.txt',data=jags.data,n.adapt=1000,n.chains=4,quiet=T)
#     samples <- coda.samples(jags.fit,variable.names = c('p'),n.iter=100000,silent=TRUE) #Fit the model
#     powerp <- as.data.frame(samples[[1]])
#     pointests[i,] <- colMeans(powerp)
#     hypo[i,] <- as.integer(apply(powerp,2,post.prob.fun)>powerprior.delta)
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
# powerprior.delta <- 0.8608096
# run <- 10000
# pw <- 0.2
# alpha0 <- 0.5
# a <- 1
# b <- 1
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
# powerprior.sim1_9 <- mclapply(list.scenarios1,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,powerprior.delta,mc.cores=10)
# save(powerprior.sim1_9,file='PowerPrior_Sim1_9.RData')
# 
# powerprior.sim10_18 <- mclapply(list.scenarios2,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,powerprior.delta,mc.cores=10)
# save(powerprior.sim10_18,file='PowerPrior_Sim10_18.RData')




#Simulation Across Scenarios-----------------------------
PowerPrior_Sim <- function(pmat,n,nh,K,q0,pw,pi,run,alpha0,a,b){
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
  for(i in 1:run){
    y <- rbinom(K,n,p)
    jags.data <- list('n'=n,'K'=K,'y'=y,'pi'=pi,'nexmu'=nexmu,'nexsigma'=nexsigma,
                      'yh'=yh,'nh'=nh,'alpha0'=alpha0,'a'=a,'b'=b,'hist'=hist,'q0'=q0)
    jags.fit <- jags.model(file='PowerPriorinNEX.txt',data=jags.data,n.adapt=1000,n.chains=4,quiet=T)
    samples <- coda.samples(jags.fit,variable.names = c('p'),n.iter=100000,silent=TRUE) #Fit the model
    powerp <- as.data.frame(samples[[1]])
    pointests[i,] <- colMeans(powerp)
    post.prob[i,] <- apply(powerp,2,post.prob.fun)
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
alpha0 <- 1
a <- 1
b <- 1
pi <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))

p1.1 <- c(0.1,0.1,0.1,0.1,0.1,1,1,1)
p2.1 <- c(0.25,0.1,0.1,0.1,0.1,1,1,1)
p3.1 <- c(0.25,0.25,0.1,0.1,0.1,1,1,1)
p4.1 <- c(0.25,0.25,0.25,0.1,0.1,1,1,1)
p5.1 <- c(0.25,0.25,0.25,0.25,0.1,1,1,1)
p6.1 <- c(0.25,0.25,0.25,0.25,0.25,1,1,1)
p1.2 <- c(0.1,0.1,0.1,0.1,0.1,3,1,1)
p2.2 <- c(0.25,0.1,0.1,0.1,0.1,3,1,1)
p3.2 <- c(0.25,0.25,0.1,0.1,0.1,3,1,1)
p4.2 <- c(0.25,0.25,0.25,0.1,0.1,3,1,1)
p5.2 <- c(0.25,0.25,0.25,0.25,0.1,3,1,1)
p6.2 <- c(0.25,0.25,0.25,0.25,0.25,3,1,1)
p1.3 <- c(0.1,0.1,0.1,0.1,0.1,3,3,1)
p2.3 <- c(0.25,0.1,0.1,0.1,0.1,3,3,1)
p3.3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,1)
p4.3 <- c(0.25,0.25,0.25,0.1,0.1,3,3,1)
p5.3 <- c(0.25,0.25,0.25,0.25,0.1,3,3,1)
p6.3 <- c(0.25,0.25,0.25,0.25,0.25,3,3,1)
p1.4 <- c(0.1,0.1,0.1,0.1,0.1,3,3,3)
p2.4 <- c(0.25,0.1,0.1,0.1,0.1,3,3,3)
p3.4 <- c(0.25,0.25,0.1,0.1,0.1,3,3,3)
p4.4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,3)
p5.4 <- c(0.25,0.25,0.25,0.25,0.1,3,3,3)
p6.4 <- c(0.25,0.25,0.25,0.25,0.25,3,3,3)
p7.1 <- c(0.1,0.1,0.1,0.25,0.1,1,1,1)
p7.2 <- c(0.1,0.1,0.1,0.25,0.1,3,1,1)
p7.3 <- c(0.1,0.1,0.1,0.25,0.1,3,3,1)
p7.4 <- c(0.1,0.1,0.1,0.25,0.1,3,3,3)
p8.1 <- c(0.25,0.1,0.1,0.25,0.1,1,1,1)
p8.2 <- c(0.25,0.1,0.1,0.25,0.1,3,1,1)
p8.3 <- c(0.25,0.1,0.1,0.25,0.1,3,3,1)
p8.4 <- c(0.25,0.1,0.1,0.25,0.1,3,3,3)


list.scenarios1 <- list(p1.1,p2.1,p3.1,p4.1,p5.1,p6.1,p7.1,p8.1,
                        p1.2,p2.2,p3.2,p4.2,p5.2,p6.2,p7.2,p8.2,
                        p1.3,p2.3,p3.3,p4.3,p5.3,p6.3,p7.3,p8.3,
                        p1.4,p2.4,p3.4,p4.4,p5.4,p6.4,p7.4,p8.4)

powerprior.sim1 <- mclapply(list.scenarios1,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=33)
save(powerprior.sim1,file='PowerPrior1_Sim.RData')
# 
# p1 <- c(0.1,0.1,0.1,0.1,0.1,1,1,1)
# p2 <- c(0.25,0.1,0.1,0.1,0.1,1,1,1)
# p3 <- c(0.25,0.25,0.1,0.1,0.1,1,1,1)
# p4 <- c(0.25,0.25,0.25,0.1,0.1,1,1,1)
# p5 <- c(0.25,0.25,0.25,0.25,0.1,1,1,1)
# p6 <- c(0.25,0.25,0.25,0.25,0.25,1,1,1)
# p7 <- c(0.1,0.1,0.1,0.25,0.1,1,1,1)
# p8 <- c(0.25,0.1,0.1,0.25,0.1,1,1,1)
# 
# list.scenarios1 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# powerprior.sim0 <- mclapply(list.scenarios1,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=9)
# save(powerprior.sim0,file='PowerPrior1_Sim0.RData')
# 
# p1 <- c(0.1,0.1,0.1,0.1,0.1,3,1,1)
# p2 <- c(0.25,0.1,0.1,0.1,0.1,3,1,1)
# p3 <- c(0.25,0.25,0.1,0.1,0.1,3,1,1)
# p4 <- c(0.25,0.25,0.25,0.1,0.1,3,1,1)
# p5 <- c(0.25,0.25,0.25,0.25,0.1,3,1,1)
# p6 <- c(0.25,0.25,0.25,0.25,0.25,3,1,1)
# p7 <- c(0.1,0.1,0.1,0.25,0.1,3,1,1)
# p8 <- c(0.25,0.1,0.1,0.25,0.1,3,1,1)
# 
# list.scenarios2 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# powerprior.sim1 <- mclapply(list.scenarios2,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=9)
# save(powerprior.sim1,file='PowerPrior1_Sim1.RData')
# 
# p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,1)
# p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,1)
# p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,1)
# p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,1)
# p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,1)
# p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,1)
# p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,1)
# p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,1)
# 
# list.scenarios3 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# powerprior.sim2 <- mclapply(list.scenarios3,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=9)
# save(powerprior.sim2,file='PowerPrior1_Sim2.RData')
# 
# p1 <- c(0.1,0.1,0.1,0.1,0.1,3,3,3)
# p2 <- c(0.25,0.1,0.1,0.1,0.1,3,3,3)
# p3 <- c(0.25,0.25,0.1,0.1,0.1,3,3,3)
# p4 <- c(0.25,0.25,0.25,0.1,0.1,3,3,3)
# p5 <- c(0.25,0.25,0.25,0.25,0.1,3,3,3)
# p6 <- c(0.25,0.25,0.25,0.25,0.25,3,3,3)
# p7 <- c(0.1,0.1,0.1,0.25,0.1,3,3,3)
# p8 <- c(0.25,0.1,0.1,0.25,0.1,3,3,3)
# 
# list.scenarios4 <- list(p1,p2,p3,p4,p5,p6,p7,p8)
# powerprior.sim3 <- mclapply(list.scenarios4,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=9)
# save(powerprior.sim3,file='PowerPrior1_Sim3.RData')

# 
# #Extra Scenarios
# p9 <- c(0.25,0.25,0.1,0.25,0.1,1,1,1)
# p10 <- c(0.1,0.1,0.1,0.25,0.25,1,1,1)
# p11 <- c(0.25,0.1,0.1,0.25,0.25,1,1,1)
# p12 <- c(0.25,0.25,0.1,0.25,0.25,1,1,1)
# 
# extra.list.scenarios1 <- list(p9,p10,p11,p12)
# extra.powerprior.sim0 <- mclapply(extra.list.scenarios1,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=5)
# save(extra.powerprior.sim0,file='PowerPrior_Sim0_Extra.RData')
# 
# p9 <- c(0.25,0.25,0.1,0.25,0.1,3,1,1)
# p10 <- c(0.1,0.1,0.1,0.25,0.25,3,1,1)
# p11 <- c(0.25,0.1,0.1,0.25,0.25,3,1,1)
# p12 <- c(0.25,0.25,0.1,0.25,0.25,3,1,1)
# 
# extra.list.scenarios2 <- list(p9,p10,p11,p12)
# extra.powerprior.sim1 <- mclapply(extra.list.scenarios2,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=5)
# save(extra.powerprior.sim1,file='PowerPrior_Sim1_Extra.RData')
# 
# p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,1)
# p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,1)
# p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,1)
# p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,1)
# 
# extra.list.scenarios3 <- list(p9,p10,p11,p12)
# extra.powerprior.sim2 <- mclapply(extra.list.scenarios3,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=5)
# save(extra.powerprior.sim2,file='PowerPrior_Sim2_Extra.RData')
# 
# 
# p9 <- c(0.25,0.25,0.1,0.25,0.1,3,3,3)
# p10 <- c(0.1,0.1,0.1,0.25,0.25,3,3,3)
# p11 <- c(0.25,0.1,0.1,0.25,0.25,3,3,3)
# p12 <- c(0.25,0.25,0.1,0.25,0.25,3,3,3)
# 
# extra.list.scenarios4 <- list(p9,p10,p11,p12)
# extra.powerprior.sim3 <- mclapply(extra.list.scenarios4,PowerPrior_Sim,n,nh,K,q0,pw,pi,run,alpha0,a,b,mc.cores=5)
# save(extra.powerprior.sim3,file='PowerPrior_Sim3_Extra.RData')

# 
# 
# 
#Results Analysis---------------------------------------------------
# calmat <- 0
# for(i in 1:5){
#   calmat <- c(calmat,powerprior.sim1[[i]]$`Posterior Probabilities`[,5])
# }
# cut.off <- quantile(calmat,0.9)

#cut.off <- mean(colQuantiles(exnex.sim[[1]]$`Posterior Probabilities`,probs=0.9))

#Sim 0------------------------------------------------
calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,powerprior.sim25[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,powerprior.sim25[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)

cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- powerprior.sim25[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(list.scenarios1[[j]][1:5]>q0)
  post.prob <- powerprior.sim25[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- powerprior.sim25[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- powerprior.sim25[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- powerprior.sim25[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


powerpriorAlpha25.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(powerpriorAlpha25.results0,file='PowerPriorAlpha25_Results0.RData')
# 
# ####Extra
# ErrorMat <- matrix(NA,ncol=K,nrow=4)
# for(i in 1:4){
#   post.probs <- extra.powerprior.sim0[[i]]$`Posterior Probabilities`
#   for(k in 1:5){
#     ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
#   }
# }
# 
# perfect <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
#   post.prob <- extra.powerprior.sim0[[j]]$`Posterior Probabilities`[,1:5]
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:5){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
#     }
#     if(all(hypo==true)){
#       perfect[j] <- perfect[j]+1
#     }
#   }
# }
# Perfect <- 100*(perfect/run)
# 
# fwer <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios1[[j]][1:5]>q0)
#   fwer.true <- which(true==0)
#   post.prob <- extra.powerprior.sim0[[j]]$`Posterior Probabilities`[,1:5]
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:5){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
#     }
#     if(sum(hypo[fwer.true])!=0){
#       fwer[j] <- fwer[j]+1
#     }
#   }
# }
# FWER <- 100*(fwer/run)
# 
# point.mean <- matrix(NA,nrow=4,ncol=K)
# point.sd <- matrix(NA,nrow=4,ncol=K)
# for(i in 1:4){
#   point.mean[i,] <- extra.powerprior.sim0[[i]]$`Point Estimate Means`[1:5]
#   point.sd[i,] <- extra.powerprior.sim0[[i]]$`Point Estimate Sds`[1:5]
# }
# point.mean
# point.sd
# 
# 
# extra.powerprior.results0 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)
# 
# save(extra.powerprior.results0,file='PowerPrior_Results0_Extra.RData')

#Sim 1------------------------------------------------
powerprior1 <- list()
listscenarios1 <- list()
for(i in 9:16){
  powerprior1[[i-8]] <- powerprior.sim25[[i]]
  listscenarios1[[i-8]] <- list.scenarios1[[i]]
}

calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,powerprior1[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat23 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat23 <- c(calmat23,powerprior1[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.23 <- quantile(calmat23,0.9)

calmat1 <- 0
for(i in 1:8){
  if(i==1|i==7){
    calmat1 <- c(calmat1,powerprior1[[i]]$`Posterior Probabilities`[,1])}
}
cut.off.1 <- quantile(calmat1,0.9)

cut.off <- c(cut.off.1,rep(cut.off.23,2),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- powerprior1[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(listscenarios1[[j]][1:5]>q0)
  post.prob <- powerprior1[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(listscenarios1[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- powerprior1[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- powerprior1[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- powerprior1[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


powerpriorAlpha25.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(powerpriorAlpha25.results1,file='PowerPriorAlpha25_Results1.RData')


####Extra
# ErrorMat <- matrix(NA,ncol=K,nrow=4)
# for(i in 1:4){
#   post.probs <- extra.powerprior.sim1[[i]]$`Posterior Probabilities`
#   for(k in 1:5){
#     ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
#   }
# }
# 
# perfect <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
#   post.prob <- extra.powerprior.sim1[[j]]$`Posterior Probabilities`[,1:5]
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:5){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
#     }
#     if(all(hypo==true)){
#       perfect[j] <- perfect[j]+1
#     }
#   }
# }
# Perfect <- 100*(perfect/run)
# 
# fwer <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios2[[j]][1:5]>q0)
#   fwer.true <- which(true==0)
#   post.prob <- extra.powerprior.sim1[[j]]$`Posterior Probabilities`[,1:5]
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:5){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
#     }
#     if(sum(hypo[fwer.true])!=0){
#       fwer[j] <- fwer[j]+1
#     }
#   }
# }
# FWER <- 100*(fwer/run)
# 
# point.mean <- matrix(NA,nrow=4,ncol=K)
# point.sd <- matrix(NA,nrow=4,ncol=K)
# for(i in 1:4){
#   point.mean[i,] <- extra.powerprior.sim1[[i]]$`Point Estimate Means`[1:5]
#   point.sd[i,] <- extra.powerprior.sim1[[i]]$`Point Estimate Sds`[1:5]
# }
# point.mean
# point.sd
# 
# 
# extra.powerprior.results1 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)
# 
# save(extra.powerprior.results1,file='PowerPrior_Results1_Extra.RData')
# 

#Sim 2------------------------------------------------
powerprior2 <- list()
listscenarios2 <- list()
for(i in 17:24){
  powerprior2[[i-16]] <- powerprior.sim25[[i]]
  listscenarios2[[i-16]] <- list.scenarios1[[i]]
}

calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,powerprior2[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat3 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat3 <- c(calmat3,powerprior2[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.3 <- quantile(calmat3,0.9)

calmat12 <- 0
for(i in 1:8){
  if(i==1|i==7|i==2|i==8){
    calmat12 <- c(calmat12,powerprior2[[i]]$`Posterior Probabilities`[,2])}
}
cut.off.12 <- quantile(calmat12,0.9)


cut.off <- c(rep(cut.off.12,2),cut.off.3,rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- powerprior2[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(listscenarios2[[j]][1:5]>q0)
  post.prob <- powerprior2[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(listscenarios2[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- powerprior2[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- powerprior2[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- powerprior2[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


powerpriorAlpha25.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(powerpriorAlpha25.results2,file='PowerPriorAlpha25_Results2.RData')

# 
# ####Extra
# ErrorMat <- matrix(NA,ncol=K,nrow=4)
# for(i in 1:4){
#   post.probs <- extra.powerprior.sim2[[i]]$`Posterior Probabilities`
#   for(k in 1:5){
#     ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
#   }
# }
# 
# perfect <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
#   post.prob <- extra.powerprior.sim2[[j]]$`Posterior Probabilities`[,1:5]
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:5){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
#     }
#     if(all(hypo==true)){
#       perfect[j] <- perfect[j]+1
#     }
#   }
# }
# Perfect <- 100*(perfect/run)
# 
# fwer <- rep(0,4)
# for(j in 1:4){
#   true <- as.numeric(extra.list.scenarios3[[j]][1:5]>q0)
#   fwer.true <- which(true==0)
#   post.prob <- extra.powerprior.sim2[[j]]$`Posterior Probabilities`[,1:5]
#   for(i in 1:run){
#     hypo <- rep(NA,K)
#     for(k in 1:5){
#       hypo[k] <- as.numeric(post.prob[i,k]>cut.off[k])
#     }
#     if(sum(hypo[fwer.true])!=0){
#       fwer[j] <- fwer[j]+1
#     }
#   }
# }
# FWER <- 100*(fwer/run)
# 
# point.mean <- matrix(NA,nrow=4,ncol=K)
# point.sd <- matrix(NA,nrow=4,ncol=K)
# for(i in 1:4){
#   point.mean[i,] <- extra.powerprior.sim2[[i]]$`Point Estimate Means`[1:5]
#   point.sd[i,] <- extra.powerprior.sim2[[i]]$`Point Estimate Sds`[1:5]
# }
# point.mean
# point.sd
# 
# 
# extra.powerprior.results2 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)
# 
# save(extra.powerprior.results2,file='PowerPrior_Results2_Extra.RData')


#Sim 3------------------------------------------------
powerprior3 <- list()
listscenarios3 <- list()
for(i in 25:32){
  powerprior3[[i-24]] <- powerprior.sim25[[i]]
  listscenarios3[[i-24]] <- list.scenarios1[[i]]
}

calmat45 <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==4|i==5|i==7|i==8){
    calmat45 <- c(calmat45,powerprior3[[i]]$`Posterior Probabilities`[,5])}
}
cut.off.45 <- quantile(calmat45,0.9)

calmat <- 0
for(i in 1:8){
  if(i==1|i==2|i==3|i==7|i==8){
    calmat <- c(calmat,powerprior3[[i]]$`Posterior Probabilities`[,3])}
}
cut.off.123 <- quantile(calmat,0.9)


cut.off <- c(rep(cut.off.123,3),rep(cut.off.45,2))

ErrorMat <- matrix(NA,ncol=K,nrow=8)
for(i in 1:8){
  post.probs <- powerprior3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,8)
for(j in 1:8){
  true <- as.numeric(listscenarios3[[j]][1:5]>q0)
  post.prob <- powerprior3[[j]]$`Posterior Probabilities`[,1:5]
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
  true <- as.numeric(listscenarios3[[j]][1:5]>q0)
  fwer.true <- which(true==0)
  post.prob <- powerprior3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- powerprior3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- powerprior3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


powerpriorAlpha25.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(powerpriorAlpha25.results3,file='PowerPriorAlpha25_Results3.RData')


####Extra
ErrorMat <- matrix(NA,ncol=K,nrow=4)
for(i in 1:4){
  post.probs <- extra.powerprior.sim3[[i]]$`Posterior Probabilities`
  for(k in 1:5){
    ErrorMat[i,k] <- mean(as.integer(post.probs[,k]>cut.off[k]))*100
  }
}

perfect <- rep(0,4)
for(j in 1:4){
  true <- as.numeric(extra.list.scenarios4[[j]][1:5]>q0)
  post.prob <- extra.powerprior.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  post.prob <- extra.powerprior.sim3[[j]]$`Posterior Probabilities`[,1:5]
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
  point.mean[i,] <- extra.powerprior.sim3[[i]]$`Point Estimate Means`[1:5]
  point.sd[i,] <- extra.powerprior.sim3[[i]]$`Point Estimate Sds`[1:5]
}
point.mean
point.sd


extra.powerprior.results3 <- list('Reject'=ErrorMat,'FWER'=FWER,'Perfect'=Perfect,'Point Estimates'=point.mean,'Point Estimates (sd)'=point.sd)

save(extra.powerprior.results3,file='PowerPrior_Results3_Extra.RData')






