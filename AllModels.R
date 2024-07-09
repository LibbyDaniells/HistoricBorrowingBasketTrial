library(rjags)
library(textmineR)
library(matrixStats)
library(RColorBrewer)
library(philentropy)
library(LaplacesDemon)

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

#EXppNEX--------------------------------------------------
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

#EXsamNEX--------------------------------------------------
SAMPrior <- function(K,y,n,q0,yh,nh,a,b,clin.diff,pi){
  hist <- as.numeric(nh>0)
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
  jags.fit.sam <- jags.model(file='SAMPriorinNEX.txt',data=jags.data.sam,n.adapt=100000,n.chains=4)
  samples.sam <- coda.samples(jags.fit.sam,variable.names=c('p','DeltaC','Epsilon'),n.iter=100000,silent=T)
  sam <- as.data.frame(samples.sam[[1]])
  return(sam)
}


#histFujikawas---------------------------------------------------
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



#mEXNEXhist--------------------------------------------------
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



#MLMixture-------------------------------------------------------
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





#Data Example------------------------------------------------------------
K <- 4
H <- 2
y <- c(3,3,8,3)
yh <- c(1,1,1,0)
n <- rep(34,K)
nh <- c(13,13,13,0)
hist <- as.numeric(nh>0)
q0 <- 0.2
a <- 1
b <- 1
clin.diff <- 0.2

#For EXNEX
pw <- 0.3
pi.exnex <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
pi.exnex.all <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))

#For EXppNEX
alpha0 <- 0.5
a <- 1
b <- 1

#For histFujikawa's
tau <- 0.2
hist.scalar <- 0.8
epsilon <- 2


#For MLMixture
pi.delta <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
pi.epsilon1 <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))
pi.epsilon2 <- rbind(c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5),c(0.5,0.5))


bhm.current <- BHM(K,y,n,q0)
bhm.all <- BHM(K+H,c(y,yh[which(nh!=0)]),c(n,nh[which(nh!=0)]),q0)
exnex.current <- EXNEX(K,y,n,q0,pw,pi.exnex)
exnex.all <- EXNEX(K+H,c(y,yh[which(nh!=0)]),c(n,nh[which(nh!=0)]),q0,pw,pi.exnex.all)
power.prior <- PowerPrior(K,y,n,q0,pi.exnex,pw,yh,nh,alpha0,a,b)
samprior <- SAMPrior(K,y,n,q0,yh,nh,a,b,clin.diff,pi.exnex)
fujikawa <- adapt.fujikawa(y,n,yh,nh,tau,hist.scalar,epsilon)$Parameters           
mexnex <- mEXNEXhell(K,y,yh,n,nh,q0,pw,H,hist.scalar,pi.exnex)
mmixture <- MLMixture(K,y,yh,n,nh,q0,H,pi.delta,pi.epsilon1,pi.epsilon2,a,b)


