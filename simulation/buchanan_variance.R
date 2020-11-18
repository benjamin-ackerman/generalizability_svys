ipsw_var <- function(data,selvar,pscore,svywt, trt,outcome) {
  cohort=data[which(data$S == 0),]
  m <- nrow(cohort)
  
  trial=data[which(data$S == 1),]
  n <- nrow(trial)
  
  b <- nrow(data)
  both <- data
  
  #define weights for selection propensity score model
  # trial$pw<-1
  # cohort$pw<-m/(nt-n)
  #combine trial and cohort (S,Z) to estimate propensity scores
  #both<-rbind(cohort,trial)
  #b<-dim(both)[1]
  #estimate selection propensity scores using logistic regression
  #using data from cohort and trial stacked
  # mylogit <- glm(s ~ selvar, data = both, family = "binomial",weights=(both$pw)^(-1))
  # both$p <- predict(mylogit,type = "response")
  # trial$p<-both$p[which(both$s==1)]
  #Compute the IPSW estimator when S=1
  #IPSW estimator: denominator (using inverse odds weights)
  w1<-sum(trial[,trt]*(1-trial[,pscore])/trial[,pscore])
  w2<-sum((1-trial[,trt]*(1-trial[,pscore]))/trial[,pscore])
  #get mu(1)
  trial$mu1i<- (trial[,trt]*trial[,outcome]*(1-trial[,pscore]))/trial[,pscore]
  mu1hat <- (w1)^(-1)*sum(trial$mu1i, na.rm=TRUE)
  #get mu(0)
  trial$mu0i <- ((1-trial[,trt])*trial[,outcome]*(1-trial[,pscore]))/trial[,pscore]
  mu0hat<- (w2)^(-1)*sum(trial$mu0i, na.rm=TRUE)
  PATEhat=mu1hat-mu0hat
  #Estimated Variance for PATE when weights estimated
  #Use sandwich estimator (Derived in Appendix A)
  #estimate of derivative of weights
  #first derivatives
  #wrt beta1,beta10
  #vector of covariates
  
  both$int=rep(1,b)
  trial$int=rep(1,n)
  
  covb=as.matrix(cbind(both$int,both[,selvar]))
  covt<-as.matrix(cbind(trial$int,trial[,selvar]))
  
  # dmy <- caret::dummyVars(" ~ .", data = both[,selvar], sep="_")
  # expanded.data = data.frame(int = both$int,
  #                            predict(dmy, newdata = both[,selvar]))
  # names(expanded.data)=stringr::str_replace_all(names(expanded.data),"\\.","-")
  # covb <- as.matrix(expanded.data)
  # 
  # dmy <- caret::dummyVars(" ~ .", data = trial[,selvar], sep="_")
  # expanded.data = data.frame(int = trial$int,
  #                            predict(dmy, newdata = trial[,selvar]))
  # names(expanded.data)=stringr::str_replace_all(names(expanded.data),"\\.","-")
  # covt <- as.matrix(expanded.data)
  
  both$wbeta <- covb*both[,pscore]*(1-both[,pscore])
  trial$wbeta <- covt*trial[,pscore]*(1-trial[,pscore])
  size=dim(covb)[2]
  size2=size+2
  bwbeta2<-list()
  twbeta2<-list()
  for(i in 1:size){
    bwbeta2[[i]]<-covb*(both$wbeta[,i]-2*both$wbeta[,i]*both[,pscore])
    twbeta2[[i]]<-covt*(trial$wbeta[,i]-2*trial$wbeta[,i]*trial[,pscore])
  }
  Aw <- matrix(rep(NA,size*size),size,size)
  for(i in 1:size){
    for(k in 1:size){
      Aw[i,k] <- sum(both[,svywt]^(-1)*(-both$wbeta[,i]*both$wbeta[,k]*(both[,pscore]*(1-both[,pscore])
                                                                   + (both$S-both[,pscore])*(1-2*both[,pscore]))/(both[,pscore]^2*(1-both[,pscore])^2)))
      + sum(both[,svywt]^(-1)*(bwbeta2[[k]][,i]*(both$S-both[,pscore]))/(both[,pscore]*(1-both[,pscore])))
    }
  }
  wAhat <- matrix(rep(NA,size2*size2),size2,size2)
  wAhat[1,1]<-(b)^(-1)*sum((trial[,trt]/trial[,pscore]))
  wAhat[1,2]<-0
  wAhat[2,1]<-0
  wAhat[2,2]<-(b)^(-1)*sum(((1-trial[,trt]))/(trial[,pscore]))
  for(k in 1:size){
    wAhat[1,k+2] <-(b)^(-1)*sum((trial[,trt]*(trial[,outcome]-mu1hat)*trial$wbeta[,k])/trial[,pscore]^2, na.rm=TRUE)
  }
  for(k in 1:size){
    wAhat[2,k+2] <-(b)^(-1)*sum(((1-trial[,trt])*(trial[,outcome]-mu0hat)*trial$wbeta[,k])/trial[,pscore]^2, na.rm=TRUE)
  }
  for(k in 3:size2){
    wAhat[k,1] <-0
    wAhat[k,2] <-0
  }
  for(i in 3:size2){
    for(k in 3:size2){
      wAhat[i,k]<-(b)^(-1)*(-Aw[i-2,k-2])
    }
  }
  wBhat <- matrix(rep(NA,size2*size2),size2,size2)
  wBhat[1,1]<-(b)^(-1)*sum((trial[,trt]*(trial[,outcome]-mu1hat)^2)/(trial[,pscore])^2, na.rm=TRUE)
  wBhat[1,2]<-wBhat[2,1]<-(b)^(-1)*sum((trial[,trt]*(trial[,outcome]-mu1hat))/trial[,pscore]*sum(((1-trial[,trt])*(trial[,outcome]-mu0hat))/trial[,pscore], na.rm=TRUE), na.rm=TRUE)
  wBhat[2,2]<-(b)^(-1)*sum(((1-trial[,trt])*(trial[,outcome]-mu0hat)^2)/(trial[,pscore])^2,na.rm=TRUE)
  for(k in 1:size){
    wBhat[1,k+2]<-(b)^(-1)*sum((trial[,trt]*(trial[,outcome]-mu1hat)*trial$wbeta[,k])/(trial[,pscore])^2, na.rm=TRUE)
    wBhat[2,k+2]<-(b)^(-1)*sum(((1-trial[,trt])*(trial[,outcome]-mu0hat)*trial$wbeta[,k])/(trial[,pscore])^2,na.rm=TRUE)
    wBhat[k+2,1]<-wBhat[1,k+2]
    wBhat[k+2,2]<-wBhat[2,k+2]
  }
  for(i in 1:size){
    for(k in 1:size){
      wBhat[i+2,k+2]<-(b)^(-1)*sum((((both$S-both[,pscore])^2)/(both[,pscore]^2*(1-both[,pscore])^2))
                                   *both$wbeta[,i]*both$wbeta[,k]*both[,svywt]^(-2), na.rm=TRUE)
    }
  }
  wsigmaM<-solve(wAhat)%*%wBhat%*%t(solve(wAhat))
  wsigma<-wsigmaM[1,1]+wsigmaM[2,2]-2*wsigmaM[2,1]
  wse<-sqrt(wsigma)/sqrt(b)
  PATE_LCL <-PATEhat-1.96*wse
  PATE_UCL <-PATEhat+1.96*wse
  #result <-c(PATEhat, mu1hat,mu0hat,wsigma,wse, PATE_LCL, PATE_UCL)
  #names(result)<-c("PATEhat", "mu1hat","mu0hat","wsigma","wse", "PATE_LCL", "PATE_UCL")
  result <- list(PATE2 = PATEhat,
                 se2 = wse)
  
  return(result)
}
