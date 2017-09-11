######Implementation of SCRDA for colon cancer data####
library(rda)
library(plsgenomics)

data(colon)
colon.x <- t(colon.x)
genenames <- genelist.rda(colon.x, colon.y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
colon.x<-colon.x[yy,]
n1<-22
n2<-40
m1<-11  
m2<-20  
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(1:22,m1, replace=FALSE)
  S2<-sample(23:62,m2, replace=FALSE)
  tr.index<-c(S1,S2)
  fit <- rda(colon.x[, tr.index], colon.y[tr.index])
  ynew <- predict(fit, x=colon.x[, tr.index], y=colon.y[tr.index],xnew=colon.x[, -tr.index], alpha=0.1, delta=0.5)
  ## calculate the prediction error
  prob<-sum(ynew == colon.y[-tr.index])/(m1+m2)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SCRDA for Leukemia cancer data####
library(rda)
library(plsgenomics)

data(leukemia)
genenames <- genelist.rda(t(leukemia$X), leukemia$Y, alpha=0.1, delta=0.9)
yy<-as.numeric(genenames)
leukemia$X<-t(leukemia$X[,yy])
n1<-27
n2<-11
m1<-15  
m2<-7  
corprob<-NULL
for(k in 1:1000){                       
  S1<-sample(1:27,m1, replace=FALSE)
  S2<-sample(28:38,m2, replace=FALSE)
  tr.index<-c(S1,S2)
  fit <- rda(leukemia$X[, tr.index], leukemia$Y[tr.index])
  ynew <- predict(fit, x=leukemia$X[, tr.index], y=leukemia$Y[tr.index],xnew=leukemia$X[, -tr.index], alpha=0.1, delta=0.5)
  ## calculate the prediction error
  prob<-sum(ynew == leukemia$Y[-tr.index])/16
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SCRDA for lymphoma cancer data####
library(rda)
library(spls)

data(lymphoma)
genenames <- genelist.rda(t(lymphoma$x), lymphoma$y, alpha=0.1, delta=1.5)
yy<-as.numeric(genenames)
lymphoma$x<-t(lymphoma$x[,yy])
y<-lymphoma$y
n1<-42
n2<-9
n3<-11
m1<-30 
m2<-6
m3<-7
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(1:42,m1, replace=FALSE)
  S2<-sample(43:51,m2, replace=FALSE)
  S3<-sample(52:62,m3, replace=FALSE)
  tr.index<-c(S1,S2,S3)
  fit <- rda(lymphoma$x[, tr.index], lymphoma$y[tr.index])
  ynew <- predict(fit, x=lymphoma$x[, tr.index], y=lymphoma$y[tr.index],xnew=lymphoma$x[, -tr.index], alpha=0.2, delta=0.5)
  ## calculate the prediction error
  prob<-sum((ynew-1) == lymphoma$y[-tr.index])/19
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SCRDA for postrate cancer data####
library(rda)
library(spls)

data(prostate)
genenames <- genelist.rda(t(prostate$x), prostate$y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
prostate$x<-t(prostate$x[,yy])
prostate$y<-prostate$y
n1<-50
n2<-52
m<-30  
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(1:50,m, replace=FALSE)
  S2<-sample(51:102,m, replace=FALSE)
  tr.index<-c(S1,S2)
  fit <- rda(prostate$x[, tr.index], prostate$y[tr.index])
  ynew <- predict(fit, x=prostate$x[, tr.index], y=prostate$y[tr.index],xnew=prostate$x[, -tr.index], alpha=0.2, delta=0.5)
  ## calculate the prediction error
  prob<-sum((ynew-1) == prostate$y[-tr.index])/42
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SCRDA for SRBCT cancer data####
library(rda)
library(plsgenomics)

data(SRBCT)
SRBCT$X1 <- t(SRBCT$X)
genenames <- genelist.rda(SRBCT$X1, SRBCT$Y, alpha=0.1, delta=0.9)
yy<-as.numeric(genenames)
SRBCT$X0<-t(SRBCT$X1[yy,])
xdata0<-cbind(SRBCT$X0,SRBCT$Y)
xdata1<-xdata0[order(SRBCT$Y),]
n0<-dim(xdata1)[2]
SRBCT$X<-t(xdata1[,1:(n0-1)])
SRBCT$Y<-xdata1[,n0]
n1<-29
n2<-11
n3<-18
n4<-25
m1<-15  
m2<-7 
m3<-9  
m4<-15  
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(1:29,m1, replace=FALSE)
  S2<-sample(30:40,m2, replace=FALSE)
  S3<-sample(41:58,m3, replace=FALSE)
  S4<-sample(59:83,m4, replace=FALSE)
  tr.index<-c(S1,S2,S3,S4)
  fit <- rda(SRBCT$X[, tr.index], SRBCT$Y[tr.index])
  ynew <- predict(fit, x=SRBCT$X[, tr.index], y=SRBCT$Y[tr.index],xnew=SRBCT$X[, -tr.index], alpha=0.1, delta=0.9)
  ## calculate the prediction error
  prob<-sum(ynew == SRBCT$Y[-tr.index])/37
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)