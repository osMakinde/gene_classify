######Implementation of NSC for colon cancer data####
library(pamr)
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
  ## calculate the prediction accuracy
  mydata <- list(x=colon.x[, tr.index],y=colon.y[tr.index])
  mytrain <- pamr.train(mydata)
  ls<-pamr.predict(mytrain, colon.x[, -tr.index], threshold=0.5)
  prob <- sum(ls==colon.y[-tr.index])/length(colon.y[-tr.index])
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of NSC for Leukemia cancer data####
library(rda)
library(plsgenomics)
library(pamr)
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
  ## calculate the prediction accuracy
  mydata <- list(x=leukemia$X[, tr.index],y=leukemia$Y[tr.index])
  mytrain <- pamr.train(mydata)
  ls<-pamr.predict(mytrain, leukemia$X[, -tr.index], threshold=0.5)
  prob <- sum(ls==leukemia$Y[-tr.index])/length(leukemia$Y[-tr.index])
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of NSC for lymphoma cancer data####
library(rda)
library(pamr)
library(spls)
data(lymphoma)
lymphoma$x0<-t(lymphoma$x)
genenames <- genelist.rda(lymphoma$x0, lymphoma$y, alpha=0.1, delta=1.5)
yy<-as.numeric(genenames)
lymphoma$x<-lymphoma$x0[yy,]
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
  ## calculate the prediction accuracy
  mydata <- list(x=lymphoma$x[, tr.index],y=lymphoma$y[tr.index])
  mytrain <- pamr.train(mydata)
  ls<-pamr.predict(mytrain, lymphoma$x[, -tr.index], threshold=0.5)
  prob <- sum(ls==lymphoma$y[-tr.index])/length(lymphoma$y[-tr.index])
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of NSC for postrate cancer data####
library(pamr)
library(rda)
library(spls)
data(prostate)
prostate$x1 <- t(prostate$x)
genenames <- genelist.rda(prostate$x1, prostate$y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
prostate$x<-prostate$x1[yy,]
n1<-50
n2<-52
m<-30  
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(1:50,m, replace=FALSE)
  S2<-sample(51:102,m, replace=FALSE)
  tr.index<-c(S1,S2)
  ## calculate the prediction accuracy
  mydata <- list(x=prostate$x[,tr.index],y=prostate$y[tr.index])
  mytrain <- pamr.train(mydata)
  ls<-pamr.predict(mytrain, prostate$x[, -tr.index], threshold=0.5)
  prob <- sum(ls==prostate$y[-tr.index])/length(prostate$y[-tr.index])
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of NSC for SRBCT cancer data####
library(rda)
library(plsgenomics)
library(pamr)
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
  ## calculate the prediction accuracy
  mydata <- list(x=SRBCT$X[, tr.index],y=SRBCT$Y[tr.index])
  mytrain <- pamr.train(mydata)
  ls<-pamr.predict(mytrain, SRBCT$X[, -tr.index], threshold=0.5)
  prob <- sum(ls==SRBCT$Y[-tr.index])/length(SRBCT$Y[-tr.index])
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)