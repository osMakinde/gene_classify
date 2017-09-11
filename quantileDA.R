######Implementation of QuantileDA for colon cancer data####
library(quantileDA)
library(spls)
library(rda)
library(plsgenomics)

data(colon)
colon.x1 <- t(colon.x)
genenames <- genelist.rda(colon.x1, colon.y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
colon.x<-t(colon.x1[yy,])
n1<-22
n2<-40
m1<-11  
m2<-20   
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(1:22,m1, replace=FALSE)
  S2<-sample(23:62,m2, replace=FALSE)
  tr.index<-c(S1,S2)
  train<-colon.x[tr.index,]
  cl.train<-colon.y[tr.index]
  test<-colon.x[-tr.index,]
  cl.test<-colon.y[-tr.index]
  out.q=quantilecldiff(train,test,cl.train,cl.test=cl.test)
  prob<-out.q$me.test
  corprob<-c(corprob,1-prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of QuantileDA for Leukemia cancer data####
library(quantileDA)
library(spls)
library(rda)
library(plsgenomics)

data(leukemia)
genenames <- genelist.rda(t(leukemia$X), leukemia$Y, alpha=0.1, delta=0.9)
yy<-as.numeric(genenames)
leukemia$X<-leukemia$X[,yy]
n1<-27
n2<-11
m1<-15  
m2<-7  
corprob<-NULL
for(k in 1:1000){                
  S1<-sample(1:27,m1, replace=FALSE)
  S2<-sample(28:38,m2, replace=FALSE)
  tr.index<-c(S1,S2)
  train<-leukemia$X[tr.index,]
  cl.train<-leukemia$Y[tr.index]
  test<-leukemia$X[-tr.index,]
  cl.test<-leukemia$Y[-tr.index]
  out.q=quantilecldiff(train,test,cl.train,cl.test=cl.test)
  prob<-out.q$me.test
  corprob<-c(corprob,1-prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of QuantileDA for lymphoma cancer data####
library(quantileDA)
library(spls)
library(rda)
data(lymphoma)
genenames <- genelist.rda(t(lymphoma$x), lymphoma$y, alpha=0.1, delta=1.5)
yy<-as.numeric(genenames)
xdata<-lymphoma$x[,yy]
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
  train<-lymphoma$x[tr.index,]
  cl.train<-lymphoma$y[tr.index]
  test<-lymphoma$x[-tr.index,]
  cl.test<-lymphoma$y[-tr.index]
  out.q=quantilecldiff(train,test,cl.train,cl.test=cl.test)
  prob<-out.q$me.test
  corprob<-c(corprob,1-prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of QuantileDA for postrate cancer data####
library(quantileDA)
library(spls)
library(rda)
library(plsgenomics)

data(prostate)
genenames <- genelist.rda(t(prostate$x), prostate$y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
prostate$x<-prostate$x[,yy]
m<-30   
corprob<-NULL
for(k in 1:1000){                
  S1<-sample(1:50,m, replace=FALSE)
  S2<-sample(51:102,m, replace=FALSE)
  tr.index<-c(S1,S2)
  train<-prostate$x[tr.index,]
  cl.train<-prostate$y[tr.index]
  test<-prostate$x[-tr.index,]
  cl.test<-prostate$y[-tr.index]
  out.q=quantilecldiff(train,test,cl.train,cl.test=cl.test)
  prob<-out.q$me.test
  corprob<-c(corprob,1-prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of QuantileDA for SRBCT cancer data####
library(quantileDA)
library(spls)
library(rda)
library(plsgenomics)

data(SRBCT)
SRBCT$X1 <- t(SRBCT$X)
genenames <- genelist.rda(SRBCT$X1, SRBCT$Y, alpha=0.1, delta=0.9)
yy<-as.numeric(genenames)
SRBCT$X0<-t(SRBCT$X1[yy,])
xdata0<-cbind(SRBCT$X0,SRBCT$Y)
SRBCT$X<-xdata0[order(SRBCT$Y),]              
n1<-29
n2<-11
n3<-18
n4<-25
m1<-15  
m2<-7 
m3<-9  
m4<-15  
corprob<-NULL
for(k in 1:294){            
  S1<-sample(1:29,m1, replace=FALSE)
  S2<-sample(30:40,m2, replace=FALSE)
  S3<-sample(41:58,m3, replace=FALSE)
  S4<-sample(59:83,m4, replace=FALSE)
  tr.index<-c(S1,S2,S3,S4)
  train<-SRBCT$X[tr.index,]
  cl.train<-SRBCT$Y[tr.index]
  test<-SRBCT$X[-tr.index,]
  cl.test<-SRBCT$Y[-tr.index]
  out.q=quantilecldiff(train,test,cl.train,cl.test=cl.test)
  prob<-out.q$me.test
  corprob<-c(corprob,1-prob)
}
quantile(corprob); mean(corprob); sd(corprob)