######Implementation of SVM for colon cancer data####
library(kernlab)
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
  model <-ksvm(train,cl.train,type="C-svc", kernel ="rbfdot", cross = 5)
  kas<-predict(model, test)
  prob<- sum(kas == cl.test)/length(cl.test)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SVM for Leukemia cancer data####
library(kernlab)
library(rda)
library(plsgenomics)

data(leukemia)
genenames <- genelist.rda(t(leukemia$X), leukemia$Y, alpha=0.1, delta=0.9)
yy<-as.numeric(genenames)
xdata<-leukemia$X[,yy]
x1a <- xdata[1:27,]
x2a <- xdata[28:38,]
n1<-27
n2<-11
m1<-15  
m2<-7  
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(n1,m1, replace=FALSE)
  S2<-sample(n2,m2, replace=FALSE)
  x1<-x1a[S1,]
  x2<-x2a[S2,]
  x<-rbind(x1,x2)
  y<-c(rep(-1,m1), rep(1,m2))
  miscount<-0
  z1<-x1a[-S1,]
  z2<-x2a[-S2,]
  xtest<-rbind(z1,z2)
  y_test<-c(rep(-1,12),rep(1,4))
  model <-ksvm(x,y,type="C-svc", kernel ="rbfdot", cross = 5)
  kas<-predict(model, xtest)
  prob<- sum(kas == y_test)/length(cl.test)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SVM for lymphoma cancer data####
library(kernlab)
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
  model <-ksvm(train,cl.train,type="C-svc", kernel ="rbfdot", cross = 5)
  kas<-predict(model, test)
  prob<- sum(kas == cl.test)/length(cl.test)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SVM for postrate cancer data####
library(kernlab)
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
  model <-ksvm(train,cl.train,type="C-svc", kernel ="rbfdot", cross = 5)
  kas<-predict(model, test)
  prob<- sum(kas == cl.test)/length(cl.test)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of SVM for SRBCT cancer data####
library(kernlab)
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
for(k in 1:1000){            
  S1<-sample(1:29,m1, replace=FALSE)
  S2<-sample(30:40,m2, replace=FALSE)
  S3<-sample(41:58,m3, replace=FALSE)
  S4<-sample(59:83,m4, replace=FALSE)
  tr.index<-c(S1,S2,S3,S4)
  train<-SRBCT$X[tr.index,]
  test<-SRBCT$X[-tr.index,]
  cl.train<-SRBCT$Y[tr.index]
  cl.test<-SRBCT$Y[-tr.index]
  model <-ksvm(train,cl.train,type="C-svc", kernel ="rbfdot", cross = 5)
  kas<-predict(model, test)
  prob<- sum(kas == cl.test)/length(cl.test)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)