######Implementation of splsda for colon cancer data####
library(rda)
library(spls)
library(plsgenomics)
data(colon)
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
  train_label<-colon.y[tr.index]
  test<-colon.x[-tr.index,]
  test_label<-colon.y[-tr.index]
  f <- splsda(train, train_label, K=3, eta=0.8, scale.x=FALSE )
  (pred.f<-predict(f, test, type="fit" ))
  prob<-sum(pred.f == test_label)/(m1+m2)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of splsda for Leukemia cancer data####
library(spls)
library(plsgenomics)

data(leukemia)
xdata<-leukemia$X
x1a <- xdata[1:27,]
x2a <- xdata[28:38,]
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
  train_label<-leukemia$Y[tr.index]
  test<-leukemia$X[-tr.index,]
  test_label<-leukemia$Y[-tr.index]
  f <- splsda(train, train_label, K=3, eta=0.8, scale.x=FALSE )
  (pred.f<-predict(f, test, type="fit" ))
  prob<-sum(pred.f == test_label)/16
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of splsda for lymphoma cancer data####
library(spls)
data(lymphoma)
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
  train_label<-lymphoma$y[tr.index]
  test<-lymphoma$x[-tr.index,]
  test_label<-lymphoma$y[-tr.index]
  f <- splsda(train, train_label, K=3, eta=0.8, scale.x=FALSE )
  (pred.f<-predict(f, test, type="fit" ))
  prob<-sum(pred.f == test_label)/19
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of splsda for postrate cancer data####
library(spls)
data(prostate)
xdata<-prostate$x
y<-prostate$y
n0<-dim(xdata)[2]
x1a <- xdata[1:50,]
x2a <- xdata[51:102,]
n1<-50
n2<-52
m<-30  
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(n1,m, replace=FALSE)
  S2<-sample(n2,m, replace=FALSE)
  x1<-x1a[S1,]
  x2<-x2a[S2,] 
  train<-rbind(x1,x2)
  train_label<-c(rep(0,m),rep(1,m))
  test<-rbind(x1a[-S1,],x2a[-S1,])
  test_label<-c(rep(0,20),rep(1,22))  
  f <- splsda(train, train_label, K=3, eta=0.8, scale.x=FALSE )
  (pred.f<-predict(f, test, type="fit" ))
  prob<-sum(pred.f == test_label)/42
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of splsda for SRBCT cancer data####
library(spls)
library(plsgenomics)

data(SRBCT)
xdata0<-cbind(SRBCT$X,SRBCT$Y)
xdata1<-xdata0[order(SRBCT$Y),]
n0<-dim(xdata1)[2]
SRBCT$X<-xdata1[,1:(n0-1)]
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
  train<-SRBCT$X[tr.index,]
  train_label<-SRBCT$Y[tr.index]
  test<-SRBCT$X[-tr.index,]
  test_label<-SRBCT$Y[-tr.index]
  f <- splsda(train, train_label, K=3, eta=0.8, scale.x=FALSE )
  (pred.f<-predict(f, test, type="fit" ))
  prob<-sum(pred.f == test_label)/37
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)