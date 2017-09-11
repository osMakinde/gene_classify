######Implementation of D3 for colon cancer data####
mrank<-function(x,xdata){
  temp<-x-xdata
  normtemp<-sum(abs(temp))
}
library(rda)
library(plsgenomics)

data(colon)
colon.x1 <- t(colon.x)
genenames <- genelist.rda(colon.x1, colon.y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
xdata<-t(colon.x1[yy,])
x1a <- xdata[1:22,]
x2a <- xdata[23:62,]
n1<-22
n2<-40
m1<-11  
m2<-20   
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(n1,m1, replace=FALSE)
  S2<-sample(n2,m2, replace=FALSE)
  x<-x1a[S1,]
  y<-x2a[S2,]
  xm<-apply(x,2,median)
  ym<-apply(y,2,median)
  rightcount<-0
  z<-x1a[-S1,]
  miscount<-0
  for(i in 1:m1){
    rx<-mrank(z[i,],xm);     ry<-mrank(z[i,],ym)
    if(rx<ry) rightcount<-rightcount+1
  }
  z<-x2a[-S2,]
  for(i in 1:m2){
    rx<-mrank(z[i,],xm);     ry<-mrank(z[i,],ym)
    if(rx>ry) rightcount<-rightcount+1
  }
  prob<-rightcount/(m1+m2)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of D3 for Leukemia cancer data####

mrank<-function(x,xdata){
  temp<-x-xdata
  normtemp<-sum(abs(temp))
}


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
  x<-x1a[S1,]
  y<-x2a[S2,]
  xm<-apply(x, 2, median)
  ym<-apply(y, 2, median)
  rightcount<-0
  z<-x1a[-S1,]
  miscount<-0
  for(i in 1:12){
    rx<-mrank(z[i,],xm);    ry<-mrank(z[i,],ym)
    if(rx<ry) rightcount<-rightcount+1
  }
  z<-x2a[-S2,]
  for(i in 1:4){
    rx<-mrank(z[i,],xm);    ry<-mrank(z[i,],ym)
    if(rx>ry) rightcount<-rightcount+1
  }
  prob<-rightcount/16
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of D3 for lymphoma cancer data####
mrank<-function(x,xdata){
  temp<-x-xdata
  normtemp<-sum(abs(temp))
}

library(spls)
library(rda)
data(lymphoma)
genenames <- genelist.rda(t(lymphoma$x), lymphoma$y, alpha=0.1, delta=1.5)
yy<-as.numeric(genenames)
xdata<-lymphoma$x[,yy]
y<-lymphoma$y
n0<-dim(xdata)[2]
x1a <- xdata[1:42,]
x2a <- xdata[43:51,]
x3a <- xdata[52:62,]
n1<-42
n2<-9
n3<-11
m1<-30 
m2<-6
m3<-7
corprob<-NULL
for(k in 1:1000){            
  S1<-sample(n1,m1, replace=FALSE)
  S2<-sample(n2,m2, replace=FALSE)
  S3<-sample(n3,m3, replace=FALSE)
  x1<-x1a[S1,]
  x2<-x2a[S2,]
  x3<-x3a[S3,]
  xm<-apply(x1,2,median)
  ym<-apply(x2,2,median)
  zm<-apply(x3,2,median)
  rightcount<-0
  z1<-x1a[-S1,]
  for(i in 1:12){
    rx<-mrank(z1[i,],xm);    ry<-mrank(z1[i,],ym);    rz<-mrank(z1[i,],zm)
    if(rx<ry && rx<rz) rightcount<-rightcount+1
  }
  z2<-x2a[-S2,]
  for(i in 1:3){
    rx<-mrank(z2[i,],xm);    ry<-mrank(z2[i,],ym);    rz<-mrank(z2[i,],zm)
    if(ry<rx && ry<rz) rightcount<-rightcount+1
  }
  z3<-x3a[-S3,]
  for(i in 1:4){
    rx<-mrank(z3[i,],xm);    ry<-mrank(z3[i,],ym);    rz<-mrank(z3[i,],zm)
    if(rz<rx && rz<ry) rightcount<-rightcount+1
  }
  prob<-rightcount/19
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of D3 for postrate cancer data####
mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  temp<-x-xdata
  normtemp<-sum(abs(temp))
}

library(spls)
library(rda)
data(prostate)
genenames <- genelist.rda(t(prostate$x), prostate$y, alpha=0.2, delta=0.3)
yy<-as.numeric(genenames)
xdata<-prostate$x[,yy]
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
  x1<-x1a[S1,];  x2<-x2a[S2,]
  xm<-apply(x1, 2, median);  ym<-apply(x2, 2, median)
  rightcount<-0
  z1<-x1a[-S1,]
  for(i in 1:20){
    rx<-mrank(z1[i,],xm);    ry<-mrank(z1[i,],ym)
    if(rx<ry) rightcount<-rightcount+1
  }
  z2<-x2a[-S2,]
  for(i in 1:22){
    rx<-mrank(z2[i,],xm);    ry<-mrank(z2[i,],ym)
    if(ry<rx) rightcount<-rightcount+1
  }
  prob<-rightcount/42
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of D3 for SRBCT cancer data####
mrank<-function(x,xdata){
  temp<-x-xdata
  normtemp<-sum(abs(temp))
}


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
xdata<-xdata1[,1:(n0-1)]
x1a <- xdata[1:29,]
x2a <- xdata[30:40,]
x3a <- xdata[41:58,]
x4a <- xdata[59:83,]
n1<-29
n2<-11
n3<-18
n4<-25
m1<-15  
m2<-7 
m3<-9  
m4<-15 
corprob<-NULL
for(k in 1:100){            
  S1<-sample(n1,m1, replace=FALSE)
  S2<-sample(n2,m2, replace=FALSE)
  S3<-sample(n3,m3, replace=FALSE)
  S4<-sample(n4,m4, replace=FALSE)
  x1<-x1a[S1,]
  x2<-x2a[S2,]
  x3<-x3a[S3,]
  x4<-x4a[S4,]
  xm<-apply(x1, 2, median)
  ym<-apply(x2, 2, median)
  zm<-apply(x3, 2, median)
  wm<-apply(x4, 2, median)
  rightcount<-0
  z1<-x1a[-S1,]
  for(i in 1:14){
    rx<-mrank(z1[i,],xm)
    ry<-mrank(z1[i,],ym)
    rz<-mrank(z1[i,],zm)
    rw<-mrank(z1[i,],wm)
    if(rx<ry && rx<rz && rx<rw) rightcount<-rightcount+1
  }
  z2<-x2a[-S2,]
  for(i in 1:4){
    rx<-mrank(z2[i,],xm)
    ry<-mrank(z2[i,],ym)
    rz<-mrank(z2[i,],zm)
    rw<-mrank(z2[i,],wm)
    if(ry<rx && ry<rz && ry<rw) rightcount<-rightcount+1
  }
  z3<-x3a[-S3,]
  for(i in 1:9){
    rx<-mrank(z3[i,],xm)
    ry<-mrank(z3[i,],ym)
    rz<-mrank(z3[i,],zm)
    rw<-mrank(z3[i,],wm)
    if(rz<rx && rz<ry && rz<rw) rightcount<-rightcount+1
  }
  z4<-x4a[-S4,]
  for(i in 1:10){
    rx<-mrank(z4[i,],xm)
    ry<-mrank(z4[i,],ym)
    rz<-mrank(z4[i,],zm)
    rw<-mrank(z4[i,],wm)
    if(rw<rx && rw<ry && rw<rz) rightcount<-rightcount+1
  }
  prob<-rightcount/37
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)