######Implementation of DD1 for colon cancer data####
mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-0
  for(i in 1:n){
    temp<-x-xdata[i,]
    normtemp<-sqrt(sum(temp^2))
    s1<-s1+normtemp
  }
  s1/n
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
  r1<-NULL
  for(j in 1:m1){ r1<-c(r1, mrank(x[j,],x[-j,])) }
  r2<-NULL
  for(j in 1:m2){ r2<-c(r2, mrank(y[j,],y[-j,])) }
  rightcount<-0
  z<-x1a[-S1,]
  miscount<-0
  for(i in 1:m1){
    rx<-mrank(z[i,],x);    ry<-mrank(z[i,],y)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2
    if(Fx<Fy) rightcount<-rightcount+1
  }
  z<-x2a[-S2,]
  for(i in 1:m2){
    rx<-mrank(z[i,],x);    ry<-mrank(z[i,],y)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2
    if(Fx>Fy) rightcount<-rightcount+1
  }
  prob<-rightcount/(m1+m2)
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of DD1 for Leukemia cancer data####
mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-0
  for(i in 1:n){
    temp<-x-xdata[i,]
    normtemp<-sqrt(sum(temp^2))
    s1<-s1+normtemp
  }
  s1/n
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
  r1<-NULL
  for(j in 1:m1){ r1<-c(r1, mrank(x[j,],x[-j,])) }
  r2<-NULL
  for(j in 1:m2){ r2<-c(r2, mrank(y[j,],y[-j,])) }
  rightcount<-0
  z<-x1a[-S1,]
  for(i in 1:12){
    rx<-mrank(z[i,],x);    ry<-mrank(z[i,],y)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2
    if(Fx<Fy) rightcount<-rightcount+1
  }
  z<-x2a[-S2,]
  for(i in 1:4){
    rx<-mrank(z[i,],x);    ry<-mrank(z[i,],y)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2
    if(Fx>Fy) rightcount<-rightcount+1
  }
  prob<-rightcount/16
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of DD1 for lymphoma cancer data####
mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-0
  for(i in 1:n){
    temp<-x-xdata[i,]
    normtemp<-sqrt(sum(temp^2))
    s1<-s1+normtemp
  }
  s1/n
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
  r1<-NULL
  for(j in 1:m1){ r1<-c(r1, mrank(x1[j,],x1[-j,])) }
  r2<-NULL
  for(j in 1:m2){ r2<-c(r2, mrank(x2[j,],x2[-j,])) }
  r3<-NULL
  for(j in 1:m3){ r3<-c(r2, mrank(x3[j,],x3[-j,])) }
  rightcount<-0
  z1<-x1a[-S1,]
  for(i in 1:12){
    rx<-mrank(z1[i,],x1);    ry<-mrank(z1[i,],x2);    rz<-mrank(z1[i,],x3)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2;    Fz<-sum(r3<rz)/m3
    if(Fx<Fy && Fx<Fz) rightcount<-rightcount+1
  }
  z2<-x2a[-S2,]
  for(i in 1:3){
    rx<-mrank(z2[i,],x1);    ry<-mrank(z2[i,],x2);    rz<-mrank(z2[i,],x3)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2;    Fz<-sum(r3<rz)/m3
    if(Fy<Fx && Fy<Fz) rightcount<-rightcount+1
  }
  z3<-x3a[-S3,]
  for(i in 1:4){
    rx<-mrank(z3[i,],x1);    ry<-mrank(z3[i,],x2);    rz<-mrank(z3[i,],x3)
    Fx<-sum(r1<rx)/m1;    Fy<-sum(r2<ry)/m2;    Fz<-sum(r3<rz)/m3
    if(Fz<Fx && Fz<Fy) rightcount<-rightcount+1
  }
  prob<-rightcount/19
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of DD1 for postrate cancer data####
mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-0
  for(i in 1:n){
    temp<-x-xdata[i,]
    normtemp<-sqrt(sum(temp^2))
    s1<-s1+normtemp
  }
  s1/n
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
  x1<-x1a[S1,]
  x2<-x2a[S2,]
  r1<-NULL
  for(j in 1:m){ r1<-c(r1, mrank(x1[j,],x1[-j,])) }
  r2<-NULL
  for(j in 1:m){ r2<-c(r2, mrank(x2[j,],x2[-j,])) }
  rightcount<-0
  z1<-x1a[-S1,]
  for(i in 1:20){
    rx<-mrank(z1[i,],x1);    ry<-mrank(z1[i,],x2)
    Fx<-sum(r1<rx)/m;    Fy<-sum(r2<ry)/m
    if(Fx<Fy) rightcount<-rightcount+1
  }
  z2<-x2a[-S2,]
  for(i in 1:22){
    rx<-mrank(z2[i,],x1);    ry<-mrank(z2[i,],x2)
    Fx<-sum(r1<rx)/m;    Fy<-sum(r2<ry)/m
    if(Fy<Fx) rightcount<-rightcount+1
  }
  prob<-rightcount/42
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)

######Implementation of DD1 for SRBCT cancer data####
mrank<-function(x,xdata){
  n<-dim(xdata)[1]
  d<-dim(xdata)[2]
  s1<-0
  for(i in 1:n){
    temp<-x-xdata[i,]
    normtemp<-sqrt(sum(temp^2))
    s1<-s1+normtemp
  }
  s1/n
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
for(k in 1:1000){            
  S1<-sample(n1,m1, replace=FALSE)
  S2<-sample(n2,m2, replace=FALSE)
  S3<-sample(n3,m3, replace=FALSE)
  S4<-sample(n4,m4, replace=FALSE)
  x1<-x1a[S1,]
  x2<-x2a[S2,]
  x3<-x3a[S3,]
  x4<-x4a[S4,]
  r1<-NULL
  for(j in 1:m1){ r1<-c(r1, mrank(x1[j,],x1[-j,])) }
  r2<-NULL
  for(j in 1:m2){ r2<-c(r2, mrank(x2[j,],x2[-j,])) }
  r3<-NULL
  for(j in 1:m3){ r3<-c(r3, mrank(x3[j,],x3[-j,])) }
  r4<-NULL
  for(j in 1:m4){ r4<-c(r4, mrank(x4[j,],x4[-j,])) }
  rightcount<-0
  z1<-x1a[-S1,]
  miscount<-0
  for(i in 1:14){
    rx<-mrank(z1[i,],x1)
    ry<-mrank(z1[i,],x2)
    rz<-mrank(z1[i,],x3)
    rw<-mrank(z1[i,],x4)
    Fx<-sum(r1<rx)/m1
    Fy<-sum(r2<ry)/m2
    Fz<-sum(r3<rz)/m3
    Fw<-sum(r4<rw)/m4
    if(Fx<Fy && Fx<Fz && Fx<Fw) rightcount<-rightcount+1
  }
  z2<-x2a[-S2,]
  for(i in 1:4){
    rx<-mrank(z2[i,],x1)
    ry<-mrank(z2[i,],x2)
    rz<-mrank(z2[i,],x3)
    rw<-mrank(z2[i,],x4)
    Fx<-sum(r1<rx)/m1
    Fy<-sum(r2<ry)/m2
    Fz<-sum(r3<rz)/m3
    Fw<-sum(r4<rw)/m4
    if(Fy<Fx && Fy<Fz && Fy<Fw) rightcount<-rightcount+1
  }
  z3<-x3a[-S3,]
  for(i in 1:9){
    rx<-mrank(z3[i,],x1)
    ry<-mrank(z3[i,],x2)
    rz<-mrank(z3[i,],x3)
    rw<-mrank(z3[i,],x4)
    Fx<-sum(r1<rx)/m1
    Fy<-sum(r2<ry)/m2
    Fz<-sum(r3<rz)/m3
    Fw<-sum(r4<rw)/m4
    if(Fz<Fx && Fz<Fy && Fz<Fw) rightcount<-rightcount+1
  }
  z4<-x4a[-S4,]
  for(i in 1:10){
    rx<-mrank(z4[i,],x1)
    ry<-mrank(z4[i,],x2)
    rz<-mrank(z4[i,],x3)
    rw<-mrank(z4[i,],x4)
    Fx<-sum(r1<rx)/m1
    Fy<-sum(r2<ry)/m2
    Fz<-sum(r3<rz)/m3
    Fw<-sum(r4<rw)/m4
    if(Fw<Fx && Fw<Fy && Fw<Fz) rightcount<-rightcount+1
  }
  prob<-rightcount/37
  corprob<-c(corprob,prob)
}
quantile(corprob); mean(corprob); sd(corprob)