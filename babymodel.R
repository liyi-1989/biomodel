library(synchrony)
data(pisco.data)
# Target Variable(log of species number)
y=subset(pisco.data,year>2000,select=c(mussel_abund))
y[y[,1]==0,]=min(y[y[,1]!=0,])
y=log(y)
# (past) species number
N=subset(pisco.data,year<2003,select=c(mussel_abund))
# site distance
D=coord2dist(pisco.data[1:48,1:2],lower.tri = F)
# Environment variable
E=subset(pisco.data,year<2003,select=c(chl,sst,upwelling))

# negative log-likelihood
logl=function(t,y,N,D,E){
  n=dim(y)[1]
  f=rep(0,n)
  yr=rep(0:2,each=48)
  for(j in 1:n){
    Ker=exp(-(D[j%%48,]-t[2])^2/(t[3]^2)) # dispersal kernal 
    Nj=N[(yr[j]*48+1):((yr[j]+1)*48),]
    f[j]=t[1]*sum(Ker*Nj)+t[4]*E[j,1]+t[5]*E[j,2]+t[6]*E[j,3]
  }
  
  return(sum((y[,1]-f)^2))
}

#logl(c(1,1,1,1,1,1),y,N,D,E)
t0=c(8,-0.2,1,0.1,0.1,0.001)
#t0=c(3.560124968,  0.052076860  , 0.128648187  , 0.083371125  , 0.026199117  , 0.006296133)
res=nlm(logl,t0,hessian=T,print.level=2,y=y,N=N,D=D,E=E,iterlim=1e4,steptol=1e-5)

t=res$estimate

n=dim(y)[1]
f=rep(0,n)
yr=rep(0:2,each=48)
for(j in 1:n){
  f[j]=t[1]*sum(exp(-(D[j%%48,]-t[2])^2/(t[3]^2))*N[(yr[j]*48+1):((yr[j]+1)*48),])+sum(t[4:6]*E[j,])
}

mean((y[,1]-f)^2)

plot(y[,1],f)

mean((exp(y[,1])-exp(f))^2)
plot(exp(y[,1]),exp(f))

#=============================================================
# negative log-likelihood
logl=function(t,y,N,D,E){
  n=dim(y)[1]
  f=rep(0,n)
  yr=rep(0:2,each=48)
  for(j in 1:n){
    Ker=t[48+2]*exp(-(D[j%%48,]-t[48+1])^2/(t[j%%48]^2)) # dispersal kernal 
    Nj=N[(yr[j]*48+1):((yr[j]+1)*48),]
    f[j]=sum(Ker*Nj)+t[48+3]*E[j,1]+t[48+4]*E[j,2]+t[48+5]*E[j,3]
  }
  
  return(sum((y[,1]-f)^2))
}

t0=c(rep(1,48),0,1,0.1,0.1,1)
res=nlm(logl,t0,hessian=T,print.level=2,y=y,N=N,D=D,E=E,iterlim=1e4)


t=res$estimate


mean((y[,1]-f)^2)

#=============================================================
# negative log-likelihood
logl=function(t,y,N,D,E){
  n=dim(y)[1]
  f=rep(0,n)
  yr=rep(0:2,each=48)
  for(j in 1:n){
    Ker=t[j%%48]*exp(-(D[j%%48,]-t[96+1])^2/(t[48+j%%48]^2)) # dispersal kernal 
    Nj=N[(yr[j]*48+1):((yr[j]+1)*48),]
    f[j]=sum(Ker*Nj)+t[96+2]*E[j,1]+t[96+4]*E[j,2]+t[96+4]*E[j,3]
  }
  
  return(sum((y[,1]-f)^2))
}

t0=c(rep(1,96),0,0.1,0.1,0.1)
res=nlm(logl,t0,hessian=T,print.level=2,y=y,N=N,D=D,E=E,iterlim=1e4)


t=res$estimate


mean((y[,1]-f)^2)






















x=rnorm(100)
y=(x-0.5)^2+0.3*rnorm(100)

fun=function(x){
  z=(x[1]-50)^2+(x[2]-1)^2+40
  print(x)
  return(z)
}

nlm(fun,c(100,50))




#load("pisco.RData")
#nsite=48
#newdata <- subset(mydata, age >= 20 | age < 10, select=c(ID, Weight)) 

X=subset(data.all.zones,speciesnum==1)

y=pisco.data[pisco.data$year>2000,"mussel_abund"]
N=pisco.data[pisco.data$year<2003,"mussel_abund"]
D=coord2dist()
(pisco.data[1:nsite,1:2])

