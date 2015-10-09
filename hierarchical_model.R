#setwd("/home/liyi/Dropbox/model")
#------------------------------------------------------------------------
# 1. Load the Data
X=read.csv("~/Dropbox/model/adults.recruits.merged.csv")
feature_name=c("year","lat","lon","sitenum","cover","filter_chla_mean","filter_sst_mean","filter_upw_mean","mytilus.rec")
X=X[feature_name]


# plot(X[,"filter_upw_mean"],X[,"mytilus.rec"])
# plot(X[,"filter_upw_mean"],X[,"cover"])
# 
# plot(X[,"filter_chla_mean"],X[,"mytilus.rec"])
# plot(X[,"filter_chla_mean"],X[,"cover"])
# 
# plot(log(X[,"filter_upw_mean"]),log(X[,"mytilus.rec"]))
# plot(log(X[,"filter_upw_mean"]),log(X[,"cover"]))

# X[X[,"year"]==1999+1&X[,"sitenum"]==24,"cover"]
# tail(X,10)

# ------------------------------------------------------------------------
# 2. derive target variable
nsample=dim(X)[1]
X["ycover"]=NA
for(i in 1:nsample){
  idx=(X[,"sitenum"]==X[i,"sitenum"])&(X[,"year"]==X[i,"year"]+1)
  
  if(sum(idx)){
    X[i,"ycover"]=X[idx,"cover"]
  }
}

X["new_upw"]=NA
X["new_chla"]=NA
X["new_sst"]=NA
for(i in 1:nsample){
  idx=(X[,"sitenum"]==X[i,"sitenum"])&(X[,"year"]==X[i,"year"]+1)
  
  if(sum(idx)){
    X[i,"new_upw"]=X[idx,"filter_upw_mean"]
    X[i,"new_chla"]=X[idx,"filter_chla_mean"]
    X[i,"new_sst"]=X[idx,"filter_sst_mean"]
    
  }
}

X1=X[complete.cases(X),]
# tail(X1)
# dim(X1)
# unique(X1["year"])


lats=subset(X1,year==2001,select = c(lat))[,1]
lons=subset(X1,year==2001,select = c(lon))[,1]
lats=floor(lats*1000)
lons=floor(lons*1000)
latlon=paste0(lats,lons) # unique identifier of each site, sorted by lats

for(i in 1:dim(X1)[1]){
  X1[i,"sitenum"]=which(latlon==paste0(floor(X1[i,"lat"]*1000),floor(X1[i,"lon"]*1000)))
}

# ------------------------------------------------------------------------
# 3. Set up the analysis data
# 3.1 Target variable, next year abundance
y=X1[,"ycover"]
y[y==0]=min(y[y!=0])/2 # replace 0 to the minimum positive number
y=log(y) # take log of the species number

# 3.1.5 recruitment data
R=X1[,"mytilus.rec"]
R[R==0]=min(R[R!=0])/2
R=log(R)

# 3.2 Year for each observation
Y=X1[,"year"]

S=X1[,"sitenum"]

# 3.3 Abundance (for current year)
A=X1["cover"]

# 3.4 Distance between sites
D=dist(X[X[,"year"]==2001,"lat"],diag = T,upper = T)
D=as.matrix(D)
for(i in 1:48){
  for(j in 1:48){
    if(i>j){
      D[i,j]=-D[i,j]
    }
  }
}

D=as.matrix(D) #D[1:6,1:6]

# 3.5 Environment variables
E=X1[,c("filter_chla_mean","filter_sst_mean","filter_upw_mean")]

Enew=X1[,c("new_upw","new_chla","new_sst")]

mean_upw=mean(X1[,"filter_upw_mean"])
sd_upw=sd(X1[,"filter_upw_mean"])
X1[,"filter_upw_mean"]=(X1[,"filter_upw_mean"]-mean_upw)/sd_upw
#------------------------------------------------------------------------
# 4. Model 1 (larval+recruitment)

m1logl=function(p,X1,R,D){
  n=dim(X1)[1]
  f=rep(0,n)
  for(iobs in 1:n){
    this_ker=0
    this_site=X1[iobs,"sitenum"]
    this_year=X1[iobs,"year"]
    Xtemp=subset(X1,year==this_year)# all the sites that have the same year as the current observation
    for(j in dim(Xtemp)[1]){
      sitej=Xtemp[j,"sitenum"]
      this_ker=this_ker+Xtemp[j,"cover"]*
        exp(p[4]+p[5]*Xtemp[j,"filter_upw_mean"]-p[6]*Xtemp[j,"filter_upw_mean"]^2
            +p[7]*D[this_site,sitej]-p[8]*D[this_site,sitej]^2)
    }
    #print(this_ker)
    f[iobs]=log(this_ker)+p[1]+p[2]*X1[iobs,"filter_upw_mean"]-p[3]*X1[iobs,"filter_upw_mean"]^2
  }
  return(sum((R-f)^2))
}

#p0=c(0,1,1,0,1,1,1,1)/10
#p0=c(0,10,1,0,1,1,1,1)/10 #210.257
#p0=c(1,1,2,0,1,1,1,1) #210.0695
p0=c(-2,1,2,-2,1,1,1,1)
m1logl(p0,X1,R,D)


resm1=nlm(m1logl,p0,hessian=T,print.level=1,X1=X1,R=R,D=D,iterlim=1e4,steptol=1e-5)
p1=resm1$estimate
p1

mu_hat=p1[7]/(2*p1[8])
sigma_hat=1/sqrt(p1[8])
E_in=p1[5]/(2*p1[6])
E_out=p1[2]/(2*p1[3])
mu_hat
sigma_hat
E_in*sd_upw+mean_upw
E_out*sd_upw+mean_upw

m1pred=function(p,X1,R,D){
  n=dim(X1)[1]
  f=rep(0,n)
  for(iobs in 1:n){
    this_ker=0
    this_site=X1[iobs,"sitenum"]
    this_year=X1[iobs,"year"]
    Xtemp=subset(X1,year==this_year)# all the sites that have the same year as the current observation
    for(j in dim(Xtemp)[1]){
      sitej=Xtemp[j,"sitenum"]
      this_ker=this_ker+Xtemp[j,"cover"]*
        exp(p[4]+p[5]*Xtemp[j,"filter_upw_mean"]-p[6]*Xtemp[j,"filter_upw_mean"]^2
            +p[7]*D[this_site,sitej]-p[8]*D[this_site,sitej]^2)
    }
    #print(this_ker)
    f[iobs]=log(this_ker)+p[1]+p[2]*X1[iobs,"filter_upw_mean"]-p[3]*X1[iobs,"filter_upw_mean"]^2
  }
  return(f)
}

plot(R,m1pred(p1,X1,R,D))
abline(a=0,b=1,col=2)
#------------------------------------------------------------------------
# 5. Model 2 (Abundance)

m2logl=function(p,X1,y){
  n=dim(X1)[1]
  f=rep(0,n)
  for(iobs in 1:n){
    f[iobs]=X1[iobs,"mytilus.rec"]*exp(p[2]+p[3]*X1[iobs,"new_upw"]-p[4]*X1[iobs,"new_upw"]^2)+p[1]*X1[iobs,"cover"]
  }
  return(sum((exp(y)-f)^2))
}

p0=c(1,0,1,1)
m2logl(p0,X1,y)

resm2=nlm(m2logl,p0,hessian=T,print.level=1,X1=X1,y=y,iterlim=1e4,steptol=1e-5)
p2=resm2$estimate
p2

E_hat=p2[3]/(2*p2[4])
E_hat*sd_upw+mean_upw


m2pred=function(p,X1){
  n=dim(X1)[1]
  f=rep(0,n)
  for(iobs in 1:n){
    f[iobs]=X1[iobs,"mytilus.rec"]*exp(p[2]+p[3]*X1[iobs,"new_upw"]-p[4]*X1[iobs,"new_upw"]^2)+p[1]*X1[iobs,"cover"]
  }
  return(f)
}

plot(m2pred(p2,X1),X1[,"ycover"])
abline(a=0,b=1,col=2)


#------------------------------------------------------------------------
plot(X1[,"filter_upw_mean"],X1[,"mytilus.rec"])
plot(X1[,"filter_upw_mean"],exp(y))
plot(X1[,"new_upw"],X1[,"ycover"])

#------------------------------------------------------------------------
# 4.1 log-likelihood for model 1(larval+recruitment)
# input: 
# 1. parameters to optimize
# 2. log of recruitment
# 3. year for each observation
# 3.5 site number
# 4. Abundance
# 5. distance between each sites
# 6. environment variables
#------------------------------------------------------------------------




