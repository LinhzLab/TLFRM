##This is the code of simulation code for example 3
##The code for remain examples are generated just by changing the data generating settings.
install.packages("fda")
install.packages("MASS")
install.packages("Matrix")
install.packages("pracma")

library(fda)
library(MASS)
library(Matrix)
library(pracma)

#E-M algorithm
EM_EST<-function(YY,cof_mu0,cof_phi10,cof_phi20,beta0,sigma0,rho0,aphla0,lamda0)
{
  N=200                    #maximum times of EM iteration
  COF_MU=matrix(NA,N,nk+3)
  COF_PHI1=matrix(NA,N,nk+3)
  COF_PHI2=matrix(NA,N,nk+3)
  
  BETA=array(NA,dim=c(N,3,2))
  
  LAMDA=matrix(NA,N,3)
  LAMDA[,1]=1
  SIGMA=matrix(NA,N,3)
  RHO=matrix(NA,N,2)
  ALPHA=rep(NA,N)
  
  COF_MU[1,]=cof_mu0
  COF_PHI1[1,]=cof_phi10
  COF_PHI2[1,]=cof_phi20
  
  BETA[1,,]=beta0
  
  LAMDA[1,]=lamda0  
  SIGMA[1,]=sigma0  
  RHO[1,]=rho0     
  ALPHA[1]=aphla0   
  object<-rep(NA,N)
  EMbeark=0
  for(em in 1:N)          # EM iteration
  { 
    #E-step
    cond_mean=matrix(NA,n,2)
    cond_var=vector(mode="list",length=n)
    hathatmu<-function(t) spline(t)%*%COF_MU[em,]
    hathatphi1<-function(t) spline(t)%*%COF_PHI1[em,]
    hathatphi2<-function(t) spline(t)%*%COF_PHI2[em,]
    for(i in 1:n) { 
      PSI=matrix(NA,p*time_numb[i],K)
      PHI=matrix(NA,time_numb[i],K)               
      for(s in 1:time_numb[i]) { PHI[s,]=c(hathatphi1(time_poin[[i]][s]),hathatphi2(time_poin[[i]][s])) }
      PSI=kronecker(PHI,LAMDA[em,])
      A=diag(RHO[em,])
      D2=diag(rep(SIGMA[em,],time_numb[i]))
      BB=sapply(time_poin[[i]],hathatmu)
      ERROR=as.vector(YY[[i]]-BETA[em,,]%*%matrix(rep(X[i,],time_numb[i]),nrow=2)-ALPHA[em]*LAMDA[em,]%o%rep(Z[i],time_numb[i])-LAMDA[em,]%o%BB)
      cond_mean[i,]=diag(RHO[em,])%*%t(PSI)%*%solve(PSI%*%diag(RHO[em,])%*%t(PSI)+D2)%*%ERROR
      cond_var[[i]]=diag(RHO[em,])-diag(RHO[em,])%*%t(PSI)%*%solve(PSI%*%diag(RHO[em,])%*%t(PSI)+D2)%*%PSI%*%diag(RHO[em,])
    }
    COF_MUmax=matrix(NA,N,nk+3)
    COF_PHI1max=matrix(NA,N,nk+3)
    COF_PHI2max=matrix(NA,N,nk+3)
    
    BETAmax=array(NA,dim=c(N,3,2))
    
    ALPHAmax=rep(NA,N)
    LAMDAmax=matrix(NA,N,3)
    SIGMAmax=matrix(NA,N,3)
    RHOmax=matrix(NA,N,2) 
    
    LAMDAmax[,1]=1
    COF_MUmax[1,]=COF_MU[em,]
    COF_PHI1max[1,]=COF_PHI1[em,]
    COF_PHI2max[1,]=COF_PHI2[em,]
    
    BETAmax[1,,]=BETA[em,,]
    
    LAMDAmax[1,]=LAMDA[em,]
    SIGMAmax[1,]=SIGMA[em,]
    RHOmax[1,]=RHO[em,]
    ALPHAmax[1]=ALPHA[em]
    Mstep_break=0
    #M-step
    for(max in 1:200) { 
      rho1=0
      rho2=0
      for(i in 1:n){
        rho1=rho1+(cond_var[[i]][1,1]+cond_mean[i,1]^2)
        rho2=rho2+(cond_var[[i]][2,2]+cond_mean[i,2]^2)
      }
      RHOmax[max+1,]=c(rho1/n,rho2/n)
         
      BBMM=t(spline(unlist(time_poin)))
      
      YYY=NULL
      for(i in 1:n){YYY=cbind(YYY,YY[[i]])}
      YYYY=t(YYY)
      
      XXX=NULL
      for(i in 1:n){
        replicat=matrix(rep(X[i,],time_numb[i]),nrow=2)
        XXX=cbind(XXX,replicat)    
      }
      
      
      
      
      ZZZ=rep(Z,time_numb)
      VVV=NULL
      for(i in 1:n){VVV=cbind(VVV,cond_var[[i]])}
      
      DDD=cbind(COF_PHI1max[max,],COF_PHI2max[max,])   
      AAA=as.matrix(as.data.frame(spline(time_poin[[1]])))%*%DDD
      for(i in 2:n){
        BBB=as.matrix(as.data.frame(spline(time_poin[[i]])))%*%DDD
        AAA=as.matrix(bdiag(AAA,BBB))
      } 
      TTT=AAA%*%as.vector(t(cond_mean))
      
      CCC=t(cond_mean[1,])
      for(i in 2:n){
        CCC=as.matrix(bdiag(CCC,t(cond_mean[i,])))
      }
      
      HAIHAI=rbind(COF_PHI1max[max,],COF_PHI2max[max,])
      MMM=HAIHAI
      for(i in 2:n){
        MMM=as.matrix(bdiag(MMM,HAIHAI))
      }
      
      HHH=t(as.matrix(as.data.frame(spline(time_poin[[1]]))))
      for(i in 2:n){
        HHH=as.matrix(bdiag(HHH,t(as.matrix(as.data.frame(spline(time_poin[[i]]))))))
      }
      ###########################################################################################
      FF1=t(YYYY[,1]-BETAmax[max,1,]%*%XXX-1*ALPHAmax[max]*ZZZ-1*(COF_MUmax[max,]%*%BBMM))
      FF11=t(FF1)%*%FF1 
      FF12=1*(YYYY[,1]-BETAmax[max,1,]%*%XXX-1*ZZZ-1*COF_MUmax[max,]%*%BBMM)%*%TTT
      FF13=1^2*sum(diag(DDD%*%(VVV+t(cond_mean)%*%CCC)%*%MMM%*%HHH%*%t(BBMM)))
      SIGMAmax[max+1,1]=(FF11-2*FF12+FF13)/NN
      
      FF2=t(YYYY[,2]-BETAmax[max,2,]%*%XXX-LAMDAmax[max,2]*ALPHAmax[max]*ZZZ-LAMDAmax[max,2]*(COF_MUmax[max,]%*%BBMM))
      FF21=t(FF2)%*%FF2
      FF22=LAMDAmax[max,2]*(YYYY[,2]-BETAmax[max,2,]%*%XXX-LAMDAmax[max,2]*ZZZ-LAMDAmax[max,2]*(COF_MUmax[max,]%*%BBMM))%*%TTT
      FF23=LAMDAmax[max,2]^2*sum(diag(DDD%*%(VVV+t(cond_mean)%*%CCC)%*%MMM%*%HHH%*%t(BBMM)))
      SIGMAmax[max+1,2]=(FF21-2*FF22+FF23)/NN
      
      FF3=t(YYYY[,3]-BETAmax[max,3,]%*%XXX-LAMDAmax[max,3]*ALPHAmax[max]*ZZZ-LAMDAmax[max,3]*(COF_MUmax[max,]%*%BBMM))
      FF31=t(FF3)%*%FF3
      FF32=LAMDAmax[max,3]*(YYYY[,3]-BETAmax[max,3,]%*%XXX-LAMDAmax[max,3]*ZZZ-LAMDAmax[max,3]*(COF_MUmax[max,]%*%BBMM))%*%TTT
      FF33=LAMDAmax[max,3]^2*sum(diag(DDD%*%(VVV+t(cond_mean)%*%CCC)%*%MMM%*%HHH%*%t(BBMM)))
      SIGMAmax[max+1,3]=(FF31-2*FF32+FF33)/NN
      #############################################################################################
      
      BETAmax[max+1,1,]=solve(XXX%*%t(XXX))%*%(XXX%*%(YYYY[,1]-1*ZZZ*ALPHAmax[max]-1*t(BBMM)%*%COF_MUmax[max,]-1*TTT))
      BETAmax[max+1,2,]=solve(XXX%*%t(XXX))%*%(XXX%*%(YYYY[,2]-LAMDAmax[max,2]*ZZZ*ALPHAmax[max]-LAMDAmax[max,2]*t(BBMM)%*%COF_MUmax[max,]-LAMDAmax[max,2]*TTT))
      BETAmax[max+1,3,]=solve(XXX%*%t(XXX))%*%(XXX%*%(YYYY[,3]-LAMDAmax[max,3]*ZZZ*ALPHAmax[max]-LAMDAmax[max,3]*t(BBMM)%*%COF_MUmax[max,]-LAMDAmax[max,3]*TTT))
      
      ##############################################################################################      
      FZ2=(YYYY[,2]-BETAmax[max+1,2,]%*%XXX)%*%(ZZZ*ALPHAmax[max]+t(BBMM)%*%COF_MUmax[max,]+TTT)
      FU1=(t(ZZZ*ALPHAmax[max])+COF_MUmax[max,]%*%BBMM)%*%t(t(ZZZ*ALPHAmax[max])+COF_MUmax[max,]%*%BBMM)
      FU2=2*(t(ZZZ*ALPHAmax[max])+COF_MUmax[max,]%*%BBMM)%*%TTT
      FU3=sum(diag(DDD%*%(VVV+t(cond_mean)%*%CCC)%*%MMM%*%HHH%*%t(BBMM)))
      FZ3=(YYYY[,3]-BETAmax[max+1,3,]%*%XXX)%*%(ZZZ*ALPHAmax[max]+t(BBMM)%*%COF_MUmax[max,]+TTT)
      LAMDAmax[max+1,2]=FZ2/(FU1+FU2+FU3)
      LAMDAmax[max+1,3]=FZ3/(FU1+FU2+FU3)
      ##############################################################################################      
      GOOD=t(ZZZ)%*%(YYYY-t(XXX)%*%t(BETAmax[max+1,,])-t(BBMM)%*%COF_MUmax[max,]%*%LAMDAmax[max+1,]-TTT%*%LAMDAmax[max+1,])%*%(LAMDAmax[max+1,]/SIGMAmax[max+1,])
      ALPHAmax[max+1]=solve(t(ZZZ)%*%ZZZ*sum(LAMDAmax[max+1,]^2/SIGMAmax[max+1,]))%*%GOOD
      
      ##############################################################      
      GOOD=BBMM%*%(YYYY-t(XXX)%*%t(BETAmax[max+1,,])-ALPHAmax[max+1]*ZZZ%o%LAMDAmax[max+1,]-TTT%*%LAMDAmax[max+1,])%*%(LAMDAmax[max+1,]/SIGMAmax[max+1,])
      COF_MUmax[max+1,]=solve(BBMM%*%t(BBMM)*sum(LAMDAmax[max+1,]^2/SIGMAmax[max+1,]))%*%GOOD
      #########################################################################
      HAHA=unlist(lapply(cond_var,function(x) x[1,1]))+cond_mean[,1]^2
      HIHI=solve(BBMM%*%diag(rep(HAHA,time_numb))%*%t(BBMM)*sum(LAMDAmax[max+1,]^2/SIGMAmax[max+1,]))
      
      BBB1=diag(rep(cond_mean[,1],time_numb))
      AA=unlist(lapply(cond_var,function(x) x[1,2]))+apply(cond_mean,1,function(x) x[1]*x[2])
      AAA=diag(rep(AA,time_numb))
      sssolut1=BBMM%*%(BBB1%*%(YYYY-t(XXX)%*%t(BETAmax[max+1,,])-ALPHAmax[max+1]*ZZZ%o%LAMDAmax[max+1,]-(t(BBMM)%*%COF_MUmax[max+1,])%*%LAMDAmax[max+1,])-AAA%*%t(BBMM)%*%COF_PHI2max[max,]%*%LAMDAmax[max+1,])%*%(LAMDAmax[max+1,]/SIGMAmax[max+1,])
      solut1=HIHI%*%sssolut1
      
      HAHA=unlist(lapply(cond_var,function(x) x[2,2]))+cond_mean[,2]^2
      HIHI=solve(BBMM%*%diag(rep(HAHA,time_numb))%*%t(BBMM)*sum(LAMDAmax[max+1,]^2/SIGMAmax[max+1,]))
      BBB2=diag(rep(cond_mean[,2],time_numb))
      sssolut2=BBMM%*%(BBB2%*%(YYYY-t(XXX)%*%t(BETAmax[max+1,,])-ALPHAmax[max+1]*ZZZ%o%LAMDAmax[max+1,]-(t(BBMM)%*%COF_MUmax[max+1,])%*%LAMDAmax[max+1,])-AAA%*%t(BBMM)%*%COF_PHI1max[max,]%*%LAMDAmax[max+1,])%*%(LAMDAmax[max+1,]/SIGMAmax[max+1,])
      solut2=HIHI%*%sssolut2
      ############################################################################################
      BTHETA=matrix(c(solut1,solut2),nrow=2,byrow=TRUE)
      
      VV=BTHETA%*%WORKING%*%t(BTHETA)
      UU=solve(sqrtm(VV)$B)%*%BTHETA
      COF_PHI1max[max+1,]=UU[1,]
      COF_PHI2max[max+1,]=UU[2,]
      
      testphi1<-function(t) bs(t,df=NULL,kont,degree=3,intercept=T,Boundary.knots=c(0,1))%*%COF_PHI1max[max+1,]
      testphi2<-function(t) bs(t,df=NULL,kont,degree=3,intercept=T,Boundary.knots=c(0,1))%*%COF_PHI2max[max+1,]
      if(testphi1(0.5)<0){COF_PHI1max[max+1,]=-COF_PHI1max[max+1,]}
      if(testphi2(0.25)<0){COF_PHI2max[max+1,]=-COF_PHI2max[max+1,]}
      AA=c(COF_MUmax[max+1,],COF_PHI1max[max+1,],COF_PHI2max[max+1,],as.vector(BETAmax[max+1,,]),LAMDAmax[max+1,],SIGMAmax[max+1,],RHOmax[max+1,],ALPHAmax[max+1])
      BB=c(COF_MUmax[max,],COF_PHI1max[max,],COF_PHI2max[max,],as.vector(BETAmax[max,,]),LAMDAmax[max,],SIGMAmax[max,],RHOmax[max,],ALPHAmax[max])
      diff=mean(abs(AA-BB))
      cat("maximum diff=",diff,"\n")
      if(diff<1e-4) { break }
      if(max==199){Mstep_break=1}
      if(max==199){break}
    }
    COF_MU[em+1,]=COF_MUmax[max+1,]
    COF_PHI1[em+1,]=COF_PHI1max[max+1,]
    COF_PHI2[em+1,]=COF_PHI2max[max+1,]
    
    BETA[em+1,,]=BETAmax[max+1,,]
    
    LAMDA[em+1,]=LAMDAmax[max+1,]
    SIGMA[em+1,]=SIGMAmax[max+1,]
    RHO[em+1,]=RHOmax[max+1,]
    ALPHA[em+1]=ALPHAmax[max+1]
    
    r1=mean(abs(COF_MU[em+1,]-COF_MU[em,]))
    r2=mean(abs(COF_PHI1[em+1,]-COF_PHI1[em,]))
    r3=mean(abs(COF_PHI2[em+1,]-COF_PHI2[em,]))
    r4=mean(abs(BETA[em+1,,]-BETA[em,,]))
    r5=mean(abs(LAMDA[em+1,]-LAMDA[em,]))
    r6=mean(abs(SIGMA[em+1,]-SIGMA[em,]))
    r7=mean(abs(RHO[em+1,]-RHO[em,]))
    r8=mean(abs(ALPHA[em+1]-ALPHA[em]))
    
    error=mean(c(r1,r2,r3,r4,r5,r6,r7,r8))
    cat("EM diff=",error,"\n")
    if(error<5e-4) {break}
    if(em==N-1){EMbeark=1}
    if(em==N-1){break}
  }
  Updated_Est<-list(COF_MU=COF_MU[em+1,],COF_PHI1=COF_PHI1[em+1,],COF_PHI2=COF_PHI2[em+1,],BETA=BETA[em+1,,],SIGMA=SIGMA[em+1,],RHO=RHO[em+1,],ALPHA=ALPHA[em+1],LAMDA=LAMDA[em+1,],EMbreak=EMbeark)
  return(Updated_Est)
}

#objective function
orig<-function(hyp,COF_MU,COF_PHI1,COF_PHI2,BETAp,SIGMAp,RHO,ALPHA,LAMDAp)
{
  out=0
  meanf<-function(t) spline(t)%*%COF_MU
  engf1<-function(t) spline(t)%*%COF_PHI1
  engf2<-function(t) spline(t)%*%COF_PHI2
  for(i in 1:n){
    for(j in 1:time_numb[i]){
      jjkzhy=c(engf1(time_poin[[i]][j]),engf2(time_poin[[i]][j]))
      AAAA=jjkzhy%*%diag(RHO)%*%jjkzhy
      sigma=sqrt(AAAA*LAMDAp^2+SIGMAp)
      out=out+pnorm((hyp-X[i,]%*%BETAp-LAMDAp%*%Z[i]%*%ALPHA-LAMDAp*meanf(time_poin[[i]][j]))/sigma)
      
    }
  }
  return(out)
}

#first derivative function
first_order<-function(hyp,COF_MU,COF_PHI1,COF_PHI2,BETAp,SIGMAp,RHO,ALPHA,LAMDAp)
{
  out=0
  meanf<-function(t) spline(t)%*%COF_MU
  engf1<-function(t) spline(t)%*%COF_PHI1
  engf2<-function(t) spline(t)%*%COF_PHI2
  for(i in 1:n){
    for(j in 1:time_numb[i]){
      jjkzhy=c(engf1(time_poin[[i]][j]),engf2(time_poin[[i]][j]))
      AAAA=jjkzhy%*%diag(RHO)%*%jjkzhy
      sigma=sqrt(AAAA*LAMDAp^2+SIGMAp)
      out=out+dnorm((hyp-X[i,]%*%BETAp-LAMDAp%*%Z[i]%*%ALPHA-LAMDAp*meanf(time_poin[[i]][j]))/sigma)/sigma
      
    }
  }
  return(out)
}

#transformation function estimation
transform_est<-function(pp,inpuevalue,COF_MU,COF_PHI1,COF_PHI2,BETAp,SIGMAp,RHO,ALPHA,LAMDAp,lowbound,upbound,quantile025,quantile975)
{
  geshu=50
  point=seq(quantile025,quantile975,length=geshu)
  HY_MATRIX<-rep(NA,geshu)
  YP=NULL
  for(i in 1:n){YP=c(YP,Y[[i]][pp,])}
  for(j in 1:geshu){
    if(j==1){initial=inpuevalue}else{initial=HY_MATRIX[j-1]}
    before=0
    after=1
    s=0
    s=sum(YP<point[j])
    while(abs(after-before)>1e-3){
      orig=orig(initial,COF_MU,COF_PHI1,COF_PHI2,BETAp,SIGMAp,RHO,ALPHA,LAMDAp)
      first_or=first_order(initial,COF_MU,COF_PHI1,COF_PHI2,BETAp,SIGMAp,RHO,ALPHA,LAMDAp)
      after=initial+(s-orig)/(first_or+1e-4)
      before=initial
      initial=after
    }
    HY_MATRIX[j]<-after
    if(j>1){if(HY_MATRIX[j]<HY_MATRIX[j-1]){HY_MATRIX[j]=HY_MATRIX[j-1]}}
    cat("j=",j,"\n")
  }
  median=0
  scale=0
  if(pp==1){
    median=H1(point[25])
    scale=mean(abs(H1(point)))
  }
  if(pp==2){
    median=H2(point[25])
    scale=mean(abs(H2(point)))
  }
  if(pp==3){
    median=H3(point[25])
    scale=mean(abs(H3(point)))
  }
  vvcc=HY_MATRIX[25]
  HY_MATRIX=HY_MATRIX-(vvcc-median)
  HY_MATRIX=HY_MATRIX*(scale/mean(abs(HY_MATRIX)))
  vvcc=HY_MATRIX[25]
  HY_MATRIX=HY_MATRIX-(vvcc-median)
  
  point=seq(quantile025,quantile975,length=geshu)
  DESIMATRX_low=cbind(rep(1,10),point[1:10])
  DESIMATRX_up=cbind(rep(1,10),point[(geshu-9):geshu])
  textpoint1=seq(quantile975,upbound,length=geshu)
  textpoint2=seq(lowbound,quantile025,length=geshu)
  coe_up<-solve(t(DESIMATRX_up)%*%DESIMATRX_up)%*%t(DESIMATRX_up)%*%HY_MATRIX[(geshu-9):geshu]
  coe_low<-solve(t(DESIMATRX_low)%*%DESIMATRX_low)%*%t(DESIMATRX_low)%*%HY_MATRIX[1:10]
  tranfun_up<-function(t) c(1,t)%*%coe_up
  tranfun_low<-function(t) c(1,t)%*%coe_low
  
  jjk1=sapply(textpoint1,tranfun_up)
  jjk2=sapply(textpoint2,tranfun_low)
  transfunc<-c(jjk2,HY_MATRIX,jjk1)
  return(transfunc)
}

identify<-function(x,vect){
  location=rep(NA,length(x))
  for(i in 1:length(x)){
    location[i]=which.min(abs(rep(x[i],length(vect))-vect))
  }
  return(location)
}


n=200    # sample size
p=3      # p outcomes
K=2      # K eigenfunctions
indict<-function(j,k) {if(j==k){1} else{0}}
square<-function(t) t^2
nk=3
kont<-c(0.35,0.75)   # spline konts
spline<-function(t) bs(t,df=NULL,kont,degree=3,intercept=T,Boundary.knots=c(0,1))

rho0=c(1,0.25)   # score variance
lamda0=c(1,1,1)  # factor loadings
beta0=matrix(c(0.25,0.25,-0.5,0.25,0.25,-0.5),3,2)  # regression coefficient,two covariate X 3 by 2
sigma0=c(0.5,0.5,0.5)                               # variance of error
aphla0=1                                            # latent process regression coefficient,just one covariate two?
mu<-function(t) t^2-t 
phi1<-function(t) sqrt(2)*sin(pi*t)   # first eigenfunction 
phi2<-function(t) sqrt(2)*sin(2*pi*t) # second eigenfunction

samp=seq(0,1,length=100)
desig<-spline(samp)
y_mu=mu(samp)
phi1_mu<-phi1(samp)
phi2_mu<-phi2(samp)
project=solve(t(desig)%*%desig)%*%t(desig)
cof_mu0<-project%*%y_mu 
cof_phi10<-project%*%phi1_mu
cof_phi20<-project%*%phi2_mu


WORK<-matrix(0,nk+3,nk+3)
for(i in 1:10000){ WORK=WORK+spline(i/10000)[1,]%o%spline(i/10000)[1,] }
WORKING=WORK/10000

SQ1<-function(t) 0.2*t
SQ2<-function(t) 0.2*t
SQ3<-function(t) 0.2*t
H1<-function(t) 5*t
H2<-function(t) 5*t
H3<-function(t) 5*t

outcome=matrix(NA,200,33)      #estimates of all parameters in 200 simulations
outH1=matrix(NA,200,150)       #estimates of transformation function  H_1 in 200 simulations
outH2=matrix(NA,200,150)       #estimates of transformation function  H_2 in 200 simulations
outH3=matrix(NA,200,150)       #estimates of transformation function  H_3 in 200 simulations

for(iter in 1:200){ # 200 times simulation
  set.seed(iter)
  ################data generate####################
  time_numb=ceiling(10*runif(n))         
  time_poin=vector(mode="list",length=n) 
  for (i in 1:n)
  { time=sort(runif(time_numb[i]))
  time_poin[[i]]=time
  }
  X=matrix(rnorm(2*n),n,2)            
  Z=rbinom(n,1,0.5)                   
  
  XI=matrix(NA,n,K)
  XI[,1]=sqrt(rho0[1])*rnorm(n)
  XI[,2]=sqrt(rho0[2])*rnorm(n)
  
  
  Y=vector(mode="list",length=n)
  for(i in 1:n)
  { 
    ZZHY=rbind(phi1(time_poin[[i]]),phi2(time_poin[[i]]))
    JJJ=beta0%*%matrix(rep(X[i,],time_numb[i]),nrow=2)+lamda0%o%(aphla0*rep(Z[i],time_numb[i])+mu(time_poin[[i]])+as.vector(XI[i,]%*%ZZHY))+matrix(rnorm(3*time_numb[i]),3,time_numb[i])*sqrt(sigma0[1])
    KKK=matrix(NA,3,time_numb[i])
    KKK[1,]=SQ1(JJJ[1,])
    KKK[2,]=SQ2(JJJ[2,])
    KKK[3,]=SQ3(JJJ[3,])
    Y[[i]]=KKK
  }
  NN=sum(time_numb)
  ################above is data generate process####################
  
  Y1=NULL
  Y2=NULL
  Y3=NULL
  for(i in 1:n){Y1=c(Y1,Y[[i]][1,])}
  for(i in 1:n){Y2=c(Y2,Y[[i]][2,])}
  for(i in 1:n){Y3=c(Y3,Y[[i]][3,])}
  
  
  lowbound1=min(Y1)
  upbound1=max(Y1)
  quantile0251=quantile(Y1,0.015)
  quantile9751=quantile(Y1,0.985)
  allpoint1=c(seq(lowbound1,quantile0251,length=50),seq(quantile0251,quantile9751,length=50),seq(quantile9751,upbound1,length=50))
  
  lowbound2=min(Y2)
  upbound2=max(Y2)
  quantile0252=quantile(Y2,0.015)
  quantile9752=quantile(Y2,0.985)
  allpoint2=c(seq(lowbound2,quantile0252,length=50),seq(quantile0252,quantile9752,length=50),seq(quantile9752,upbound2,length=50))
  
  
  lowbound3=min(Y3)
  upbound3=max(Y3)
  quantile0253=quantile(Y3,0.015)
  quantile9753=quantile(Y3,0.985)
  allpoint3=c(seq(lowbound3,quantile0253,length=50),seq(quantile0253,quantile9753,length=50),seq(quantile9753,upbound3,length=50))
  
  ##################################################################################following is EM
  N=200
  HHHY1<-list()
  HHHY2<-list()
  HHHY3<-list()
  HHHY1[[1]]=rep(1,150)
  HHHY2[[1]]=rep(1,150)
  HHHY3[[1]]=rep(1,150)
  
  cofmu=matrix(NA,200,6)
  cof_phi1=matrix(NA,200,6)
  cof_phi2=matrix(NA,200,6)
  
  beta=array(NA,dim=c(200,3,2))
  
  sigma=matrix(NA,200,3)
  rho=matrix(NA,200,2)
  aphla=rep(NA,200)
  lamda=matrix(NA,200,3)
  
  cofmu[1,]=cof_mu0
  cof_phi1[1,]=cof_phi10
  cof_phi2[1,]=cof_phi20
  
  beta[1,,]=beta0
  
  sigma[1,]=sigma0
  rho[1,]=rho0
  aphla[1]=aphla0
  lamda[1,]=lamda0
  par(mfrow=c(1,3))
  all_break=0
  for(glob in 1:200){ # transformation function estimation
    HHHY1[[glob+1]]=transform_est(1,H1(quantile0251),cofmu[glob,],cof_phi1[glob,],cof_phi2[glob,],beta[glob,1,],sigma[glob,1],rho[glob,],aphla[glob],lamda[glob,1],lowbound1,upbound1,quantile0251,quantile9751)
    HHHY2[[glob+1]]=transform_est(2,H2(quantile0252),cofmu[glob,],cof_phi1[glob,],cof_phi2[glob,],beta[glob,2,],sigma[glob,2],rho[glob,],aphla[glob],lamda[glob,2],lowbound2,upbound2,quantile0252,quantile9752)
    HHHY3[[glob+1]]=transform_est(3,H3(quantile0253),cofmu[glob,],cof_phi1[glob,],cof_phi2[glob,],beta[glob,3,],sigma[glob,3],rho[glob,],aphla[glob],lamda[glob,3],lowbound3,upbound3,quantile0253,quantile9753)
    plot(allpoint1,HHHY1[[glob+1]],ylim=c(-10,10))
    curve(H1,add=T)
    plot(allpoint2,HHHY2[[glob+1]],ylim=c(-10,10))
    curve(H2,add=T)
    plot(allpoint3,HHHY3[[glob+1]],ylim=c(-10,10))
    curve(H3,add=T)
    ##############################################################################Above is transformation estimation
    YY=vector(mode="list",length=n)
    for(i in 1:n){
      yuanshi=Y[[i]]
      respon=matrix(NA,3,time_numb[i])
      for(k in 1:3){
        if(k==1){respon[k,]=HHHY1[[glob+1]][identify(yuanshi[k,],allpoint1)]}
        if(k==2){respon[k,]=HHHY2[[glob+1]][identify(yuanshi[k,],allpoint2)]}
        if(k==3){respon[k,]=HHHY3[[glob+1]][identify(yuanshi[k,],allpoint3)]}
      }
      YY[[i]]=respon
    }
    ##############################################################################
    EMout<-EM_EST(YY,cofmu[glob,],cof_phi1[glob,],cof_phi2[glob,],beta[glob,,],sigma[glob,],rho[glob,],aphla[glob],lamda[glob,])
    if(EMout$EMbreak==1){
      all_break=1
      break}
    cofmu[glob+1,]=EMout$COF_MU
    cof_phi1[glob+1,]=EMout$COF_PHI1
    cof_phi2[glob+1,]=EMout$COF_PHI2
    
    beta[glob+1,,]=EMout$BETA
    
    sigma[glob+1,]=EMout$SIGMA
    rho[glob+1,]=EMout$RHO
    aphla[glob+1]=EMout$ALPHA
    lamda[glob+1,]=EMout$LAMDA 
    #################################################################################Above is EM
    H_diff=mean(c(abs(HHHY1[[glob+1]][50:100]-HHHY1[[glob]][50:100]),abs(HHHY2[[glob+1]][50:100]-HHHY2[[glob]][50:100]),abs(HHHY3[[glob+1]][50:100]-HHHY3[[glob]][50:100])))
    aft=c(cofmu[glob+1,],cof_phi1[glob+1,],cof_phi2[glob+1,],as.vector(beta[glob+1,,]),sigma[glob+1,],rho[glob+1,],aphla[glob+1],lamda[glob+1,])
    bef=c(cofmu[glob,],cof_phi1[glob,],cof_phi2[glob,],as.vector(beta[glob,,]),sigma[glob,],rho[glob,],aphla[glob],lamda[glob,])
    par_diff=mean(abs(aft-bef))
    cat("iter=",iter)
    cat("glob=",glob)
    cat("H_diff=",H_diff)
    cat("par_diff=",par_diff)
    if(H_diff<1e-3&par_diff<1e-3){break}
    if(glob==199){
      all_break=1
      break}
  }
  if(all_break==0){
  outH1[iter,]=HHHY1[[glob+1]]
  outH2[iter,]=HHHY2[[glob+1]]
  outH3[iter,]=HHHY3[[glob+1]]
  outcome[iter,]=c(cofmu[glob+1,],cof_phi1[glob+1,],cof_phi2[glob+1,],as.vector(beta[glob+1,,]),sigma[glob+1,],rho[glob+1,],aphla[glob+1],lamda[glob+1,])
}
  }

