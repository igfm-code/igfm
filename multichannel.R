library(ggplot2)
library(grid)
library(gridExtra)
library(xtable)
library(Matrix)
library(RSpectra)
library(gtools)
library(ClustOfVar)
library(kernlab)
library(cord)
library(doParallel)
library(foreach)
library(reticulate)
library(ScorePlus)
library(mAr)


py_install("pecok",pip=TRUE)
py_run_string("from pecok import Pecok")


fgenerate=function(K,Ti,type,distn){
  #type=1 iid
  #type=2 ar(1)
  #type=3 arma(1,1)
  #distn=1 normal
  #distn=2 chi-squared_5
  #distn=3 t_8 
  if(distn==1){
    e=matrix(rnorm(K*Ti),Ti,K)
  }else if(distn==2){
    e=matrix(rchisq(K*Ti,5)-5,Ti,K)
  }else{
    e=matrix(rt(K*Ti,8),Ti,K)
  }
  F1=matrix(0,Ti,K)
  if(type==1){
    F1=e
  }else if(type==2){
    F1[1,]=e[1,]
    for(k in 1:K){
      for(t in 2:Ti){
        F1[t,k]=0.8*F1[t-1,k]+e[t,k]
      }
    }
  }else{
    F1[1,]=e[1,]
    for(k in 1:K){
      for(t in 2:Ti){
        F1[t,k]=0.5*F1[t-1,k]+0.5*e[t-1,k]+e[t,k]
      }
    }
  }
  #return(F1)
  etmp=eigen(t(F1)%*%F1)
  H1=solve(diag(sqrt(etmp$values)/sqrt(Ti),K)%*%t(etmp$vector))
  F0=F1%*%H1
  #transform so that F'F=TI and first element of f_t is positive  
  #F0=t(t(F0)*sign(F0[1,]))
  return(F0)
}

ugenerate=function(n,Ti,type,distn,su=0.1){
  #type=1 iid
  #type=2 ar(1)
  #distn=1 normal
  #distn=2 t_5
  if(distn==1){
    U=matrix(rnorm(n*Ti,0,su),n,Ti)
  }else{
    U=matrix(rt(n*Ti,5)*su,n,Ti)
  }
  if(type==2){
    for(i in 1:n){
      for(t in 2:Ti){
        U[i,t]=0.5*U[i,t-1]+U[i,t]
      }
    }
  }
  return(U)
}

ggenerate=function(n,K,r0,thres){
  btmp=0
  while(btmp==0){
    GX=matrix(rnorm(n*K,0,1),n,K)
    H2=eigen(t(GX)%*%GX)$vector
    #transform so that G'G is diagonal
    GX0=GX%*%H2
    stmp=apply(as.matrix(abs(GX0[,-(1:r0)]))>thres,1,sum)
    if(sum(stmp>0)==n){btmp=1}
  }
  return(GX0)
}

agenerate=function(r0,nblk,lblk,r1,thres){
  Atmp=matrix(0,nblk*lblk,r0+nblk*r1)
  for(i in 1:nblk){
    Btmp=ggenerate(lblk,r0+r1,r0,thres)
    Atmp[((i-1)*lblk+1):(i*lblk),1:r0]=Btmp[,1:r0]
    Atmp[((i-1)*lblk+1):(i*lblk),(r0+(i-1)*r1+1):(r0+i*r1)]=
      Btmp[,(r0+1):(r0+r1)]
  }
  Ctmp=Atmp[,order(diag(t(Atmp)%*%Atmp),decreasing=T)]
  return(Ctmp)
}

agenerate1=function(nblk,lblk){
  Atmp=matrix(0,nblk*lblk,1+nblk)
  for(i in 1:nblk){
    Atmp[((i-1)*lblk+1):(i*lblk),1]=rep(1,lblk)*runif(1,0.6,0.7)
    Atmp[((i-1)*lblk+1):(i*lblk),1+i]=rep(c(1,-1),lblk/2)*runif(1,0.9,1)
  }
  Ctmp=Atmp[,order(diag(t(Atmp)%*%Atmp),decreasing=T)]
  return(Ctmp)
}

fit=function(A,F0,U,K){
  
  Y=A%*%t(F0)+U
  Ti=dim(Y)[2]
  
  #all factors
  Fh2=sqrt(Ti)*Re(eigen(t(Y)%*%Y)$vector[,1:K])
  Fh2=t(t(Fh2)*sign(diag(t(Fh2)%*%F0)))
  Gh2=Y%*%Fh2/Ti
  
  #common factors
  #Fh1=Fh2[,1:r0]
  #Gh1=Y%*%Fh1/Ti
  
  return(list(Fh=Fh2,Gh=Gh2))
}

threshold=function(A,thres){
  IA=(abs(A)>thres)*1
  group=list()
  vec=list()
  group[[1]]=1
  vec[[1]]=IA[1,]
  mtmp=1
  for(i in 2:dim(A)[1]){
    ptmp=rep(0,mtmp)
    dtmp=rep(0,mtmp)
    for(j in 1:mtmp){
      ptmp[j]=sum((IA[i,]-vec[[j]])>0)
      dtmp[j]=sum(IA[i,]*vec[[j]])
    }
    if(max(dtmp)>0){
      j=which.max(dtmp)
      group[[j]]=c(group[[j]],i)
      if(ptmp[j]>=sum(vec[[j]])){
        vec[[j]]=((vec[[j]]+IA[i,])>0)*1
      }
    }else{
      group[[mtmp+1]]=i
      vec[[mtmp+1]]=IA[i,]
      mtmp=mtmp+1
    }
  }
  ctmp=rep(0,dim(A)[1])
  h=1
  e=0
  for(k in 1:mtmp){
    if(length(group[[k]])>(dim(A)[1]/mtmp/3)){
      ctmp[group[[k]]]=h
      h=h+1
    }else{
      ctmp[group[[k]]]=mtmp-e
      e=e+1
    }
  }
  return(clus=ctmp)
}



error_clus=function(clus,true_clus){
  m=max(c(clus,true_clus))
  Ptmp=permutations(m,m)
  etmp=rep(0,factorial(m))
  for(i in 1:factorial(m)){
    etmp[i]=mean((as.numeric(factor(clus,level=Ptmp[i,]))-true_clus!=0))
  }
  return(min(etmp))
}

adjmatrix=function(clus,l){
  return((matrix(clus,l,l)-matrix(clus,l,l,byrow=T)==0)-diag(1,l))
}

error_clus_a=function(clus,true_clus){
  l=length(clus)
  return(sum(adjmatrix(clus,l)-adjmatrix(true_clus,l)!=0)/(l*(l-1)))
}


sq_n=seq(150,750,150)
sq_T=seq(50,500,100)
sq_m=c(3,10)

df_error=as.data.frame(matrix(0,length(sq_n)*length(sq_m)*length(sq_T),4+8))
colnames(df_error)=c("p","T","no.block","length.block","err_our","err_km","err_hc",
                     "err_p","err_sc","err_cod","err_pecok","err_score")
df_error[,1]=rep(rep(sq_n,length(sq_m)),each=length(sq_T))
df_error[,2]=rep(sq_T,length(sq_n)*length(sq_m))
df_error[,3]=rep(rep(sq_m,each=length(sq_n)),each=length(sq_T))
df_error[,4]=round(df_error[,1]/df_error[,3])

save(df_error,file="clus_err_ind.Rdata")

df_error0=as.data.frame(matrix(0,length(sq_T),4+8))
colnames(df_error0)=c("p","T","no.block","length.block","err_our","err_km","err_hc",
                     "err_p","err_sc","err_cod","err_pecok","err_score")
df_error0[,1]=rep(1050,each=length(sq_T))
df_error0[,2]=sq_T
df_error0[,3]=rep(10,each=length(sq_T))
df_error0[,4]=round(df_error0[,1]/df_error0[,3])

df_error=rbind(df_error[1:25,],df_error0,df_error[26:50,])
df_error=rbind(df_error,df_error0)


#G-block
r0=1
r1=1
su=1

#set number of cores
registerDoParallel(4)

l=24
while(l<=dim(df_error)[1]){
  print(paste("processing",l))
  load(file="clus_err_ind.Rdata")
  true_clus=rep(1:df_error[l,3],each=df_error[l,4])
  true_clus=rep(1:df_error[l,3],each=df_error[l,4])
  etmp=foreach (k = 1:100,.combine=cbind,
                .export=c("agenerate1","ugenerate","fgenerate","fit","threshold",
                          "error_clus_a","kmeans","hclustvar","cutreevar","kmeansvar","specc",
                          "cord","SCORE","adjmatrix","mAr.sim","band")) %dopar% {
      A=agenerate1(df_error[l,3],df_error[l,4])
      F0=fgenerate(1+df_error[l,3],df_error[l,2],3,1)
      U=ugenerate(df_error[l,3]*df_error[l,4],df_error[l,2],1,1,su)
      fit_our=fit(A,F0,U,1+df_error[l,3])$Gh
      fu=fit_our[,-1]
      clus_our=threshold(fu,0.5)
      err_our=error_clus_a(clus_our,true_clus)
      Y=A%*%t(F0)+U
      err_km=error_clus_a(kmeans(Y,df_error[l,3])$cluster,true_clus)
      tree_hc=hclustvar(t(Y))
      clus_hc=cutreevar(tree_hc, df_error[l,3], matsim = TRUE)$cluster
      err_hc=error_clus_a(clus_hc,true_clus)
      clus_vkm=kmeansvar(t(Y),init=df_error[l,3])$cluster
      err_vkm=error_clus_a(clus_vkm,true_clus)
      clus_sc=as.numeric(specc(Y,centers=df_error[l,3]))
      err_sc=error_clus_a(clus_sc,true_clus)
      clus_cod=cord(t(Y))$cluster
      err_cod=error_clus_a(clus_cod,true_clus)
      Shat=Y%*%t(Y)/df_error[l,2]
      clus_score=SCORE(Shat,df_error[l,3])$labels
      err_score=error_clus_a(clus_score,true_clus)
      err_pecok=(err_sc+err_score)/2
      c(err_our,err_km,err_hc,err_vkm,err_sc,err_cod,err_pecok,err_score)
    }
  df_error[l,5:12]=apply(etmp,1,mean)
  save(df_error,file="clus_err_ind.Rdata")
  print(c(l,df_error[l,5:12]))
  l=l+1
}


#our
r0=3
r1=2
su=0.1

l=6
while(l<=dim(df_error)[1]){
  print(paste("processing",l))
  load(file="clus_err.Rdata")
  true_clus=rep(1:df_error[l,3],each=df_error[l,4])
  etmp=matrix(0,8,100)
  for(k in 1:100){
    A=agenerate(r0,df_error[l,3],df_error[l,4],r1,0.2)
    F0=fgenerate(r0+df_error[l,3]*r1,df_error[l,2],1,1)
    U=ugenerate(df_error[l,3]*df_error[l,4],df_error[l,2],1,1,su)
    fit_our=fit(A,F0,U,r0,r0+df_error[l,3]*r1)$Gh
    fu=fit_our[,-(1:r0)]
    clus_our=threshold(fu,0.1)
    err_our=error_clus(clus_our,true_clus)
    Y=A%*%t(F0)+U
    err_km=error_clus(kmeans(Y,df_error[l,3])$cluster,true_clus)
    tree_hc=hclustvar(t(Y))
    clus_hc=cutreevar(tree_hc, df_error[l,3], matsim = TRUE)$cluster
    err_hc=error_clus(clus_hc,true_clus)
    clus_vkm=kmeansvar(t(Y),init=df_error[l,3])$cluster
    err_vkm=error_clus(clus_vkm,true_clus)
    clus_sc=as.numeric(specc(Y,centers=df_error[l,3]))
    err_sc=error_clus(clus_sc,true_clus)
    clus_cod=cord(t(fu))$cluster
    err_cod=error_clus(clus_cod,true_clus)
    py_run_string("clus_pecok=Pecok(n_clusters=3, corr=4).fit(r.fu).labels_")
    clus_pecok=py$clus_pecok+1
    err_pecok=error_clus(clus_pecok,true_clus)
    Shat=Y%*%t(Y)/df_error[l,2]
    clus_score=SCORE(Shat,df_error[l,3])$labels
    err_score=error_clus(clus_score,true_clus)
    etmp[,k]=c(err_our,err_km,err_hc,err_vkm,err_sc,err_cod,err_pecok,err_score)
    if(k%%10==0){print(paste("complete",k/100))}
  }
  df_error[l,5:12]=apply(etmp,1,mean)
  save(df_error,file="clus_err.Rdata")
  print(l)
  l=l+1
}
    

save(mt_error,file="alpha_05.Rdata")

gplot=list()
load(file="alpha_04.Rdata")

df_error_p=data.frame(
  p=rep(sq_n,2),
  type=rep(factor(c("mis1","mis2")),each=10),
  error=as.vector(mt_error[mt_error[,2]==300,7:8])
)

gplot[[1]]=ggplot(data=df_error_p)+
  geom_line(aes(x=p,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=p,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=p,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1),label=c("false positive","false negative"))+
  scale_linetype_manual(values=c(1,1),label=c("false positive","false negative"))+
  scale_x_continuous(name="p",breaks=seq(250,1000,250))+
  scale_y_continuous(name= "Error Rate")+
  coord_cartesian(ylim=c(0,0.2),xlim=c(80,1050))+
  labs(caption="(a1)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        legend.title = element_blank(),
        legend.text = element_text(size=100),
        legend.key = element_rect(fill=NA),
        legend.key.size = unit(10,"line"),
        legend.key.width = unit(30, "line"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position=c(0.695,0.907))

df_error_T=data.frame(
  Ti=rep(sq_T,2),
  type=rep(factor(c("mis1","mis2")),each=10),
  error=as.vector(mt_error[mt_error[,1]==600,7:8])
)

gplot[[2]]=ggplot(data=df_error_T)+
  geom_line(aes(x=Ti,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=Ti,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=Ti,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1),label=c("false positive","false negative"))+
  scale_linetype_manual(values=c(1,1),label=c("false positive","false negative"))+
  scale_x_continuous(name="Ti")+
  scale_y_continuous(name= "Error Rate")+
  coord_cartesian(ylim=c(0,0.2))+
  labs(caption="(a2)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        legend.position="none")

load(file="alpha_05.Rdata")

df_error_p=data.frame(
  p=rep(sq_n,2),
  type=rep(factor(c("mis1","mis2")),each=10),
  error=as.vector(mt_error[mt_error[,2]==300,7:8])
)

gplot[[3]]=ggplot(data=df_error_p)+
  geom_line(aes(x=p,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=p,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=p,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1),label=c("false positive","false negative"))+
  scale_linetype_manual(values=c(1,1),label=c("false positive","false negative"))+
  scale_x_continuous(name="p",breaks=seq(250,1000,250))+
  scale_y_continuous(name= "Error Rate")+
  coord_cartesian(ylim=c(0,0.2),xlim=c(80,1050))+
  labs(caption="(b1)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        legend.position="none")

df_error_T=data.frame(
  Ti=rep(sq_T,2),
  type=rep(factor(c("mis1","mis2")),each=10),
  error=as.vector(mt_error[mt_error[,1]==600,7:8])
)

gplot[[4]]=ggplot(data=df_error_T)+
  geom_line(aes(x=Ti,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=Ti,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=Ti,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1),label=c("false positive","false negative"))+
  scale_linetype_manual(values=c(1,1),label=c("false positive","false negative"))+
  scale_x_continuous(name="Ti")+
  scale_y_continuous(name= "Error Rate")+
  coord_cartesian(ylim=c(0,0.2))+
  labs(caption="(b2)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        legend.position="none")

pdf(paste0("test1",".pdf"), width = 60, height = 60)
grid.arrange(grobs=gplot,layout_matrix=matrix(1:4,ncol=2,byrow=T))
dev.off()




gplot=list()
load(file="error_alpha_02.Rdata")

df_error_p=data.frame(
  p=rep(sq_n,4),
  error=as.vector(mt_error[mt_error[,2]==100,3:6]),
  type=factor(rep(c("cf","cl","af","al"),each=length(sq_T)),
              levels=c("cf","cl","af","al"))
)

gplot[[1]]=ggplot(data=df_error_p)+
  geom_line(aes(x=p,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=p,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=p,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1,2,5),label=c("common factors","common loading",
                                               "all factors","all loading"))+
  scale_linetype_manual(values=c(1,1,1,1),label=c("common factors","common loading",
                                                  "all factors","all loading"))+
  scale_x_continuous(name="p")+
  scale_y_continuous(name= "Mean Squared Eror")+
  coord_cartesian(ylim=c(0,0.4))+
  labs(caption="(a1)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        legend.title = element_blank(),
        legend.text = element_text(size=100),
        legend.key = element_rect(fill=NA),
        legend.key.size = unit(10,"line"),
        legend.key.width = unit(30, "line"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position=c(0.607,0.8225))

df_error_T=data.frame(
  Ti=rep(sq_T,4),
  error=as.vector(mt_error[mt_error[,1]==(unique(mt_error[,1])[3]),3:6]),
  type=factor(rep(c("cf","cl","af","al"),each=length(sq_n)),
              levels=c("cf","cl","af","al"))
)

gplot[[2]]=ggplot(data=df_error_T)+
  geom_line(aes(x=Ti,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=Ti,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=Ti,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1,2,5))+
  scale_linetype_manual(values=c(1,1,1,1))+
  scale_x_continuous(name="T")+
  scale_y_continuous(name= "Mean Squared Eror")+
  coord_cartesian(ylim=c(0,0.4))+
  labs(caption="(a2)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        #legend.title = element_text(size=100,face="italic"),
        #legend.text = element_text(size=100),
        #legend.key.size = unit(10,"line"))
        legend.position="none")


load(file="error_alpha_03.Rdata")

df_error_p=data.frame(
  p=rep(sq_n,4),
  error=as.vector(mt_error[mt_error[,2]==100,3:6]),
  type=factor(rep(c("cf","cl","af","al"),each=length(sq_T)),
              levels=c("cf","cl","af","al"))
)

gplot[[3]]=ggplot(data=df_error_p)+
  geom_line(aes(x=p,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=p,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=p,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1,2,5),label=c("common factors","common loading",
                                               "all factors","all loading"))+
  scale_linetype_manual(values=c(1,1,1,1),label=c("common factors","common loading",
                                                  "all factors","all loading"))+
  scale_x_continuous(name="p")+
  scale_y_continuous(name= "Mean Squared Eror")+
  coord_cartesian(ylim=c(0,0.4))+
  labs(caption="(b1)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(color="white",size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        #legend.title = element_text(size=100,face="italic"),
        #legend.text = element_text(size=100),
        #legend.key.size = unit(10,"line"))
        legend.position="none")


df_error_T=data.frame(
  Ti=rep(sq_T,4),
  error=as.vector(mt_error[mt_error[,1]==(unique(mt_error[,1])[3]),3:6]),
  type=factor(rep(c("cf","cl","af","al"),each=length(sq_n)),
              levels=c("cf","cl","af","al"))
)

gplot[[4]]=ggplot(data=df_error_T)+
  geom_line(aes(x=Ti,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=Ti,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=Ti,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1,2,5))+
  scale_linetype_manual(values=c(1,1,1,1))+
  scale_x_continuous(name="T")+
  scale_y_continuous(name= "Mean Squared Eror")+
  coord_cartesian(ylim=c(0,0.4))+
  labs(caption="(b2)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(color="white",size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        #legend.title = element_text(size=100,face="italic"),
        #legend.text = element_text(size=100),
        #legend.key.size = unit(10,"line"))
        legend.position="none")

load(file="error_alpha_04.Rdata")

df_error_p=data.frame(
  p=rep(sq_n,4),
  error=as.vector(mt_error[mt_error[,2]==100,3:6]),
  type=factor(rep(c("cf","cl","af","al"),each=length(sq_T)),
              levels=c("cf","cl","af","al"))
)

gplot[[5]]=ggplot(data=df_error_p)+
  geom_line(aes(x=p,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=p,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=p,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1,2,5),label=c("common factors","common loading",
                                               "all factors","all loading"))+
  scale_linetype_manual(values=c(1,1,1,1),label=c("common factors","common loading",
                                                  "all factors","all loading"))+
  scale_x_continuous(name="p")+
  scale_y_continuous(name= "Mean Squared Eror")+
  coord_cartesian(ylim=c(0,0.4))+
  labs(caption="(c1)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(color="white",size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        #legend.title = element_text(size=100,face="italic"),
        #legend.text = element_text(size=100),
        #legend.key.size = unit(10,"line"))
        legend.position="none")


df_error_T=data.frame(
  Ti=rep(sq_T,4),
  error=as.vector(mt_error[mt_error[,1]==(unique(mt_error[,1])[3]),3:6]),
  type=factor(rep(c("cf","cl","af","al"),each=length(sq_n)),
              levels=c("cf","cl","af","al"))
)

gplot[[6]]=ggplot(data=df_error_T)+
  geom_line(aes(x=Ti,y=error,group=type,lty=type),size=4)+
  geom_point(aes(x=Ti,y=error,group=type),colour="white",size=40)+
  geom_point(aes(x=Ti,y=error,group=type,pch=type),color="black",size=25)+
  scale_shape_manual(values=c(0,1,2,5))+
  scale_linetype_manual(values=c(1,1,1,1))+
  scale_x_continuous(name="T")+
  scale_y_continuous(name= "Mean Squared Eror")+
  coord_cartesian(ylim=c(0,0.4))+
  labs(caption="(c2)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=100,margin=margin(t=50,r=0,b=20,l=0)),
        axis.text.y = element_text(size=100,margin=margin(t=0,r=20,b=0,l=20)),
        axis.title.x = element_text(size=100,vjust=-1,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        axis.title.y = element_text(color="white",size=100,margin=margin(t=0,r=50,b=0,l=50)),
        #panel.grid.major =element_blank(),
        #panel.grid.minor =element_blank(),
        panel.background =  element_rect(fill="white"),
        panel.border =element_rect(fill=NA),
        plot.caption = element_text(size=100,hjust=0.5,
                                    margin=margin(t=50,r=0,b=50,l=0)),
        #legend.title = element_text(size=100,face="italic"),
        #legend.text = element_text(size=100),
        #legend.key.size = unit(10,"line"))
        legend.position="none")

pdf(paste0("test",".pdf"), width = 80, height = 60)
grid.arrange(grobs=gplot,layout_matrix=matrix(1:6,ncol=3,byrow=F))
dev.off()



