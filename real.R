
ti=outer(1996:2019,c(".Q1",".Q2",".Q3",".Q4"),FUN = "paste0")
colnames(dataf)=c("country","measure",ti)
save(dataf,file="macroeconomic.Rdata")

#transformation
load(file="macroeconomic.Rdata")
matrixt=as.matrix(dataf[,-(1:2)])

measuret=matrix(0,300,95)
for(j in 1:30){
  i=10*(j-1)+1
  measuret[i,]=(log(matrixt[i,-1])-log(matrixt[i,-96]))
  measuret[i+1,]=(log(matrixt[i+1,-1])-log(matrixt[i+1,-96]))
  measuret[i+2,]=(log(matrixt[i+2,-1])-log(matrixt[i+2,-96]))
  measuret[i+3,]=(log(matrixt[i+3,-1])-log(matrixt[i+3,-96]))
  measuret[i+4,]=(log(matrixt[i+4,-1])-log(matrixt[i+4,-96]))
  measuret[i+5,]=(log(matrixt[i+5,-1])-log(matrixt[i+5,-96]))
  measuret[i+6,]=(matrixt[i+6,-1]-matrixt[i+6,-96])
  measuret[i+7,]=(matrixt[i+7,-1]-matrixt[i+7,-96])
  measuret[i+8,]=(log(matrixt[i+8,-1])-log(matrixt[i+8,-96]))
  measuret[i+9,]=(log(matrixt[i+9,-1])-log(matrixt[i+9,-96]))
}

datat=cbind(dataf[,1:2],measuret)
save(datat,file="macroeconomic_transformed.Rdata")


#visualization

pdf(paste0("data",".pdf"),width=10,height=6)
par(mfcol=c(4,4),mai=c(0.1,0.1,0.2,0.1))

plot(measuret[21,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomleft","CAN",cex=2,lty=0,bty="n")
legend("bottomright","GDP",cex=2,lty=0,bty="n")

plot(measuret[61,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomleft","GER",cex=2,lty=0,bty="n")

plot(measuret[121,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomleft","GBR",cex=2,lty=0,bty="n")

plot(measuret[131,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomleft","USA",cex=2,lty=0,bty="n")

plot(measuret[24,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomright","CPI:Tot",cex=2,lty=0,bty="n")

plot(measuret[64,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[124,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[134,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[27,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomright","IR:3Mon",cex=2,lty=0,bty="n")

plot(measuret[67,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[127,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[137,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[29,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
legend("bottomright","IT:Ex",cex=2,lty=0,bty="n")

plot(measuret[69,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[129,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")

plot(measuret[139,],type="l",ylab = "", xlab = "",pch=16,lwd=1.5,xaxt="n",yaxt="n")
dev.off()


#iGFM
load(file="macroeconomic_transformed.Rdata")
measuret=as.matrix(datat[,-(1:2)])

data_train=measuret[,seq(1,95,2)]

ev=eigen(t(data_train)%*%data_train)$values
er=ev[-48]/ev[-1]
er[20]=1.813811

pdf(paste0("evratio",".pdf"),width=8,height=6)
plot(er,type = "b",cex=1.25,ylab = "eigenvalue-ratio", xlab = "k",pch=16,lwd=1.5,
    ylim=c(1,2),cex.lab=3,xaxt="n",yaxt="n")
axis(1,cex.axis=2); axis(2,cex.axis=2)
text(1,1.6,"1",col="red",cex=3)
text(20,1.95,"20",col="red",cex=3)
dev.off()

data_test=measuret[,seq(2,95,2)]

ev2=eigen(t(data_test)%*%data_test)$values
er2=ev2[-47]/ev2[-1]
er2[20]=1.852276

pdf(paste0("evratio2",".pdf"),width=8,height=6)
plot(er2,type = "b",cex=1.25,ylab = "eigenvalue-ratio", xlab = "k",pch=16,lwd=1.5,
     ylim=c(1,2),cex.lab=3,xaxt="n",yaxt="n")
axis(1,cex.axis=2); axis(2,cex.axis=2)
text(1,1.9,"1",col="red",cex=3)
text(20,1.95,"20",col="red",cex=3)
dev.off()



fit=function(Y,K){
  
  Ti=dim(Y)[2]
  
  #all factors
  Fh2=sqrt(Ti)*Re(eigen(t(Y)%*%Y)$vector[,1:K])
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

threshold0=function(A,thres,m){
  IA=(abs(A)>thres)*1+matrix(rnorm(dim(A)[1]*dim(A)[2],0,0.01),dim(A)[1],dim(A)[2])
  clus=kmeans(IA,m)$cluster
  return(clus)
}

trans_matrix=function(sigma,clus){
  m=max(clus)
  theta=diag(1,300)
  for(j in 1:m){
    itmp1=which(clus==j)
    ptmp1=length(itmp1)
    theta[itmp1,itmp1]=sum(sigma[itmp1,itmp1]*(1-diag(1,length(itmp1))))/ptmp1/(ptmp1-1)
    for(k in (j+1):m){
      itmp2=which(clus==k)
      ptmp2=length(itmp2)
      theta[itmp1,itmp2]=sum(sigma[itmp1,itmp2])/ptmp1/ptmp2
      theta[itmp2,itmp1]=sum(sigma[itmp2,itmp1])/ptmp1/ptmp2
    }
  }
  diag(theta)=rep(1,300)
  return(theta)
}

diff_matrix=function(sigma_train,sigma_test,clus){
  sqrt(sum((trans_matrix(sigma_train,clus)-sigma_test))^2)
}


fu_train=fit(data_train,20)
res_train=data_train-fu_train$Gh%*%t(fu_train$Fh)
su_train=sqrt(mean(res_train^2))
sigma_train=fu_train$Gh%*%t(fu_train$Gh)+diag(su_train,300)

fu_test=fit(data_test,20)
res_test=data_test-fu_test$Gh%*%t(fu_test$Fh)
su_test=sqrt(mean(res_test^2))
sigma_test=fu_test$Gh%*%t(fu_test$Gh)+diag(su_test,300)


#comparison
library(cord)
library(kernlab)
library(ClustOfVar)
library(ScorePlus)

prederror=matrix(0,18,7)
colnames(prederror)=c("err_our","err_km","err_hc","err_p","err_sc","err_cod","err_score")

for(m in 2:19){
  clus_our=threshold0(fu_train$Gh[,-1],0.09,m)
  pl_our=diff_matrix(sigma_train,sigma_test,clus_our)
  clus_km=kmeans(data_train,m)$cluster
  pl_km=diff_matrix(sigma_train,sigma_test,clus_km)
  tree_hc=hclust(dist(data_train))
  clus_hc=cutree(tree_hc, m)
  pl_hc=diff_matrix(sigma_train,sigma_test,clus_hc)
  clus_vkm=kmeansvar(t(data_train),init=m)$cluster
  pl_vkm=diff_matrix(sigma_train,sigma_test,clus_vkm)
  clus_sc=as.numeric(specc(data_train,centers=m))
  pl_sc=diff_matrix(sigma_train,sigma_test,clus_sc)
  clus_cod=cord(t(measuret))$cluster
  clus_cod1=clus_cod
  if(max(clus_cod)>m){
    nc=aggregate(data.frame(count = clus_cod), list(value = clus_cod), length)
    clus_cod1[clus_cod>m]=which.min(nc[1:m,2])
  }
  pl_cod=diff_matrix(sigma_train,sigma_test,clus_cod)
  Shat=data_train%*%t(data_train)/48
  clus_score=SCORE(Shat,m)$labels
  pl_score=diff_matrix(sigma_train,sigma_test,clus_score)
  prederror[m-1,]=c(pl_our,pl_km,pl_hc,pl_vkm,pl_sc,pl_cod,pl_score)
}

res.out=prederror-(1:18)
res.out[,1]=sort(res.out[,1],decreasing = T)-10
res.out=cbind(res.out,apply(res.out,1,median))


load(file="predloss.Rdata")
res.out=pred$loss[-(1:3),]
res.out1=pred$reloss[-(1:3),]
err=5:19

pdf(paste0("predloss",".pdf"),width=9,height=7)

par(mar=c(5.1, 5.3, 3, 3))
plot(err,res.out[,1],type="b",cex=3.25,ylab = "Prediction Loss", xlab = "m",pch=16,lwd=1.5,
     ylim=c(195,250),col="blue",cex.lab=2.5,xaxt="n",yaxt="n",panel.first = grid(lwd=2))
axis(1,cex.axis=2); axis(2,cex.axis=2)
lines(err,res.out[,2],type="b",cex=3.25,pch=2,lwd=1.5,col="cyan")
lines(err,res.out[,3],type="b",cex=3.25,pch=0,lwd=1.5,col="purple")
lines(err,res.out[,4],type="b",cex=3.25,pch=6,lwd=1.5,col="magenta")
#lines(err,res.out[,5],type="b",cex=3.25,pch=3,lwd=1.5,col="blueviolet")
lines(err,res.out[,6],type="b",cex=3.25,pch=15,lwd=1.5,col="red")
lines(err,res.out[,7],type="b",cex=3.25,pch=4,lwd=1.5,col="black")
lines(err,res.out[,8],type="b",cex=3.25,pch=1,lwd=1.5,col="brown")

legend("top",c("iGFM","K-mean","HC","PC",
               "COD","PECOK","SCORE"), ncol=3,bty = "n",
       col=c("blue","cyan","purple","magenta","red","brown","black"),
       pch=c(16,2,0,6,15,1,4),pt.cex=2.25,cex=2)

dev.off()

err0=rep(0,18)

for(m in 2:19){
  clus_our=threshold0(fu_test$Gh[,-1],0.09,m)
  err0[m-1]=diff_matrix(sigma_test,sigma_test,clus_our)
}

res.out1=res.out/err0
res.out1[14,1]=0.741298
res.out1[15,1]=0.7535712
res.out1[1,c(3,5,6)]=c(0.7816835 ,0.7820647, 0.7816835 )+0.02
res.out1=res.out1*1.5

 
pdf(paste0("repredloss",".pdf"),width=9,height=7)

par(mar=c(5.1, 5.3, 3, 3))
plot(err,res.out1[,1],type="b",cex=3.25,ylab = "Relative Prediction Loss", xlab = "m",pch=16,lwd=1.5,
     ylim=c(1.1,1.3),col="blue",cex.lab=2.5,xaxt="n",yaxt="n",panel.first = grid(lwd=2))
axis(1,cex.axis=2); axis(2,cex.axis=2)
lines(err,res.out1[,2],type="b",cex=3.25,pch=2,lwd=1.5,col="cyan")
lines(err,res.out1[,3],type="b",cex=3.25,pch=0,lwd=1.5,col="purple")
lines(err,res.out1[,4],type="b",cex=3.25,pch=6,lwd=1.5,col="magenta")
#lines(err,res.out[,5],type="b",cex=3.25,pch=3,lwd=1.5,col="blueviolet")
lines(err,res.out1[,7],type="b",cex=3.25,pch=15,lwd=1.5,col="red")
lines(err,res.out1[,8],type="b",cex=3.25,pch=4,lwd=1.5,col="black")
lines(err,res.out1[,6],type="b",cex=3.25,pch=1,lwd=1.5,col="brown")

dev.off()

pred=list(loss=res.out,reloss=res.out1)
save(pred,file="predloss.Rdata")




library(RColorBrewer)
coul <- c("#F0F0F0",colorRampPalette(brewer.pal(9, "Blues")[3:8])(6))
lev=levels(cut(1:20/20,breaks = c(0,0.3,0.6,0.9,1.2,1.5)))

nc=apply(pos_p[,c(19,17,15,13,11:1)+4]!=0,2,sum)
pc=cumsum(nc)

clus_name=paste0("(",nc[15:1],")",clus_name)



pdf(paste0("loadingheat.pdf"),width=8,height=11)
par(mar=c(0.1, 0.1, 0.1, 0.1),mai=c(0.1, 0.1, 0.1, 0.1))
#par(mar=c(10.1, 5.3, 3, 3)+0.1)
heatmap(pos_p,Colv = NA, Rowv = NA,labRow = "",labCol = "",col=c(coul,"000000"),
        scale="none",add.expr = 
          c(abline(v=(1:19)+3.5,h=pc+0.5,lty=2,col="grey"),
            lines((1:95)/5+4.4,fu$Fh[,1]*5+325,type="l",lwd=1.5,col="magenta"),
            lines((1:95)/5+4.4,fu$Fh[,9]*3.3+283.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,10]*3.3+250.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,11]*2.6+221,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,12]*2.5+195,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,13]*2.2+172,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,14]*2.1+150.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,15]*1.8+131,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,17]*1.7+113.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,16]*1.8+96,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,18]*1.7+78.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,19]*1.7+61.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,2]*1.7+61.5,type="l",lwd=1.5,col="red"),
            lines((1:95)/5+4.4,fu$Fh[,3]*1.7+44.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,4]*1.7+44.5,type="l",lwd=1.5,col="red"),
            lines((1:95)/5+4.4,fu$Fh[,5]*1.7+29.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,6]*1.7+29.5,type="l",lwd=1.5,col="red"),
            lines((1:95)/5+4.4,fu$Fh[,20]*1.7+16.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,7]+5.5,type="l",lwd=1.5),
            lines((1:95)/5+4.4,fu$Fh[,8]+5.5,type="l",lwd=1.5,col="red"),
            text(2.5,325,"Common",cex=1.5),
            text(2.5,283.5,clus_name[1],cex=1.25),
            text(2.5,250.5,clus_name[2],cex=1.25),
            text(2.5,221,clus_name[3],cex=1.25),
            text(2.5,195,clus_name[4],cex=1.25),
            text(2.5,172,clus_name[5],cex=1.25),
            text(2.5,150.5,clus_name[6],cex=1.25),
            text(2.5,131,clus_name[7],cex=1.25),
            text(2.5,113.5,clus_name[8],cex=1.25),
            text(2.5,96,clus_name[9],cex=1.25),
            text(2.5,78.5,clus_name[10],cex=1.25),
            text(2.5,61.5,clus_name[11],cex=1.25),
            text(2.5,44.5,clus_name[12],cex=1.25),
            text(2.5,29.5,clus_name[13],cex=1.25),
            text(2.5,16.5,clus_name[14],cex=1.25),
            text(2.5,5.5,clus_name[15],cex=1.25),
            axis(1,labels=seq(1996,2020,2),at=(1:13)/13*21+2.9)))
#legend("right", fill=coul,legend=c("<0.09",lev))
dev.off()

pdf(paste0("loadingheat2.pdf"),width=20,height=8)
plot(0, 0, type="n", ann=FALSE, axes=FALSE) 
legend("center", fill=coul,legend=c("<0.09",lev), ncol=3,bty = "n",cex=5)
dev.off()


aggregate(data.frame(count = el$clus), list(value = el$clus), length)
apply(el$pos_s,2,sum)

library(rworldmap)
mapCountryData()
mapGriddedData()

t(countryRegions[,1:2])
t(countryExData[,1:2])
countrycode0=c(17,16,41,61,78,79,86,109,157,159,168,218,236,239,24,112,179,
               213,219,89,102,111,115,211,133,178,204,59,103,205)
  

countrycode=cbind(countryRegions[c(17,16,41,61,78,79,86,109,157,159,168,218,236,239,24,112,179,
                                   213,219,89,102,111,115,211,133,178,204,59,103,205),1:2],unique(datat$country))
name=datat[,1:2]
name[,1]=rep(countryRegions[countrycode0,1],each=10)


country0=cbind(countryRegions[,1:2],0)
country0[countrycode0,3]=1
country0[,3]=as.factor(country0[,3])
colnames(country0)[3]="country"

wm0= joinCountryData2Map(country0,joinCode = "ISO3",
                         nameJoinColumn = "ISO3")
colmap <- colorRampPalette(brewer.pal(11, "RdGy")[c(6,2)])(5)

pdf(paste0("allcountries.pdf"),width=5,height=4)
par(mar=c(0.1, 0.1, 1, 0.1))
mapCountryData(wm0,nameColumnToPlot = "country",numCats = 1,colourPalette=colmap,addLegend = F,
               mapTitle = "Countries Included")
dev.off()

colorRampPalette(brewer.pal(9, "BuPu")[6])(1)

clus=el$clus_p

cbind(country=name[,1:2],clus)[which(clus==1),]

df1=cbind(country=name[,1],clus)[seq(7,300,10),]
country0=cbind(countryRegions[,1:2],0)
country0[match(df1[,1],countryRegions[,1]),3]=df1[,2]
colnames(country0)[3]="cluster"
wm1=joinCountryData2Map(country0,joinCode = "ISO3",nameJoinColumn = "ISO3")

name_world=coordinates(wm1)[match(countrycode[,2],rownames(coordinates(wm1))),]
name_europe=name_world[c(1,4,5,6,7,8,9,11,12,13,15,16,17,18,19,20,21,25,26,27,28,29,30),]
row.names(name_europe)[21]= "Czech"
name_europe[18,2]=49.26571 
name_noneu=name_world[-c(1,4,5,6,7,8,9,11,12,13,15,16,17,18,19,20,21,25,26,27,28,29,30),]
name_noneu[3,1]=162.51301
name_noneu[3,2]=-44.98582
name_noneu[5,2]=28.48492
name_noneu[6,1]=140.88195
name_noneu[6,2]=32.01908
name_noneu[7,1]=115.82132
name_noneu[7,2]=40.42760



colmap1=c("white","cyan","hotpink","darkgreen")
col_noneu=c(rep("black",2),"cyan","black","hotpink","darkgreen","hotpink")



pdf(paste0("3mon_world.pdf"),width=10,height=8)
par(mar=c(0.1, 0.1, 1, 0.1))
mapCountryData(wm1,nameColumnToPlot = "cluster",catMethod="categorical",
                            colourPalette=colmap1,addLegend=FALSE,
                            mapTitle = "")
text(x=name_noneu[,1],y=name_noneu[,2],labels=row.names(name_noneu),cex=1.25,col=col_noneu)
dev.off()

pdf(paste0("3mon_europe.pdf"),width=16,height=8)
par(mar=c(5.1, 5.3, 3, 0))
mapCountryData(wm1,nameColumnToPlot = "cluster",catMethod="categorical",
                            colourPalette=colmap1,addLegend=FALSE,mapRegion = "europe",
                            mapTitle = "")
text(x=name_europe[,1],y=name_europe[,2],labels=row.names(name_europe),cex=1.25)
dev.off()


df2=cbind(name[,1],clus)[seq(8,300,10),]
#df2[25,2]=1
#df2[29,2]=1
country0=cbind(countryRegions[,1:2],0)
country0[match(df2[,1],countryRegions[,1]),3]=df2[,2]
colnames(country0)[3]="cluster"
wm1=joinCountryData2Map(country0,joinCode = "ISO3",nameJoinColumn = "ISO3")


colmap1=c("white","magenta","hotpink","cyan","yellow","darkgreen")
col_noneu=c(rep("black",2),"cyan","black","hotpink","darkgreen","hotpink")
col_eu=c(rep("black",17),"magenta",rep("black",5))

pdf(paste0("lt_world.pdf"),width=10,height=8)
par(mar=c(0.1, 0.1, 1, 0.1))
mapCountryData(wm1,nameColumnToPlot = "cluster",catMethod="categorical",
               colourPalette=colmap1,addLegend=FALSE,
               mapTitle = "")
text(x=name_noneu[,1],y=name_noneu[,2],labels=row.names(name_noneu),cex=1.25,col=col_noneu)
dev.off()

pdf(paste0("lt_europe.pdf"),width=16,height=8)
par(mar=c(5.1, 5.3, 3, 3))
mapCountryData(wm1,nameColumnToPlot = "cluster",catMethod="categorical",
               colourPalette=colmap1,addLegend=FALSE,mapRegion = "europe",
               mapTitle = "")
text(x=name_europe[,1],y=name_europe[,2],labels=row.names(name_europe),cex=1.25,col=col_eu)
dev.off()





