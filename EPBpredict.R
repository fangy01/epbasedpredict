EPBpredict=function(emat,cls,scorecut=0.8,npairs=10,sepat=0,doplot=TRUE) {
# emat -- expression matrix
# cls -- annotate charactor factor or vector, present class label of samples. 
# scorecut -- criterion for discarding low score pairs.
# napirs-- top npairs will be evalated and further selected to built marker profile
# sepat -- cutoff point for prediction decision in correlation iffrence base prediction  
# doplot -- a logical value indicates whether to draw profile plot and resamplling results plot
if(!is.matrix(emat)) emat <- as.matrix(emat)
if(!is.numeric(emat)) stop("expression matrix should be numeric") 
if(any(is.null(names(cls)))) stop("cls elements are not annotated") 
if(any(!(names(cls) %in% colnames(emat)))) stop("Samples of expression matrix and samples of cls are not same set")
if(length(unlist(cls)) != ncol(emat)) stop("the size of samples in emat and in cls are different")

if(!is.factor(cls)) {
ucls=unique((cls))
cls=factor(cls,level=ucls)
}
if(length(levels(cls)) !=2) stop("categories diffined by cls is not 2, exact 2 catogories required")
idx1=which(cls==levels(cls)[1])
idx2=which(cls==levels(cls)[2])
clsorig=cls
cls=cls[c(idx1,idx2)]
idx=match(names(cls),colnames(emat))
mat=emat[,idx]
naidx=which(rowSums(is.na(mat) | is.infinite(mat)) >0)
if(length(naidx)!=0) mat=mat[-naidx,]
i1=which(cls==levels(cls)[1])
i2=which(cls==levels(cls)[2])
g1=mat[,i1]
g2=mat[,i2]
##################################################### scoring ------------------------------------------------------
sc=apply(mat,1,FUN=function(x,y,i1,i2) {
a=rep(1,nrow(y)) %o% x
yk=a[,i1];
zk=a[,i2];
out1=(rowSums(yk < y[,i1])/length(i1) + rowSums(zk > y[,i2])/length(i2)) / 2
out2=(rowSums(yk > y[,i1])/length(i1) + rowSums(zk < y[,i2])/length(i2)) / 2
out=pmax(out1,out2)
names(out) <- 1:nrow(y)
if(sum(out > scorecut)>0) out=out[out > scorecut] else out=NULL
out
},mat,i1,i2)

names(sc) <- 1:nrow(mat)
X=unlist(sc)
o=order(X,decreasing=T)
sc=X[o]
#####################################################  pick top pairs and zero centralisation ---------------------------
if(length(sc) == 0) {message(paste0("No gene pair pass the scorecut > ", scorecut, " threshold, please try lower threshold"))
prof=NULL
return(prof=NULL)
}
idxx=matrix(as.numeric(unlist(strsplit(names(sc),split="\\."))),nrow=2)   # get index of gene pairs
keeps=which(idxx[1,] < idxx[2,])
idxx=idxx[,keeps]
idxx=idxx[1:length(idxx)]
sco=sc[keeps]
selectedMat=mat[idxx,]
pair=rep(0,min(length(sco),npairs))
ip=length(pair)
g=mat[idxx[1:(2*ip)],]

gnorm=matrix(g,nrow=2)
gnorm=c(1,1) %o% colMeans(gnorm)
gnorm=matrix(gnorm,nrow=2*ip)

g=g-gnorm
m1=rowMeans(g[,i1])
m2=rowMeans(g[,i2])
prof0=list(g1=g[,i1],g2=g[,i2],mat=g,sc=sc[1:ip],cls=cls)

#print(sc)
###################################################### maximize predict accuracy by removing neative contribtion pairs

accuEstimate <- function(obj,cutat=sepat) {
# obj is the prof list generated within the function
obj$g1prof=rowMeans(obj$g1);
obj$g2prof=rowMeans(obj$g2);
obj$col=c(rep(1,ncol(obj$g1)),rep(3,ncol(obj$g2)));
obj$genes=rownames(obj$mat);
obj$group=levels(obj$cls);
obj$corr1= cor(obj$mat,matrix(obj$g1prof,ncol=1))
obj$corr2= cor(obj$mat,matrix(obj$g2prof,ncol=1))
obj$corrdiff=obj$corr1 - obj$corr2
#x11();par(mfcol=c(1,3)); plot(obj$corr1);plot(obj$corr2);plot(obj$corrdiff)
## leave one out cross validation  (LOOCV)
obj$corr=NULL;
for(k in 1:ncol(obj$g1)) {a=cor(obj$mat,matrix(rowMeans(obj$g1[,-k]),ncol=1))-obj$corr2;obj$corr=cbind(obj$corr,a)}
for(k in 1:ncol(obj$g2)) {a=obj$corr1-cor(obj$mat,matrix(rowMeans(obj$g2[,-k]),ncol=1));obj$corr=cbind(obj$corr,a)}
colnames(obj$corr) <- rownames(obj$corr) <- colnames(obj$mat)
##
predlab = rep(obj$group[2],length(obj$cls))
predlab[obj$corrdiff > cutat] = obj$group[1]
accuval= sum(as.character(predlab) == as.character(obj$cls))/length(obj$cls)
out=list(accuval=accuval,diff=obj$corrdiff,obj=obj)
}
gs=g
prof=prof0
while(nrow(gs) > 4) {
accu0 = accuEstimate(prof)
accu=NULL
for(k in 1:(nrow(gs)/2)) {
profk=list(g1=gs[-(1:2+2*(k-1)),i1],g2=gs[-(1:2+2*(k-1)),i2],mat=gs[-(1:2+2*(k-1)),],sc=prof$sc[-(1:2+2*(k-1))],cls=cls)
accu =c(accu,accuEstimate(profk)$accuval)
}
accudiff = accu - unlist(accu0$accuval);
#print(accudiff)
if(any(accudiff >= 0)) {
dis=which(accudiff == max(accudiff))
dis=dis[length(dis)]
gs=gs[-c(2*dis-1,2*dis),]
prof=list(g1=gs[,i1],g2=gs[,i2],mat=gs,sc=prof$sc[-dis],cls=cls)
} else break 
}


#print("ok2")
############################################## bootstrapping the prediction accuracy based on generated "prof"
boot.epbPredict <- function(ds,indices,cutat=0) {
d=ds[indices,];
cls=d[,1]
i1=which(cls %in% levels(cls)[1])
i2=which(cls %in% levels(cls)[2])
mat=t(d[,-1])
prof=list(g1=mat[,i1],g2=mat[,i2],mat=mat[,c(i1,i2)],sc=NA,cls=cls)
accu=accuEstimate(prof,cutat=cutat)$accuval
}

profmat=cbind(prof$g1,prof$g2)
df=data.frame(cls=cls,t(profmat)) 
strata=rep(1,length(cls))
strata[cls==levels(cls)[2]]=2

require(boot)
bootout=boot(df,statistic=boot.epbPredict,R=1000,strata=strata,cutat=sepat)
profileAndBootRes=list(profile=prof,bootres=bootout,outmat=mat,outcls=cls,score=sco,selectedMat=selectedMat)
#print("ok3")
##################################################### produce plot
plot.profile <- function(prof,cutat=0) {
# prof is a list of prof generated in the maon function
obj=accuEstimate(prof,cutat=sepat)$obj
x11()
par(mfcol=c(1,2))
par(mar=c(8,4.5,2,1))
matplot(obj$mat,type="l",col=obj$col,xlab="",xaxt="n",ylab=expression(log[2]~expression),main="");
lines(1:nrow(obj$g1),rowMeans(obj$g1),col=4,lwd=4)
lines(1:nrow(obj$g1),rowMeans(obj$g2),col=6,lwd=4)
axis(side=1,at=1:length(obj$genes),labels=obj$genes,las=2,cex.axis=0.8)
legd=c(paste0(obj$group,"s"),paste0(obj$group,"_profile"))
legend("topright",legend=legd,lty=1,lwd=c(1,1,4,4),col=c(1,3,4,6))
mtext(side=3,line=0.5,cex=1.4,text="A",adj=0)

plot(1:ncol(obj$mat),obj$corrdiff,col=obj$col,main="",xlab="",ylab="Correlation difference",xaxt="n");
abline(h=0,lty=2)
lines(c(1,ncol(obj$g1)),c(par()$usr[3],par()$usr[3]),lwd=8,col=1)
lines(c(ncol(obj$g1)+1,ncol(obj$mat)),c(par()$usr[3],par()$usr[3]),lwd=8,col=3)
abline(h=cutat,col="red")
at=c(ncol(obj$g1)/2,ncol(obj$g1)+ncol(obj$g2)/2)
axis(side=1,at=at,labels=obj$group,las=2,cex.axis=0.8)
#legend("bottomleft",pch=1,legend=obj$group,col=c(1,3))
#axis(side=1,at=1:length(obj$corrdiff),labels=rownames(obj$corrdiff),las=2,cex.axis=0.8)
mtext(side=3,line=0.5,cex=1.4,text="B",adj=0)
}

plot.bootout <- function(bootout) {
## bootout is a list object geberated by "boot" and it run within EPBpredict function
x=bootout
a0=x$t0
a=x$t
b=c(a,a0)
x11()
par(mfcol=c(1,2))
hist(a,n=31,xlab="Accuracy",ylab="Frequency",main="Histogram")
abline(v=mean(a),lty=2,col=3)
ci=mean(a)+c(-1,1)*qnorm(0.025,0,1)*sd(a)
lines(ci,c(0,0),lwd=6,col=3)
points(a0,0,cex=2,pch=4,col=2)
mtext(side=3,line=2,cex=1.4,text="A",adj=0)

xy=qqnorm(b,xlab="Standard Normal", ylab="Accuracy",main="qqplot",col=c(rep(1,length(a)),2),cex=c(rep(0.8,length(a)),2))
abline(lm(xy$y ~ xy$x),lty=1,col=2)
abline(h=mean(xy$y),col=3,lty=3)
lines(c(par()$usr[1],par()$usr[1]),ci,lwd=6,col=3)
mtext(side=3,line=2,cex=1.4,text="B",adj=0)

#print(intcept)
}

if(doplot){
if(!is.null(prof)) plot.profile(prof,cutat=sepat) else warning("profile does not exist, no figure generated2")
if(!is.null(bootout)) plot.bootout(bootout) else warning("bootout does not exist, no figure generated2")
}
profileAndBootRes

}

## end 





