dindex_calc<-function(data){
  dind<-pnorm(data$est,0,1)
  upper<-pnorm(data$est+1.96*data$se)
  lower<-pnorm(data$est-1.96*data$se)
  return(data.frame(Sample=data$Sample,est=dind,upper_ci=upper,lower_ci=lower))
}

dscore_mkr<-function(data,main="",colr="black"){
  len<-dim(data)[1]
  ord=order(data$est)
  plot(data$est[ord],1:len,yaxt="n",ylab ="",xlab="",pch=21,bg="gray",xlim=c(-4,4),
       main=main,bty='n',col=colr)
  arrows(data$est[ord],1:len,data$est[ord]+1.96*data$se[ord],length = 0.1,angle = 90,col=colr)
  arrows(data$est[ord],1:len,data$est[ord]-1.96*data$se[ord],length = 0.1,angle = 90,col=colr)
  abline(v=0,lwd=2)
  points(data$est[ord],1:len,pch=21,bg=colr,cex=1.25)
  axis(side=2,1:len,labels = data$Sample[ord],las=T,cex.axis=.75)
  mtext("Station",side = 2,cex=1.5,padj=-3)
  mtext("Disturbance Score",side = 1,cex=1.5,padj=3)
}

dindex_mkr<-function(data,main=""){
  len<-dim(data)[1]
  ord=order(data$est)
  plot(pnorm(data$est[ord],0,1),1:len,yaxt="n",ylab ="",xlab="",pch=21,bg="gray",xlim=c(0,1),
       main=main,bty='n')
  arrows(pnorm(data$est[ord]),1:len,pnorm(data$est[ord]+1.96*data$se[ord]),length = 0.1,angle = 90)
  arrows(pnorm(data$est[ord]),1:len,pnorm(data$est[ord]-1.96*data$se[ord]),length = 0.1,angle = 90)
  abline(v=.5,lwd=2)
  points(pnorm(data$est[ord]),1:len,pch=21,bg="gray")
  axis(side=2,1:len,labels = data$Sample[ord],las=T,cex.axis=.5)
  mtext("Station",side = 2,cex=1.5,padj=-3)
  mtext("Disturbance Index",side = 1,cex=1.5,padj=3)
}

alpha_mkr<-function(data,main=""){
  nsp=dim(data)[1]
  ord=order(data$est)
  par(mar=c(5,8,2,2))
  plot(data$est[ord],1:nsp,yaxt="n",ylab = "",xlab="",pch=21,bg="gray",xlim=c(-.3,.4))
  arrows(data$est[ord],1:nsp,data$est[ord]+2*data$se[ord],length = 0.1,angle = 90)
  arrows(data$est[ord],1:nsp,data$est[ord]-2*data$se[ord],length = 0.1,angle = 90)
  abline(v=0,lwd=2)
  axis(side=2,1:nsp,labels = data$names[ord],las=T,cex.axis=.35)
  mtext("Taxon",side = 2,cex=1.5,padj=-6)
  mtext(bquote("Sensitivity Index ("*alpha*")"),side = 1,cex=1.5,padj=3)
}
  