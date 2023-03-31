#This code updates D scores and D index for test data
#fill in your path and the file names of the habitat and benthos data
#the updated D scores, with and without prior will appear in the /Output folder

#set path to folders
#This is be done automatically by starting RStudio by double clicking the Puget_Sound_MBI.Rproj
#otherwise uncomment, update and run the following 2 lines.
#path<-"your/path/to/Puget_Sound_MBI/"
#setwd(path)


#read in Habitat Data make sure the file names match those
# in the /Data folder 
hab_dat<-read.csv("Data/habitat data for Long-term 2019.csv")
#read in Benthos Data
benthos_dat<-read.csv("Data/Long-term 2019 abundances standardized 129.csv")

##############HIGHLIGHT THE REST OF THE CODE AND RUN IT AS ONE CHUNK ###################
#read in model values 
scale_vars<-read.csv("Data/Model/E_calibration_means.csv")
alphas<-read.csv("Data/Model/alphas_calibration.csv")
sp_coefs<-read.csv("Data/Model/TaxaCoefficients.csv")
priors<-read.csv("Data/Model/D_ScoresPriors.csv")
nms<-alphas$names

####read in new data###
#turn to pres/absence
benthos_dat_Sample<-benthos_dat$Sample
benthos_dat[benthos_dat>0]<-1
benthos_dat<-data.frame(Sample=benthos_dat_Sample,benthos_dat)


#merge them
df_test<-merge(hab_dat,benthos_dat,by="Sample")
names(df_test)
#transform Fines data to remove the effect of TOC
logit_trans<-function(x)log((x/100)/(1-x/100))
df_test$Fines<-logit_trans(df_test$Fines)-(as.numeric(scale_vars[5])+as.numeric(scale_vars[6])*sqrt(df_test$TOC_Sed))


###scale test data by calibration data mean and sd
mod_vars<-df_test[,c("Depth","Penetration","Salinity","Temperature","Fines","Gravel","InorgC")]
mod_vars<-sweep(mod_vars,MARGIN = 2,STATS = as.numeric(c(scale_vars[1:4],0,scale_vars[7:8])),FUN = "-")

names(mod_vars)<-paste0(names(mod_vars),"_x")
df_test<-data.frame(df_test,mod_vars)
#######
n_sites<-dim(df_test)[1]
res_new<-matrix(NA,nrow=n_sites,ncol=length(nms))
colnames(res_new)<-nms
rownames(res_new)<-df_test$Sample


#setup data for calc
i<-nms[1]
fmls<-paste0(i,'~Depth_x+Penetration_x+Salinity_x+Temperature_x+Fines_x+Gravel_x+InorgC_x+Salinity_x:Temperature_x')
x<-model.matrix(formula(fmls),data = df_test)[,-1]

#test data structure
all(colnames(x)==sp_coefs[-1,1])
#calculate residuals
for(i in nms){
  pdct<-1/(1+exp(-cbind(rep(1,dim(df_test)[1]),x)%*%sp_coefs[,i]))
  pdct[which(pdct>9.999999e-01)]<-9.999999e-01
  a<-pbinom(df_test[,i]-1,size = 1,pdct)
  b<-dbinom(df_test[,i],size = 1,pdct)
  tmp_scr<-(2*a+b)/2
  res_new[,i]<-qnorm(tmp_scr,0,1)
  
}



D_ests_new<-list()
for(i in 1:n_sites){
  ll<-function(x)sum(-log(dnorm(res_new[i,],alphas$int+alphas$est*x[1],x[2])),na.rm = T)
  f_tem<-optim(c(0,.3),ll,hessian = T)
  D_ests_new[[i]]<-c(D=f_tem$par[1],se=sqrt(diag(solve(f_tem$hessian)))[1])
}

d_new<-data.frame(Sample=rownames(res_new),est=sapply(D_ests_new,"[",1),se=sapply(D_ests_new,"[",2))
nm_bits<-strsplit(d_new$Sample,"_")
d_new_short_name<-paste(sapply(nm_bits,"[",1),sapply(nm_bits,"[",3),sapply(nm_bits,"[",4),sep="_")
d_new_long_names<-d_new$Sample
#write out results with no prior
write.csv(d_new,paste0("Output/D_Score_update_no_prior",Sys.Date(),".csv"),row.names = F)

######stop here if no prior is being used#####
#priors
nm_bits<-strsplit(priors$Sample,"_")
prior_short_name<-paste(sapply(nm_bits,"[",1),sapply(nm_bits,"[",3),sapply(nm_bits,"[",4),sep="_")
t1<-tapply(priors$est,INDEX = prior_short_name,FUN = mean)
t2<-tapply(priors$est,INDEX = prior_short_name,FUN = function(x)sqrt(sum(x^2)/length(x)))
new_prior=data.frame(Sample=names(t1),est=t1,se=t2)

#update based on priors

new_dat = d_new
old_dat = priors
update_prior<-function(old_dat,new_dat){
  who<-new_dat$Sample%in%old_dat$Sample
  cf<-merge(old_dat,new_dat,by="Sample",all.x = T,all.y = T)
  new_dat_update<-data.frame(Sample=cf$Sample,est=NA,se=NA)
  who2<-is.na(cf$est.x)
  who3<-is.na(cf$est.y)
  #copy over new not in old
  new_dat_update[who2,c("est","se")]<-cf[who2,c("est.y","se.y")]
  #copy over old not in new
  new_dat_update[who3,c("est","se")]<-cf[who3,c("est.x","se.x")]
  #update the rest
  !(who2|who3)
  new_dat_update[!(who2|who3),"est"]<-(cf$est.x[!(who2|who3)]*cf$se.y[!(who2|who3)]^2+cf$est.y[!(who2|who3)]*cf$se.x[!(who2|who3)]^2)/(cf$se.x[!(who2|who3)]^2+cf$se.y[!(who2|who3)]^2)
  new_dat_update[!(who2|who3),"se"]<-sqrt((cf$se.y[!(who2|who3)]^2*cf$se.x[!(who2|who3)]^2)/(cf$se.y[!(who2|who3)]^2+cf$se.x[!(who2|who3)]^2))
  return(new_dat_update)}

d_new$Sample<-d_new_short_name
dat_out<-update_prior(new_dat = d_new,old_dat = new_prior)

#write out only estimates from latest input
d_out<-dat_out[dat_out$Sample%in%d_new$Sample,]
d_out$Sample<-d_new_long_names
write.csv(d_out,paste0("Output/D_Score_update",Sys.Date(),".csv"),row.names = F)
