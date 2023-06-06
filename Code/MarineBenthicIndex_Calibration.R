####This code reads in training data from the PSSMP and fits taxa models to environmental variables####
#this version separates the signal between fines and TOC

#load packages
library(glmnet)
library(lavaan)
library(reshape2)


##load data##

#This is be done automatically by starting RStudio by double clicking the Puget_Sound_MBI.Rproj
#otherwise uncomment, update and run the following 2 lines.
path<-"your/path/to/marine-benthic-index-master/"
setwd(path)


#load species data#
df<-read.csv("Data/habitat and TOC and species abundances for Long-term 2017 and 2018 samples.csv")
#grab species data
names(df)
#Look at column names for data set and find where the first taxon data column is
sp_start<-14
#grab just the taxa and turn them into pres/abs
sp_dat<-df[,sp_start:dim(df)[2]]
names(sp_dat)
sp_dat[sp_dat>0]<-1
#check to make sure
apply(sp_dat,2,max)

####assemble training data set###
#set weights for replicated stations within year
who<-which(table(df$Station,df$Year)[,1]>1)
train_weights<-rep(1,dim(df)[1])
train_weights[df$Station%in%names(who)&df$Year==2017]<-1/3

#check them
data.frame(df[,1:4],train_weights)

#remove any taxa with 10 or fewer occurrences across all observations
names(df)
dim(df)
occs<-apply(sp_dat,2,sum)
ubiq<-names(occs[occs/dim(sp_dat)[1]>=.3])
max(occs/dim(sp_dat)[1])
tail(occs[order(occs)])
lost_em<-which(occs<=10)
#look at names of rare taxa
length(lost_em)
names(occs)[lost_em]
#remove them from the data.frame
sp_dat<-sp_dat[,-(which(occs<=10))]

dim(sp_dat)
#get the names of those remaining
sp_nms<-names(sp_dat)

#look at the column names of df and grab the id vars and all the var we want in E
#in this case these are in columns 1-12
df_test<-data.frame(df[,1:12],sp_dat)

#linearize the relationshiop between sed_TOC and Fines
#plot(log((df_test$Fines/100)/(1-df_test$Fines/100))~sqrt(df$Sed_TOC))
#cor(log((df_test$Fines/100)/(1-df_test$Fines/100)),sqrt(df$Sed_TOC))
f<-lm((log((df_test$Fines/100)/(1-df_test$Fines/100))~sqrt(df$Sed_TOC)),weights = train_weights)

#use this to replace fines in the E data
df$Fines<-resid(f)


#re-scale continuous environmental var to have zero mean
#grab the vars for E
evars<-c("Depth","Penetration","Salinity","Temperature","Fines","Gravel")
env_vars<-df[,evars]
#center them and keep the value used to center each

col_means_mod_vars<-apply(env_vars,2,sum)*1/sum(train_weights)
col_sd_mod_vars<-apply(env_vars,2,sd)
env_vars<-sweep(env_vars,2,col_means_mod_vars,FUN = "-")
#env_vars<-sweep(env_vars,2,col_sd_mod_vars,FUN = "/")
apply(env_vars,2,mean)


df_test<-data.frame(df[,"Sample"],env_vars,sp_dat)
head(df_test)

# fit taxa models using lasso regularization and cross-validation to select "lambda" parameter
fts<-list()
lambda_min<-rep(NA,length(sp_nms))
names(lambda_min)<-sp_nms

for(i in sp_nms){
  fmls<-paste0(i,paste0("~",paste0(names(env_vars),collapse = "+")))
  fmls<-paste0(fmls,'+Salinity:Temperature')
  x<-model.matrix(formula(fmls),data = df_test)[,-1]
  y<-df_test[,i]
  tryCatch({cv.lasso<-glmnet::cv.glmnet(x,y,alpha=1,nfolds = dim(x)[1], family="binomial",grouped = F,weights = train_weights)},warning=function(w)print(i))
  lambda_min[i]<-cv.lasso$lambda.min
  tryCatch({fts[[i]]<-glmnet::glmnet(x,y,alpha=1,family = "binomial", lambda = cv.lasso$lambda,type.logistic = "modified.Newton",weights = train_weights)},warning=function(w)print(i))
}

#calc percent deviance
pct_deviance<-sapply(sp_nms,function(x)fts[[x]]$dev.ratio[fts[[x]]$lambda==lambda_min[x]])
#look at the distribution of deviance explained

hist(pct_deviance,xlab="Percent Devianace Explained",main="")
abline(v=.2,lwd=2,lty=2)

#grap those with pct_deviance .20 plus some special considerations
nms<-c(names(pct_deviance[pct_deviance>=.2]),ubiq,
c("Leitoscoloplos_pugettensis",
"Paraprionospio",
"Micronephthys_cornuta",
"Prionospio_lighti"))
#colony_spp%in%nms

#look at the distribution of deviance explained for this subset of taxa
hist(pct_deviance[nms])
#how many did we end up with?
nms<-unique(nms)
length(nms)

#cleanup column name for first column
names(df_test)[1]<-"Sample"

#calculate quantile residuals#
#this calculates the residuals as related to Dunn & Smyth (1996)
#if we think of the prediction p(y=1|x) as dividing the interval [0,1], then a
#and b are used to select the sub-interval based on whether the event was observed
#or not. if Y_obs=0, the a=0 and b= p(y=1|x), if Y_obs=1 the a=p(y=1|x) and b=1
# the residual is then calculated as the probability int_x(pdf(norm(0,1),x))|(0,mid), where mid is the midpoint of the 
#interval set by a and b. This residual does not have the weird continuity issues of the other log_reg residuals

res<-matrix(NA,nrow=dim(df)[1],ncol=length(nms))
colnames(res)<-nms
rownames(res)<-df_test$Sample
for(i in nms){
  a<-pbinom(df_test[,i]-1,size = 1,predict(fts[[i]],lambda_min[i],newx=x,type="respo"))
  b<-dbinom(df_test[,i],size = 1,predict(fts[[i]],lambda_min[i],newx=x,type="respo"))
  res[,i]<-qlogis((2*a+b)/2,0,1)
}

####set up SEM to calculate D and alphas####
#check our the distribution of residuals if any ar crazy high (or low) with
#might deal with those separately to preserve our assumption of species correlation 
# only through common effects of E
res_cor<-cor(res)
hist(res_cor)
#set up the model of lavaan
#we will set the var(D)=1 to allow identification
model<-paste("D=~NA*",paste(colnames(res),sep=" ",collapse ="+"),"\nD~~1*D",sep="")

#fit the model
fit<-sem(model,data = res,meanstructure = FALSE)

#grab the alphas
#get the estimates
alpha_est<-fit@ParTable$est[which(fit@ParTable$op=="=~")]

#get the se 
alpha_se<-fit@ParTable$se[which(fit@ParTable$op=="=~")]
nsp<-length(nms)
#grab the alphas
alphas<-data.frame(names=fit@ParTable$rhs[which(fit@ParTable$op=="=~")],
                   est=alpha_est,se=alpha_se)
#look at them
alphas

#calculate D using maximum likelihood 
wts<-1/alphas$se^2/min(1/alphas$se^2)

D_ests<-D_ests_new<-list()
len<-dim(df_test)[1]
for(i in 1:len){
  ll<-function(x)sum(-log(dlogis(res[i,],location = alphas$est*x[1],scale = x[2])),na.rm = T)
  f_tem<-optim(c(0,.3),ll,hessian = T)
  D_ests[[i]]<-c(D=f_tem$par[1],se=sqrt(diag(solve(f_tem$hessian)))[1])
  }

#collect the estimates into a single data.frame
d_ans<-data.frame(Sample=rownames(res),est=sapply(D_ests,"[",1),se=sapply(D_ests,"[",2))

#stuff needed for updating test sets

#write out the coefs estimates
sp_coefs<-sapply(nms,function(x)as.vector(coef(fts[[x]],s=lambda_min[x])))
rownames(sp_coefs)<-(coef(fts[[1]])@Dimnames[[1]])

write.csv(x = data.frame(Parameter=rownames(sp_coefs),sp_coefs),file = "Data/Model/TaxaCoefficients.csv",row.names = F)

#alphas
write.csv(alphas,"Data/Model/alphas_calibration.csv",row.names = F)

#Means for the carbonation continuous E var
col_mod_vars<-c(col_means_mod_vars[c("Depth","Penetration","Salinity","Temperature")],Fines=coef(f),col_means_mod_vars["Gravel"])
col_mod_vars<-as.matrix(t(col_mod_vars))

write.csv(col_mod_vars,"Data/Model/E_calibration_means.csv",row.names = F)

#prior to output folder and model data/model folder for update
write.csv(d_ans,"Data/Model/D_ScoresPriors.csv",row.names = F)

write.csv(d_ans,paste0("Output/D_Scores",Sys.Date(),".csv"),row.names = F)

#create example figures and toss it in the Output folder
nm_bits<-strsplit(d_ans$Sample,"_")
nm_bits
d_ans_short_name<-paste(sapply(nm_bits,"[",3),sapply(nm_bits,"[",2),sapply(nm_bits,"[",4),sep="_")
d_ans_short_name
d_ans$Sample<-d_ans_short_name
source("Code/Fig_makers.R")
pdf(paste0("Output/D_Scores",Sys.Date(),".pdf"))
dscore_mkr(d_ans)
dev.off()
pdf(paste0("Output/alpha_values",Sys.Date(),".pdf"))
alpha_mkr(alphas)
dev.off()

