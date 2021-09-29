library(cplm)
library(TMB)
library(glmmTMB)
library(lme4)
library(mnormt)
# 
# data structure
Nsubject<-300
robs<-7
DF <- data.frame(id = rep(seq_len(Nsubject), each = robs),
                 bmi = rep(gl(2, Nsubject/2, labels = c("normal", "obese")), each = robs))
# design matrices for the fixed and random effects non-zero part
X <- model.matrix(~ bmi, data = DF)
Z <- model.matrix(~ 1, data = DF)
# effects on sleep time
thetas<-c(5.0,-0.1)
# effects on number of wake-ups
deltas<-c(-2.5,0.2)
# effects on waking up time
betas<-c(0.5,0.15)
# all true coefficients
truepara<-c(deltas[2],betas[2],thetas[2],deltas[2]+betas[2],deltas[2]+thetas[2],deltas[2]+betas[2]+thetas[2])
# random effects
rho12<-0.2
rho13<--0.4 
rho23<--0.20
sigma1<-0.35
sigma2<-0.30
sigma3<-0.18
sigmapara<-c(sigma1,sigma2,sigma3,sqrt(sigma1^2+sigma2^2+2*rho12*sigma1*sigma2),
sqrt(sigma1^2+sigma3^2+2*rho13*sigma1*sigma3),sqrt(sigma1^2+sigma2^2+sigma3^2
+2*rho12*sigma1*sigma2+2*rho13*sigma1*sigma3+2*rho23*sigma2*sigma3))
# covariance matrix for three random effects
covmatrix<-matrix(c(sigma1^2,rho12*sigma1*sigma2,sigma1*sigma3*rho13, 
rho12*sigma1*sigma2,sigma2^2,sigma2*sigma3*rho23, 
sigma1*sigma3*rho13, sigma2*sigma3*rho23,sigma3^2),nrow=3,ncol=3)
corpara<-c(rho12,rho13,rho23)
#generate data
set.seed(123)
# generate the data
# total sleep time
# random effects
b<-rmnorm(Nsubject,mean=c(0,0,0),covmatrix)
#
# generate total sleep time
DF$Sleep<-rgamma(Nsubject*robs,shape=2.0,scale=exp(as.vector(X%*%thetas+rowSums(Z * b[DF$id,3,drop=FALSE]))))
#
# mean on number of wake ups conditional on total sleep time
eta_m<-as.vector(log(DF$Sleep)+X %*% deltas + rowSums(Z * b[DF$id, 1, drop = FALSE]))
# number of wake ups per day 
DF$Mobs<-rpois(Nsubject * robs, exp(eta_m))
#
# total number of wake time per day
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, 2, drop = FALSE]))
DF$y <- rgamma(Nsubject * robs, shape=(2.0*DF$Mobs),scale=exp(eta_y))
#
# regression models
# model on number of wake bouts with sleep time as offset
result1<-glmmTMB(Mobs ~as.factor(bmi)  +(1 | id),family = poisson(link="log"), offset=log(Sleep), data = DF)
#
# model on number of wake bouts without sleep time as offset (v1+v3)
result1c<-glmmTMB(Mobs ~as.factor(bmi)  +(1 | id),family = nbinom2(link="log"),   data = DF)
#
# model on wake time conditional on wake bouts (number of wake bout>0)
result2<-glmmTMB(y ~ as.factor(bmi)  +(1 | id),offset=log(Mobs), weights=Mobs,family=Gamma(link="log"),data=subset(DF,DF$Mobs>0)) 
#
# model on total sleep time
#result3<-lmer(log(Sleep) ~as.factor(bmi)  +(1 | id), data = DF)
result3<-glmmTMB(Sleep ~as.factor(bmi)  +(1 | id), family=Gamma(link="log"),data = DF)
#
# model on marginal wake time using mixed compound Poisson model
result4<-cpglmm(y ~as.factor(bmi) +(1 | id), link="log",offset=log(Sleep), data = DF) 
result4c<-cpglmm(y ~as.factor(bmi) +(1 | id), link="log", data = DF) 
#
#
# random effects
# note model 1 wake time random effect standard deviation was divided by its overdispersion parameters
sigmaest<-c(sqrt(VarCorr(result1)$cond$id),sqrt(VarCorr(result2)$cond$id),sqrt(VarCorr(result3)$cond$id),
as.numeric(summary(result4)$REmat[1,4]),sqrt(VarCorr(result1c)$cond$id),as.numeric(summary(result4c)$REmat[1,4]))
#
# correlation between random effects
rho12<-(sigmaest[4]^2-sigmaest[1]^2-sigmaest[2]^2)/(2*sigmaest[1]*sigmaest[2])
rho13<-(sigmaest[5]^2-sigmaest[1]^2-sigmaest[3]^2)/(2*sigmaest[1]*sigmaest[3])
rho23<-(sigmaest[6]^2-sigmaest[4]^2-sigmaest[5]^2+sigmaest[1]^2)/(2*sigmaest[2]*sigmaest[3])
#
# bootstrap CI for random effects parameters
nboot<-300
sigma_boot<-matrix(0,nboot,6)
rho_boot<-matrix(0,nboot,3)
for (k in 1:nboot) {
newid<-sample(unique(DF$id),replace=TRUE)
tempdata<-do.call(rbind,lapply(newid,function(x) DF[DF$id==x,]))
tempdata$id<- rep(seq_len(Nsubject),each=robs)
# model on number of wake bouts with sleep time as offset
result1<-glmmTMB(Mobs ~as.factor(bmi)  +(1 | id),family = poisson(link="log"), offset=log(Sleep), data = tempdata)
sigma_boot[k,1]<-sqrt(VarCorr(result1)$cond$id)
# model on number of wake bouts without sleep time as offset (v1+v3)
result1c<-glmmTMB(Mobs ~as.factor(bmi)  +(1 | id),family = nbinom2(link="log"),   data = tempdata)
sigma_boot[k,5]<-sqrt(VarCorr(result1c)$cond$id)
#
# model on wake time conditional on wake bouts (number of wake bout>0)
result2<-glmmTMB(y ~ as.factor(bmi)  +(1 | id),offset=log(Mobs), weights=Mobs,family=Gamma(link="log"),data=subset(tempdata,tempdata$Mobs>0)) 
sigma_boot[k,2]<-sqrt(VarCorr(result2)$cond$id)/exp(summary(result2)$sigma^2)
#
# model on total sleep time 
result3<-glmmTMB(Sleep ~as.factor(bmi)  +(1 | id), family=Gamma(link="log"),data = tempdata)
sigma_boot[k,3]<-sqrt(VarCorr(result3)$cond$id) 
#
# model on total number of wake time marginally as a mixed effect compound Poisson distribution conditioning on total sleep time
result4<-cpglmm(y ~as.factor(bmi) +(1 | id), link="log",offset=log(Sleep), data = tempdata) 
# model on total number of wake time marginally as a mixed effect compound Poisson distribution conditioning on total sleep time
result4c<-cpglmm(y ~as.factor(bmi) +(1 | id), link="log", data = tempdata) 
sigma_boot[k,4]<-as.numeric(summary(result4)$REmat[1,4])
sigma_boot[k,6]<-as.numeric(summary(result4c)$REmat[1,4])
#
rho_boot[k,1]<-(sigma_boot[k,4]^2-sigma_boot[k,1]^2-sigma_boot[k,2]^2)/(2*sigma_boot[k,1]*sigma_boot[k,2])
rho_boot[k,2]<-(sigma_boot[k,5]^2-sigma_boot[k,1]^2-sigma_boot[k,3]^2)/(2*sigma_boot[k,1]*sigma_boot[k,3])
rho_boot[k,3]<-(sigma_boot[k,6]^2-sigma_boot[k,4]^2-sigma_boot[k,5]^2+sigma_boot[k,1]^2)/(2*sigma_boot[k,2]*sigma_boot[k,3])
                     }
# end of bootstrap
# bootstap CI
bootCI(sigma_boot)
bootCI(rho_boot)
#
# fuction bootstrap CI
bootCI<-function(boot) {
bootCI<-apply(boot,2,function(x) quantile(x,probs=c(0.025,0.975)))
bootCI
                        }

 

