library(survival)
###
set.seed(123)
# set up data 
nsubject<-2000
#binary covariate
trt<-c(rep(0,nsubject/2),rep(1,nsubject/2))
# continuous variable
age<-rnorm(nsubject,0,1)
X<-cbind(rep(1,nsubject),trt,age)
#
# assuming X and Z are the same
Z<-X
# # probability of uncure baseline 
pbaseunc<-0.7
trueb<-c(log(pbaseunc/(1-pbaseunc)),-log(2.0),log(1.5))
shapepara<-2.0
truebeta<-c(0.05,-log(1.5),0.2)*shapepara
# proportion of observed among cured
pobscure<-0.8
#
# Step 1
# to generate probability of uncure
#
uncurepi<-exp(X%*%trueb)/(1+exp(X%*%trueb))
uncstatus<-rbinom(nsubject,1,uncurepi)
#
# step 2
# survival time
# assign all very large survival time 
stime<-rep(100,nsubject)
#
# to generate survival data for the uncure
scalepara<-1/exp(X%*%(truebeta/shapepara))   
# give time to event to uncured
stime[uncstatus==1]<-rweibull(sum(uncstatus),shape=shapepara,scale=scalepara[uncstatus==1])
#
# generate censoring time and censoring status
ctime<-runif(nsubject,0,6)
# 
eventstatus<-ifelse(stime<ctime,1,0)
#
# event subject observed status=1
obstime<-ifelse(eventstatus==0, ctime,stime)
#
# generate observed status 
obstatus<-rep(0,nsubject)
obstatus[uncstatus==0]<-rbinom(nsubject-sum(uncstatus),1,pobscure)
#
obstime[uncstatus==0 & obstatus==1]<-10 
obstatus[eventstatus==1]<-1
#
#
# start em algorithm
# intitial estimate
# w=1 if event or observed cure
w<-uncstatus
b<-glm(w~Z[,-1],family=binomial)$coef
beta<-coxph(Surv(obstime, eventstatus) ~ X[, -1] + offset(log(w)), 
            subset = w != 0, method = "breslow")$coef
offsetvar<-NULL
#
emfit<-emnewmod(obstime, eventstatus,obstatus, X, Z, offsetvar, b, beta,"logit", emmax = 50, eps = 1e-07) 
bmatrix<-emfit$b
betamat<-emfit$latencyfit
pobsvec<-emfit$Pobs 
#
# prepare for bootstrapping
#
b <- emfit$b
beta <- emfit$latencyfit
s <- emfit$Survival
logistfit <- emfit$logistfit
# 
# set up bootstrap parameters
nboot<-100
bnm <- colnames(Z)[-1]
nb <- ncol(Z)
betanm <- colnames(X)[-1]
nbeta <- ncol(X) - 1
eps<-1.0e-7
#
# start bootstrapping variance
b_boot <- matrix(rep(0, nboot * nb), nrow = nboot)
beta_boot <- matrix(rep(0, nboot * (nbeta)), nrow = nboot)
pcure_boot<-rep(0,nboot)
iter <- matrix(rep(0, nboot), ncol = 1)
        
tempdata <- cbind(obstime, eventstatus,obstatus, X)
        data1 <- subset(tempdata, eventstatus == 1)
        data0 <- subset(tempdata, eventstatus == 0)
        n1 <- nrow(data1)
        n0 <- nrow(data0)
        i <- 1
        while (i <= nboot) {
            id1 <- sample(1:n1, n1, replace = TRUE)
            id0 <- sample(1:n0, n0, replace = TRUE)
            bootdata <- rbind(data1[id1, ], data0[id0, ])
            bootZ <- as.matrix(cbind(rep(1, nsubject), bootdata[, bnm]))
            bootX <- as.matrix(cbind(rep(1, nsubject), bootdata[, betanm]))
            bootfit <- emnewmod(bootdata[, 1], bootdata[, 2], bootdata[,3],bootX, 
                bootZ, offsetvar, b, beta, "logit", emmax=50, eps=1e-08)
            b_boot[i, ] <- bootfit$b
            beta_boot[i, ] <- bootfit$latencyfit
            pcure_boot[i]<-bootfit$Pobs
            if (bootfit$tau < eps) 
            i <- i + 1
        }
# variance estimates for b, beta and pobs
        b_var <- apply(b_boot, 2, var)
        beta_var <- apply(beta_boot, 2, var)
        pcure_var<-var(pcure_boot)
        bsdmatrix <- sqrt(b_var)
        betasdmat <- sqrt(beta_var)
        pobssd<-sqrt(pcure_var)
#                    

