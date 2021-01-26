
# function to calculate baseline survival function
smsurv <-
function(Time,Status,X,beta,w){    
    death_point <- sort(unique(subset(Time, Status==1)))
    coxexp <- exp((beta)%*%t(X[,-1]))    
    lambda <- numeric()
    event <- numeric()
      for(i in 1: length(death_point)){
       event[i] <- sum(Status*as.numeric(Time==death_point[i]))
                 temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp),na.rm=T)
       	     temp1 <- event[i]
       lambda[i] <- temp1/temp
        }
    HHazard <- numeric()
    for(i in 1:length(Time)){
        HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
        if(Time[i]>max(death_point))HHazard[i] <- Inf
        if(Time[i]<min(death_point))HHazard[i] <- 0
        }
   survival <- exp(-HHazard)
   list(survival=survival)
}
#
# this EM algorithm estimates P(observed|uncure) 
emnewmod <-
function(Time,Status,observed,X,Z,offsetvar,b,beta,link,emmax,eps)
{     
# initial estimate
   	w <- Status	
      pobserve<-(sum(observed)-sum(Status))/(length(Time)-sum(Status))
	n <- length(Status)
	s <- smsurv(Time,Status,X,beta,w,model)$survival
	convergence<- 1000;i <-1
	while (convergence > eps & i < emmax){  
		uncureprob <- matrix(exp((b)%*%t(Z))/(1+exp((b)%*%t(Z))),ncol=1)
		survival<-drop(s^(exp((beta)%*%t(X[,-1]))))
		## E step 
		w <- Status+(1-observed)*(uncureprob*survival)/((1-pobserve)*(1-uncureprob)+uncureprob*survival)

	## M step
	logistfit<- eval(parse(text = paste("glm", "(", "w~Z[,-1]",",family = quasibinomial(link='", link, "'",")",")",sep = "")))
	update_cureb <- logistfit$coef
      if(!is.null(offsetvar)) update_cureb <- as.numeric(eval(parse(text = paste("glm", "(", "w~Z[,-1]+offset(offsetvar)",",family = quasibinomial(link='", link, "'",")",")",sep = "")))$coef)
	update_beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)), subset=w!=0, method="breslow")$coef
	if(!is.null(offsetvar)) update_beta <- coxph(Surv(Time, Status)~X[,-1]+offset(offsetvar+log(w)), subset=w!=0, method="breslow")$coef
	update_s <-smsurv(Time,Status,X,beta,w,model)$survival
	convergence<-sum(c(update_cureb-b,update_beta-beta)^2)+sum((s-update_s)^2)
	b <- update_cureb
     	beta <- update_beta 
	s<-update_s
	uncureprob <- matrix(exp((b)%*%t(Z))/(1+exp((b)%*%t(Z))),ncol=1)
      pobserve<-(sum(observed)-sum(Status))/(length(Time)-sum(w))
   	i <- i+1 
		}
 
