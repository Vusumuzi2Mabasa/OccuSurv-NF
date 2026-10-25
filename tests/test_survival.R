library(survival)
# Hand-computed competing events: three workers, diagnosis then death.
d<-data.frame(time=c(1,2,3),cause=factor(c(1,2,0),levels=0:2))
f<-survfit(Surv(time,cause)~1,d,timefix=FALSE)
stopifnot(max(abs(f$pstate[2,]-rep(1/3,3)))<1e-12)
stopifnot(max(abs(rowSums(f$pstate)-1))<1e-12)
# All-censored data retain unit survival and zero cumulative incidence.
d$cause<-factor(rep(0,3),levels=0:2)
f<-survfit(Surv(time,cause)~1,d,timefix=FALSE)
stopifnot(all(f$pstate[,1]==1),all(f$pstate[,2:3]==0))
# Start-stop boundary likelihood at fixed beta; only one interval per worker.
d<-data.frame(id=c(1,2,2,3,3),start=c(0,0,1,0,1),stop=c(1,1,3,1,4),event=c(1,0,1,0,0),x=c(0,1,2,2,3))
b<-.2
f<-coxph(Surv(start,stop,event)~x,d,init=b,control=coxph.control(iter.max=0,timefix=FALSE),ties='breslow')
manual<-sum(sapply(which(d$event==1),function(j){risk<-d$start<d$stop[j]&d$stop>=d$stop[j];d$x[j]*b-log(sum(exp(d$x[risk]*b)))}))
stopifnot(abs(f$loglik[2]-manual)<1e-12)
cat('PASS: hand Aalen-Johansen, probability mass, all-censored states, start-stop boundary likelihood\n')
