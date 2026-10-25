#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(survival);library(jsonlite);library(ggplot2)})
args <- commandArgs(TRUE)
clean <- args[1]; raw <- args[2]; lng <- args[3]
for (d in c('analysis','longitudinal','figures','longitudinal_figures')) dir.create(d,showWarnings=FALSE)
read <- function(p) read.csv(p,stringsAsFactors=FALSE)
write <- function(x,p) write.csv(x,p,row.names=FALSE,na='')
js <- function(x,p) write_json(x,p,pretty=TRUE,auto_unbox=TRUE,digits=16,na='null')
w <- read(file.path(clean,'workers.csv')); c <- read(file.path(clean,'cohort.csv')); m <- read(file.path(clean,'measurements.csv'))
groups <- function(d,col) c(list(Overall=d),split(d,d[[col]]))
wilson <- function(k,n) {z<-qnorm(.975); p<-k/n; den<-1+z*z/n; mid<-(p+z*z/(2*n))/den; half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den;c(p,mid-half,mid+half)}
rate <- function(k,t) c(1000*k/t,if(k==0) 0 else 500*qchisq(.025,2*k)/t,500*qchisq(.975,2*(k+1))/t)
prev<-inc<-tab<-list()
for (name in names(groups(w,'industry'))) {
 ww<-groups(w,'industry')[[name]];cc<-groups(c,'industry')[[name]]
 p<-wilson(sum(ww$baseline_disease),nrow(ww));r<-rate(sum(cc$event),sum(cc$time_years))
 prev[[name]]<-data.frame(industry=name,n=nrow(ww),cases=sum(ww$baseline_disease),prevalence=p[1],lower=p[2],upper=p[3])
 inc[[name]]<-data.frame(industry=name,n=nrow(cc),events=sum(cc$event),person_years=sum(cc$time_years),rate_per_1000=r[1],lower=r[2],upper=r[3])
 tab[[name]]<-data.frame(industry=name,n=nrow(ww),age_mean=mean(ww$age),age_sd=sd(ww$age),smoking_observed=sum(!is.na(ww$smoker)),smoking_missing=sum(is.na(ww$smoker)),smokers=sum(ww$smoker,na.rm=TRUE),cumulative_median=median(ww$cumulative_mg_year_m3))
}
write(do.call(rbind,prev),'analysis/prevalence.csv');write(do.call(rbind,inc),'analysis/incidence.csv');write(do.call(rbind,tab),'analysis/table1.csv')
means<-aggregate(concentration_mg_m3~worker_id,m,mean);names(means)[2]<-'mean_measured_mg_m3';means<-merge(means,w[c('worker_id','industry')]);write(means,'analysis/exposure_means.csv')
curves<-risk<-list()
for (name in names(split(c,c$industry))) {
 d<-c[c$industry==name,]; f<-survfit(Surv(time_years,event)~1,d,conf.type='log-log',timefix=FALSE)
 f$lower[f$surv==1|f$surv==0]<-f$surv[f$surv==1|f$surv==0];f$upper[f$surv==1|f$surv==0]<-f$surv[f$surv==1|f$surv==0]
 curves[[name]]<-data.frame(time=c(0,f$time),survival=c(1,f$surv),lower=c(1,f$lower),upper=c(1,f$upper),at_risk=c(nrow(d),f$n.risk),events=c(0,f$n.event),industry=name)
 risk[[name]]<-data.frame(industry=name,time=c(0,2,4,6,8),at_risk=sapply(c(0,2,4,6,8),function(t)sum(d$time_years>=t)))
}
km<-do.call(rbind,curves);write(km,'analysis/km.csv');write(do.call(rbind,risk),'analysis/risk_table.csv')
fit_model <- function(d,terms,cause=NULL,path,assumed) {
 d<-d[complete.cases(d[terms]),];event<-if(is.null(cause))d$event else as.integer(d$cause==cause)
 n<-if(is.null(cause))nrow(d) else length(unique(d$worker_id))
 status<-list(status='skipped',events=sum(event),n=n,subjects=n,intervals=nrow(d),minimum_events=30,reason='fewer than 30 cause-specific events')
 if(sum(event)<30)return(list(status=status))
 d$endpoint<-event
 response<-if(is.null(cause))'Surv(time_years, endpoint)' else 'Surv(start, stop, endpoint)'
 form<-as.formula(paste(response,'~',paste(terms,collapse='+')))
 warns<-character()
 fit<-withCallingHandlers(coxph(form,d,ties='breslow',robust=FALSE,x=TRUE,y=TRUE,control=coxph.control(eps=1e-10,iter.max=100,timefix=FALSE)),warning=function(e){warns<<-c(warns,conditionMessage(e));invokeRestart('muffleWarning')})
 b<-coef(fit);se<-sqrt(diag(vcov(fit)));stopifnot(all(is.finite(b)),all(is.finite(se)),all(se>0))
 score<-max(abs(colSums(residuals(fit,type='score'))));stopifnot(score<1e-4)
 tbl<-data.frame(term=terms,beta=b,se=se,hr=exp(b),lower=exp(b-qnorm(.975)*se),upper=exp(b+qnorm(.975)*se),p=2*pnorm(-abs(b/se)),assumed_hr=assumed,row.names=NULL)
 if(!is.null(cause))tbl$cause<-cause
 write(tbl,path)
 status$status<-'fitted';status$reason<-NULL;status$score_max<-score;status$warnings<-warns;status$ties<-'breslow';status$covariance<-'model-based counting-process information; independent workers';status$log_likelihood<-unname(fit$loglik[2])
 list(status=status,fit=fit,data=d,table=tbl)
}
base<-fit_model(c,c('log2_1plus_cumulative','age_per_decade','smoker'),path='analysis/cox.csv',assumed=c(1.7,1.3,1.5))
if(base$status$status=='fitted') {
 sr<-residuals(base$fit,type='schoenfeld'); tt<-as.numeric(rownames(sr));rr<-dg<-list()
 for(j in seq_len(ncol(sr))) {rr[[j]]<-data.frame(term=colnames(sr)[j],time=tt,residual=sr[,j]);dg[[j]]<-data.frame(term=colnames(sr)[j],spearman_rho=cor(log(tt),sr[,j],method='spearman'),n_events=length(tt))}
 write(do.call(rbind,rr),'analysis/schoenfeld.csv');write(do.call(rbind,dg),'analysis/ph_diagnostics.csv')
 z<-cox.zph(base$fit,transform='log');write(data.frame(term=rownames(z$table),z$table,row.names=NULL),'analysis/cox_zph.csv')
}
js(list(synthetic=TRUE,workers=nrow(w),prevalent_cases=sum(w$baseline_disease),incident_cohort=nrow(c),events=sum(c$event),person_years=sum(c$time_years),model=base$status,truth=fromJSON(file.path(raw,'truth.json')),versions=list(R=R.version.string,survival=as.character(packageVersion('survival')))),'analysis/analysis.json')
file.copy(file.path(lng,'raw'),'longitudinal',recursive=TRUE);file.copy(file.path(lng,'qc.json'),'longitudinal')
i<-read(file.path(lng,'raw/intervals.csv'));s<-read(file.path(lng,'raw/subjects.csv'));h<-read(file.path(lng,'raw/history.csv'))
aj <- function(d) {
 d$state<-factor(d$cause,levels=0:2,labels=c('censored','disease','death'))
 f<-survfit(Surv(time,state)~1,d,timefix=FALSE)
 net<-survfit(Surv(time,cause==1)~1,d,timefix=FALSE)
 data.frame(time=c(0,f$time),survival=c(1,f$pstate[,1]),disease=c(0,f$pstate[,2]),death=c(0,f$pstate[,3]),naive_one_minus_km=c(0,1-net$surv))
}
at_grid<-function(f,grid,col) f[[col]][pmax(1,findInterval(grid,f$time))]
set.seed(2028);grid<-seq(0,8,.1);B<-200;ac<-ci<-boot_indices<-list()
for(name in names(groups(s,'industry'))) {
 d<-groups(s,'industry')[[name]];f<-aj(d);f$group<-name;ac[[name]]<-f
 indices<-replicate(B,sample.int(nrow(d),nrow(d),replace=TRUE));boot_indices[[name]]<-indices
 draws<-sapply(seq_len(B),function(b){bd<-d[indices[,b],];bc<-aj(bd);cbind(at_grid(bc,grid,'disease'),at_grid(bc,grid,'death'))},simplify='array')
 for(j in 1:2) {lim<-apply(draws[,j,],1,quantile,probs=c(.025,.975));cause<-c('disease','death')[j];ci[[paste(name,j)]]<-data.frame(group=name,cause=cause,time=grid,cif=at_grid(f,grid,cause),lower=lim[1,],upper=lim[2,],bootstrap_replicates=B)}
}
ac<-do.call(rbind,ac);ci<-do.call(rbind,ci);write(ac,'longitudinal/aj_curves.csv');write(ci,'longitudinal/cif_intervals.csv')
js(boot_indices,'longitudinal/bootstrap_indices.json')
mods<-lapply(1:2,function(k)fit_model(i,c('log2_1plus_lagged','age_per_decade','smoker'),k,paste0('longitudinal/cox_cause',k,'.csv'),if(k==1)c(1.7,1.3,1.5)else c(1,1.5,1.2)))
statuses<-lapply(1:2,function(k)c(list(cause=k),mods[[k]]$status))
js(list(scenario='longitudinal competing events',seed=fromJSON(file.path(lng,'raw/truth.json'))$seed,subjects=nrow(s),intervals=nrow(i),disease=sum(s$cause==1),deaths=sum(s$cause==2),censored=sum(s$cause==0),person_years=sum(s$time),model_complete_subjects=sum(!is.na(s$smoker)),models=statuses,bootstrap_replicates=B,bootstrap_seed=2028,worker_sampling='iid within each displayed group; no residual workplace frailty assumed'),'longitudinal/summary.json')
theme_set(theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),legend.position='bottom',plot.title=element_text(face='bold')))
saveplot<-function(p,path,height=4.3) ggsave(path,p,width=8,height=height,dpi=180,bg='white')
saveplot(ggplot(means,aes(industry,mean_measured_mg_m3,colour=industry))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=.12,alpha=.45,size=1)+labs(x=NULL,y=expression('Mean measured silica (mg/'*m^3*')'),title='Synthetic baseline exposure')+guides(colour='none'),'figures/exposure.png')
incd<-do.call(rbind,inc);incd<-incd[incd$industry!='Overall',]
saveplot(ggplot(incd,aes(rate_per_1000,industry))+geom_errorbarh(aes(xmin=lower,xmax=upper),height=.15)+geom_point(size=3,colour='#176b85')+labs(x='Diagnoses per 1000 person years',y=NULL,title='Crude incidence with 95% Poisson intervals'),'figures/incidence.png')
saveplot(ggplot(km,aes(time,survival,colour=industry,fill=industry))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.12,colour=NA)+geom_step()+labs(x='Years since baseline',y='Disease free survival',title='Kaplan Meier estimates'),'figures/km.png')
forest<-function(tbl,title) ggplot(tbl,aes(hr,term))+geom_vline(xintercept=1,linetype=2)+geom_errorbarh(aes(xmin=lower,xmax=upper),height=.15)+geom_point(colour='#176b85',size=2.5)+geom_point(aes(x=assumed_hr),shape=4,colour='#cb761c',size=3)+scale_x_log10()+labs(x='Hazard ratio and 95% Wald interval',y=NULL,title=title)
empty<-function()ggplot()+annotate('text',x=0,y=0,label='Regression skipped\nFewer than 30 complete case events')+theme_void()
saveplot(if(base$status$status=='fitted')forest(base$table,'Baseline adjusted associations')else empty(),'figures/cox_forest.png')
saveplot(if(base$status$status=='fitted')ggplot(do.call(rbind,rr),aes(log(time),residual))+geom_point(alpha=.5,size=1)+geom_smooth(method='loess',se=FALSE,span=.75)+facet_wrap(~term,scales='free_y')+labs(x='Log event time',y='Schoenfeld residual',title='Exploratory proportional hazards diagnostics')else empty(),'figures/ph_diagnostics.png')
observed<-merge(h,s[c('worker_id','time')]);observed<-observed[observed$year<observed$time,]
eh<-rbind(data.frame(observed[c('industry','year')],value=observed$concentration,metric='Annual concentration mg per cubic metre'),data.frame(observed[c('industry','year')],value=observed$lagged_cumulative,metric='Five year lagged cumulative exposure'))
eh<-aggregate(value~industry+year+metric,eh,median)
saveplot(ggplot(eh,aes(year,value,colour=industry))+geom_line()+geom_point()+facet_wrap(~metric,scales='free_y')+labs(x='Year since baseline',y='Median among observed workers',title='Time varying exposure histories'),'longitudinal_figures/exposure_history.png')
saveplot(ggplot(ci[ci$group=='Overall',],aes(time,cif,colour=cause,fill=cause))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.12,colour=NA)+geom_step()+geom_step(data=ac[ac$group=='Overall',],aes(time,naive_one_minus_km),inherit.aes=FALSE,linetype=2)+labs(x='Years since baseline',y='Cumulative incidence',title='Disease and competing death'), 'longitudinal_figures/cumulative_incidence.png')
saveplot(ggplot(ci[ci$group!='Overall' & ci$cause=='disease',],aes(time,cif,colour=group,fill=group))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.12,colour=NA)+geom_step()+labs(x='Years since baseline',y='Disease cumulative incidence',title='Industry specific cumulative incidence'),'longitudinal_figures/cif_by_industry.png')
tables<-lapply(1:2,function(k)if(mods[[k]]$status$status=='fitted')transform(mods[[k]]$table,endpoint=c('Diagnosis','Competing death')[k])else NULL)
ft<-do.call(rbind,tables)
saveplot(if(!is.null(ft))forest(ft,'Time varying cause specific associations')+facet_wrap(~endpoint)else empty(),'longitudinal_figures/time_varying_hr.png')
capture.output(sessionInfo(),file='analysis/R-sessionInfo.txt')
