# Summarize complete primary checkpoints; never refit or modify a checkpoint.
suppressPackageStartupMessages(library(jsonlite))
out <- normalizePath(Sys.getenv("BRM_MULTIVERSE_OUT","generated-artifacts/brm_multiverse"),mustWork=TRUE)
root <- normalizePath("BRM impact",mustWork=TRUE)
registry <- fromJSON(file.path(out,"datasets.json"),simplifyVector=FALSE)
exclusions <- fromJSON(file.path(root,"analysis","excluded_datasets.json"),simplifyVector=FALSE)
excluded <- vapply(exclusions,`[[`,character(1),"dataset")
registry <- Filter(function(m) m$tier=="primary" && !m$dataset%in%excluded,registry)
files <- unlist(lapply(registry,function(m) file.path(out,"fits",
  paste0(m$dataset,"__",names(m$region_counts),".rds"))))
if(any(!file.exists(files))) stop("Missing primary checkpoints")
objects <- lapply(files,readRDS)
for(x in objects) {
  m <- Filter(function(m)m$dataset==x$metadata$dataset,registry)[[1]]
  stopifnot(x$version=="multiverse-4",isTRUE(x$complete),
    identical(x$source_sha,m$sha256),
    identical(x$input_fingerprint,m$region_fingerprints[[x$region]]))
}
decorate <- function(x,component) {
  d <- x[[component]]
  if(is.null(d)||!nrow(d)) return(NULL)
  d$dataset<-x$metadata$dataset; d$citation<-x$metadata$citation
  d$task<-x$metadata$task; d$region<-x$region
  d
}
components <- c(coefficients="coefficients",status="pipeline_status",
  attempts="model_attempts",random_variances="random_variances",
  conditional_followups="conditional_followups")
tables<-lapply(names(components),function(k)do.call(rbind,lapply(objects,decorate,component=k)))
names(tables)<-names(components)
# Restore the native-intercept null without altering archived fits.
cc<-tables$coefficients; ix<-cc$term=="(Intercept)"
cc$t[ix]<-cc$estimate[ix]/cc$SE[ix]
cc$p[ix]<-2*pt(-abs(cc$t[ix]),cc$df[ix]); tables$coefficients<-cc
for(k in names(components)) if(!is.null(tables[[k]]))
  write.csv(tables[[k]],gzfile(file.path(out,paste0(components[[k]],".csv.gz"))),row.names=FALSE)
stopifnot(nrow(tables$status)==length(files)*78L)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
policies<-c("none","sd1.5","sd2","sd2.5","sd3","fixed100_2000","fixed100_3000",
 "fixed200_2000","fixed200_3000","median","mad2","mad2.5","mad3")
labels<-c("None","1.5 SD","2 SD","2.5 SD","3 SD","100-2000","100-3000",
 "200-2000","200-3000","Median","2 MAD","2.5 MAD","3 MAD")
types<-c("identity","ln","inverse","sqrt","boxcox")
type_labels<-c("Raw ms","Natural log","Inverse","Square root","Box-Cox")
colors<-c("#1F74AD","#E98E2E","#61A9BD","#9267C6","#8D6C5B")
status<-tables$status[tables$status$transformation!="log10",]
coef<-cc[cc$term!="(Intercept)" & cc$transformation!="log10",]
coef$key<-paste(coef$dataset,coef$region,coef$term,sep="__")
base<-coef[coef$policy=="none"&coef$transformation=="ln",]
stopifnot(!anyDuplicated(base$key))
reference<-base[,c("key","p","estimate","SE","standardized_estimate","formula")]
names(reference)<-c("key","reference_p","reference_estimate","reference_SE",
                   "reference_standardized_estimate","reference_formula")
coef<-merge(coef,reference,by="key",all.x=TRUE)
coef$significant<-coef$p<.05
coef$disagreement<-coef$significant!=(coef$reference_p<.05)
coef$same_structure<-coef$formula==coef$reference_formula
coef$standardized_shift<-coef$standardized_estimate-coef$reference_standardized_estimate
coef$order<-lengths(strsplit(coef$term,":",fixed=TRUE))
coef$effect_kind<-ifelse(coef$order>1,"Interaction","Main effect")
write.csv(coef,gzfile(file.path(out,"coefficient_comparisons.csv.gz")),row.names=FALSE)
finite_mean<-function(x)if(any(is.finite(x)))mean(x[is.finite(x)]) else NA_real_
effects<-do.call(rbind,lapply(split(coef,coef$key),function(d) {
  ref<-is.finite(d$reference_p[1])
  same<-d$same_structure%in%TRUE
  se<-d$standardized_SE[is.finite(d$standardized_SE)&d$standardized_SE>0]
  data.frame(dataset=d$dataset[1],region=d$region[1],term=d$term[1],
    citation=d$citation[1],task=d$task[1],effect_kind=d$effect_kind[1],
    successful=nrow(d),reference_available=ref,
    significance_rate=mean(d$significant),inference_changes=length(unique(d$significant))>1,
    direction_changes=length(unique(sign(d$estimate)))>1,
    standardized_min=min(d$standardized_estimate),standardized_max=max(d$standardized_estimate),
    standardized_range=diff(range(d$standardized_estimate)),
    SE_ratio=if(length(se))max(se)/min(se) else NA_real_,
    reference_standardized_estimate=d$reference_standardized_estimate[1],
    disagreement=finite_mean(as.numeric(d$disagreement)),
    same_structure_n=sum(same),
    same_structure_disagreement=finite_mean(as.numeric(d$disagreement[same])))
}))
write.csv(effects,file.path(out,"effect_stability.csv"),row.names=FALSE)
dataset_effects<-do.call(rbind,lapply(split(effects,effects$dataset),function(d)
  data.frame(dataset=d$dataset[1],citation=d$citation[1],task=d$task[1],
    disagreement=finite_mean(d$disagreement),
    inference_changes=mean(d$inference_changes),direction_changes=mean(d$direction_changes),
    standardized_range=median(d$standardized_range),SE_ratio=median(d$SE_ratio))))
metadata<-do.call(rbind,lapply(registry,function(m)data.frame(
 dataset=m$dataset,participants=m$participants,items=m$items,skewness=m$skewness,
 design=m$design)))
dataset_effects<-merge(dataset_effects,metadata,by="dataset")
write.csv(dataset_effects,file.path(out,"dataset_stability.csv"),row.names=FALSE)
# A publication receives one vote, after equal coefficient/region and dataset weights.
source_values<-function(d,value) {
  by_dataset<-aggregate(d[[value]],list(dataset=d$dataset,citation=d$citation),
                       finite_mean)
  names(by_dataset)[3]<-"value"
  v<-tapply(by_dataset$value,by_dataset$citation,finite_mean)
  v[is.finite(v)]
}
set.seed(20261118)
cluster_mean<-function(d,value) {
  v<-source_values(d,value)
  if(!length(v))return(list(mean=NA_real_,lower=NA_real_,upper=NA_real_,sources=0))
  bs<-replicate(1000,mean(sample(v,length(v),replace=TRUE)))
  list(mean=mean(v),lower=unname(quantile(bs,.025)),
       upper=unname(quantile(bs,.975)),sources=length(v))
}
quantiles<-function(x)as.list(setNames(as.numeric(quantile(x,c(.25,.5,.75),na.rm=TRUE)),
                                      c("q25","median","q75")))
retention<-do.call(rbind,lapply(policies,function(p){
  d<-status[status$policy==p & status$transformation=="ln",]
  b<-cluster_mean(d,"proportion_excluded")
  data.frame(policy=p,as.data.frame(b))
}))
write.csv(retention,file.path(out,"retention_summary.csv"),row.names=FALSE)
inference<-do.call(rbind,lapply(policies,function(p)do.call(rbind,lapply(types,function(t){
  d<-coef[coef$policy==p & coef$transformation==t,]
  b<-cluster_mean(d,"disagreement")
  st<-status[status$policy==p & status$transformation==t,]
  data.frame(policy=p,transformation=t,as.data.frame(b),
    accepted=sum(st$success),attempted=nrow(st))
}))))
write.csv(inference,file.path(out,"inference_summary.csv"),row.names=FALSE)
subgroups<-list()
for(variable in c("task","effect_kind","region","design")) {
  d<-effects
  if(variable=="region")d$region<-ifelse(grepl("^postcritical",d$region),"Post-critical","Critical")
  if(variable=="design")d$design<-metadata$design[match(d$dataset,metadata$dataset)]
  for(level in sort(unique(d[[variable]]))) {
    z<-d[d[[variable]]==level,]; b<-cluster_mean(z,"disagreement")
    subgroups[[length(subgroups)+1]]<-data.frame(group=variable,level=level,
      coefficient_regions=nrow(z),unstable=sum(z$inference_changes),as.data.frame(b))
  }
}
subgroups<-do.call(rbind,subgroups)
write.csv(subgroups,file.path(out,"subgroup_summary.csv"),row.names=FALSE)
moderators<-list()
for(variable in c("skewness","participants","items","reference_magnitude")) {
  d<-dataset_effects
  if(variable=="reference_magnitude") {
    mag<-tapply(abs(effects$reference_standardized_estimate),effects$dataset,finite_mean)
    d$reference_magnitude<-mag[d$dataset]
  }
  d<-d[is.finite(d[[variable]])&is.finite(d$disagreement),]
  sources<-unique(d$citation)
  rho<-suppressWarnings(cor(d[[variable]],d$disagreement,method="spearman"))
  bs<-replicate(1000,{
    sampled<-sample(sources,length(sources),replace=TRUE)
    b<-do.call(rbind,lapply(sampled,function(s)d[d$citation==s,]))
    suppressWarnings(cor(b[[variable]],b$disagreement,method="spearman"))
  })
  moderators[[variable]]<-list(rho=rho,lower=unname(quantile(bs,.025,na.rm=TRUE)),
    upper=unname(quantile(bs,.975,na.rm=TRUE)),datasets=nrow(d),sources=length(sources))
}
ss<-tables$status; att<-tables$attempts
accepted<-att[att$accepted%in%TRUE,]
lambda<-status[status$transformation=="boxcox" & is.finite(status$lambda),]
summary<-list(
 corpus=list(datasets=length(registry),publications=length(unique(dataset_effects$citation)),
   region_jobs=length(files),observations=sum(vapply(registry,`[[`,numeric(1),"observations")),
   SPR=sum(dataset_effects$task=="SPR"),Maze=sum(dataset_effects$task=="Maze"),
   critical_jobs=sum(!grepl("^postcritical",vapply(objects,`[[`,character(1),"region"))),
   postcritical_jobs=sum(grepl("^postcritical",vapply(objects,`[[`,character(1),"region"))),
   exclusions=exclusions),
 pipelines=list(nominal=nrow(ss),accepted_nominal=sum(ss$success),
   distinct=nrow(status),accepted_distinct=sum(status$success),
   failed_distinct=sum(!status$success),fit_attempts=nrow(att),
   accepted_first=sum(accepted$step==1),accepted_reduced=sum(accepted$step>1),
   reference_coefficients=nrow(base),
   coefficient_regions=nrow(effects),unstable=sum(effects$inference_changes),
   direction_changes=sum(effects$direction_changes),
   same_structure_disagreement=cluster_mean(coef[coef$same_structure%in%TRUE,],"disagreement"),
   disagreement=cluster_mean(coef,"disagreement")),
 effect_range=quantiles(effects$standardized_range),precision_ratio=quantiles(effects$SE_ratio),
 lambda=list(min=min(lambda$lambda),max=max(lambda$lambda),
   boundaries=sum(lambda$lambda_at_boundary)),
 moderators=moderators,
 computational_seconds=sum(vapply(objects,`[[`,numeric(1),"elapsed_seconds")))
source(file.path(root,"analysis","decision_sensitivity.R"),local=TRUE)
summary$decision_sensitivity <- decision_summary
write_json(summary,file.path(out,"summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA,na="null")
cat(toJSON(summary,pretty=TRUE,auto_unbox=TRUE),"\n")

# Native vector figures for the print manuscript; shapes/line types survive grayscale.
figures<-file.path(out,"figures"); dir.create(figures,showWarnings=FALSE)
eps<-function(name,width=7.1,height=4.5) {
 postscript(file.path(figures,paste0(name,".eps")),onefile=FALSE,horizontal=FALSE,
   paper="special",width=width,height=height,family="Helvetica",colormodel="srgb",pointsize=10)
 par(mar=c(6.5,4.3,1,1),las=1,mgp=c(2.6,.65,0),tcl=-.25,bty="l")
}
xlabels<-function(){axis(1,at=1:13,labels=FALSE)
 text(1:13,par("usr")[3]-.04*diff(par("usr")[3:4]),labels,srt=42,adj=1,xpd=TRUE,cex=.79)}
eps("retention")
ylim<-c(0,ceiling(max(retention$upper)*100/5)*5)
plot(1:13,100*retention$mean,type="n",xaxt="n",ylim=ylim,xlab="",ylab="Additional RTs excluded (%)")
arrows(1:13,100*retention$lower,1:13,100*retention$upper,angle=90,code=3,length=.04,col=colors[1])
points(1:13,100*retention$mean,pch=16,col=colors[1]);xlabels();dev.off()
eps("inference_stability",height=5)
par(mar=c(2,7.5,3.4,5.8))
mat <- matrix(NA_real_,13,5)
for(i in 1:13)for(j in 1:5)
 mat[i,j]<-100*inference$mean[inference$policy==policies[i]&
                             inference$transformation==types[j]]
shades<-colorRampPalette(c("#FFFFFF","#8CCEDE","#154C75"))(100)
image(1:5,1:13,t(mat[13:1,]),col=shades,zlim=c(0,25),
 xlim=c(.5,5.5),ylim=c(.5,13.5),
 axes=FALSE,xlab="",ylab="")
axis(3,at=1:5,labels=type_labels,tick=FALSE,cex.axis=.87)
axis(2,at=1:13,labels=rev(labels),tick=FALSE,las=1,cex.axis=.85)
for(i in 1:13)for(j in 1:5)
 text(j,14-i,sprintf("%.1f",mat[i,j]),cex=.83,
      col=if(mat[i,j]>17)"white" else "#154C75")
for(x in seq(.5,5.5,1))abline(v=x,col="white",lwd=.7)
for(y in seq(.5,13.5,1))abline(h=y,col="white",lwd=.7)
par(xpd=TRUE)
for(i in 0:99)rect(5.65,.5+i*.1,5.8,.5+(i+1)*.1,
                   col=shades[i+1],border=NA)
text(5.87,c(.5,5.5,10.5),c("0","12.5","25%"),pos=4,cex=.8)
dev.off()
eps("decision_sensitivity",width=7.1,height=6)
layout(matrix(c(1,1,2,3),nrow=2,byrow=TRUE),heights=c(1,3))
upper<-ceiling(max(decision_marginals$upper)*100/5)*5
dot_panel<-function(z,names,title,margin=7) {
 yy<-rev(seq_len(nrow(z)));par(mar=c(3.6,margin,2.2,1),mgp=c(2.3,.5,0))
 plot(c(0,upper),c(.5,nrow(z)+.5),type="n",yaxt="n",xlab="Changed classification (%)",ylab="")
 abline(v=seq(0,upper,5),col="gray90",lty=3)
 arrows(100*z$lower,yy,100*z$upper,yy,angle=90,code=3,length=.035,col=colors[1])
 points(100*z$mean,yy,pch=16,col=colors[1],cex=.8)
 axis(2,at=yy,labels=names,las=1,tick=FALSE,cex.axis=.86)
 mtext(title,3,line=.6,adj=0,cex=.93)
}
z<-decision_marginals[decision_marginals$level=="ALL",]
dot_panel(z,c("Exclusion rule","Transformation"),"A  Change one dimension only",margin=10)
z<-decision_marginals[decision_marginals$dimension=="exclusion"&
                       decision_marginals$level!="ALL",]
z<-z[order(z$mean,decreasing=TRUE),]
dot_panel(z,labels[match(z$level,policies)],"B  Exclusion rules",margin=7.5)
z<-decision_marginals[decision_marginals$dimension=="transformation"&
                       decision_marginals$level!="ALL",]
z<-z[order(z$mean,decreasing=TRUE),]
dot_panel(z,type_labels[match(z$level,types)],"C  Transformations",margin=7.5)
dev.off()
eps("lambda")
boxplot(lapply(policies,function(p)lambda$lambda[lambda$policy==p]),xaxt="n",
 col="#8CCEDE",ylab="Box-Cox lambda",outline=TRUE)
abline(h=0,lty=2,col="gray40");xlabels();dev.off()
eps("precision",height=4)
par(mfrow=c(1,2),mar=c(4.3,4,1,1))
hist(effects$standardized_range,breaks=25,col="#8CCEDE",border="white",main="",
 xlab="Range of standardized estimates",ylab="Coefficient-region combinations")
hist(log2(effects$SE_ratio),breaks=25,col="#FAB770",border="white",main="",
 xlab="log2(maximum / minimum standardized SE)",ylab="Coefficient-region combinations")
dev.off()
eps("subgroups",height=4.2)
z<-subgroups[subgroups$group%in%c("task","effect_kind","region"),]
yy<-rev(seq_len(nrow(z)));par(mar=c(4,9,1,1))
plot(c(0,max(z$upper,na.rm=TRUE)*105),range(yy)+c(-.5,.5),type="n",yaxt="n",
 xlab="Classifications differing from reference (%)",ylab="")
axis(2,at=yy,labels=paste0(z$level," (",z$sources," sources)"),las=1,cex.axis=.8)
segments(100*z$lower,yy,100*z$upper,yy,col=colors[1],lwd=2)
points(100*z$mean,yy,pch=16,col=colors[1]);dev.off()
