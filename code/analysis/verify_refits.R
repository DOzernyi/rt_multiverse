# Refit a small deterministic sample of accepted final formulas, not a new multiverse.
source("BRM impact/analysis/multiverse.R")
checks<-list()
for(id in c("fujita2021_e1","winkowski2024_e2","brockmeyer2026_prediction")) {
  m<-Filter(function(m)m$dataset==id,registry)[[1]]
  x<-readRDS(file.path(out,"fits",paste0(id,"__critical.rds")))
  rows<-x$coefficients[x$coefficients$policy=="median" &
    x$coefficients$transformation=="ln",]
  input<-readRDS(file.path(out,m$data))$data
  input<-input[input$region=="critical",]
  d<-input[screen(input,"median"),]
  design<-prepare_design(d,m);d<-design$d
  y<-log(d$RT);center<-mean(y);unit<-sd(y);d$y<-(y-center)/unit
  engine<-if(grepl("^bobyqa",x$optimizer))"bobyqa" else "nloptwrap"
  controls<-if(engine=="bobyqa")list(maxfun=20000) else list(maxeval=20000)
  fit<-lme4::lmer(as.formula(rows$formula[1]),data=d,REML=TRUE,
    control=lmerControl(optimizer=engine,optCtrl=controls))
  fresh<-coef(summary(lmerTest::as_lmerModLmerTest(fit)))
  rows<-rows[rows$term!="(Intercept)",]
  fresh<-fresh[match(rows$term,rownames(fresh)),,drop=FALSE]
  error<-max(abs(fresh[,"Estimate"]*unit-rows$estimate)/rows$SE)
  se_error<-max(abs(fresh[,"Std. Error"]*unit/rows$SE-1))
  same<-all((fresh[,"Pr(>|t|)"]<.05)==(rows$p<.05))
  stopifnot(error<.005,se_error<.005,same,!isSingular(fit,tol=1e-4))
  checks[[id]]<-list(max_estimate_discrepancy_in_SE=error,
    max_relative_SE_discrepancy=se_error,same_significance=same)
}
write_json(list(passed=TRUE,models=checks),
  file.path(out,"independent_refit_verification.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Independent final-model refits passed for",length(checks),"datasets\n")
