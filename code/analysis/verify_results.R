# Independent release gate: check sources, inputs and inference arithmetic.
suppressPackageStartupMessages(library(jsonlite))
root<-"generated-artifacts/brm_multiverse"
registry<-fromJSON(file.path(root,"datasets.json"),simplifyVector=FALSE)
excluded<-fromJSON("BRM impact/analysis/excluded_datasets.json",simplifyVector=FALSE)
excluded<-vapply(excluded,`[[`,character(1),"dataset")
registry<-Filter(function(m)m$tier=="primary"&&!m$dataset%in%excluded,registry)
jobs<-0L;tested<-0L;all_status<-0L
for(m in registry) {
  input<-readRDS(file.path(root,m$data))
  stopifnot(identical(unname(tools::md5sum(file.path(root,m$data))),m$input_md5))
  stopifnot(nrow(input$data)==m$observations,
    all(is.finite(input$data$RT)),all(input$data$RT>0),
    !anyNA(input$data[,c("subject","item",names(m$factors))]))
  for(roi in names(m$region_counts)) {
    local<-droplevels(input$data[input$data$region==roi,,drop=FALSE])
    rownames(local)<-NULL
    tmp<-tempfile()
    saveRDS(list(data=local,factors=m$factors,levels=m$levels),tmp,compress=FALSE)
    stopifnot(identical(unname(tools::md5sum(tmp)),m$region_fingerprints[[roi]]))
    unlink(tmp)
    x<-readRDS(file.path(root,"fits",paste0(m$dataset,"__",roi,".rds")))
    stopifnot(isTRUE(x$complete),identical(x$input_fingerprint,m$region_fingerprints[[roi]]),
      identical(x$source_sha,m$sha256),x$version=="multiverse-4",
      nrow(x$status)==78L,length(unique(x$status$policy))==13L)
    stopifnot(sum(input$data$region==roi)==m$region_counts[[roi]])
    a<-x$attempts[x$attempts$accepted%in%TRUE,]
    stopifnot(!any(a$singular%in%TRUE))
    d<-x$coefficients
    # Native intercept correction is an export-only operation.
    d<-d[d$term!="(Intercept)",]
    stopifnot(all(is.finite(d$estimate)),all(d$SE>0),all(d$df>0),
      max(abs(d$t-d$estimate/d$SE))<1e-6,
      max(abs(d$p-2*pt(-abs(d$t),d$df)))<1e-7,
      max(abs(d$CI_low-(d$estimate-qt(.975,d$df)*d$SE)))<1e-5,
      max(abs(d$CI_high-(d$estimate+qt(.975,d$df)*d$SE)))<1e-5)
    ln<-d[d$transformation=="ln",];lg<-d[d$transformation=="log10",]
    key<-function(d)paste(d$policy,d$term)
    lg<-lg[match(key(ln),key(lg)),]
    stopifnot(nrow(ln)==nrow(lg),max(abs(lg$estimate-ln$estimate/log(10)))<1e-7,
      max(abs(lg$p-ln$p))<1e-10)
    jobs<-jobs+1L;tested<-tested+nrow(d);all_status<-all_status+nrow(x$status)
  }
}
summary<-fromJSON(file.path(root,"summary.json"))
stopifnot(jobs==summary$corpus$region_jobs,all_status==summary$pipelines$nominal)
write_json(list(passed=TRUE,datasets=length(registry),region_jobs=jobs,
  coefficient_rows_verified=tested,nominal_pipeline_rows=all_status,
  checks=c("source/input checkpoint identity","sample sizes and valid responses",
    "complete policy coverage","accepted fits nonsingular","t and p arithmetic",
    "confidence interval arithmetic","exact log10 rescaling")),
 file.path(root,"verification_results.json"),pretty=TRUE,auto_unbox=TRUE)
cat("Verified",jobs,"completed region jobs and",tested,"non-intercept coefficient rows\n")
