# Read-only audit of the exact compact inputs and completed checkpoint metadata.
suppressPackageStartupMessages(library(jsonlite))
root <- "generated-artifacts/brm_multiverse"
out <- file.path(root,"release_audit")
dir.create(file.path(out,"rds_exports"),recursive=TRUE,showWarnings=FALSE)
registry <- fromJSON(file.path(root,"datasets.json"),simplifyVector=FALSE)
excluded <- fromJSON("BRM impact/analysis/excluded_datasets.json",simplifyVector=FALSE)
excluded <- vapply(excluded,`[[`,character(1),"dataset")
registry <- Filter(function(m)m$tier=="primary"&&!m$dataset%in%excluded,registry)
normalise <- function(v) {
 s <- as.character(v)
 hit <- grepl("^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$",s)
 s[hit] <- vapply(as.numeric(s[hit]),function(n)format(n,scientific=FALSE,trim=TRUE,digits=16),character(1))
 s
}
inventory <- list(); levels <- list(); cells <- list(); jobs <- list(); privacy <- list()
coding <- list(); structures <- list()
# Evaluate only the design-construction function, not the runner or any fit.
design_function <- Filter(function(e)is.call(e)&&identical(e[[1]],as.name("<-"))&&
 identical(e[[2]],as.name("prepare_design")),parse("BRM impact/analysis/multiverse.R"))
stopifnot(length(design_function)==1)
eval(design_function[[1]])
for(m in registry) {
 input <- readRDS(file.path(root,m$data));d <- input$data
 stopifnot(identical(unname(tools::md5sum(file.path(root,m$data))),m$input_md5),
   identical(input$metadata$dataset,m$dataset),
   identical(input$metadata$factors,m$factors),
   identical(input$metadata$levels,m$levels),
   nrow(d)==m$observations,length(unique(d$subject))==m$participants,
   length(unique(d$item))==m$items,setequal(unique(as.character(d$region)),names(m$region_counts)))
 actual_levels <- lapply(names(m$factors),function(f)sort(unique(normalise(d[[f]]))))
 names(actual_levels) <- names(m$factors)
 design <- paste(lengths(actual_levels),collapse="x")
 stopifnot(identical(design,m$design))
 for(f in names(m$factors)) {
  stopifnot(setequal(actual_levels[[f]],normalise(unlist(m$levels[[f]]))))
  levels[[length(levels)+1]] <- data.frame(dataset=m$dataset,citation=m$citation,
   factor=f,source_factor=m$factors[[f]],design=m$design,
   registry_levels=as.character(toJSON(m$levels[[f]],auto_unbox=FALSE)),
   observed_levels=as.character(toJSON(sort(unique(as.character(d[[f]]))),auto_unbox=FALSE)),
   numeric_representation_normalised=TRUE)
 }
 for(region in names(m$region_counts)) {
  z <- droplevels(d[d$region==region,])
  stopifnot(nrow(z)==m$region_counts[[region]])
  specification <- prepare_design(z,m)
  stopifnot(qr(specification$X)$rank==ncol(specification$X))
  for(f in names(m$factors)) {
   contrast <- contrasts(specification$d[[f]])
   for(j in seq_len(nrow(contrast)))for(k in seq_len(ncol(contrast)))
    coding[[length(coding)+1]] <- data.frame(dataset=m$dataset,region=region,
     factor=f,source_factor=m$factors[[f]],level=levels(specification$d[[f]])[j],
     contrast_column=k,value=contrast[j,k])
  }
  for(group in names(specification$active))structures[[length(structures)+1]] <- data.frame(
   dataset=m$dataset,region=region,group=group,
   supported_random_terms=paste(specification$active[[group]],collapse="+"),
   slope_mapping=as.character(toJSON(setNames(as.list(colnames(specification$X)[-1]),
                                             specification$zs),auto_unbox=TRUE)),
   fixed_design_columns=ncol(specification$X),fixed_design_rank=qr(specification$X)$rank)
  tab <- as.data.frame(table(z[,names(m$factors),drop=FALSE]),stringsAsFactors=FALSE)
  stopifnot(all(tab$Freq>0))
  for(i in seq_len(nrow(tab)))cells[[length(cells)+1]] <- data.frame(
   dataset=m$dataset,region=region,
   condition=paste(paste(names(m$factors),tab[i,names(m$factors)],sep="="),collapse="; "),
   observations=tab$Freq[i])
  fit <- readRDS(file.path(root,"fits",paste0(m$dataset,"__",region,".rds")))
  stopifnot(isTRUE(fit$complete),identical(fit$metadata$factors,m$factors),
   identical(fit$metadata$levels,m$levels),
   identical(fit$source_sha,m$sha256),
   identical(fit$input_fingerprint,m$region_fingerprints[[region]]))
  jobs[[length(jobs)+1]] <- data.frame(dataset=m$dataset,region=region,
   elapsed_seconds=fit$elapsed_seconds,optimizer=fit$optimizer,
   fit_attempts=nrow(fit$attempts),nominal_outputs=nrow(fit$status))
 }
 subjects <- unique(as.character(d$subject))
 privacy[[length(privacy)+1]] <- data.frame(dataset=m$dataset,
  email_like_subject_ids=sum(grepl("[[:alnum:]._%+-]+@[[:alnum:].-]+\\.[A-Za-z]{2,}",subjects)),
  human_name_like_subject_ids=sum(grepl("^[[:alpha:]]{2,}[[:space:]]+[[:alpha:]]{2,}",subjects)),
  subject_id_screen="Heuristic only; original source IDs retained, values not printed.")
 inventory[[length(inventory)+1]] <- data.frame(dataset=m$dataset,citation=m$citation,
  task=m$task,design=design,participants=length(subjects),items=length(unique(d$item)),
  observations=nrow(d),regions=length(m$region_counts),input_md5=m$input_md5,
  source=m$source,source_sha256=m$sha256)
 # Lossless decimal export for independent reconstruction comparison.
 for(f in names(d))if(is.factor(d[[f]]))d[[f]] <- as.character(d[[f]])
 d$RT <- sprintf("%.17g",d$RT)
 write.csv(d,file.path(out,"rds_exports",paste0(m$dataset,".csv")),row.names=FALSE)
}
for(name in c("inventory","levels","cells","jobs","privacy","coding","structures"))
 write.csv(do.call(rbind,get(name)),file.path(out,paste0(name,".csv")),row.names=FALSE)
jobs <- do.call(rbind,jobs)
timing <- list(
 complete_region_jobs=nrow(jobs),recorded_job_elapsed_seconds=sum(jobs$elapsed_seconds),
 recorded_job_elapsed_hours=sum(jobs$elapsed_seconds)/3600,
 longest_job_seconds=max(jobs$elapsed_seconds),
 longest_job=as.list(jobs[which.max(jobs$elapsed_seconds),c("dataset","region")]),
 interpretation=paste("Sum of recorded job elapsed times, including recorded resumptions.",
 "Not CPU time or end-to-end wall time; jobs overlapped in parallel.",
 "Excluded/abandoned runs and unrecorded setup time are not included.",
 "The study is an empirical multiverse reanalysis, not a Monte Carlo simulation."))
write_json(timing,file.path(out,"runtime_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat("Verified inventory and factor/region structure for",length(registry),"inputs and",nrow(jobs),"jobs.\n")
cat("Privacy heuristic flags:",sum(do.call(rbind,privacy)$email_like_subject_ids),
    "email-like;",sum(do.call(rbind,privacy)$human_name_like_subject_ids),"name-like identifiers.\n")
