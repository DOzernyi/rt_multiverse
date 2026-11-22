# Compact, reusable analysis inputs. No transformations/outlier rules here.
suppressPackageStartupMessages(library(jsonlite))
out <- normalizePath("generated-artifacts/brm_multiverse", mustWork=TRUE)
registry <- fromJSON(file.path(out,"datasets.json"),simplifyVector=FALSE)
dir.create(file.path(out,"primary_inputs"),showWarnings=FALSE)
for (i in seq_along(registry)) {
  m <- registry[[i]]
  if (!identical(m$tier,"primary")) next
  d <- if(grepl("\\.rds$",m$data)) readRDS(file.path(out,m$data))$data else
    read.csv(file.path(out,m$data),stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(nrow(d)==m$observations,all(is.finite(d$RT)),all(d$RT>0))
  for (name in c("subject","item","region",names(m$factors))) d[[name]]<-factor(d[[name]])
  path <- file.path("primary_inputs",paste0(m$dataset,".rds"))
  metadata <- m
  metadata$data <- path
  metadata$input_md5 <- NULL
  metadata$region_fingerprints <- NULL
  metadata$input_definition <- paste(
    "Original-scale positive recorded RTs; source-defined regions and common",
    "task eligibility only. No multiverse outlier removal or transformation.")
  fingerprints <- list()
  tmp <- tempfile(tmpdir=out)
  for(roi in names(m$region_counts)) {
    local <- droplevels(d[d$region==roi,,drop=FALSE])
    rownames(local) <- NULL
    saveRDS(list(data=local,factors=m$factors,levels=m$levels),tmp,compress=FALSE)
    fingerprints[[roi]] <- unname(tools::md5sum(tmp))
  }
  unlink(tmp)
  connection <- xzfile(file.path(out,path),"wb",compression=3)
  saveRDS(list(data=d,metadata=metadata),connection,version=3)
  close(connection)
  restored <- readRDS(file.path(out,path))
  stopifnot(identical(restored$data,d))
  registry[[i]]$data <- path
  registry[[i]]$input_md5 <- unname(tools::md5sum(file.path(out,path)))
  registry[[i]]$region_fingerprints <- fingerprints
  cat(m$dataset,nrow(d),file.info(file.path(out,path))$size,"bytes\n")
}
write_json(registry,file.path(out,"datasets.json"),auto_unbox=TRUE,pretty=TRUE,null="null")
