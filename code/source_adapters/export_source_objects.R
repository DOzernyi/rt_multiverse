out <- "generated-artifacts/brm_multiverse/source_exports"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
base <- "BRM impact/preprocessing_robustness/f8qwh-paape_local_2025"
env <- new.env(parent=emptyenv())
load(file.path(base,"Experiment 1 (SPR)/lc_spr_trialinfo2.Rda"),envir=env)
stopifnot(is.data.frame(env$lc_spr_trialinfo2))
write.csv(env$lc_spr_trialinfo2,file.path(out,"paape_spr2.csv"),row.names=FALSE)
