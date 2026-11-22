# Reconstruct the relative-clause ROI using the supplied authors' label joins.
root <- "BRM impact/preprocessing_robustness/github-boyce_maze_2020/experiment"
out <- "generated-artifacts/brm_multiverse/boyce"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
labels <- as.data.frame(readRDS(file.path(root,"Labelled/rel_clause_labelled.rds")))
names(labels)[names(labels)=="item_num"] <- "group"
keys <- c("word_num","sentence","type","word","group")
stopifnot(!anyDuplicated(labels[keys]))
for(edition in c("spr","l_maze","g_maze","one_b_maze","gulo_maze")) {
  d <- as.data.frame(readRDS(file.path(root,"Data/Processed",paste0(edition,".rds"))))
  d <- d[d$citizen=="yes" & d$native=="yes" & d$resident=="yes",]
  d <- d[!is.na(d$type) & d$type%in%c("relative_high","relative_low"),]
  if(edition=="spr") {
    d <- d[!is.na(d$accuracy)&d$accuracy>.8,]
    d$word_num <- d$word_num-1
  } else {
    d <- d[!is.na(d$correct)&d$correct=="yes",]
  }
  d <- merge(d,labels,by=keys,all.x=TRUE)
  d <- d[!is.na(d$mod_word_number)&d$mod_word_number%in%c(0,1),]
  d <- d[is.finite(d$rt)&d$rt>0,]
  d <- d[,c("subject","group","type","rt","mod_word_number")]
  write.csv(d,file.path(out,paste0(edition,".csv")),row.names=FALSE)
  cat(edition,nrow(d),"ROI observations\n")
}
