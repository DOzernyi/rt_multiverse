# Independent base-R reconstruction of all 468 plotted Python summary values.
root <- "generated-artifacts/brm_multiverse"
out <- file.path(root,"effect_sensitivity")
c <- read.csv(gzfile(file.path(root,"coefficients.csv.gz")),stringsAsFactors=FALSE)
c <- c[c$term!="(Intercept)",]
keys <- c("dataset","region","term")
workflow <- c("policy","transformation")
base <- c[c$policy=="none"&c$transformation=="identity",]
c$width <- (c$CI_high-c$CI_low)*c$standardized_SE/c$SE
base$width <- (base$CI_high-base$CI_low)*base$standardized_SE/base$SE
d <- merge(c,base[,c(keys,"standardized_estimate","standardized_SE","width")],
           by=keys,suffixes=c("",".reference"))
d$effect_shift <- abs(d$standardized_estimate-d$standardized_estimate.reference)
d$SE_change <- abs(100*(d$standardized_SE/d$standardized_SE.reference-1))
d$CI_width_change <- abs(100*(d$width/d$width.reference-1))
d$interaction_shift <- abs(abs(d$standardized_estimate)-abs(d$standardized_estimate.reference))

publication_average <- function(z,metric) {
 z$value <- z[[metric]]
 input <- aggregate(value~policy+transformation+dataset+citation,z,mean)
 source <- aggregate(value~policy+transformation+citation,input,mean)
 result <- aggregate(value~policy+transformation,source,mean)
 result$metric <- metric
 result
}
expected <- list()
for(metric in c("effect_shift","SE_change","CI_width_change","interaction_shift")) {
 z <- if(metric=="interaction_shift")d[grepl(":",d$term,fixed=TRUE),] else d
 expected[[metric]] <- publication_average(z,metric)
}
v <- read.csv(gzfile(file.path(root,"random_variances.csv.gz")),stringsAsFactors=FALSE)
v <- v[v$grp!="Residual",]
v$group <- sub("\\.[0-9]+$","",v$grp)
b <- ifelse(is.na(v$var2),v$var1,v$var2)
v$component <- paste(v$group,ifelse(v$var1<=b,v$var1,b),
                    ifelse(v$var1<=b,b,v$var1),sep="|")
lv <- v[v$transformation=="ln",];lv$transformation <- "log10"
v <- rbind(v,lv)
reference <- v[v$policy=="none"&v$transformation=="identity",]
z <- merge(v,reference[,c("dataset","region","component","vcov")],
           by=c("dataset","region","component"),suffixes=c("",".reference"))
z$random_covariance_shift <- abs(z$vcov-z$vcov.reference)
jobs <- aggregate(random_covariance_shift~dataset+region+policy+transformation+citation,z,mean)
expected$random_covariance_shift <- publication_average(jobs,"random_covariance_shift")

id <- function(z)paste(z$dataset,z$region,z$policy,z$transformation,sep="||")
sets <- split(v$component,id(v))
unique_jobs <- unique(v[,c("dataset","region","policy","transformation","citation")])
structural <- vapply(seq_len(nrow(unique_jobs)),function(i){
 z <- unique_jobs[i,]
 current <- sets[[id(z)]]
 z$policy <- "none";z$transformation <- "identity"
 ref <- sets[[id(z)]]
 100*(1-length(intersect(current,ref))/length(union(current,ref)))
},numeric(1))
unique_jobs$random_structure_change <- structural
expected$random_structure_change <- publication_average(unique_jobs,"random_structure_change")
expected <- do.call(rbind,expected)
actual <- read.csv(file.path(out,"figure_values.csv"),stringsAsFactors=FALSE)
check <- merge(expected,actual,by=c("metric","policy","transformation"))
stopifnot(nrow(check)==468L,all(is.finite(check$mean)),
          max(abs(check$value-check$mean))<1e-10)
zero <- actual[actual$policy=="none"&actual$transformation=="identity",]
stopifnot(max(abs(zero$mean))<1e-12)
ln <- actual[actual$transformation=="ln",]
lg <- actual[actual$transformation=="log10",]
eq <- merge(ln,lg,by=c("metric","policy"),suffixes=c(".ln",".log10"))
stopifnot(max(abs(eq$mean.ln-eq$mean.log10))<1e-10)
writeLines(c("PASSED: independently reconstructed all 468 plotted summary values.",
 "PASSED: reference cells are zero and log10 equals ln after standardization.",
 "Random-effect comparison canonicalized participant/item and named terms;",
 "no missing random-effect entry was imputed as zero.",
 paste("Maximum absolute R/Python discrepancy:",format(max(abs(check$value-check$mean)),digits=6))),
 file.path(out,"verification.txt"))
cat(readLines(file.path(out,"verification.txt")),sep="\n")
