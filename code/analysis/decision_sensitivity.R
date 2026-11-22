# Called by summarize.R with the verified, nonredundant coefficient table.
# No refits: compare workflows differing in one decision dimension only.
set.seed(20261121)
decision_sources <- list()
metrics <- c("disagreement","significant_gain","significant_loss",
             "absolute_t_shift","same_formula_fraction","same_formula_disagreement")
for(dimension in c("exclusion","transformation")) {
  varied <- if(dimension=="exclusion") "policy" else "transformation"
  held <- if(dimension=="exclusion") "transformation" else "policy"
  levels <- if(dimension=="exclusion") policies else types
  pairs <- combn(levels,2,simplify=FALSE)
  for(pair in pairs) {
    a <- coef[coef[[varied]]==pair[1],]
    b <- coef[coef[[varied]]==pair[2],]
    keep <- c("key","dataset","citation",held,"p","t","formula")
    d <- merge(a[,keep],b[,keep],by=c("key","dataset","citation",held),
               suffixes=c("_a","_b"))
    stopifnot(nrow(d)==nrow(effects)*if(dimension=="exclusion")5L else 13L)
    d$disagreement <- as.numeric((d$p_a<.05)!=(d$p_b<.05))
    d$significant_gain <- as.numeric(d$p_a>=.05 & d$p_b<.05)
    d$significant_loss <- as.numeric(d$p_a<.05 & d$p_b>=.05)
    d$absolute_t_shift <- abs(d$t_b-d$t_a)
    d$same_formula_fraction <- as.numeric(d$formula_a==d$formula_b)
    d$same_formula_disagreement <- ifelse(d$same_formula_fraction==1,
                                          d$disagreement,NA_real_)
    values <- lapply(metrics,function(m)source_values(d,m))
    names(values) <- metrics
    z <- data.frame(dimension=dimension,level_a=pair[1],level_b=pair[2],
                    citation=names(values$disagreement))
    for(m in metrics)z[[m]] <- values[[m]][z$citation]
    decision_sources[[length(decision_sources)+1]] <- z
  }
}
decision_sources <- do.call(rbind,decision_sources)
write.csv(decision_sources,file.path(out,"decision_pair_sources.csv"),row.names=FALSE)
vector_interval <- function(v) {
  v <- v[is.finite(v)]
  bs <- replicate(1000,mean(sample(v,length(v),replace=TRUE)))
  data.frame(mean=mean(v),lower=unname(quantile(bs,.025)),
             upper=unname(quantile(bs,.975)),sources=length(v))
}
summarize_decision <- function(z) {
  d <- vector_interval(z$disagreement)
  for(m in setdiff(metrics,"disagreement")) d[[m]] <- finite_mean(z[[m]])
  d
}
pair_groups <- split(decision_sources,
  interaction(decision_sources$dimension,decision_sources$level_a,
              decision_sources$level_b,drop=TRUE))
decision_pairs <- do.call(rbind,lapply(pair_groups,function(z)
  cbind(z[1,c("dimension","level_a","level_b")],summarize_decision(z))))
rownames(decision_pairs) <- NULL
write.csv(decision_pairs,file.path(out,"decision_pair_summary.csv"),row.names=FALSE)
decision_marginals <- list()
for(dimension in c("exclusion","transformation")) {
  z <- decision_sources[decision_sources$dimension==dimension,]
  # Average each distinct pair equally within source, then sources equally.
  for(level in c("ALL",if(dimension=="exclusion")policies else types)) {
    a <- if(level=="ALL")z else z[z$level_a==level|z$level_b==level,]
    b <- aggregate(a[,metrics],list(citation=a$citation),finite_mean)
    decision_marginals[[length(decision_marginals)+1]] <-
      cbind(data.frame(dimension=dimension,level=level),summarize_decision(b))
  }
}
decision_marginals <- do.call(rbind,decision_marginals)
decision_marginals$significant_gain <- NULL
decision_marginals$significant_loss <- NULL
write.csv(decision_marginals,file.path(out,"decision_marginals.csv"),row.names=FALSE)
source_dimension <- aggregate(decision_sources[,metrics],
  decision_sources[,c("citation","dimension")],finite_mean)
ex <- source_dimension[source_dimension$dimension=="exclusion",]
tr <- source_dimension[source_dimension$dimension=="transformation",]
delta <- vector_interval(ex$disagreement-
                          tr$disagreement[match(ex$citation,tr$citation)])
decision_summary <- list(
 exclusion=as.list(decision_marginals[decision_marginals$dimension=="exclusion" &
                                     decision_marginals$level=="ALL",]),
 transformation=as.list(decision_marginals[decision_marginals$dimension=="transformation" &
                                          decision_marginals$level=="ALL",]),
 exclusion_minus_transformation=as.list(delta),
 exclusion_pairs=choose(length(policies),2),
 transformation_pairs=choose(length(types),2),
 comparisons=as.integer(nrow(effects)*(choose(13,2)*5+choose(5,2)*13)))

# Explicit effect/interval summaries requested in the original paper outline.
width <- (coef$CI_high-coef$CI_low)*coef$standardized_SE/coef$SE
width_ratio <- tapply(width,coef$key,function(v)max(v)/min(v))
effects$CI_width_ratio <- width_ratio[paste(effects$dataset,effects$region,
                                           effects$term,sep="__")]
write.csv(effects,file.path(out,"effect_stability.csv"),row.names=FALSE)
effect_kind_summary <- do.call(rbind,lapply(split(effects,effects$effect_kind),function(z)
 data.frame(effect_kind=z$effect_kind[1],coefficient_regions=nrow(z),
  range_median=median(z$standardized_range),range_q25=quantile(z$standardized_range,.25),
  range_q75=quantile(z$standardized_range,.75),SE_ratio_median=median(z$SE_ratio),
  CI_width_ratio_median=median(z$CI_width_ratio))))
write.csv(effect_kind_summary,file.path(out,"effect_kind_summary.csv"),row.names=FALSE)
