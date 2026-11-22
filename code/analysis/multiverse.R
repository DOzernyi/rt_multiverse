# Reproducible exclusion x transformation multiverse. No source data are edited.
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(jsonlite)
})
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath("BRM impact", mustWork=TRUE)
out <- normalizePath("generated-artifacts/brm_multiverse", mustWork=TRUE)
registry <- fromJSON(file.path(out, "datasets.json"), simplifyVector=FALSE)
exclusion_file <- file.path(root, "analysis", "excluded_datasets.json")
exclusion_records <- if (file.exists(exclusion_file))
  fromJSON(exclusion_file, simplifyVector=FALSE) else list()
excluded_datasets <- if (length(exclusion_records))
  vapply(exclusion_records, `[[`, character(1), "dataset") else character()
registry <- Filter(function(m) !(m$dataset %in% excluded_datasets), registry)
workers <- if (length(args)) as.integer(args[1]) else 2L
only <- if (length(args)>1) args[2] else "all"
dir.create(file.path(out, "fits"), showWarnings=FALSE)
policies <- c("none", paste0("sd", c(1.5,2,2.5,3)),
              "fixed100_2000","fixed100_3000","fixed200_2000","fixed200_3000",
              "median", "mad2", "mad2.5", "mad3")
transforms <- c("identity","ln","log10","inverse","sqrt","boxcox")
version <- "multiverse-4"
optimizer <- "nloptwrap"

screen <- function(d, policy) {
  keep <- rep(TRUE,nrow(d))
  groups <- interaction(d[,grep("^F[0-9]",names(d)),drop=FALSE], drop=TRUE)
  for (ix in split(seq_len(nrow(d)), groups)) {
    x <- d$RT[ix]
    bounds <- c(-Inf,Inf)
    if (startsWith(policy,"fixed")) {
      bounds <- as.numeric(strsplit(sub("fixed","",policy),"_",fixed=TRUE)[[1]])
    } else if (startsWith(policy,"sd")) {
      k <- as.numeric(sub("sd","",policy))
      bounds <- mean(x) + c(-1,1)*k*sd(x)
    } else if (startsWith(policy,"mad")) {
      k <- as.numeric(sub("mad","",policy))
      bounds <- median(x) + c(-1,1)*k*mad(x,constant=1.4826)
    } else if (policy=="median") {
      # rt_process a=1.8 is an additive median-distance parameter:
      # M - .5*M <= RT <= M + 1.8*M, i.e. [.5*M, 2.8*M].
      bounds <- median(x)*c(.5,2.8)
    }
    keep[ix] <- x >= bounds[1] & x <= bounds[2]
  }
  keep
}
bc <- function(x, lambda) {
  if (abs(lambda)<1e-7) log(x) else expm1(lambda*log(x))/lambda
}
transform_rt <- function(x, type, X) {
  lambda <- NA_real_
  if(type=="boxcox") {
    # Fixed-design Gaussian profile with the required Jacobian. This estimates
    # one lambda per dataset/ROI/policy, not separate lambdas per condition.
    Q <- qr(X)
    objective <- function(l) {
      z <- bc(x,l)
      rss <- sum(qr.resid(Q,z)^2)
      if(!is.finite(rss) || rss<=0) return(Inf)
      length(x)/2*log(rss/length(x))-(l-1)*sum(log(x))
    }
    optimum <- optimize(objective,c(-2,2),tol=1e-5)
    lambda <- optimum$minimum
    y <- bc(x,lambda)
  } else y <- switch(type, identity=x, ln=log(x), log10=log10(x),
                      inverse=-1000/x,sqrt=sqrt(x))
  list(y=y,lambda=lambda)
}
formula_for <- function(fixed, active, correlated) {
  parts <- character()
  for(g in names(active)) {
    v <- active[[g]]
    if(!length(v)) next
    if(correlated) {
      parts <- c(parts,paste0("(",paste(v,collapse="+"),"|",g,")"))
    } else {
      parts <- c(parts,vapply(v,function(t) {
        if(t=="1") paste0("(1|",g,")") else paste0("(0+",t,"|",g,")")
      },character(1)))
    }
  }
  if(!length(parts)) return(NULL)
  as.formula(paste("y ~",fixed,"+",paste(parts,collapse="+")))
}
variance_candidates <- function(fit) {
  v <- as.data.frame(VarCorr(fit))
  v <- v[is.na(v$var2) & v$grp!="Residual",]
  v$group <- sub("\\.[0-9]+$","",v$grp)
  v$term <- ifelse(v$var1=="(Intercept)","1",v$var1)
  v[order(v$vcov,v$group,v$term),]
}
fit_sequence <- function(d, fixed, active) {
  attempts <- list()
  correlated <- TRUE
  selected <- NULL
  for(step in seq_len(40L)) {
    form <- formula_for(fixed,active,correlated)
    if(is.null(form)) break
    # lmerTest reconstructs the deviance function from the model call. Retain
    # the fitting data/formula environment after this helper returns.
    environment(form) <- environment()
    warning_text <- character()
    # Derivatives are omitted for rejecting preliminary singular fits, then
    # computed for every provisionally acceptable final fit.
    fit <- tryCatch(withCallingHandlers(
      lme4::lmer(form,d,REML=TRUE,
        control=lmerControl(optimizer=optimizer,calc.derivs=FALSE,
                            optCtrl=list(maxeval=20000))),
      warning=function(w) {warning_text<<-c(warning_text,conditionMessage(w));invokeRestart("muffleWarning")}),
      error=function(e)e)
    error <- inherits(fit,"error")
    singular <- if(error) NA else isSingular(fit,tol=1e-4)
    good <- !error && fit@optinfo$conv$opt==0 && !singular
    if(good) {
      fit <- tryCatch(withCallingHandlers(
        update(fit,start=getME(fit,"theta"),
          control=lmerControl(optimizer=optimizer,calc.derivs=TRUE,
                              optCtrl=list(maxeval=20000))),
        warning=function(w) {warning_text<<-c(warning_text,conditionMessage(w));invokeRestart("muffleWarning")}),
        error=function(e)e)
      error <- inherits(fit,"error")
      singular <- if(error) NA else isSingular(fit,tol=1e-4)
      messages <- if(error) NULL else fit@optinfo$conv$lme4$messages
      good <- !error && fit@optinfo$conv$opt==0 && is.null(messages) && !singular
    }
    attempts[[step]] <- data.frame(step=step,
      formula=paste(deparse(form,width.cutoff=500L),collapse=""),
      correlated=correlated, accepted=good, singular=singular,
      error=if(error) conditionMessage(fit) else "",
      warnings=paste(unique(warning_text),collapse=" | "),
      removed_group="",removed_term="",removed_variance=NA_real_)
    if(good) {selected<-fit;break}
    if(correlated) {correlated<-FALSE;next}
    if(error) break # No estimable variances: do not invent a reduction ordering.
    v <- variance_candidates(fit)
    v <- v[mapply(function(g,t) g%in%names(active) && t%in%active[[g]],v$group,v$term),]
    if(!nrow(v)) break
    g <- v$group[1]; t <- v$term[1]
    attempts[[step]]$removed_group <- g
    attempts[[step]]$removed_term <- t
    attempts[[step]]$removed_variance <- v$vcov[1]
    active[[g]] <- setdiff(active[[g]],t)
  }
  list(model=selected,attempts=do.call(rbind,attempts))
}
prepare_design <- function(d, metadata) {
  fs <- names(metadata$factors)
  for(f in fs) {
    raw <- as.character(d[[f]])
    lev <- as.character(unlist(metadata$levels[[f]]))
    raw_num <- suppressWarnings(as.numeric(raw))
    lev_num <- suppressWarnings(as.numeric(lev))
    numeric_levels <- length(lev_num)>0 && all(is.finite(lev_num)) &&
      all(is.na(raw) | is.finite(raw_num))
    if(numeric_levels) {
      d[[f]] <- factor(raw_num,levels=lev_num,labels=lev)
    } else {
      d[[f]] <- factor(raw,levels=lev)
    }
    # Binary coefficients represent the full level2-level1 difference.
    # Higher-level factors use explicit sum contrasts /2, retaining K-1 dfs.
    contrasts(d[[f]]) <- if(nlevels(d[[f]])==2) matrix(c(-.5,.5),ncol=1) else contr.sum(nlevels(d[[f]]))/2
  }
  fixed <- paste(fs,collapse="*")
  X <- model.matrix(as.formula(paste("~",fixed)),d)
  zs <- paste0("z",seq_len(ncol(X)-1))
  for(i in seq_along(zs)) d[[zs[i]]] <- X[,i+1]
  terms <- attr(terms(as.formula(paste("~",fixed))),"term.labels")
  assignment <- attr(X,"assign")[-1]
  slope_factors <- lapply(assignment,function(i) strsplit(terms[i],":",fixed=TRUE)[[1]])
  active <- list()
  for(g in c("subject","item")) {
    d[[g]] <- factor(d[[g]])
    within <- vapply(fs,function(f) {
      any(tapply(d[[f]],d[[g]],function(z)length(unique(z)))>1)
    },logical(1))
    active[[g]] <- c("1", zs[vapply(slope_factors,function(ff)all(within[ff]),logical(1))])
  }
  list(d=d,X=X,fixed=fixed,active=active,zs=zs)
}
run_one <- function(job) {
  m <- job$metadata; roi <- job$region
  file <- file.path(out,"fits",paste0(m$dataset,"__",roi,".rds"))
  fingerprint <- m$region_fingerprints[[roi]]
  resumed <- NULL
  if(file.exists(file)) {
    old <- readRDS(file)
    compatible <- identical(old$version,version) &&
      identical(old$source_sha,m$sha256) &&
      !is.null(fingerprint) && identical(old$input_fingerprint,fingerprint)
    # Bootstrap only the source-verified, numerically unchanged critical
    # adapters from the earlier run. No ellipsis or postcritical cache migration.
    legacy <- identical(old$version,"multiverse-3") &&
      isTRUE(m$legacy_critical_unchanged) && identical(roi,"critical") &&
      identical(old$source_sha,m$sha256) &&
      identical(old$metadata$factors,m$factors) &&
      identical(old$metadata$levels,m$levels) &&
      identical(as.integer(old$metadata$region_counts[[roi]]),
                as.integer(m$region_counts[[roi]]))
    if(legacy) {
      old$version<-version;old$input_fingerprint<-fingerprint
      old$metadata<-m;old$complete<-TRUE
      old$optimizer<-"bobyqa (unchanged-input legacy cache)"
      saveRDS(old,file,compress="gzip")
      compatible<-TRUE
    }
    if(compatible && isTRUE(old$complete)) {
      cat("CACHED",m$dataset,roi,"\n");return(TRUE)
    }
    if(compatible) resumed<-old
  }
  input <- if(grepl("\\.rds$",m$data)) readRDS(file.path(out,m$data))$data else
    read.csv(file.path(out,m$data),stringsAsFactors=FALSE,check.names=FALSE)
  input <- input[input$region==roi,]
  results <- list(); traces <- list(); statuses <- list()
  random_variances <- list(); nested <- list()
  if(!is.null(resumed)) {
    results<-resumed$partial$results;traces<-resumed$partial$traces
    statuses<-resumed$partial$statuses
    random_variances<-resumed$partial$random_variances
    nested<-resumed$partial$nested
  }
  completed_policies<-if(is.null(resumed)) character() else resumed$completed_policies
  start <- Sys.time()
  elapsed_before<-if(is.null(resumed)) 0 else resumed$elapsed_seconds
  persist <- function(complete=FALSE) {
    answer <- list(version=version,source_sha=m$sha256,metadata=m,region=roi,
      input_fingerprint=fingerprint,optimizer=optimizer,complete=complete,
      completed_policies=completed_policies,
      coefficients=do.call(rbind,results),status=do.call(rbind,statuses),
      attempts=do.call(rbind,traces),random_variances=do.call(rbind,random_variances),
      conditional_followups=do.call(rbind,nested),
      elapsed_seconds=elapsed_before+as.numeric(difftime(Sys.time(),start,units="secs")))
    if(!complete) answer$partial<-list(results=results,traces=traces,statuses=statuses,
      random_variances=random_variances,nested=nested)
    tmp<-paste0(file,".tmp-",Sys.getpid())
    saveRDS(answer,tmp,compress="gzip")
    if(!file.rename(tmp,file)) stop("Could not atomically save job checkpoint")
    answer
  }
  for(policy in setdiff(policies,completed_policies)) {
    keep <- screen(input,policy)
    d <- input[keep,]
    if(nrow(d)<30) {
      for(type in transforms) statuses[[paste(policy,type,sep="__")]]<-data.frame(
        policy=policy,transformation=type,success=FALSE,n=nrow(d),
        excluded=sum(!keep),proportion_excluded=mean(!keep),lambda=NA_real_,
        lambda_at_boundary=FALSE,attempts=0)
      completed_policies<-c(completed_policies,policy);persist()
      next
    }
    design <- prepare_design(d,m)
    d <- design$d
    if(qr(design$X)$rank<ncol(design$X)) {
      for(type in transforms) statuses[[paste(policy,type,sep="__")]]<-data.frame(
        policy=policy,transformation=type,success=FALSE,n=nrow(d),
        excluded=sum(!keep),proportion_excluded=mean(!keep),lambda=NA_real_,
        lambda_at_boundary=FALSE,attempts=0)
      completed_policies<-c(completed_policies,policy);persist()
      next
    }
    # Only ln is numerically fitted. log10 outputs follow the exact analytical
    # unit conversion: identical p/t/df and standardized effect, scaled beta/SE.
    for(type in setdiff(transforms,"log10")) {
      label <- paste(policy,type,sep="__")
      tr <- transform_rt(d$RT,type,design$X)
      center <- mean(tr$y); unit <- sd(tr$y)
      if(!is.finite(unit)||unit==0) {
        statuses[[label]]<-data.frame(policy=policy,transformation=type,
          success=FALSE,n=nrow(d),excluded=sum(!keep),
          proportion_excluded=mean(!keep),lambda=tr$lambda,
          lambda_at_boundary=is.finite(tr$lambda)&&abs(tr$lambda)>1.99,attempts=0)
        if(type=="ln") {
          ls<-statuses[[label]];ls$transformation<-"log10"
          statuses[[paste(policy,"log10",sep="__")]]<-ls
        }
        next
      }
      d$y <- (tr$y-center)/unit
      fitted <- fit_sequence(d,design$fixed,design$active)
      trace <- fitted$attempts
      trace$policy <- policy; trace$transformation <- type
      traces[[label]] <- trace
      statuses[[label]] <- data.frame(policy=policy,transformation=type,
        success=!is.null(fitted$model), n=nrow(d), excluded=sum(!keep),
        proportion_excluded=mean(!keep), lambda=tr$lambda,
        lambda_at_boundary=is.finite(tr$lambda) && abs(tr$lambda)>1.99,
        attempts=nrow(trace))
      if(!is.null(fitted$model)) {
        model <- fitted$model
        tested <- tryCatch(as_lmerModLmerTest(model),error=function(e)NULL)
        stats <- if(is.null(tested)) NULL else coef(summary(tested))
        if(is.null(stats)) {
          statuses[[label]]$success<-FALSE
          if(type=="ln") {
            ls<-statuses[[label]];ls$transformation<-"log10"
            statuses[[paste(policy,"log10",sep="__")]]<-ls
          }
          next
        }
        native <- stats
        native[,"Estimate"] <- stats[,"Estimate"]*unit
        native[1,"Estimate"] <- native[1,"Estimate"]+center
        native[,"Std. Error"] <- stats[,"Std. Error"]*unit
        rr <- data.frame(term=rownames(stats),estimate=native[,"Estimate"],
          SE=native[,"Std. Error"],df=stats[,"df"],t=stats[,"t value"],
          p=stats[,"Pr(>|t|)"],standardized_estimate=stats[,"Estimate"],
          standardized_SE=stats[,"Std. Error"], stringsAsFactors=FALSE)
        # Centering changes the intercept's zero hypothesis. Test the restored
        # native intercept against native zero, not the centered mean.
        rr$t[1] <- rr$estimate[1]/rr$SE[1]
        rr$p[1] <- 2*pt(-abs(rr$t[1]),rr$df[1])
        rr$CI_low <- rr$estimate-qt(.975,rr$df)*rr$SE
        rr$CI_high <- rr$estimate+qt(.975,rr$df)*rr$SE
        rr$direction <- sign(rr$estimate)
        rr$policy<-policy;rr$transformation<-type
        rr$n<-nrow(d);rr$excluded<-sum(!keep);rr$lambda<-tr$lambda
        rr$residual_SD <- sigma(model)*unit
        rr$formula <- paste(deparse(formula(model),width.cutoff=500),collapse="")
        results[[label]] <- rr
        vv <- as.data.frame(VarCorr(model))
        vv$policy <- policy;vv$transformation<-type
        vv$vcov_native <- vv$vcov*unit^2
        random_variances[[label]] <- vv
        # Conditional follow-ups are computed as equivalent reparameterized
        # contrasts of the accepted model, not a second arbitrary random model.
        interactions <- grepl(":",rr$term)
        if(any(interactions & rr$p<.05,na.rm=TRUE)) {
          fs <- names(m$factors)
          grid <- expand.grid(lapply(fs,function(f)levels(d[[f]])),stringsAsFactors=FALSE)
          names(grid)<-fs
          for(f in fs) {
            grid[[f]]<-factor(grid[[f]],levels=levels(d[[f]]))
            contrasts(grid[[f]])<-contrasts(d[[f]])
          }
          G <- model.matrix(as.formula(paste("~",design$fixed)),grid)
          followups<-list()
          for(f in fs) {
            rest<-setdiff(fs,f)
            combinations<-if(length(rest)) unique(grid[,rest,drop=FALSE]) else data.frame(dummy=1)
            for(ci in seq_len(nrow(combinations))) {
              in_cell<-rep(TRUE,nrow(grid))
              for(other in rest) in_cell<-in_cell & grid[[other]]==combinations[[other]][ci]
              indices<-which(in_cell)
              for(pair in combn(indices,2,simplify=FALSE)) {
                L<-G[pair[2],]-G[pair[1],]
                ct<-tryCatch(contest1D(tested,L),error=function(e)NULL)
                if(is.null(ct)) next
                followups[[length(followups)+1]]<-data.frame(
                  contrast=paste(f,as.character(grid[[f]][pair[2]]),"-",
                    as.character(grid[[f]][pair[1]]),"given",
                    paste(paste(rest,combinations[ci,rest,drop=TRUE],sep="="),collapse=";")),
                  estimate=ct[1,"Estimate"]*unit,SE=ct[1,"Std. Error"]*unit,
                  df=ct[1,"df"],t=ct[1,"t value"],p=ct[1,"Pr(>|t|)"],
                  policy=policy,transformation=type)
              }
            }
          }
          if(length(followups)) nested[[label]]<-do.call(rbind,followups)
        }
        if(type=="ln") {
          lr <- rr
          for(col in c("estimate","SE","CI_low","CI_high","residual_SD")) lr[[col]]<-lr[[col]]/log(10)
          lr$transformation <- "log10"
          results[[paste(policy,"log10",sep="__")]] <- lr
          ls <- statuses[[label]];ls$transformation<-"log10"
          statuses[[paste(policy,"log10",sep="__")]]<-ls
        }
      } else if(type=="ln") {
        ls <- statuses[[label]];ls$transformation<-"log10"
        statuses[[paste(policy,"log10",sep="__")]]<-ls
      }
    }
    cat(format(Sys.time()),m$dataset,roi,policy,"completed\n")
    flush.console()
    completed_policies<-c(completed_policies,policy)
    persist()
  }
  answer<-persist(complete=TRUE)
  cat("DONE",m$dataset,roi,round(answer$elapsed_seconds,1),"seconds\n")
  TRUE
}
if(sys.nframe()==0) {
  jobs<-list()
  for(m in registry) {
    if(only!="all" && m$tier!=only && m$dataset!=only) next
    for(roi in names(m$region_counts)) jobs[[length(jobs)+1]]<-list(metadata=m,region=roi)
  }
  # Small jobs finish first; source strata remain separate in every summary.
  jobs<-jobs[order(vapply(jobs,function(j)j$metadata$observations,numeric(1)))]
  cat("Starting",length(jobs),"dataset/ROI jobs with",workers,"workers\n")
  if(!length(jobs)) stop("No dataset/ROI jobs matched the requested stratum")
  answers<-parallel::mclapply(jobs,run_one,mc.cores=workers,mc.preschedule=FALSE)
  failed<-which(!vapply(answers,isTRUE,logical(1)))
  if(length(failed)) stop(length(failed)," dataset/ROI jobs failed; inspect explicit worker errors and saved checkpoints")
  cat("FINISHED",length(jobs),"complete dataset/ROI jobs\n")
}
