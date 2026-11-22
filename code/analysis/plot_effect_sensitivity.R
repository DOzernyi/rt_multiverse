# Print-ready single-page A4 landscape figure; no additional R packages.
out <- "generated-artifacts/brm_multiverse/effect_sensitivity"
d <- read.csv(file.path(out,"figure_values.csv"),stringsAsFactors=FALSE)
policies <- c("none","sd1.5","sd2","sd2.5","sd3","fixed100_2000",
 "fixed100_3000","fixed200_2000","fixed200_3000","median","mad2","mad2.5","mad3")
scales <- c("identity","ln","log10","inverse","sqrt","boxcox")
labels <- c("No removal","1.5 SD","2 SD","2.5 SD","3 SD","100-2000",
 "100-3000","200-2000","200-3000","Median","2 MAD","2.5 MAD","3 MAD")
metrics <- c("effect_shift","SE_change","CI_width_change",
 "random_covariance_shift","interaction_shift","random_structure_change")
titles <- c("A  Fixed-effect shift (outcome SD)",
 "B  SE change (%; standardized)","C  95% CI-width change (%; standardized)",
 "D  Random-effect estimate shift (SD squared)",
 "E  Interaction-magnitude shift (outcome SD)",
 "F  Random-effect entries not shared (%)")
sequential <- colorRampPalette(c("#eff1f3","#bdd1df","#8ab1cc","#5891b9","#2672a6"))(101)
render <- function() {
 par(mfrow=c(2,3),oma=c(4.8,0,.2,.4),family="Helvetica",
     mgp=c(1.8,.5,0),tcl=-.2,xaxs="i",yaxs="i")
 for(k in seq_along(metrics)) {
  col <- (k-1)%%3+1
  par(mar=c(4.7,if(col==1)7 else 1.2,2.3,.9))
  z <- d[d$metric==metrics[k],]
  m <- matrix(NA_real_,13,6)
  for(i in 1:13)for(j in 1:6)
   m[i,j] <- z$mean[z$policy==policies[i]&z$transformation==scales[j]]
  if(k%in%c(2,3)) {
   # Identical scales in the SE and interval-width panels.
   lim <- ceiling(max(abs(d$mean[d$metric%in%c("SE_change","CI_width_change")]))/5)*5
   limits <- c(0,lim); palette <- sequential
  } else {
   lim <- if(k==6)100 else if(k%in%c(1,5))
    max(d$mean[d$metric%in%c("effect_shift","interaction_shift")])*1.015 else max(m)*1.015
   limits <- c(0,lim);palette <- sequential
  }
  image(1:6,1:13,t(m[13:1,]),col=palette,zlim=limits,
    xlim=c(.5,6.5),ylim=c(.5,13.5),axes=FALSE,xlab="",ylab="")
  abline(v=seq(.5,6.5,1),h=seq(.5,13.5,1),col="white",lwd=.55)
  for(i in 1:13)for(j in 1:6) {
   value <- m[i,j]
   fmt <- if(k%in%c(1,4,5))"%.3f" else if(k%in%c(2,3))"%.1f" else "%.0f"
   label <- if(abs(value)<1e-11)if(k%in%c(1,4,5))"0.000" else "0" else sprintf(fmt,value)
   position <- (value-limits[1])/diff(limits)
   ink <- if(position>.72)"white" else "#173B54"
   text(j,14-i,label,cex=.77,col=ink)
  }
  # First row = transformation alone; first column = exclusion alone.
  rect(.5,12.5,6.5,13.5,lwd=1.05)
  rect(.5,.5,1.5,13.5,lwd=1.05)
  axis(1,at=1:6,labels=c("Raw","ln","log10","Inverse","sqrt","Box-Cox"),
       tick=FALSE,line=-.25,cex.axis=.72)
  if(col==1)axis(2,at=13:1,labels=labels,tick=FALSE,las=1,cex.axis=.76)
  mtext(titles[k],3,line=.65,adj=0,cex=.82)
  par(xpd=NA)
  xx <- seq(.5,6.5,length.out=102)
  for(b in 1:101)rect(xx[b],-1.3,xx[b+1],-.95,col=palette[b],border=NA)
  ticks <- c(limits[1],mean(limits),limits[2])
  tx <- .5+6*(ticks-limits[1])/diff(limits)
  tl <- if(k%in%c(1,4,5))sprintf("%.3f",ticks) else sprintf("%.0f%%",ticks)
  text(tx,-1.65,tl,adj=c(.5,.5),cex=.65)
  par(xpd=FALSE)
 }
 mtext("Reference: no removal + untransformed RTs. Outlined top row: transformations alone; outlined first column: exclusions alone; other cells: combinations.",
       1,outer=TRUE,line=1.15,adj=.015,cex=.69)
 mtext("Mean absolute changes weight publications equally. A-C: 380 coefficient-regions; E: 150 interactions. A/E share a colour scale; B/C share a scale.",
       1,outer=TRUE,line=2.25,adj=.015,cex=.69)
 mtext("D compares shared variance/covariance entries only; missing entries are not zero-imputed (see F). Log10 duplicates ln. Source-bootstrap intervals and full definitions accompany the figure.",
       1,outer=TRUE,line=3.35,adj=.015,cex=.66)
}
pdf(file.path(out,"preprocessing_effects.pdf"),width=297/25.4,height=210/25.4,
    family="Helvetica",useDingbats=FALSE,pointsize=10)
render();dev.off()
postscript(file.path(out,"preprocessing_effects.eps"),width=297/25.4,height=210/25.4,
           onefile=FALSE,horizontal=FALSE,paper="special",family="Helvetica",
           pointsize=10,colormodel="srgb")
render();dev.off()
cat("Saved one-page vector PDF and EPS\n")
