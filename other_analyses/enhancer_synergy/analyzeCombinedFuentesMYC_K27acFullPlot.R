library(dplyr)
library(ggplot2)
library(data.table)
library(nloptr)
library("viridis") 
library(scales)
library(ggpubr)

# Always run:

myc <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Fulco2016_DHSPeaks_noDistFilter/H3K27ac.table.txt"
ltr <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Fuentes2018_DHSPeaks_noDistFilter/H3K27ac.table.txt"
ltr.regions <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Fuentes2018_DHSPeaks_noDistFilter/cas9.regions.filtered.tsv"
out <-'/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/K27ac_plots/DHS_noDistFilter_LTRFilter/'

# myc <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Fulco2016_DHSPeaksComparisonTest/H3K27ac.table.txt"
# ltr <- "/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/Fuentes2018_DHSPeaks/H3K27ac.table.txt"
# out <-'/oak/stanford/groups/engreitz/Users/sheth/EnhancerSynergy/K27ac_plots/DHS/'

dir.create(out, showWarnings = FALSE)

setwd(out)
ltr <- fread(ltr)
ltr.regions <- fread(ltr.regions)
myc <- fread(myc)

ltr$distance <- ltr$enh.krab.dist
ltr$H3K27ac.KO.RPM <- ltr$H3K27ac_KRAB.RPM
ltr$H3K27ac.RPM.ctrl <- ltr$H3K27ac.RPM
ltr$H3K27ac.RPM.ctrl.at.PerturbationSite <- ltr$H3K27ac.RPM.of.cas9.region
ltr$H3K27ac.RPM.KO.at.PerturbationSite <- ltr$H3K27ac_KRAB.RPM.of.cas9.region
ltr$Reference <- "Fuentes2018"
myc$Reference <- "Fulco2016"

print(nrow(ltr))
ltr = subset(ltr, class != "promoter" & is.na(region.name))
print(nrow(ltr))
ltr = dplyr::filter(ltr, region.name.of.cas9.region %in% ltr.regions$region.name)
message("LTR-peak pairs after filtering: ", nrow(ltr))

print(nrow(myc))
combined <- rbind(myc, ltr, fill=TRUE)
print(nrow(combined))

compute.sem.ratio <- function(mean.x, sem.x, mean.y, sem.y) {
  #Calculate sem(x/y)
  
  sem.ratio <- (mean.x/mean.y) * sqrt((sem.x/mean.x)**2 + (sem.y/mean.y)**2)
  return(sem.ratio)
}

# based on normal power law
getPlContact2 <- function(dist, x0 = 5000, GAMMA = 0.992786448095086) {
  return(exp(-1*GAMMA * log(pmax(dist, x0))) / exp(-1*GAMMA * log(x0)))
}

# normalized power law for 5kb resolution Intact Hi-C, K562, use this!!
getPlContact3 <- function(dist, x0=5000, scale=5.522049023040865, gamma = 1.0044322152129266) {
  #return(exp(-1*GAMMA * log(pmax(dist, x0))) / exp(-1*GAMMA * log(x0)))
  return(exp(scale + -1*gamma * log(pmax(dist, x0))) / exp(scale + -1*gamma * log(x0)))
}

#SCALE <- -4.80 + 11.63 * GAMMA
#combined$pl.contact <- with(combined,  exp(-1*GAMMA * log(pmax(enh.krab.dist/500, 1)))) / exp(-1*GAMMA * log(pmax(0/500, 1)))
#combined$pl.contact <- with(combined,  exp(-1*GAMMA * log(pmax(enh.krab.dist, x0)))) / exp(-1*GAMMA * log(x0))
combined$pl.contact = with(combined, getPlContact3(enh.krab.dist))

dist = ggplot(data=combined, aes(x=enh.krab.dist, y=pl.contact)) +
  geom_point(color=muted('blue'), size=0.5) + 
  geom_point(aes(x=combined$enh.krab.dist, y=getPlContact2(combined$distance)),  color=muted('red'), size=0.5) +
  xlim(c(0, 1e5)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
#ggsave(file.path(out, "pl.contact.new.norm.pdf"), dist, width = 4, height = 3)

# pdf(file.path(out, "pl.contact.pdf"))
# with(subset(combined, enh.krab.dist < 1e5), plot(log10(enh.krab.dist), pl.contact, xlab = 'Distance (log10)', ylab = 'Estimated Contact', cex.lab = 1.3, cex.axis = 1.3))
# dev.off()

combined$delta.K27ac.at.pert.site <- with(combined, H3K27ac.RPM.ctrl.at.PerturbationSite - H3K27ac.RPM.KO.at.PerturbationSite)

combined.to.plot = combined


predict.K27ac <- function(x) {
  ##############
  # Find best powerlaw parameters
  ##############
  pred.func <- function(x) {
    #pl.contact <- with(combined,  exp(-1*x[[2]] * log(pmax(enh.krab.dist, x[[1]])))) / exp(-1*x[[2]] * log(x[[1]]))
    #k27ac.pred <- x[[3]] + pmax(0, combined$H3K27ac.RPM.ctrl - x[[4]]*pl.contact*(combined$H3K27ac.RPM.ctrl.at.PerturbationSite - combined$H3K27ac.RPM.KO.at.PerturbationSite))
    
    ## standard power law
    scale = 5.522049023040865
    gamma = 1.0044322152129266
    x0 = 5000
    pl.contact = with(combined, exp(scale + -1*gamma * log(pmax(x0, enh.krab.dist))) / exp(scale + -1*gamma * log(x0)))
    k27ac.pred <- x[[1]] + pmax(0, combined$H3K27ac.RPM.ctrl - x[[2]]*pl.contact*(combined$H3K27ac.RPM.ctrl.at.PerturbationSite - combined$H3K27ac.RPM.KO.at.PerturbationSite))
    
    # max delta is H3K27ac.RPM.ctrl because you cannot have negative acetylation!
    delta.k27ac.pred <- x[[1]] + pmin(combined$H3K27ac.RPM.ctrl, x[[2]]*pl.contact*(combined$H3K27ac.RPM.ctrl.at.PerturbationSite - combined$H3K27ac.RPM.KO.at.PerturbationSite))
    #return(k27ac.pred)
    return(delta.k27ac.pred)
  }
  
  eval.func <- function(x) {
    idx <- with(combined, which(H3K27ac.RPM.ctrl > 5 & H3K27ac.RPM.ctrl < 100 & enh.krab.dist < 1e5 & is.na(enhancer.name)))
    #sqrt(sum(abs((log2(as.vector(combined$H3K27ac.KO.RPM) + 1) - log2(pred.func(x) + 1)))[idx], na.rm = T)) # predicting K27ac RPM
    
    actual.delta = as.vector(combined$H3K27ac.RPM.ctrl - combined$H3K27ac.KO.RPM)
    cost = sqrt(sum(abs((log2(actual.delta + 1) - log2(pred.func(x) + 1)))[idx], na.rm = T)) # predicting delta K27ac RPM
    message(cost)
    return(cost)
  }
  

  init.values <- c(0, 1)
  local_opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                local_opts=local_opts)
  # nonlinear optimization function
  res <- nloptr( x0=init.values,
                 eval_f=eval.func,
                 opts=opts,
                 lb = c(0, 0),
                 ub = c(0, 10))
  res.df <- as.data.frame(t(res$solution))
  colnames(res.df) <- c("a", "b")
  write.table(res.df, file.path(out, "optimal.params.new.norm.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  
  a <- res$solution[[1]]
  b <- res$solution[[2]]
  message("a: ", a)
  message("b: ", b)

  combined.to.plot$H3K27ac.KO.RPM.prediction <- a + with(combined, pmax(0, H3K27ac.RPM.ctrl -  b*pl.contact*(H3K27ac.RPM.ctrl.at.PerturbationSite - H3K27ac.RPM.KO.at.PerturbationSite)))
  combined.to.plot$H3K27ac.KO.RPM.residual <- with(combined.to.plot, H3K27ac.KO.RPM.prediction -  H3K27ac.KO.RPM)
  write.table(combined, file=file.path(out, 'combined.to.plot.predictions.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
  
  this.cols <- c("chr","start","end","enhancer.name","PerturbationName","enh.krab.dist","pl.contact","H3K27ac.RPM.ctrl","H3K27ac.KO.RPM","H3K27ac.KO.RPM.prediction","H3K27ac.KO.RPM.residual", "Reference")
  
  # check if this works...
  # test without filtering distance<1e5
  #combined.to.plot <- subset(combined, enh.krab.dist < 1e5 & abs(H3K27ac.KO.RPM.residual) > 10 & is.na(enhancer.name)) #[, ..this.cols]

  # filter data to analyze
  # control RPM > 5
  # change in K27ac at target > 5
  # distance > 1kb
  
  combined.to.plot = dplyr::filter(combined.to.plot, abs(H3K27ac.RPM.ctrl)>5,
                                   abs(delta.K27ac.at.pert.site)>5,
                                   enh.krab.dist>1000)
  
  # plot k27 predicted vs actual
  # delta
  k27 = ggplot(combined.to.plot, aes(x = H3K27ac.RPM.ctrl - H3K27ac.KO.RPM,
                                     y = pmin(H3K27ac.RPM.ctrl, b * pl.contact * (H3K27ac.RPM.ctrl.at.PerturbationSite - H3K27ac.RPM.KO.at.PerturbationSite)))) +
    geom_point(color=muted("blue")) +
    geom_abline(aes(slope = 1, intercept = 0), lty = 2, col = 'black') + 
    stat_regline_equation(label.y=0.1, label.x=0.3, aes(label=..rr.label..)) +
    labs(x = 'Delta K27ac - Enhancer (measured)', y = 'Delta K27ac - Enhancer (predicted)') +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
  ggsave(file.path(out, "K27ac.delta.prediction.new.norm.pdf"), k27, width = 4, height = 3)
  
  # actual RPM
  k27.2 = ggplot(combined.to.plot, aes(x = H3K27ac.KO.RPM,
                                       y = H3K27ac.KO.RPM.prediction)) +
    geom_point(color=muted("blue"), size=0.7, alpha=0.7) +
    geom_abline(aes(slope = 1, intercept = 0), lty = 2, col = 'black') + 
    stat_regline_equation(label.y=60, label.x=10, aes(label=..rr.label..)) +
    labs(x = 'H3K27ac Ci condition\n(measured at enhancer)', y = 'H3K27ac Ci condition\n(predicted at enhancer)') +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
  ggsave(file.path(out, "K27ac.linear.model.new.norm.pdf"), k27.2, width = 4, height = 3)
}




K27ac.by.distance <- function() {
  ## Define distance bins for enh.krab.dist (1-10kb, 10kb-100kb, >100kb)
  # n for each...
  
  combined.to.plot = dplyr::filter(combined.to.plot, abs(H3K27ac.RPM.ctrl)>5,
                                   abs(delta.K27ac.at.pert.site)>5,
                                   enh.krab.dist>1000)
  print(nrow(combined.to.plot))
  combined.to.plot$enh.krab.dist = as.numeric(combined.to.plot$enh.krab.dist)
  n.d1 = nrow(dplyr::filter(combined.to.plot, enh.krab.dist>=1000, enh.krab.dist<10000))
  n.d2 = nrow(dplyr::filter(combined.to.plot, enh.krab.dist>=10000, enh.krab.dist<100000))
  n.d3 = nrow(dplyr::filter(combined.to.plot, enh.krab.dist>=100000))
  message(n.d1)
  message(n.d2)
  message(n.d3)
  
  combined.to.plot$distanceBin = 'NA'
  combined.to.plot$distanceBin[combined.to.plot$enh.krab.dist<10e3&combined.to.plot$enh.krab.dist>=1e3] = paste0("1-10kb\nN=", n.d1)
  combined.to.plot$distanceBin[combined.to.plot$enh.krab.dist>=10e3&combined.to.plot$enh.krab.dist<100e3] = paste0("10-100kb\nN=", n.d2)
  combined.to.plot$distanceBin[combined.to.plot$enh.krab.dist>=100e3] = paste0(">100kb\nN=", n.d3)
  
  # order factor
  combined.to.plot$distanceBin = factor(combined.to.plot$distanceBin)
  combined.to.plot$distanceBin = ordered(combined.to.plot$distanceBin, levels = c(paste0("1-10kb\nN=", n.d1), 
                                                                                  paste0("10-100kb\nN=", n.d2), 
                                                                                  paste0(">100kb\nN=", n.d3)))
  combined.to.plot$delta.K27ac.at.enhancer = -1*(combined.to.plot$H3K27ac.RPM.ctrl - combined.to.plot$H3K27ac.KO.RPM)
  combined.to.plot$FC.K27ac.at.enhancer = log10(combined.to.plot$H3K27ac.RPM.ctrl/combined.to.plot$H3K27ac.KO.RPM)
  
  write.table(dplyr::select(combined.to.plot, -distanceBin), file=file.path(out, 'combinedToPlotData.LTRfilt.tsv'), col.names=T, row.names=F, quote=F, sep='\t')
  
  fig = ggplot(combined.to.plot, aes(x = distanceBin,
                                     y = -1*(H3K27ac.RPM.ctrl - H3K27ac.KO.RPM))) +
    geom_boxplot(width=0.5) +
    geom_point() +
    geom_hline(yintercept=0, lty='solid') +
    geom_hline(yintercept=0, lty='dashed', size=0.5, color='red') +
    labs(x="Distance between\nenhancer and target site", y="Change in H3K27ac RPM (Ci - WT) at enhancer") + 
    ylim(c(-30,10.1)) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
  
  fig2 = ggplot(combined.to.plot, aes(x = distanceBin,
                                      y = (log2(H3K27ac.KO.RPM/H3K27ac.RPM.ctrl)))) +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    geom_point(position=position_jitter(width=0.1), size=0.5, alpha=0.5) +
    geom_hline(yintercept=0, lty='solid') +
    geom_hline(yintercept=0, lty='dashed', size=0.5, color='red') +
    labs(x="Distance from\nCRISPRi element to nearby element", y="H3K27ac at nearby enhancer\nlog2 (CRISPRi / WT)") + 
    #ylim(c(-30,10.1)) +
    theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))
  
  ggsave(file.path(out, 'main_fig4d_fc_ltrfilt.pdf'), fig2, height=3, width=3)
  
  # Calculate significance values of distributions vs each other
  ## Wilcoxon rank sum test
  g1 = filter(combined.to.plot, distanceBin== paste0("1-10kb\nN=", n.d1))
  g2 = filter(combined.to.plot, distanceBin==paste0("10-100kb\nN=", n.d2))
  g3 = filter(combined.to.plot, distanceBin==paste0(">100kb\nN=", n.d3))
  
  g1.v.g3 = wilcox.test(g1$delta.K27ac.at.enhancer, g3$delta.K27ac.at.enhancer)
  g1.v.g2 = wilcox.test(g1$delta.K27ac.at.enhancer, g2$delta.K27ac.at.enhancer)
  g2.v.g3 = wilcox.test(g2$delta.K27ac.at.enhancer, g3$delta.K27ac.at.enhancer)
  
  print("Wilcoxon\n1 vs 3")
  print(g1.v.g3)
  
  print("1 vs 2")
  print(g1.v.g2)
  
  print("2 vs 3")
  print(g2.v.g3)
  
  ## Plotting
  this.breaks <- c(seq(100,900,100), seq(1000, 9000, 1000), seq(1e4, 1e5, 1e4))
  
  x=ggplot(combined.to.plot,
           aes(x = enh.krab.dist,
               y = -1*(H3K27ac.RPM.ctrl - H3K27ac.KO.RPM),
               color = -1*delta.K27ac.at.pert.site,
               shape = Reference)) + 
    geom_point(size = 3) +
    scale_x_log10(breaks = this.breaks,
                  minor_breaks = NULL,
                  labels = ifelse( log10(this.breaks) == floor(log10(this.breaks)), this.breaks, "")) + 
    coord_cartesian(ylim = c(-30,10)) + 
    scale_color_viridis(option = "plasma", direction = 1) +
    labs(x = "Enhancer - TargetSite Distance (bp)", y = "Change H3K27ac (Ci - WT) - Enhancer", color = "Change H3K27ac (Ci - WT) - TargetSite") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  ggsave(file.path(out, "K27ac.full.plot.ltrFilt.pdf"), x, width = 12, height = 6)
  
  tmp <- subset(combined, H3K27ac.RPM.ctrl > 5 & enh.krab.dist < 1e5 & (H3K27ac.RPM.ctrl.at.PerturbationSite - H3K27ac.RPM.KO.at.PerturbationSite) > 5)
  ggplot(tmp,
         aes(x = log10(enh.krab.dist),
             y = (H3K27ac.RPM.ctrl - H3K27ac.KO.RPM) / (H3K27ac.RPM.ctrl.at.PerturbationSite - H3K27ac.RPM.KO.at.PerturbationSite),
             color = Reference)) + #Reference
    geom_point(size = 2) + 
    geom_line(data = tmp,
              col = 'black',
              aes(x = log10(enh.krab.dist),
                  y = pl.contact)) + 
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) + 
    labs(x = "Enhancer - TargetSite Distance log10(bp)", y = "(Change H3K27ac Enhancer) / (Change H3K27ac TargetSite)")
  
  ggplot(subset(tmp, Reference == "Fulco2016"),
         aes(x = log10(enh.krab.dist),
             y = (H3K27ac.RPM.ctrl - H3K27ac.KO.RPM) / (H3K27ac.RPM.ctrl.at.PerturbationSite - H3K27ac.RPM.KO.at.PerturbationSite),
             color = PerturbationName)) + #Reference
    geom_point(size = 2) + 
    geom_line(data = tmp,
              col = 'black',
              aes(x = log10(enh.krab.dist),
                  y = pl.contact)) + 
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) + 
    labs(x = "Enhancer - TargetSite Distance log10(bp)", y = "(Change H3K27ac Enhancer) / (Change H3K27ac TargetSite)")
}

K27ac.by.distance()
#predict.K27ac()
#warnings()

