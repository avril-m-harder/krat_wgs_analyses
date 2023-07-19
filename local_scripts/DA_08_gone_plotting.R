##### Results for the same n=12 ~unrelated individuals born in 2005, different #s of contigs used #####
setwd('/Users/Avril/Documents/krat_genetics/data/gone_results/n12_100_replicates/')
library(scales)
library(ghibli)

OUT <- NULL
dirs <- list.files()   ## list of directories (different #s contigs sampled)
for(d in dirs){
  n.con <- unlist(strsplit(d, split = '_', fixed = TRUE))[1]  ## # of contigs sampled
  fns <- list.files(d)                                        ## results files for current contig sampling level
  r <- 1
  for(f in unique(fns)){
    temp <- read.table(paste0(d,'/',f), skip = 1, header = TRUE)
    temp$rep <- r
    temp$n.contigs <- n.con
    OUT <- rbind(OUT, temp)
    r <- r+1
  }
}

### Make plots for each contig sampling level
MINS <- NULL
MAXS <- NULL
pdf('/Users/Avril/Documents/krat_genetics/gone/figures/n12_100replicates_5contiglevels_multiplot.pdf', width = 10, height = 7)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 2.1, 2.1))
for(c in sort(as.numeric(unique(OUT$n.contigs)))){
  sub <- OUT[OUT$n.contigs == c,]
  
  MEANS <- NULL
  LO <- NULL
  HI <- NULL
  
  ## summarize data
  for(g in unique(sub$Generation)){
    reps <- max(sub$rep)
    if(nrow(sub[sub$Generation == g,]) == reps){
      m <- mean(sub[sub$Generation == g, 'Geometric_mean'])
      sd <- sd(sub[sub$Generation == g, 'Geometric_mean'])
      MEANS <- c(MEANS, m)
      LO <- c(LO, m-sd)
      HI <- c(HI, m+sd)
    }
  }
  MINS <- c(MINS, min(LO))
  MAXS <- c(MAXS, max(HI))
  
  ## plot individual line for each replicate
  plot(0, 0, xlim = c(0, max(sub$Generation)), ylim = c(min(LO), max(c(HI, sub$Geometric_mean))), col = 'transparent',
       xlab = 'Generations ago', ylab = 'Ne', main = paste0(c,' contigs'))
    for(r in unique(sub$rep)){
      temp <- sub[sub$rep == r,]
      lines(temp$Generation, temp$Geometric_mean, col = alpha('black', 0.3))
    }
  
  ## individual lines, zoomed in
  plot(0, 0, xlim = c(0, 200), ylim = c(min(LO), max(c(HI, sub$Geometric_mean))), col = 'transparent',
       xlab = 'Generations ago', ylab = 'Ne', main = paste0(c,' contigs'))
  for(r in unique(sub$rep)){
    temp <- sub[sub$rep == r,]
    lines(temp$Generation, temp$Geometric_mean, col = alpha('black', 0.3))
  }
  
  ## mean with SD polygon
  plot(c(1:length(MEANS)), MEANS, xlim = c(0, max(sub$Generation)), ylim = c(min(LO), max(HI)),
       xlab = 'Generations ago', ylab = 'Ne', main = paste0(c,' contigs'), col = 'transparent')
    polygon(c(1:length(MEANS), length(MEANS):1), c(LO, HI[length(HI):1]),
            border = NA, col = alpha('springgreen4', 0.4))
    lines(c(1:length(MEANS)), MEANS, col = 'springgreen4')
    
  ## mean with SD polygon, zoomed in 
  plot(c(1:length(MEANS)), MEANS, xlim = c(0, 200), ylim = c(min(LO), max(HI)),
       xlab = 'Generations ago', ylab = 'Ne', main = paste0(c,' contigs'), col = 'transparent')
    polygon(c(1:length(MEANS), length(MEANS):1), c(LO, HI[length(HI):1]),
            border = NA, col = alpha('springgreen4', 0.4))
    lines(c(1:length(MEANS)), MEANS, col = 'springgreen4')
}
dev.off()

### one plot to compare contig levels
cols <- ghibli_palettes$PonyoMedium[c(1,3,4,6,2)]

pdf('/Users/Avril/Documents/krat_genetics/gone/figures/n12_100replicates_5contiglevels.pdf', width = 7, height = 5)
## all contig levels, zoomed out
plot(0, 0, xlim = c(0, max(OUT$Generation)), ylim = c(min(MINS), max(MAXS)), col = 'transparent',
     xlab = 'Generations ago', ylab = 'Ne', main = 'All contig levels')
col <- 1
for(c in sort(as.numeric(unique(OUT$n.contigs)), decreasing = TRUE)){
  sub <- OUT[OUT$n.contigs == c,]
  
  MEANS <- NULL
  LO <- NULL
  HI <- NULL
  
  ## summarize data
  for(g in unique(sub$Generation)){
    reps <- max(sub$rep)
    if(nrow(sub[sub$Generation == g,]) == reps){
      m <- mean(sub[sub$Generation == g, 'Geometric_mean'])
      sd <- sd(sub[sub$Generation == g, 'Geometric_mean'])
      MEANS <- c(MEANS, m)
      LO <- c(LO, m-sd)
      HI <- c(HI, m+sd)
    }
  }
  
  ## mean with SD polygon
  polygon(c(1:length(MEANS), length(MEANS):1), c(LO, HI[length(HI):1]),
          border = NA, col = alpha(col = cols[col], 0.4))
  lines(c(1:length(MEANS)), MEANS, col = cols[col])
  
  col <- col+1
}

## all contig levels, zoomed in y-axis
par(mar = c(5.1, 4.1, 4.1, 5.1), xpd = FALSE)
plot(0, 0, xlim = c(0, max(OUT$Generation)), ylim = c(min(MINS), 1e5), col = 'transparent',
     xlab = 'Generations ago', ylab = 'Ne', main = 'All contig levels')
col <- 1
for(c in sort(as.numeric(unique(OUT$n.contigs)), decreasing = TRUE)){
  sub <- OUT[OUT$n.contigs == c,]
  
  MEANS <- NULL
  LO <- NULL
  HI <- NULL
  
  ## summarize data
  for(g in unique(sub$Generation)){
    reps <- max(sub$rep)
    if(nrow(sub[sub$Generation == g,]) == reps){
      m <- mean(sub[sub$Generation == g, 'Geometric_mean'])
      sd <- sd(sub[sub$Generation == g, 'Geometric_mean'])
      MEANS <- c(MEANS, m)
      LO <- c(LO, m-sd)
      HI <- c(HI, m+sd)
    }
  }
  
  ## mean with SD polygon
  polygon(c(1:length(MEANS), length(MEANS):1), c(LO, HI[length(HI):1]),
          border = NA, col = alpha(col = cols[col], 0.4))
  lines(c(1:length(MEANS)), MEANS, col = cols[col])
  
  col <- col+1
}
par(xpd = TRUE)
legend('right', legend = c('200','100','50','25','10'), col = cols, lty = 1, inset = -0.2)

## all contig levels, zoomed in y-axis and x-axis
par(mar = c(5.1, 4.1, 4.1, 5.1), xpd = FALSE)
plot(0, 0, xlim = c(0, 100), ylim = c(min(MINS), 6e4), col = 'transparent',
     xlab = 'Generations ago', ylab = 'Ne', main = 'All contig levels')
col <- 1
for(c in sort(as.numeric(unique(OUT$n.contigs)), decreasing = TRUE)){
  sub <- OUT[OUT$n.contigs == c,]
  
  MEANS <- NULL
  LO <- NULL
  HI <- NULL
  
  ## summarize data
  for(g in unique(sub$Generation)){
    reps <- max(sub$rep)
    if(nrow(sub[sub$Generation == g,]) == reps){
      m <- mean(sub[sub$Generation == g, 'Geometric_mean'])
      sd <- sd(sub[sub$Generation == g, 'Geometric_mean'])
      MEANS <- c(MEANS, m)
      LO <- c(LO, m-sd)
      HI <- c(HI, m+sd)
    }
  }
  
  ## mean with SD polygon
  polygon(c(1:length(MEANS), length(MEANS):1), c(LO, HI[length(HI):1]),
          border = NA, col = alpha(col = cols[col], 0.4))
  lines(c(1:length(MEANS)), MEANS, col = cols[col])
  
  col <- col+1
}
par(xpd = TRUE)
legend('right', legend = c('200','100','50','25','10'), col = cols, lty = 1, inset = -0.2)

dev.off()

##### 12 samples, 100 replicates #####
setwd('/Users/Avril/Documents/krat_genetics/data/gone_results/n12_100_replicates/')
library(scales)

OUT <- NULL
for(c in 1:100){
  print(c)
  temp <- read.table(paste0('Output_Ne_gone_subset_plink_',c), skip = 1, header = TRUE)
  temp$rep <- c
  OUT <- rbind(OUT, temp)
}

## plot all raw data, just to see variation
plot(0, 0, xlim = c(1, max(OUT$Generation)), ylim = c(min(OUT$Geometric_mean), max(OUT$Geometric_mean)), col = 'transparent',
     xlab = 'Generations ago', ylab = 'Ne')
  for(r in unique(OUT$rep)){
    sub <- OUT[OUT$rep == r,]
    lines(sub$Generation, sub$Geometric_mean, col = alpha('black', 0.4))
  }
## wow what a mess

## 100 individual plots (lol)
pdf('/Users/Avril/Desktop/test.pdf', width = 25, height = 25)
par(mfrow = c(10, 10))
for(r in unique(OUT$rep)){
  sub <- OUT[OUT$rep == r,]
  plot(sub$Generation, sub$Geometric_mean, xlim = c(1, max(OUT$Generation)), 
       ylim = c(min(OUT$Geometric_mean), max(OUT$Geometric_mean)), col = 'transparent',
       xlab = 'Generations ago', ylab = 'Ne')
    lines(sub$Generation, sub$Geometric_mean, col = 'blue', lwd = 2)
}
dev.off()

## calculate means and 95% CIs
MEANS <- NULL
LO <- NULL
HI <- NULL
for(g in unique(OUT$Generation)){
  m <- mean(OUT[OUT$Generation == g, 'Geometric_mean'])
  MEANS <- c(MEANS, m)
  LO <- c(LO, quantile(OUT[OUT$Generation == g, 'Geometric_mean'], probs = c(0.025, 0.975))[1])
  HI <- c(HI, quantile(OUT[OUT$Generation == g, 'Geometric_mean'], probs = c(0.025, 0.975))[2])
}
plot(0, 0, xlim = c(1, 100), ylim = c(min(OUT$Geometric_mean), max(OUT$Geometric_mean)), col = 'transparent',
     xlab = 'Generations ago', ylab = 'Ne')
  polygon(x = c(1:length(MEANS), length(MEANS):1), y = c(LO, HI),
          border = NA, col = alpha('blue', 0.2))
  lines(c(1:length(MEANS)), MEANS, col = 'blue')
  