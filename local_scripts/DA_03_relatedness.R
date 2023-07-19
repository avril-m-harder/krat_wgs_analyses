setwd('/Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/')

library(vcfR)
library(scales)
library(ghibli)

##### Read in data #####
## contig information
contigs <- read.table('../contigs_easley.txt', header = TRUE)
contigs <- contigs[contigs$length >= 1e6,]

## f(ROH)verlap results for non-LD-pruned data
frohverlap <- read.csv('../corrected_plink_overlap.csv')
frohverlap <- frohverlap[frohverlap$c.index %in% contigs$c.index,]  ## added on Dec. 2, 2022 to limit to contigs >= 1 Mb in length
frohverlap$pair.id <- paste0(frohverlap$id1,'-',frohverlap$id2)
OUT <- NULL
for(p in unique(frohverlap$pair.id)){
  sub <- frohverlap[frohverlap$pair.id == p,]
  save <- c(sub$id1[1], sub$id2[1], sum(sub$len.overlap))
  OUT <- rbind(OUT, save)
}
frohverlap <- as.data.frame(OUT)
colnames(frohverlap) <- c('id1','id2','len.overlap')
tot.len <- max(contigs$cumsum) ## length of filtered contigs
frohverlap$frohverlap <- frohverlap$len.overlap/tot.len

## NGSrelate2 relatedness results
ngs.LD <- read.table('ngsrelate_krat_allsamps_LDpruned.res', header = TRUE)
ngs <- read.table('../ngsrelate_krat_allsamps.res', header = TRUE)
samp.order <- read.table('../vcf_sample_order.txt')
samp.order[,2] <- c(0:47)
## replace indices with sample names in NGS output
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'a', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id1'
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'b', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id2'
ngs.LD <- ngs.LD[order(ngs.LD$a, ngs.LD$b),]

ngs <- merge(ngs, samp.order, by.x = 'a', by.y = 'V2')
colnames(ngs)[ncol(ngs)] <- 'id1'
ngs <- merge(ngs, samp.order, by.x = 'b', by.y = 'V2')
colnames(ngs)[ncol(ngs)] <- 'id2'
ngs <- ngs[order(ngs$a, ngs$b),]


## read in genotypes to calculate Rxy in related package -- ope, 'related' is a pain for this version of R on a Mac
# gts <- read.table('related_formatted_genotypes.txt', sep = ' ', header = FALSE)

## read in relationship information
relates <- read.csv('../../data/pairwise_samp_relationships.csv')

## read in pedigree data
ped.dat <- read.csv('../../data/summary_allstats2_withgen_zeroes.csv')


##### Plot R1:KING #####
## KING cutoffs
king.un <- 1/(2^(9/2))
king.3 <- 1/(2^(7/2))
king.2 <- 1/(2^(5/2))
king.1 <- 1/(2^(3/2))

k <- 1
while(k == 1){
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/KING_R1_pruned_vs_nonpruned.pdf', width = 10, height = 5)
par(mfrow = c(1,2))
plot(ngs$R1, ngs$KING, col = 'transparent', xlim = c(0.1, 0.8), ylim = c(-0.15, 0.3), 
     xlab = 'R1', ylab = 'KING', main = 'No LD pruning')
  abline(h = king.1, lty = 2)
  abline(h = king.2, lty = 2)
  abline(h = king.3, lty = 2)
  abline(h = king.un, lty = 2)
  text(0.07, king.un-((king.un+king.3)/2), labels = 'Unrelated', pos=4)
  text(0.07, (king.un+king.3)/2, labels = '3rd', pos=4)
  text(0.07, (king.3+king.2)/2, labels = '2nd', pos=4)
  text(0.07, (king.2+king.1)/2, labels = '1st', pos=4)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs[ngs$id1 %in% mates[r, c('id1','id2')] & ngs$id2 %in% mates[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs[ngs$id1 %in% half.sibs[r, c('id1','id2')] & ngs$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs[ngs$id1 %in% full.sibs[r, c('id1','id2')] & ngs$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){
    sub <- ngs[ngs$id1 %in% par.off[r, c('id1','id2')] & ngs$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))

  
plot(ngs.LD$R1, ngs.LD$KING, col = 'transparent', xlim = c(0.1, 0.8), ylim = c(-0.15, 0.3), 
       xlab = 'R1', ylab = 'KING', main = 'LD pruning')
  abline(h = king.1, lty = 2)
  abline(h = king.2, lty = 2)
  abline(h = king.3, lty = 2)
  abline(h = king.un, lty = 2)
  text(0.07, king.un-((king.un+king.3)/2), labels = 'Unrelated', pos=4)
  text(0.07, (king.un+king.3)/2, labels = '3rd', pos=4)
  text(0.07, (king.3+king.2)/2, labels = '2nd', pos=4)
  text(0.07, (king.2+king.1)/2, labels = '1st', pos=4)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){
    sub <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
dev.off()
k <- k+1
}


##### Compare pruned/non-pruned KING and Rab stats #####
k <- 1
while(k == 1){
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/pruned_vs_nonpruned_KING_Rab.pdf', width = 6, height = 6)
## KING
hist(ngs$KING, breaks = 50, col = 'dodgerblue3', main = 'KING', xlab = 'KING')
  hist(ngs.LD$KING, add = TRUE, breaks = 50, col = alpha('grey', 0.4))
  legend('topright', legend = c('No pruning', 'Pruning'), fill = c('dodgerblue3',alpha('grey', 0.4)))
  
plot(0, 0, col = 'transparent', xlim = c(min(ngs$KING), max(ngs$KING)), ylim = c(min(ngs.LD$KING), max(ngs.LD$KING)),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'KING')
  abline(0, 1, lty = 2)
  for(r in 1:nrow(ngs)){
    points(ngs$KING[r], ngs.LD$KING[r], col = alpha('darkorchid4', 0.4), pch = 19, cex = 0.4)
  }
plot(0, 0, col = 'transparent', xlim = c(-0.15, 0), ylim = c(-0.1, 0.1),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'KING - Zoomed')
  abline(0, 1, lty = 2)
  for(r in 1:nrow(ngs)){
    points(ngs$KING[r], ngs.LD$KING[r], col = alpha('darkorchid4', 0.4), pch = 19, cex = 0.4)
  }
  
OUT <- NULL
plot(0, 0, col = 'transparent', xlim = c(-0.2, 0.3), ylim = c(-0.1, 0.3),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'KING')
  abline(0, 1, lty = 2)
  lines(c(king.2, king.2), c(-1, king.2), lty = 3)
  lines(c(-1, king.2), c(king.2, king.2), lty = 3)
  lines(c(king.3, king.3), c(-1, king.3), lty = 3)
  lines(c(-1, king.3), c(king.3, king.3), lty = 3)
  lines(c(king.un, king.un), c(-1, king.un), lty = 3)
  lines(c(-1, king.un), c(king.un, king.un), lty = 3)
  text(-0.22, king.un-((king.un+king.3)/2), labels = 'Unrelated', pos=4)
  text(-0.22, (king.un+king.3)/2, labels = '3rd', pos=4)
  text(-0.22, (king.3+king.2)/2, labels = '2nd', pos=4)
  text(-0.22, (king.2+king.1)/2, labels = '1st', pos=4)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs[ngs$id1 %in% mates[r, c('id1','id2')] & ngs$id2 %in% mates[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
    points(sub$KING, sub1$KING, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs[ngs$id1 %in% half.sibs[r, c('id1','id2')] & ngs$id2 %in% half.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$KING, sub1$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
    if(sub1$KING < king.un | sub$KING > king.2){
      sub1 <- cbind(half.sibs[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs[ngs$id1 %in% full.sibs[r, c('id1','id2')] & ngs$id2 %in% full.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$KING, sub1$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
    if(sub1$KING < king.2){
      sub1 <- cbind(full.sibs[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){
    sub <- ngs[ngs$id1 %in% par.off[r, c('id1','id2')] & ngs$id2 %in% par.off[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$KING, sub1$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
    if(sub1$KING < king.2){
      sub1 <- cbind(par.off[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]), bg = 'white')
  king.mismatches <- as.data.frame(OUT)
  
## Rab
hist(ngs$rab, xlim = c(0, 0.1), breaks = 100, col = 'dodgerblue3', main = 'Rab - Zoomed', xlab = 'Rab')
  hist(ngs.LD$rab, add = TRUE, breaks = 100, col = alpha('grey', 0.4))
  legend('topright', legend = c('No pruning', 'Pruning'), fill = c('dodgerblue3',alpha('grey', 0.4)))
  
plot(0, 0, col = 'transparent', xlim = c(min(ngs$rab), max(ngs$rab)), ylim = c(min(ngs.LD$rab), max(ngs.LD$rab)),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'Rab')
  abline(0, 1, lty = 2)
  for(r in 1:nrow(ngs)){
    points(ngs$rab[r], ngs.LD$rab[r], col = alpha('darkorchid4', 0.4), pch = 19, cex = 0.4)
  }
plot(0, 0, col = 'transparent', xlim = c(0, 0.1), ylim = c(0, 0.1),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'Rab - Zoomed')
  abline(0, 1, lty = 2)
  for(r in 1:nrow(ngs)){
    points(ngs$rab[r], ngs.LD$rab[r], col = alpha('darkorchid4', 0.4), pch = 19, cex = 0.4)
  }
  
plot(0, 0, col = 'transparent', xlim = c(min(ngs$rab), max(ngs$rab)), ylim = c(min(ngs.LD$rab), max(ngs.LD$rab)),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'Rab')
  abline(0, 1, lty = 2)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs[ngs$id1 %in% mates[r, c('id1','id2')] & ngs$id2 %in% mates[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs[ngs$id1 %in% half.sibs[r, c('id1','id2')] & ngs$id2 %in% half.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs[ngs$id1 %in% full.sibs[r, c('id1','id2')] & ngs$id2 %in% full.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){
    sub <- ngs[ngs$id1 %in% par.off[r, c('id1','id2')] & ngs$id2 %in% par.off[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))

OUT <- NULL
plot(0, 0, col = 'transparent', xlim = c(0, 0.15), ylim = c(0, 0.15),
     xlab = 'Non-pruned', ylab = 'Pruned', main = 'Rab - Zoomed')
  abline(0, 1, lty = 2)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs[ngs$id1 %in% mates[r, c('id1','id2')] & ngs$id2 %in% mates[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs[ngs$id1 %in% half.sibs[r, c('id1','id2')] & ngs$id2 %in% half.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
    if(sub1$rab < 0.05){
      sub1 <- cbind(half.sibs[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs[ngs$id1 %in% full.sibs[r, c('id1','id2')] & ngs$id2 %in% full.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
    if(sub1$rab < 0.15){
      sub1 <- cbind(full.sibs[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){
    sub <- ngs[ngs$id1 %in% par.off[r, c('id1','id2')] & ngs$id2 %in% par.off[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
    if(sub1$rab < 0.15){
      sub1 <- cbind(par.off[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
  rab.mismatches <- as.data.frame(OUT)
  
dev.off()
k <- k+1
}


##### Compare all relatedness calculations (KING, Rab, f(ROH)verlap) #####
k <- 1
while(k == 1){
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/pruned_KING_Rab_vs_fROHverlap.pdf', width = 6, height = 6)
## use all pairs with f(ROH)verlap information
## KING vs. f(ROH)verlap
plot(0, 0, col = 'transparent', xlim = c(min(ngs.LD$KING), max(ngs.LD$KING)), 
     ylim = c(min(frohverlap$frohverlap), max(frohverlap$frohverlap)), xlab = 'KING', ylab = 'f(ROH)verlap')
  for(r in 1:nrow(frohverlap)){
    k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                   (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
    points(k, frohverlap$frohverlap[r], pch = 19, col = alpha('lightgrey', 0.7), cex = 0.7)
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                        relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'mates'] == 1){
        k <- ngs.LD[which(ngs.LD$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                            ngs.LD$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 19, col = 'darkgrey')
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'half.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[4])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'full.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[5])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'parent.off'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[6])
      }
    }
  }
  legend('topleft', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'), col=c('darkgrey', ghibli_palettes$YesterdayMedium[4], ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))

## Rab vs. f(ROH)verlap
plot(0, 0, col = 'transparent', xlim = c(min(ngs.LD$rab), max(ngs.LD$rab)), 
     ylim = c(min(frohverlap$frohverlap), max(frohverlap$frohverlap)), xlab = 'Rab', ylab = 'f(ROH)verlap')
  for(r in 1:nrow(frohverlap)){
    k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                        (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
    points(k, frohverlap$frohverlap[r], pch = 19, col = alpha('lightgrey', 0.7), cex = 0.7)
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                       relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'mates'] == 1){
        k <- ngs.LD[which(ngs.LD$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                            ngs.LD$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 19, col = 'darkgrey')
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'half.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[4])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'full.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[5])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'parent.off'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[6])
      }
    }
  }
  legend('topleft', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'), col=c('darkgrey', ghibli_palettes$YesterdayMedium[4], ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
  dev.off()
  
## highlight all pairs in KING mismatches in each type of plot
# pdf('../../krat_genetics_scripts/figures_output/LD_pruned/KING_vs_fROHverlap_mismatches.pdf', width = 6, height = 6)
# for(m in 1:nrow(king.mismatches)){
#   plot(0, 0, col = 'transparent', xlim = c(-0.08, max(ngs.LD$KING)), 
#        ylim = c(min(frohverlap$frohverlap), max(frohverlap$frohverlap)), xlab = 'KING', ylab = 'f(ROH)verlap',
#        main = paste0(king.mismatches$id1[m],' - ',king.mismatches$id2[m]))
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                        relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'mates'] == 1){
#         k <- ngs.LD[which(ngs.LD$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                             ngs.LD$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'KING']
#         points(k, frohverlap$frohverlap[r], pch = 19, col = 'darkgrey')
#       }
#     }
#   }
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
#                        (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'half.sibs'] == 1){
#         k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
#                             (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
#         points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[4])
#       }
#     }
#   }
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
#                        (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'full.sibs'] == 1){
#         k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
#                             (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
#         points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[5])
#       }
#     }
#   }
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
#                        (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'parent.off'] == 1){
#         k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
#                             (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
#         points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[6])
#       }
#     }
#   }
#   legend('topleft', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'), col=c('darkgrey', ghibli_palettes$YesterdayMedium[4], ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
#   
#   ## get coordinates for mismatched pair
#   k <- ngs.LD[which((ngs.LD$id1 == king.mismatches$id1[m] & ngs.LD$id2 == king.mismatches$id2[m]) |
#                       (ngs.LD$id1 == king.mismatches$id2[m] & ngs.LD$id2 == king.mismatches$id1[m])), 'KING']
#   o <- frohverlap[which((frohverlap$id1 == king.mismatches$id1[m] & frohverlap$id2 == king.mismatches$id2[m]) |
#                       (frohverlap$id1 == king.mismatches$id2[m] & frohverlap$id2 == king.mismatches$id1[m])),
#               'frohverlap']
#   points(k, o, pch = 21, col = 'red', cex = 1.7, lwd = 2)
# }
# dev.off()
# 
# pdf('../../krat_genetics_scripts/figures_output/LD_pruned/rab_vs_fROHverlap_mismatches.pdf', width = 6, height = 6)
# for(m in 1:nrow(king.mismatches)){
#   plot(0, 0, col = 'transparent', xlim = c(min(ngs.LD$rab), max(ngs.LD$rab)), 
#        ylim = c(min(frohverlap$frohverlap), max(frohverlap$frohverlap)), xlab = 'rab', ylab = 'f(ROH)verlap',
#        main = paste0(king.mismatches$id1[m],' - ',king.mismatches$id2[m]))
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                        relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'mates'] == 1){
#         k <- ngs.LD[which(ngs.LD$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                             ngs.LD$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'rab']
#         points(k, frohverlap$frohverlap[r], pch = 19, col = 'darkgrey')
#       }
#     }
#   }
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
#                        (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'half.sibs'] == 1){
#         k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
#                             (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
#         points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[4])
#       }
#     }
#   }
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
#                        (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'full.sibs'] == 1){
#         k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
#                             (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
#         points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[5])
#       }
#     }
#   }
#   for(r in 1:nrow(frohverlap)){
#     if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
#                           relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
#       if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
#                        (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'parent.off'] == 1){
#         k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
#                             (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
#         points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[6])
#       }
#     }
#   }
#   legend('topleft', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'), col=c('darkgrey', ghibli_palettes$YesterdayMedium[4], ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
#   
#   ## get coordinates for mismatched pair
#   k <- ngs.LD[which((ngs.LD$id1 == king.mismatches$id1[m] & ngs.LD$id2 == king.mismatches$id2[m]) |
#                       (ngs.LD$id1 == king.mismatches$id2[m] & ngs.LD$id2 == king.mismatches$id1[m])), 'rab']
#   o <- frohverlap[which((frohverlap$id1 == king.mismatches$id1[m] & frohverlap$id2 == king.mismatches$id2[m]) |
#                           (frohverlap$id1 == king.mismatches$id2[m] & frohverlap$id2 == king.mismatches$id1[m])),
#                   'frohverlap']
#   points(k, o, pch = 21, col = 'red', cex = 1.7, lwd = 2)
# }
# dev.off()

k <- k+1
}

##### Get IDs for pairs with apparent genetic : pedigree relationship mismatches ###dddd##
## Rab mismatches: full-sibs and P-O pairs with Rab < 0.15 & half-sibs with Rab < 0.05, just based on visual inspection
head(rab.mismatches)
head(king.mismatches)

## for each mismatch, plot separately and highlight to show where it is exactly
## Rab
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/Rab_mismatches.pdf', width = 6, height = 6)
for(n in 1:nrow(rab.mismatches)){
  plot(0, 0, col = 'transparent', xlim = c(min(ngs$rab), max(ngs$rab)), ylim = c(min(ngs.LD$rab), max(ngs.LD$rab)),
       xlab = 'Non-pruned', ylab = 'Pruned', main = paste0(rab.mismatches$id1[n],' - ',rab.mismatches$id2[n]))
    abline(0, 1, lty = 2)
    mates <- relates[relates$mate == 1,]
    for(r in 1:nrow(mates)){
      sub <- ngs[ngs$id1 %in% mates[r, c('id1','id2')] & ngs$id2 %in% mates[r, c('id1','id2')],]
      sub1 <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
      points(sub$rab, sub1$rab, pch = 19, col = 'grey', cex = 0.7)
    }
    half.sibs <- relates[relates$half.sibs == 1,]
    for(r in 1:nrow(half.sibs)){
      sub <- ngs[ngs$id1 %in% half.sibs[r, c('id1','id2')] & ngs$id2 %in% half.sibs[r, c('id1','id2')],]
      sub1 <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
      points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
      if(sub1$rab < 0.05){
        sub1 <- cbind(half.sibs[r,], sub1)
        OUT <- rbind(OUT, sub1)
      }
    }
    full.sibs <- relates[relates$full.sibs == 1,]
    for(r in 1:nrow(full.sibs)){
      sub <- ngs[ngs$id1 %in% full.sibs[r, c('id1','id2')] & ngs$id2 %in% full.sibs[r, c('id1','id2')],]
      sub1 <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
      points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
      if(sub1$rab < 0.15){
        sub1 <- cbind(full.sibs[r,], sub1)
        OUT <- rbind(OUT, sub1)
      }
    }
    par.off <- relates[relates$parent.off == 1,]
    for(r in 1:nrow(par.off)){
      sub <- ngs[ngs$id1 %in% par.off[r, c('id1','id2')] & ngs$id2 %in% par.off[r, c('id1','id2')],]
      sub1 <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
      points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
      if(sub1$rab < 0.15){
        sub1 <- cbind(par.off[r,], sub1)
        OUT <- rbind(OUT, sub1)
      }
    }
    legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
           col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
                 ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
    sub <- ngs[ngs$id1 %in% rab.mismatches[n, c('id1','id2')] & ngs$id2 %in% rab.mismatches[n, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% rab.mismatches[n, c('id1','id2')] & ngs.LD$id2 %in% rab.mismatches[n, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 1, col = 'red', cex = 2, lwd = 2)
}
dev.off()

## KING
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/KING_mismatches.pdf', width = 6, height = 6)
for(n in 1:nrow(king.mismatches)){
  plot(ngs.LD$R1, ngs.LD$KING, col = 'transparent', xlim = c(0.1, 0.8), ylim = c(-0.07, 0.28), 
       xlab = 'R1', ylab = 'KING', main = paste0(king.mismatches$id1[n],' - ',king.mismatches$id2[n],'\n',
                                                 colnames(king.mismatches[n,c(3:6)])[which(king.mismatches[n,c(3:6)] == 1)]))
    abline(h = king.1, lty = 2)
    abline(h = king.2, lty = 2)
    abline(h = king.3, lty = 2)
    abline(h = king.un, lty = 2)
    text(0.07, king.un-((king.un+king.3)/2), labels = 'Unrelated', pos=4)
    text(0.07, (king.un+king.3)/2, labels = '3rd', pos=4)
    text(0.07, (king.3+king.2)/2, labels = '2nd', pos=4)
    text(0.07, (king.2+king.1)/2, labels = '1st', pos=4)
    mates <- relates[relates$mate == 1,]
    for(r in 1:nrow(mates)){
      sub <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
      points(sub$R1, sub$KING, pch = 19, col = 'grey', cex = 0.7)
    }
    half.sibs <- relates[relates$half.sibs == 1,]
    for(r in 1:nrow(half.sibs)){
      sub <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
      points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
    }
    full.sibs <- relates[relates$full.sibs == 1,]
    for(r in 1:nrow(full.sibs)){
      sub <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
      points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
    }
    par.off <- relates[relates$parent.off == 1,]
    for(r in 1:nrow(par.off)){ 
      sub <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
      points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
    }
    legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
           col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
                 ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
    sub1 <- ngs.LD[ngs.LD$id1 %in% king.mismatches[n, c('id1','id2')] & ngs.LD$id2 %in% king.mismatches[n, c('id1','id2')],]
    points(sub1$R1, sub1$KING, pch = 1, col = 'red', cex = 2, lwd = 2)
}
dev.off()


##### Manual changes to pedigree-defined relationships based on genetic information #####
## ** see genetic_pedigree_mismatches.pptx and genetic_pedigree_mismatches_NOTES.txt for details **
## low probability of both parents correct (0.5525 for 4915)
relates[which((relates$id1 == 4915 & relates$id2 == 4795) |
              (relates$id2 == 4915 & relates$id1 == 4795)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4915 & relates$id2 == 4796) |
              (relates$id2 == 4915 & relates$id1 == 4796)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4915 & relates$id2 == 4952) |
              (relates$id2 == 4915 & relates$id1 == 4952)), 'half.sibs'] <- 0

## all 3 genetic estimates are low, even though pedigree probz aren't awful
relates[which((relates$id1 == 4952 & relates$id2 == 4795) |
              (relates$id2 == 4952 & relates$id1 == 4795)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4952 & relates$id2 == 4796) |
              (relates$id2 == 4952 & relates$id1 == 4796)), 'half.sibs'] <- 0

## low probability for shared parent (father), low genetic estimates, f(ROH)verlap inconclusive
relates[which((relates$id1 == 4910 & relates$id2 == 4897) |
              (relates$id2 == 4910 & relates$id1 == 4897)), 'half.sibs'] <- 0

## all 3 genetic estimates are low
relates[which((relates$id1 == 4901 & relates$id2 == 4943) |
              (relates$id2 == 4901 & relates$id1 == 4943)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4901 & relates$id2 == 4962) |
              (relates$id2 == 4901 & relates$id1 == 4962)), 'half.sibs'] <- 0

relates[which((relates$id1 == 4976 & relates$id2 == 5054) |
              (relates$id2 == 4976 & relates$id1 == 5054)), 'full.sibs'] <- 0


#### Questionable reassignments because marginal posterior probz are pretty high for these two sets in the pedigree
## fROHverlap low here, too
relates[which((relates$id1 == 4962 & relates$id2 == 4943) |
              (relates$id2 == 4962 & relates$id1 == 4943)), 'full.sibs'] <- 0

## below is extra strange because their parents are also full-siblings in the pedigree? 
# relates[which((relates$id1 == 4195 & relates$id2 == 4459) |
#               (relates$id2 == 4195 & relates$id1 == 4459)), 'full.sibs'] <- 0
#### these 2 have huuuuuuuge shared ROH tracts, keeping them as full-sibs.


##### Visualize genetic estimates and pedigree relationships again after manual revisions #####
k <- 1
while(k == 1){
## Rab
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/Rab_mismatches_after_relates_corrections.pdf', width = 6, height = 6)
for(n in 1:nrow(rab.mismatches)){
  plot(0, 0, col = 'transparent', xlim = c(min(ngs$rab), max(ngs$rab)), ylim = c(min(ngs.LD$rab), max(ngs.LD$rab)),
       xlab = 'Non-pruned', ylab = 'Pruned', main = paste0(rab.mismatches$id1[n],' - ',rab.mismatches$id2[n]))
  abline(0, 1, lty = 2)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs[ngs$id1 %in% mates[r, c('id1','id2')] & ngs$id2 %in% mates[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs[ngs$id1 %in% half.sibs[r, c('id1','id2')] & ngs$id2 %in% half.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
    if(sub1$rab < 0.05){
      sub1 <- cbind(half.sibs[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs[ngs$id1 %in% full.sibs[r, c('id1','id2')] & ngs$id2 %in% full.sibs[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
    if(sub1$rab < 0.15){
      sub1 <- cbind(full.sibs[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){
    sub <- ngs[ngs$id1 %in% par.off[r, c('id1','id2')] & ngs$id2 %in% par.off[r, c('id1','id2')],]
    sub1 <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$rab, sub1$rab, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
    if(sub1$rab < 0.15){
      sub1 <- cbind(par.off[r,], sub1)
      OUT <- rbind(OUT, sub1)
    }
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
  sub <- ngs[ngs$id1 %in% rab.mismatches[n, c('id1','id2')] & ngs$id2 %in% rab.mismatches[n, c('id1','id2')],]
  sub1 <- ngs.LD[ngs.LD$id1 %in% rab.mismatches[n, c('id1','id2')] & ngs.LD$id2 %in% rab.mismatches[n, c('id1','id2')],]
  points(sub$rab, sub1$rab, pch = 1, col = 'red', cex = 2, lwd = 2)
}
dev.off()

## KING
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/KING_mismatches_after_relates_corrections.pdf', width = 6, height = 6)
for(n in 1:nrow(king.mismatches)){
  plot(ngs.LD$R1, ngs.LD$KING, col = 'transparent', xlim = c(0.1, 0.8), ylim = c(-0.07, 0.28), 
       xlab = 'R1', ylab = 'KING', main = paste0(king.mismatches$id1[n],' - ',king.mismatches$id2[n],'\n',
                                                 colnames(king.mismatches[n,c(3:6)])[which(king.mismatches[n,c(3:6)] == 1)]))
  abline(h = king.1, lty = 2)
  abline(h = king.2, lty = 2)
  abline(h = king.3, lty = 2)
  abline(h = king.un, lty = 2)
  text(0.07, king.un-((king.un+king.3)/2), labels = 'Unrelated', pos=4)
  text(0.07, (king.un+king.3)/2, labels = '3rd', pos=4)
  text(0.07, (king.3+king.2)/2, labels = '2nd', pos=4)
  text(0.07, (king.2+king.1)/2, labels = '1st', pos=4)
  mates <- relates[relates$mate == 1,]
  for(r in 1:nrow(mates)){
    sub <- ngs.LD[ngs.LD$id1 %in% mates[r, c('id1','id2')] & ngs.LD$id2 %in% mates[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 19, col = 'grey', cex = 0.7)
  }
  half.sibs <- relates[relates$half.sibs == 1,]
  for(r in 1:nrow(half.sibs)){
    sub <- ngs.LD[ngs.LD$id1 %in% half.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% half.sibs[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[4])
  }
  full.sibs <- relates[relates$full.sibs == 1,]
  for(r in 1:nrow(full.sibs)){
    sub <- ngs.LD[ngs.LD$id1 %in% full.sibs[r, c('id1','id2')] & ngs.LD$id2 %in% full.sibs[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[5])
  }
  par.off <- relates[relates$parent.off == 1,]
  for(r in 1:nrow(par.off)){ 
    sub <- ngs.LD[ngs.LD$id1 %in% par.off[r, c('id1','id2')] & ngs.LD$id2 %in% par.off[r, c('id1','id2')],]
    points(sub$R1, sub$KING, pch = 17, col = ghibli_palettes$YesterdayMedium[6])
  }
  legend('bottomright', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'),
         col=c('lightgrey', ghibli_palettes$YesterdayMedium[4], 
               ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
  sub1 <- ngs.LD[ngs.LD$id1 %in% king.mismatches[n, c('id1','id2')] & ngs.LD$id2 %in% king.mismatches[n, c('id1','id2')],]
  points(sub1$R1, sub1$KING, pch = 1, col = 'red', cex = 2, lwd = 2)
}
dev.off()

pdf('../../krat_genetics_scripts/figures_output/LD_pruned/pruned_KING_Rab_vs_fROHverlap_after_relates_corrections.pdf', width = 6, height = 6)
## use all pairs with f(ROH)verlap information
## KING vs. f(ROH)verlap
plot(0, 0, col = 'transparent', xlim = c(-0.08, max(ngs.LD$KING)), 
     ylim = c(min(frohverlap$frohverlap), max(frohverlap$frohverlap)), xlab = 'KING', ylab = 'f(ROH)verlap')
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                       relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'mates'] == 1){
        k <- ngs.LD[which(ngs.LD$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                            ngs.LD$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 19, col = 'darkgrey')
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'half.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[4])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'full.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[5])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'parent.off'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'KING']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[6])
      }
    }
  }
  legend('topleft', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'), col=c('darkgrey', ghibli_palettes$YesterdayMedium[4], ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))

## Rab vs. f(ROH)verlap
plot(0, 0, col = 'transparent', xlim = c(min(ngs.LD$rab), max(ngs.LD$rab)), 
     ylim = c(min(frohverlap$frohverlap), max(frohverlap$frohverlap)), xlab = 'Rab', ylab = 'f(ROH)verlap')
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                       relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'mates'] == 1){
        k <- ngs.LD[which(ngs.LD$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                            ngs.LD$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 19, col = 'darkgrey')
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'half.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[4])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'full.sibs'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[5])
      }
    }
  }
  for(r in 1:nrow(frohverlap)){
    if(nrow(relates[which(relates$id1 %in% c(frohverlap$id1[r], frohverlap$id2[r]) & 
                          relates$id2 %in% c(frohverlap$id1[r], frohverlap$id2[r])),]) > 0){
      if(relates[which((relates$id1 == frohverlap$id1[r] & relates$id2 == frohverlap$id2[r]) |
                       (relates$id1 == frohverlap$id2[r] & relates$id2 == frohverlap$id1[r])), 'parent.off'] == 1){
        k <- ngs.LD[which((ngs.LD$id1 == frohverlap$id1[r] & ngs.LD$id2 == frohverlap$id2[r]) |
                            (ngs.LD$id1 == frohverlap$id2[r] & ngs.LD$id2 == frohverlap$id1[r])), 'rab']
        points(k, frohverlap$frohverlap[r], pch = 17, col = ghibli_palettes$YesterdayMedium[6])
      }
    }
  }
  legend('topleft', inset=0.05, pch=c(19,17,17,17), legend=c('Mate pairs','Half-siblings','Full-siblings','Parent-offspring'), col=c('darkgrey', ghibli_palettes$YesterdayMedium[4], ghibli_palettes$YesterdayMedium[5], ghibli_palettes$YesterdayMedium[6]))
  dev.off()
k <- k+1
}

##### Compare pruned and non-pruned stats for f(ROH) results (this should probably go in DA_02.R) #####
# plot(froh$prop.miss, froh.LD$prop.miss)
# plot(froh$norm.prop.het, froh.LD$norm.prop.het)
# plot(froh$pl.froh, froh.LD$pl.froh)
#   abline(0,1)
# plot(froh$gt.froh, froh.LD$gt.froh)  
#   abline(0,1)
# plot(froh$plink.froh, froh.LD$plink.froh)  
#   abline(0, 1)
  
##### Manually examine relatedness stats for purposefully-sequenced pairs of full siblings and their mates #####
