setwd('/Users/Avril/Documents/krat_genetics/data/')

library(scales)

### set colors
library(ghibli)
plink.col <- ghibli_palette('PonyoMedium')[3]

## contig information
contigs <- read.csv('contig_lengths.csv')
contigs <- contigs[order(-contigs$length),]
contigs$c.index <- c(1:nrow(contigs)) ## create contig index, where c.index = 1 == longest contig
tot.len <- sum(contigs$length) ## 2,571,239,112 bp ## length of filtered contigs

## 4 problematic samples
prob <- c(5075, 4901, 5018, 5060)

##### Effects of missingness (while het == 2) #####
## phwh = 10
plink.roh.10 <- read.table('plink_round6/round6_hom_files/krat_plink_roh_phwh_2_phwm_10_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh.10 <- plink.roh.10[,c(1, 4, 7, 8, 10)]
colnames(plink.roh.10) <- c('id','contig','start','end','n.snps')
plink.roh.10$length <- plink.roh.10$end - plink.roh.10$start + 1
plink.roh.10 <- plink.roh.10[plink.roh.10$length >= 100000,]

## phwm = 30
plink.roh.30 <- read.table('plink_round6/round6_hom_files/krat_plink_roh_phwh_2_phwm_30_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh.30 <- plink.roh.30[,c(1, 4, 7, 8, 10)]
colnames(plink.roh.30) <- c('id','contig','start','end','n.snps')
plink.roh.30$length <- plink.roh.30$end - plink.roh.30$start + 1
plink.roh.30 <- plink.roh.30[plink.roh.30$length >= 100000,]

## phwm = 50
plink.roh.50 <- read.table('plink_round6/round6_hom_files/krat_plink_roh_phwh_2_phwm_50_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh.50 <- plink.roh.50[,c(1, 4, 7, 8, 10)]
colnames(plink.roh.50) <- c('id','contig','start','end','n.snps')
plink.roh.50$length <- plink.roh.50$end - plink.roh.50$start + 1
plink.roh.50 <- plink.roh.50[plink.roh.50$length >= 100000,]

## set f(ROH) formatting for plots
froh.lab <- substitute(paste(italic('F')[ROH]))

##### >>> Calculate and compare individual f(ROH) values between methods #####
OUT <- NULL
for(i in unique(plink.roh.10$id)){
  save <- c(as.numeric(i), sum(plink.roh.10[plink.roh.10$id == i, 'length'])/tot.len, 
            sum(plink.roh.30[plink.roh.30$id == i, 'length'])/tot.len,
            sum(plink.roh.50[plink.roh.50$id == i, 'length'])/tot.len)
  OUT <- rbind(OUT, save)
}
indiv.froh <- as.data.frame(OUT)
colnames(indiv.froh) <- c('id','froh.10','froh.30','froh.50')

pdf('../krat_genetics_scripts/figures_output/missingness_fROH_comps.pdf', width = 6, height = 6)
plot(indiv.froh$froh.10, indiv.froh$froh.30, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'missingness = 10', ylab = 'missingness = 30')
  points(indiv.froh[indiv.froh$id %in% prob, 'froh.10'], indiv.froh[indiv.froh$id %in% prob, 'froh.30'], col = 'orange1', pch = 19)
  abline(0,1, lty = 2)

plot(indiv.froh$froh.10, indiv.froh$froh.50, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'missingness = 10', ylab = 'missingness = 50')
  points(indiv.froh[indiv.froh$id %in% prob, 'froh.10'], indiv.froh[indiv.froh$id %in% prob, 'froh.50'], col = 'orange1', pch = 19)
  abline(0,1, lty = 2)

plot(indiv.froh$froh.30, indiv.froh$froh.50, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'missingness = 30', ylab = 'missingness = 50')
  points(indiv.froh[indiv.froh$id %in% prob, 'froh.30'], indiv.froh[indiv.froh$id %in% prob, 'froh.50'], col = 'orange1', pch = 19)
  abline(0,1, lty = 2)

dev.off()  

##### >>> Calculate length-specific f(ROH) values #####
b1 <- 5e5
b2 <- 1e6
b3 <- 2e6
b4 <- 4e6

for(i in unique(indiv.froh$id)){
  sub.10 <- plink.roh.10[plink.roh.10$id == i,]
  sub.30 <- plink.roh.30[plink.roh.30$id == i,]
  sub.50 <- plink.roh.50[plink.roh.50$id == i,]
  
  ## bin1
  indiv.froh[indiv.froh$id == i, 'froh.10.bin1'] <- sum(sub.10[sub.10$length < b1, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin1'] <- sum(sub.30[sub.30$length < b1, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin1'] <- sum(sub.50[sub.50$length < b1, 'length'])/tot.len
  
  ## bin2
  indiv.froh[indiv.froh$id == i, 'froh.10.bin2'] <- sum(sub.10[sub.10$length >= b1 & sub.10$length < b2, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin2'] <- sum(sub.30[sub.30$length >= b1 & sub.30$length < b2, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin2'] <- sum(sub.50[sub.50$length >= b1 & sub.50$length < b2, 'length'])/tot.len
  
  ## bin3
  indiv.froh[indiv.froh$id == i, 'froh.10.bin3'] <- sum(sub.10[sub.10$length >= b2 & sub.10$length < b3, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin3'] <- sum(sub.30[sub.30$length >= b2 & sub.30$length < b3, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin3'] <- sum(sub.50[sub.50$length >= b2 & sub.50$length < b3, 'length'])/tot.len
  
  ## bin4
  indiv.froh[indiv.froh$id == i, 'froh.10.bin4'] <- sum(sub.10[sub.10$length >= b3 & sub.10$length < b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin4'] <- sum(sub.30[sub.30$length >= b3 & sub.30$length < b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin4'] <- sum(sub.50[sub.50$length >= b3 & sub.50$length < b4, 'length'])/tot.len
  
  ## bin5
  indiv.froh[indiv.froh$id == i, 'froh.10.bin5'] <- sum(sub.10[sub.10$length >= b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin5'] <- sum(sub.30[sub.30$length >= b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin5'] <- sum(sub.50[sub.50$length >= b4, 'length'])/tot.len
}

## visualize individual f(ROH) values across length bins
ln.alph <- 0.3
pt.alph <- 0.5
y.max <- max(indiv.froh[,(grep('bin', colnames(indiv.froh)))])


k <- 1
while(k == 1){
  pdf('../krat_genetics_scripts/figures_output/length_bin_froh_missingness.pdf', width = 10, height = 4)
  par(mfrow = c(1,3))
  ## m - 10
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = 'missingness = 10')
    axis(1, at = c(1:5))
    for(i in unique(indiv.froh$id)){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, temp$froh.10.bin5), col = alpha(plink.col, ln.alph))
      points(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, temp$froh.10.bin5), col = alpha(plink.col, pt.alph),
             pch = 19)
    }
    for(i in prob){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, 
                      temp$froh.10.bin5), col = alpha('orange1', orange.lines))
      points(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, 
                       temp$froh.10.bin5), col = alpha('orange1', orange.pts),
             pch = 19)
    }
  
  ## m = 30
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = 'missingness = 30')
    axis(1, at = c(1:5))
    for(i in unique(indiv.froh$id)){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, temp$froh.30.bin5), col = alpha(plink.col, ln.alph))
      points(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, temp$froh.30.bin5), col = alpha(plink.col, pt.alph),
             pch = 19)
    }
    for(i in prob){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, 
                      temp$froh.30.bin5), col = alpha('orange1', orange.lines))
      points(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, 
                       temp$froh.30.bin5), col = alpha('orange1', orange.pts),
             pch = 19)
    }
  
  ## m = 50
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = 'missingness = 50')
    axis(1, at = c(1:4))
    for(i in unique(indiv.froh$id)){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, temp$froh.50.bin5), col = alpha(plink.col, ln.alph))
      points(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, temp$froh.50.bin5), col = alpha(plink.col, pt.alph),
             pch = 19)
    }
    for(i in prob){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, 
                      temp$froh.50.bin5), col = alpha('orange1', orange.lines))
      points(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, 
                       temp$froh.50.bin5), col = alpha('orange1', orange.pts),
             pch = 19)
    }

  ## add legend with bin cutoff information
  legend('topright', legend = c(paste0('cut 1: ',b1), paste0('cut 2: ',b2), paste0('cut 3: ',b3), paste0('cut 4: ',b4)), bty = 'n')
  
  dev.off()
  k <- k+1
}


##### EFfects of het (while missingness == 20) #####
## phwh = 1
plink.roh.10 <- read.table('plink_round5/round5_hom_files/krat_plink_roh_phwh_1_phwm_20_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh.10 <- plink.roh.10[,c(1, 4, 7, 8, 10)]
colnames(plink.roh.10) <- c('id','contig','start','end','n.snps')
plink.roh.10$length <- plink.roh.10$end - plink.roh.10$start + 1
plink.roh.10 <- plink.roh.10[plink.roh.10$length >= 100000,]

## phwh = 2
plink.roh.30 <- read.table('plink_round5/round5_hom_files/krat_plink_roh_phwh_2_phwm_20_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh.30 <- plink.roh.30[,c(1, 4, 7, 8, 10)]
colnames(plink.roh.30) <- c('id','contig','start','end','n.snps')
plink.roh.30$length <- plink.roh.30$end - plink.roh.30$start + 1
plink.roh.30 <- plink.roh.30[plink.roh.30$length >= 100000,]

## phwh = 3
plink.roh.50 <- read.table('plink_round5/round5_hom_files/krat_plink_roh_phwh_3_phwm_20_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh.50 <- plink.roh.50[,c(1, 4, 7, 8, 10)]
colnames(plink.roh.50) <- c('id','contig','start','end','n.snps')
plink.roh.50$length <- plink.roh.50$end - plink.roh.50$start + 1
plink.roh.50 <- plink.roh.50[plink.roh.50$length >= 100000,]

## set f(ROH) formatting for plots
froh.lab <- substitute(paste(italic('F')[ROH]))

##### >>> Calculate and compare individual f(ROH) values between methods #####
OUT <- NULL
for(i in unique(plink.roh.10$id)){
  save <- c(as.numeric(i), sum(plink.roh.10[plink.roh.10$id == i, 'length'])/tot.len, 
            sum(plink.roh.30[plink.roh.30$id == i, 'length'])/tot.len,
            sum(plink.roh.50[plink.roh.50$id == i, 'length'])/tot.len)
  OUT <- rbind(OUT, save)
}
indiv.froh <- as.data.frame(OUT)
colnames(indiv.froh) <- c('id','froh.10','froh.30','froh.50')

pdf('../krat_genetics_scripts/figures_output/het_sites_fROH_comps.pdf', width = 6, height = 6)
plot(indiv.froh$froh.10, indiv.froh$froh.30, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = '1 het site', ylab = '2 het sites')
  points(indiv.froh[indiv.froh$id %in% prob, 'froh.10'], indiv.froh[indiv.froh$id %in% prob, 'froh.30'], col = 'orange1', pch = 19)
  abline(0,1, lty = 2)

plot(indiv.froh$froh.10, indiv.froh$froh.50, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = '1 het site', ylab = '3 het sites')
  points(indiv.froh[indiv.froh$id %in% prob, 'froh.10'], indiv.froh[indiv.froh$id %in% prob, 'froh.50'], col = 'orange1', pch = 19)
  abline(0,1, lty = 2)

plot(indiv.froh$froh.30, indiv.froh$froh.50, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = '2 het sites', ylab = '3 het sites')
  points(indiv.froh[indiv.froh$id %in% prob, 'froh.30'], indiv.froh[indiv.froh$id %in% prob, 'froh.50'], col = 'orange1', pch = 19)
  abline(0,1, lty = 2)

dev.off()  

##### >>> Calculate length-specific f(ROH) values #####
b1 <- 5e5
b2 <- 1e6
b3 <- 2e6
b4 <- 4e6

for(i in unique(indiv.froh$id)){
  sub.10 <- plink.roh.10[plink.roh.10$id == i,]
  sub.30 <- plink.roh.30[plink.roh.30$id == i,]
  sub.50 <- plink.roh.50[plink.roh.50$id == i,]
  
  ## bin1
  indiv.froh[indiv.froh$id == i, 'froh.10.bin1'] <- sum(sub.10[sub.10$length < b1, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin1'] <- sum(sub.30[sub.30$length < b1, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin1'] <- sum(sub.50[sub.50$length < b1, 'length'])/tot.len
  
  ## bin2
  indiv.froh[indiv.froh$id == i, 'froh.10.bin2'] <- sum(sub.10[sub.10$length >= b1 & sub.10$length < b2, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin2'] <- sum(sub.30[sub.30$length >= b1 & sub.30$length < b2, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin2'] <- sum(sub.50[sub.50$length >= b1 & sub.50$length < b2, 'length'])/tot.len
  
  ## bin3
  indiv.froh[indiv.froh$id == i, 'froh.10.bin3'] <- sum(sub.10[sub.10$length >= b2 & sub.10$length < b3, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin3'] <- sum(sub.30[sub.30$length >= b2 & sub.30$length < b3, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin3'] <- sum(sub.50[sub.50$length >= b2 & sub.50$length < b3, 'length'])/tot.len
  
  ## bin4
  indiv.froh[indiv.froh$id == i, 'froh.10.bin4'] <- sum(sub.10[sub.10$length >= b3 & sub.10$length < b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin4'] <- sum(sub.30[sub.30$length >= b3 & sub.30$length < b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin4'] <- sum(sub.50[sub.50$length >= b3 & sub.50$length < b4, 'length'])/tot.len
  
  ## bin5
  indiv.froh[indiv.froh$id == i, 'froh.10.bin5'] <- sum(sub.10[sub.10$length >= b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.30.bin5'] <- sum(sub.30[sub.30$length >= b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'froh.50.bin5'] <- sum(sub.50[sub.50$length >= b4, 'length'])/tot.len
}

## visualize individual f(ROH) values across length bins
ln.alph <- 0.3
pt.alph <- 0.5
orange.pts <- 0.8
orange.lines <- 0.6
y.max <- max(indiv.froh[,(grep('bin', colnames(indiv.froh)))])


k <- 1
while(k == 1){
  pdf('../krat_genetics_scripts/figures_output/length_bin_froh_het_sites.pdf', width = 10, height = 4)
  par(mfrow = c(1,3))
  ## m - 10
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = '1 het sites')
    axis(1, at = c(1:5))
    for(i in unique(indiv.froh$id)){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, 
                      temp$froh.10.bin5), col = alpha(plink.col, ln.alph))
      points(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, 
                       temp$froh.10.bin5), col = alpha(plink.col, pt.alph),
             pch = 19)
    }
    for(i in prob){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, 
                      temp$froh.10.bin5), col = alpha('orange1', orange.lines))
      points(c(1:5), c(temp$froh.10.bin1, temp$froh.10.bin2, temp$froh.10.bin3, temp$froh.10.bin4, 
                       temp$froh.10.bin5), col = alpha('orange1', orange.pts),
             pch = 19)
    }
  
  ## m = 30
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = '2 het sites')
    axis(1, at = c(1:5))
    for(i in unique(indiv.froh$id)){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, temp$froh.30.bin5), col = alpha(plink.col, ln.alph))
      points(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, temp$froh.30.bin5), col = alpha(plink.col, pt.alph),
             pch = 19)
    }
    for(i in prob){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, 
                      temp$froh.30.bin5), col = alpha('orange1', orange.lines))
      points(c(1:5), c(temp$froh.30.bin1, temp$froh.30.bin2, temp$froh.30.bin3, temp$froh.30.bin4, 
                       temp$froh.30.bin5), col = alpha('orange1', orange.pts),
             pch = 19)
    }
  
  ## m = 50
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = '3 het sites')
    axis(1, at = c(1:5))
    for(i in unique(indiv.froh$id)){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, temp$froh.50.bin5), col = alpha(plink.col, ln.alph))
      points(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, temp$froh.50.bin5), col = alpha(plink.col, pt.alph),
             pch = 19)
    }
    for(i in prob){
      temp <- indiv.froh[indiv.froh$id == i,]
      lines(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, 
                      temp$froh.50.bin5), col = alpha('orange1', orange.lines))
      points(c(1:5), c(temp$froh.50.bin1, temp$froh.50.bin2, temp$froh.50.bin3, temp$froh.50.bin4, 
                       temp$froh.50.bin5), col = alpha('orange1', orange.pts),
             pch = 19)
    }
  
  ## add legend with bin cutoff information
  legend('topright', legend = c(paste0('cut 1: ',b1), paste0('cut 2: ',b2), paste0('cut 3: ',b3), paste0('cut 4: ',b4)), bty = 'n')
  
  dev.off()
  k <- k+1
}
