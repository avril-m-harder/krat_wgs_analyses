setwd('/Users/Avril/Documents/krat_genetics/data/')

library(vcfR)
library(scales)

### set colors
library(ghibli)
gt.col <- ghibli_palette('PonyoMedium')[4]
pl.col <- ghibli_palette('PonyoMedium')[2]
plink.col <- ghibli_palette('PonyoMedium')[3]


##### Read in data #####
## SNP data (final VCF used for calling ROHs)
# vcf <- read.vcfR('krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz') ### way too big to be read in

## sample ROH data (bcftools/ROH + final PLINK settings)
pl.roh <- read.table('krat_final_allfiltcontigs_all_samps_PL_RG_ONLY.txt', header = FALSE)
colnames(pl.roh) <- c('RG','id','contig','start','end','length','n.snps','quality')
pl.roh <- pl.roh[pl.roh$length >= 100000,]
pl.roh <- pl.roh[,-1]
pl.roh$id <- do.call(rbind, strsplit(pl.roh$id, split = '_'))[,1]

gt.roh <- read.table('krat_final_allfiltcontigs_all_samps_GT_RG_ONLY.txt', header = FALSE)
colnames(gt.roh) <- c('RG','id','contig','start','end','length','n.snps','quality')
gt.roh <- gt.roh[gt.roh$length >= 100000,]
gt.roh <- gt.roh[,-1]
gt.roh$id <- do.call(rbind, strsplit(gt.roh$id, split = '_'))[,1]

## optional filtering of BCFtools results by quality (not sure what this means)
gt.roh <- gt.roh[gt.roh$quality >= 40,]
pl.roh <- pl.roh[pl.roh$quality >= 40,]

## sample ROH data (PLINK)
plink.roh <- read.table('plink_final_settings/krat_plink_roh_phwh_2_phwm_30_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh <- plink.roh[,c(1, 4, 7, 8, 10)]
colnames(plink.roh) <- c('id','contig','start','end','n.snps')
plink.roh$length <- plink.roh$end - plink.roh$start + 1
plink.roh <- plink.roh[plink.roh$length >= 100000,]

## sample information
samps <- read.csv('../preseq_sample_information/summary_allstats2_withgen_zeroes.csv')
relates <- read.csv('../data/pairwise_samp_relationships.csv')

## below calculated using bcftools stats on krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz
vcf.stats <- read.csv('../qc_stuff/bcftools_stats/just_sample_stats.txt', check.names = FALSE)
vcf.stats <- vcf.stats[,c('sample','nHets','average depth','nMissing')]
colnames(vcf.stats) <- c('id','n.hets','avg.depth','n.missing')
n.sites <- 11796181 ## # of SNPs in this VCF
vcf.stats$prop.miss <- vcf.stats$n.missing/n.sites
vcf.stats$norm.prop.het <- vcf.stats$n.hets/(n.sites - vcf.stats$n.missing)

## contig information
contigs <- read.csv('contig_lengths.csv')
contigs$length <- as.numeric(contigs$length)
contigs <- contigs[order(-contigs$length),]
contigs$c.index <- c(1:nrow(contigs)) ## create contig index, where c.index = 1 == longest contig
tot.len <- sum(contigs$length) ## 2,571,239,112 bp ## length of filtered contigs
## add this contig information to ROH data
gt.roh <- merge(gt.roh, contigs[,c(1,3)], by = 'contig')
pl.roh <- merge(pl.roh, contigs[,c(1,3)], by = 'contig')
plink.roh <- merge(plink.roh, contigs[,c(1,3)], by = 'contig')

## limit all data to contigs >= 1 Mb in length
contigs <- contigs[contigs$length >= 1e6,]
gt.roh <- gt.roh[gt.roh$c.index %in% contigs$c.index,]
pl.roh <- pl.roh[pl.roh$c.index %in% contigs$c.index,]
plink.roh <- plink.roh[plink.roh$c.index %in% contigs$c.index,]

## set f(ROH) formatting for plots
froh.lab <- substitute(paste(italic('F')[ROH]))

# ##### Examine ROH distribution across contig characteristics #####
# temp.contigs <- contigs
# for(r in 1:nrow(temp.contigs)){
#   temp.contigs$gt.tot.len[r] <- sum(gt.roh[gt.roh$contig == temp.contigs$contig[r], 'length'])
#   temp.contigs$gt.num.roh[r] <- nrow(gt.roh[gt.roh$contig == temp.contigs$contig[r],])
#   temp.contigs$pl.tot.len[r] <- sum(pl.roh[pl.roh$contig == temp.contigs$contig[r], 'length'])
#   temp.contigs$pl.num.roh[r] <- nrow(pl.roh[pl.roh$contig == temp.contigs$contig[r],])
#   temp.contigs$plink.tot.len[r] <- sum(plink.roh[plink.roh$contig == temp.contigs$contig[r], 'length'])
#   temp.contigs$plink.num.roh[r] <- nrow(plink.roh[plink.roh$contig == temp.contigs$contig[r],])
# }
# ## total ROH length across individuals
# plot(temp.contigs$length, temp.contigs$gt.tot.len, pch = 19, cex = 0.5, col = alpha(gt.col, 0.4),
#      xlab = 'Contig length', ylab = 'Total ROH length across individuals', main = 'Genotypes')
# plot(temp.contigs$length, temp.contigs$gt.tot.len, pch = 19, cex = 0.5, col = alpha(gt.col, 0.4),
#      xlim = c(0, 1e7), ylim = c(0, 5e7), xlab = 'Contig length', ylab = 'Total ROH length across individuals', main = 'Genotypes')
# plot(temp.contigs$length, temp.contigs$pl.tot.len, pch = 19, cex = 0.5, col = alpha(pl.col, 0.4),
#      xlab = 'Contig length', ylab = 'Total ROH length across individuals', main = 'Likelihoods')
# plot(temp.contigs$length, temp.contigs$pl.tot.len, pch = 19, cex = 0.5, col = alpha(pl.col, 0.4),
#      xlim = c(0, 1e7), ylim = c(0, 5e7), xlab = 'Contig length', ylab = 'Total ROH length across individuals', main = 'Likelihoods')
# plot(temp.contigs$length, temp.contigs$plink.tot.len, pch = 19, cex = 0.5, col = alpha(plink.col, 0.4),
#      xlab = 'Contig length', ylab = 'Total ROH length across individuals', main = 'PLINK')
# plot(temp.contigs$length, temp.contigs$plink.tot.len, pch = 19, cex = 0.5, col = alpha(plink.col, 0.4),
#      xlim = c(0, 1e7), ylim = c(0, 5e7), xlab = 'Contig length', ylab = 'Total ROH length across individuals', main = 'PLINK')
# 
# ## total # of ROHs across individuals
# plot(temp.contigs$length, temp.contigs$gt.num.roh, pch = 19, cex = 0.5, col = alpha(gt.col, 0.4),
#      xlab = 'Contig length', ylab = 'Total number of ROHs across individuals', main = 'Genotypes')
# plot(temp.contigs$length, temp.contigs$gt.num.roh, pch = 19, cex = 0.5, col = alpha(gt.col, 0.4),
#      xlim = c(0, 1e7), ylim = c(0, 250), xlab = 'Contig length', ylab = 'Total number of ROHs across individuals', main = 'Genotypes')
# plot(temp.contigs$length, temp.contigs$pl.num.roh, pch = 19, cex = 0.5, col = alpha(pl.col, 0.4),
#      xlab = 'Contig length', ylab = 'Total number of ROHs across individuals', main = 'Likelihoods')
# plot(temp.contigs$length, temp.contigs$pl.num.roh, pch = 19, cex = 0.5, col = alpha(pl.col, 0.4),
#      xlim = c(0, 1e7), ylim = c(0, 250), xlab = 'Contig length', ylab = 'Total number of ROHs across individuals', main = 'Likelihoods')
# plot(temp.contigs$length, temp.contigs$plink.num.roh, pch = 19, cex = 0.5, col = alpha(plink.col, 0.4),
#      xlab = 'Contig length', ylab = 'Total number of ROHs across individuals', main = 'PLINK')
# plot(temp.contigs$length, temp.contigs$plink.num.roh, pch = 19, cex = 0.5, col = alpha(plink.col, 0.4),
#      xlim = c(0, 1e7), ylim = c(0, 250), xlab = 'Contig length', ylab = 'Total number of ROHs across individuals', main = 'PLINK')
# 
# rm(temp.contigs)
## nothing super concerning here. main note is that PLINK appears to call ROHs more frequently on short contigs than 
## either BCFtools method (short ~ <= 5 Mb)

##### Calculate and compare individual f(ROH) values between methods #####
OUT <- NULL
for(i in unique(gt.roh$id)){
  save <- c(as.numeric(i), sum(gt.roh[gt.roh$id == i, 'length'])/tot.len, 
            sum(pl.roh[pl.roh$id == i, 'length'])/tot.len,
            sum(plink.roh[plink.roh$id == i, 'length'])/tot.len)
  OUT <- rbind(OUT, save)
}
indiv.froh <- as.data.frame(OUT)
colnames(indiv.froh) <- c('id','gt.froh','pl.froh','plink.froh')
indiv.froh <- merge(indiv.froh, vcf.stats, by = 'id')

plot(indiv.froh$plink.froh, indiv.froh$pl.froh, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'PLINK', ylab = 'Likelihoods')
abline(0,1, lty = 2)

plot(indiv.froh$plink.froh, indiv.froh$gt.froh, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'PLINK', ylab = 'Genotypes')
abline(0,1, lty = 2)

plot(indiv.froh$gt.froh, indiv.froh$pl.froh, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'Genotypes', ylab = 'Likelihoods')
abline(0,1, lty = 2)

### plot cumulative f(ROH) for all individuals (PLINK only)
plot(0, 0, col = 'transparent', xlim = c(1e5, max(plink.roh$length)+1e6), ylim = c(0, max(indiv.froh$plink.froh)),
     xlab = 'ROH length', ylab = froh.lab)
for(i in unique(plink.roh$id)){
  sub <- plink.roh[plink.roh$id == i,]
  sub <- sub[order(sub$length),]
  sub$cumsum <- cumsum(sub$length)
  sub$cumfroh <- sub$cumsum/tot.len
  lines(sub$length, sub$cumfroh)
  points(8.5e6, indiv.froh[indiv.froh$id == i, 'plink.froh'], col = alpha(plink.col, 0.2), pch = 19)
}



##### Calculate length-specific f(ROH) values #####
b1 <- 5e5
b2 <- 1e6
b3 <- 2e6
b4 <- 4e6

for(i in unique(indiv.froh$id)){
  gt.sub <- gt.roh[gt.roh$id == i,]
  pl.sub <- pl.roh[pl.roh$id == i,]
  plink.sub <- plink.roh[plink.roh$id == i,]
  
  ## bin1
  indiv.froh[indiv.froh$id == i, 'gt.bin1'] <- sum(gt.sub[gt.sub$length < b1, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'pl.bin1'] <- sum(pl.sub[pl.sub$length < b1, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'plink.bin1'] <- sum(plink.sub[plink.sub$length < b1, 'length'])/tot.len
  
  ## bin2
  indiv.froh[indiv.froh$id == i, 'gt.bin2'] <- sum(gt.sub[gt.sub$length >= b1 & gt.sub$length < b2, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'pl.bin2'] <- sum(pl.sub[pl.sub$length >= b1 & pl.sub$length < b2, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'plink.bin2'] <- sum(plink.sub[plink.sub$length >= b1 & plink.sub$length < b2, 'length'])/tot.len
  
  ## bin3
  indiv.froh[indiv.froh$id == i, 'gt.bin3'] <- sum(gt.sub[gt.sub$length >= b2 & gt.sub$length < b3, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'pl.bin3'] <- sum(pl.sub[pl.sub$length >= b2 & pl.sub$length < b3, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'plink.bin3'] <- sum(plink.sub[plink.sub$length >= b2 & plink.sub$length < b3, 'length'])/tot.len
  
  ## bin4
  indiv.froh[indiv.froh$id == i, 'gt.bin4'] <- sum(gt.sub[gt.sub$length >= b3 & gt.sub$length < b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'pl.bin4'] <- sum(pl.sub[pl.sub$length >= b3 & pl.sub$length < b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'plink.bin4'] <- sum(plink.sub[plink.sub$length >= b3 & plink.sub$length < b4, 'length'])/tot.len
  
  ## bin5
  indiv.froh[indiv.froh$id == i, 'gt.bin5'] <- sum(gt.sub[gt.sub$length >= b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'pl.bin5'] <- sum(pl.sub[pl.sub$length >= b4, 'length'])/tot.len
  indiv.froh[indiv.froh$id == i, 'plink.bin5'] <- sum(plink.sub[plink.sub$length >= b4, 'length'])/tot.len
}

## visualize individual f(ROH) values across length bins
ln.alph <- 0.3
pt.alph <- 0.5
y.max <- max(indiv.froh[,(grep('bin', colnames(indiv.froh)))])


k <- 1
while(k == 1){
  par(mfrow = c(1,3))
  ## Genotypes
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = 'BCFtools Genotypes')
  axis(1, at = c(1:4))
  for(i in unique(indiv.froh$id)){
    temp <- indiv.froh[indiv.froh$id == i,]
    lines(c(1:5), c(temp$gt.bin1, temp$gt.bin2, temp$gt.bin3, temp$gt.bin4, temp$gt.bin5), col = alpha(gt.col, ln.alph))
    points(c(1:5), c(temp$gt.bin1, temp$gt.bin2, temp$gt.bin3, temp$gt.bin4, temp$gt.bin5), col = alpha(gt.col, pt.alph),
           pch = 19)
  }
  # lines(c(4.8, 5.2), c(median(indiv.froh$gt.froh), median(indiv.froh$gt.froh)), col = gt.col, lwd = 2)
  
  ## Likelihoods
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = 'BCFtools Likelihoods')
  axis(1, at = c(1:4))
  for(i in unique(indiv.froh$id)){
    temp <- indiv.froh[indiv.froh$id == i,]
    lines(c(1:5), c(temp$pl.bin1, temp$pl.bin2, temp$pl.bin3, temp$pl.bin4, temp$pl.bin5), col = alpha(pl.col, ln.alph))
    points(c(1:5), c(temp$pl.bin1, temp$pl.bin2, temp$pl.bin3, temp$pl.bin4, temp$pl.bin5), col = alpha(pl.col, pt.alph),
           pch = 19)
  }
  # lines(c(4.8, 5.2), c(median(indiv.froh$pl.froh), median(indiv.froh$pl.froh)), col = pl.col, lwd = 2)
  
  ## PLINK
  plot(0, 0, xlim = c(1,5), ylim = c(0, y.max), ylab = froh.lab, xaxt = 'n', xlab = 'Length bin',
       main = 'PLINK')
  axis(1, at = c(1:4))
  for(i in unique(indiv.froh$id)){
    temp <- indiv.froh[indiv.froh$id == i,]
    lines(c(1:5), c(temp$plink.bin1, temp$plink.bin2, temp$plink.bin3, temp$plink.bin4, temp$plink.bin5), col = alpha(plink.col, ln.alph))
    points(c(1:5), c(temp$plink.bin1, temp$plink.bin2, temp$plink.bin3, temp$plink.bin4, temp$plink.bin5), col = alpha(plink.col, pt.alph),
           pch = 19)
  }
  # lines(c(4.8, 5.2), c(median(indiv.froh$plink.froh), median(indiv.froh$plink.froh)), col = plink.col, lwd = 2)
  
  ## add legend with bin cutoff information
  legend('topright', legend = c(paste0('cut 1: ',b1), paste0('cut 2: ',b2), paste0('cut 3: ',b3), paste0('cut 4: ',b4)), bty = 'n')
  
  k <- k+1
}

##### Look for relationships between individual f(ROH) values and fitness characteristics #####
## get information from pedigree results
samps$age.at.death <- samps$deathyear - samps$birthyear
temp.samps <- samps[,c('id','inb','off','off_survive','nmates','age.at.death','genback','sex')]
indiv.froh <- merge(indiv.froh, temp.samps, by = 'id')
colnames(indiv.froh)[25] <- 'ped.inb'
rm(temp.samps)
write.csv(indiv.froh, 'individual_frohs_1Mbcontigs.csv', row.names = FALSE)
