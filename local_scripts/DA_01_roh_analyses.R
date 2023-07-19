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


pdf('../krat_genetics_scripts/figures_output/method_fROH_comps.pdf', width = 6, height = 6)
plot(indiv.froh$plink.froh, indiv.froh$pl.froh, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'PLINK', ylab = 'Likelihoods')
  abline(0,1, lty = 2)
  
plot(indiv.froh$plink.froh, indiv.froh$gt.froh, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'PLINK', ylab = 'Genotypes')
  abline(0,1, lty = 2)
  
plot(indiv.froh$gt.froh, indiv.froh$pl.froh, pch = 19, col = alpha('darkorchid4', 0.3),
     xlab = 'Genotypes', ylab = 'Likelihoods')
  abline(0,1, lty = 2)

dev.off()

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
pdf('../krat_genetics_scripts/figures_output/length_bin_froh_indiv_lines_3methods.pdf', width = 10, height = 4)
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
  
dev.off()
k <- k+1
}

# ##### Look for relationships between individual f(ROH) values and fitness characteristics #####
# ## get information from pedigree results
# samps$age.at.death <- samps$deathyear - samps$birthyear
# temp.samps <- samps[,c('id','inb','off','off_survive','nmates','age.at.death','genback','sex')]
# indiv.froh <- merge(indiv.froh, temp.samps, by = 'id')
# colnames(indiv.froh)[25] <- 'ped.inb'
# rm(temp.samps)
# write.csv(indiv.froh, 'individual_frohs.csv', row.names = FALSE)
# 
# ## overall f(ROH) relationships
# ## # offspring
# pdf('../krat_genetics_scripts/figures_output/froh_and_bins_vs_offspring.pdf', width = 9, height = 3.5)
# par(mfrow = c(1,3))
# for(f in c('froh','bin1','bin2','bin3','bin4','bin5')){
#   plot(indiv.froh[,grep(paste0('gt.',f), colnames(indiv.froh))], indiv.froh$off, pch = 19, col = gt.col, 
#        xlab = froh.lab, ylab = 'Number of offspring', main = paste0('Genotypes - ',f))
#   plot(indiv.froh[,grep(paste0('pl.',f), colnames(indiv.froh))], indiv.froh$off, pch = 19, col = pl.col, 
#        xlab = froh.lab, ylab = 'Number of offspring', main = paste0('Likelihoods - ',f))
#   plot(indiv.froh[,grep(paste0('plink.',f), colnames(indiv.froh))], indiv.froh$off, pch = 19, col = plink.col, 
#        xlab = froh.lab, ylab = 'Number of offspring', main = paste0('PLINK - ',f))
# }
# dev.off()
# 
# ## # surviving offspring 
# pdf('../krat_genetics_scripts/figures_output/froh_and_bins_vs_surviving_offspring.pdf', width = 9, height = 3.5)
# par(mfrow = c(1,3))
# for(f in c('froh','bin1','bin2','bin3','bin4','bin5')){
#   plot(indiv.froh[,grep(paste0('gt.',f), colnames(indiv.froh))], indiv.froh$off_survive, pch = 19, col = gt.col, 
#        xlab = froh.lab, ylab = 'Number of surviving offspring', main = paste0('Genotypes - ',f))
#   plot(indiv.froh[,grep(paste0('pl.',f), colnames(indiv.froh))], indiv.froh$off_survive, pch = 19, col = pl.col, 
#        xlab = froh.lab, ylab = 'Number of surviving offspring', main = paste0('Likelihoods - ',f))
#   plot(indiv.froh[,grep(paste0('plink.',f), colnames(indiv.froh))], indiv.froh$off_survive, pch = 19, col = plink.col, 
#        xlab = froh.lab, ylab = 'Number of surviving offspring', main = paste0('PLINK - ',f))
# }
# dev.off()
# 
# ## age at death
# pdf('../krat_genetics_scripts/figures_output/froh_and_bins_vs_ageatdeath.pdf', width = 9, height = 3.5)
# par(mfrow = c(1,3))
# for(f in c('froh','bin1','bin2','bin3','bin4','bin5')){
#   plot(indiv.froh[,grep(paste0('gt.',f), colnames(indiv.froh))], indiv.froh$age.at.death, pch = 19, col = gt.col, 
#        xlab = froh.lab, ylab = 'Age at death', main = paste0('Genotypes - ',f))
#   plot(indiv.froh[,grep(paste0('pl.',f), colnames(indiv.froh))], indiv.froh$age.at.death, pch = 19, col = pl.col, 
#        xlab = froh.lab, ylab = 'Age at death', main = paste0('Likelihoods - ',f))
#   plot(indiv.froh[,grep(paste0('plink.',f), colnames(indiv.froh))], indiv.froh$age.at.death, pch = 19, col = plink.col, 
#        xlab = froh.lab, ylab = 'Age at death', main = paste0('PLINK - ',f))
# }
# dev.off()
# 
# ## Nothing obvious or super compelling going on here, just because there are so many points clustered to the left in each plot.
 

##### Identify ROHverlap between all pairs of individuals #####
##### >>> PLINK #####
# PLINK.OVERLAP <- NULL
# plink.roh$roh.id <- c(1:nrow(plink.roh))
# 
## should be 1,128 combos
pwise <- expand.grid(unique(plink.roh$id), unique(plink.roh$id))
pwise <- pwise[-which(pwise$Var1 == pwise$Var2),]
KEEP.PWISE <- NULL
for(r in 1:nrow(pwise)){
  if(r == 1){
    KEEP.PWISE <- rbind(KEEP.PWISE, pwise[r,])
  }
  if(r > 1 &
     nrow(KEEP.PWISE[KEEP.PWISE[,1] == pwise[r,1] & KEEP.PWISE[,2] == pwise[r,2],]) < 1 &
     nrow(KEEP.PWISE[KEEP.PWISE[,2] == pwise[r,1] & KEEP.PWISE[,1] == pwise[r,2],]) < 1){
    KEEP.PWISE <- rbind(KEEP.PWISE, pwise[r,])
  }
}
pwise <- KEEP.PWISE
rm(KEEP.PWISE)
colnames(pwise) <- c('id1','id2')
# 
# ct <- 1
# for(i in unique(plink.roh$id)){              ## for each individual,
#   # print(paste0('PLINK - ',i))
#   ## and for each individual the focal individual has a relationship with,
#   comps <- unique(pwise[pwise$id1 == i, 'id2'])
#   print(paste0(ct,'/48 - ',length(comps),' comps'))
#   for(c in comps){
#     for(ind in contigs$c.index){
#       sub.i.plink <- plink.roh[plink.roh$id == i & plink.roh$c.index == ind,]
#       sub.c.plink <- plink.roh[plink.roh$id == c & plink.roh$c.index == ind,]
#       if(nrow(sub.i.plink) > 0 & nrow(sub.c.plink) > 0){
#         for(r in 1:nrow(sub.i.plink)){                     ## loop over ROHs in focal individual,
#           c.index <- sub.i.plink$c.index[r]                ## save c.index ID,
#           s <- sub.i.plink$start[r]                        ## save start,
#           e <- sub.i.plink$end[r]                          ## and end.
#           ## check for overlaps with individual being compared against
#           ## comp ROHs beginning outside of focal ROH, ending inside
#           if(nrow(sub.c.plink[sub.c.plink$start < s & sub.c.plink$end >= s & sub.c.plink$end <= e,]) > 0){
#             temp <- sub.c.plink[sub.c.plink$start < s & sub.c.plink$end >= s & sub.c.plink$end <= e,]
#             for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#               rohverlap.len <- temp$end[t] - s + 1                             ## length of overlap with the focal ROH
#               save <- as.numeric(c(i, c, c.index, s, temp$end[t], rohverlap.len, temp$roh.id[t], sub.i.plink$roh.id[r]))
#               PLINK.OVERLAP <- rbind(PLINK.OVERLAP, save)
#             }
#           }
#           ## comp ROHs beginning inside of a focal ROH, ending outside
#           if(nrow(sub.c.plink[sub.c.plink$start >= s & sub.c.plink$start <= e & sub.c.plink$end > e,]) > 0){
#             temp <- sub.c.plink[sub.c.plink$start >= s & sub.c.plink$start <= e & sub.c.plink$end > e,]
#             for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#               rohverlap.len <- e - temp$start[t] + 1                           ## length of overlap with the focal ROH
#               save <- as.numeric(c(i, c, c.index, temp$start[t], e, rohverlap.len, temp$roh.id[t], sub.i.plink$roh.id[r]))
#               PLINK.OVERLAP <- rbind(PLINK.OVERLAP, save)
#             }
#           }
#           ## comp ROHs completely covering a focal ROH
#           if(nrow(sub.c.plink[sub.c.plink$start < s & sub.c.plink$end > e,]) > 0){
#             temp <- sub.c.plink[sub.c.plink$start < s & sub.c.plink$end > e,]
#             for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#               rohverlap.len <- e - s + 1                                       ## length of overlap with the focal ROH,
#               save <- as.numeric(c(i, c, c.index, s, e, rohverlap.len, temp$roh.id[t], sub.i.plink$roh.id[r]))
#               PLINK.OVERLAP <- rbind(PLINK.OVERLAP, save)
#             }
#           }
#           ## comp ROHs completely within a focal ROH
#           if(nrow(sub.c.plink[sub.c.plink$start >= s & sub.c.plink$end <= e,]) > 0){
#             temp <- sub.c.plink[sub.c.plink$start >= s & sub.c.plink$end <= e,]
#             for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#               rohverlap.len <- temp$length[t]                                  ## length of overlap with the focal ROH,
#               save <- as.numeric(c(i, c, c.index, temp$start[t], temp$end[t], rohverlap.len, temp$roh.id[t], sub.i.plink$roh.id[r]))
#               PLINK.OVERLAP <- rbind(PLINK.OVERLAP, save)
#             }
#           }
#         }
#       }
#     }
#   }
#   ct <- ct+1
# }
# colnames(PLINK.OVERLAP) <- c('id1','id2','c.index','start.overlap','end.overlap','len.overlap','roh.id.2','roh.id.1')
# 
# plink.overlap <- as.data.frame(PLINK.OVERLAP)
# write.csv(plink.overlap, 'plink_rohverlap_ALL_SAMPLES.csv', row.names = FALSE)

##### >>> read in ROHverlap results #####
plink.overlap <- read.csv('plink_rohverlap_ALL_SAMPLES.csv')

##### >>> for each pair, calculate f(ROHverlap) #####
# for(r in 1:nrow(relates)){
#   relates$plink.frohverlap[r] <- sum(plink.overlap[plink.overlap$id1 == relates$id1[r] & 
#                                                      plink.overlap$id2 == relates$id2[r], 'len.overlap'])/tot.len
# }
# 
# ## plot some stuff for different relationship types
# plot(relates$mates, relates$plink.frohverlap, col = alpha(plink.col, 0.5), pch = 19, xlab = 'Mates or no', ylab = 'f(ROH)verlap')
# plot(relates$parent.off, relates$plink.frohverlap, col = alpha(plink.col, 0.5), pch = 19, xlab = 'Parent-offspring or no', ylab = 'f(ROH)verlap')
# plot(relates$half.sibs, relates$plink.frohverlap, col = alpha(plink.col, 0.5), pch = 19, xlab = 'Half-sibs or no', ylab = 'f(ROH)verlap')
# plot(relates$full.sibs, relates$plink.frohverlap, col = alpha(plink.col, 0.5), pch = 19, xlab = 'Full-sibs or no', ylab = 'f(ROH)verlap')
# 
# ## I think there's too much weird stuff going on in this little incestuous population for ROHverlap to turn out like expected (e.g., super high between 2 full-sibs because their parents were also full-sibs to one another)
# 
# plot(relates[relates$mates == 1, 'plink.frohverlap'], relates[relates$mates == 1, 'n.offspring'], 
#      col = alpha(plink.col, 0.5), pch = 19, ylab = '# offspring', xlab = 'f(ROH)verlap')
# 
# plot(relates[relates$mates == 1, 'plink.frohverlap'], relates[relates$mates == 1, 'n.surv.offspring'], 
#      col = alpha(plink.col, 0.5), pch = 19, ylab = '# surving offspring', xlab = 'f(ROH)verlap')

##### *Going with just PLINK from now on, conclusions are pretty conserved across the 3 methods (see comp figures) #####

##### Collapse PLINK ROHverlap regions to build minimized VCF from #####
# OUT <- NULL
# for(c in unique(plink.overlap$c.index)){
#   sub <- plink.overlap[plink.overlap$c.index == c,]
#   summary <- rep(0, contigs[contigs$c.index == c, 'length'])
#   for(r in 1:nrow(sub)){
#     print(paste0(c,' - ',r/nrow(sub)))
#     summary[c(sub$start.overlap[r]:sub$end.overlap[r])] <- summary[c(sub$start.overlap[r]:sub$end.overlap[r])] + 1
#   }
# 
#   ## collapse loc list into ranges --> bcftools regions file for filtering full VCF (just need chrom / start / end)
#   summary[which(summary > 0)] <- 1 ## convert to just 1's (ROHverlap present) and 0's (ROHverlap absent)
#   lens <- rle(summary)$lengths     ## lengths of consecutive value runs
#   p.a <- rle(summary)$values       ## values of those runs
#   for(p in 1:length(p.a)){
#     if(p.a[p] == 0){
#       if(p != length(p.a)){
#         s <- sum(lens[c(1:p)]) + 1
#         e <- sum(lens[c(1:(p + 1))])
#         save <- c(c, s, e)
#         OUT <- rbind(OUT, save)
#       } else{
#         next
#       }
#     } else{
#       next
#     }
#   }
# }

## diagnostic plot to make sure the above worked correctly, trying with largest contig (use different x lim values to zoom in and examine)
# xmin <- 12e6
# xmax <- 13e6 ## max for 1st contig = 60933957
# plot(0, 0, xlim = c(xmin, xmax), col = 'transparent', ylim = c(0.5, 2.5))
#   sub <- plink.overlap[plink.overlap$c.index == 1,]
#   for(r in 1:nrow(sub)){
#     lines(c(sub$start.overlap[r], sub$end.overlap[r]), c(2, 2))
#   }
#   temp <- OUT[OUT[,1] == 1,]
#   for(r in 1:nrow(temp)){
#     lines(c(temp[r,2], temp[r,3]), c(1,1))
#   }
## (looks good)
# 
# coords <- as.data.frame(OUT)
# colnames(coords) <- c('c.index','start','end')
# coords <- coords[order(coords$contig, coords$start),]
# sum(coords$end - coords$start + 1)/tot.len ## cuts it down to about 50% of total contig lengths (75% for all samples)
# coords <- merge(coords, contigs, by = 'c.index')
# coords <- coords[,c('contig','start','end')]
# coords <- coords[order(coords$contig, coords$start),]
# write.table(coords, 'rohverlap_coordinates_ALL_SAMPLES.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

##### Use genotypes from filtered VCF to check for shared alleles in ROHverlap regions #####
## read in VCF sample order and genotype information (collected using bcftools query)
samp.order <- read.table('vcf_sample_order.txt')
samp.order <- unlist(samp.order)

### genotypes_to_numeric.sh script used to convert cluster output (krat_plink_rohverlap_GTs_ALL_SAMPLES.txt) to an easier format
# gts <- read.table('krat_plink_rohverlap_ALL_SAMPLES_numericGTs.txt', sep = '\t')
# gts <- gts[,-3]
# colnames(gts) <- c('contig','pos',samp.order)
# gts <- merge(gts, contigs[,c(1,3)], by = 'contig')
# gts <- gts[,c(51,2:50)]
# for(c in 1:ncol(gts)){
#   gts[,c] <- as.numeric(gts[,c])
# }

## for each pair of individuals with ROHverlap, check that alleles are the same throughout each ROH
# OUT <- NULL    ## base-position mismatches
# OUT1 <- NULL   ## enclosing ROHverlap with at least 1 mismatch
# OUT2 <- NULL   ## ROHverlap with zero shared SNPs included

## don't analyze ROHverlaps < 5 kb in length; some may include zero SNPs
plink.overlap <- plink.overlap[plink.overlap$len.overlap >= 5000, ]

## genotypes from characters --> numeric
## missing or heterozygous == 0
## homozygous reference == 1
## homozygous alternate == 2
# keep.gts <- c(1,2)

# for(r in 1:nrow(pwise)){
# for(r in c(970:nrow(pwise))){
#   print(r/nrow(pwise))
#   sub.overlap <- plink.overlap[plink.overlap$id1 == pwise$id1[r] & plink.overlap$id2 == pwise$id2[r],]
#   sub.gt <- gts[,c(1, 2, grep(pwise$id1[r], colnames(gts)), grep(pwise$id2[r], colnames(gts)))]
#   for(n in 1:nrow(sub.overlap)){
#     # print(paste0(r/nrow(pwise),' - ',n/nrow(sub.overlap)))
#     temp.gt <- sub.gt[sub.gt$pos >= sub.overlap$start.overlap[n] & sub.gt$pos <= sub.overlap$end.overlap[n] & sub.gt$c.index == sub.overlap$c.index[n],]
#     ## some heterozygous and missing sites are allowed, so only compare genotypes that are homozygous for both individuals
#     temp.gt <- temp.gt[which(temp.gt[,3] %in% keep.gts & temp.gt[,4] %in% keep.gts),]
#     ## this could be way more efficient by making all the data numeric only (e.g., 0/0 --> 0)
#     if(nrow(temp.gt) == 0){
#       print(paste0('EMPTY INTERVAL: ',pwise$id1[r],' & ',pwise$id2[r]))
#       OUT2 <- rbind(OUT2, sub.overlap[n,])
#     }
#     if(nrow(temp.gt) > 0 & all(temp.gt[,3] != temp.gt[,4])){ ## if a set of homozygous genotypes doesn't match
#       print(paste0('MISMATCH: ',pwise$id1[r],' & ',pwise$id2[r]))
#       save <- temp.gt[which(temp.gt[,3] != temp.gt[,4]),]
#       save[,5] <- pwise$id1[r]
#       save[,6] <- pwise$id2[r]
#       colnames(save) <- c('contig','pos','gt.1','gt.2','id.1','id.2')
#       OUT <- rbind(OUT, save)
#       OUT1 <- rbind(OUT1, sub.overlap[n,])
#     }
#   }
# }
# 
# write.csv(OUT, 'mismatches_within_plink_fROHverlaps_ALL_SAMPLES.csv', row.names = FALSE, quote = FALSE)
# write.csv(OUT1, 'rohverlap_segments_with_mismatches_ALL_SAMPLES.csv', row.names = FALSE, quote = FALSE)
# write.csv(OUT2, 'rohverlap_segments_with_no_shared_SNPs_ALL_SAMPLES.csv', row.names = FALSE, quote = FALSE)

# pos.mismatch <- read.csv('mismatches_within_plink_fROHverlaps_ALL_SAMPLES.csv')
# roh.mismatch <- read.csv('rohverlap_segments_with_mismatches_ALL_SAMPLES.csv')
# roh.mismatch <- merge(roh.mismatch, contigs, by = 'c.index')
# roh.no.snps <- read.csv('rohverlap_segments_with_no_shared_SNPs_ALL_SAMPLES.csv')

# save.image('DA_01_roh_analyses_ALL_SAMPLES.RData') ## after saving this, delete the temporary .RData file saved above
##### START HERE #####
load('DA_01_roh_analyses_ALL_SAMPLES.RData') ## can read everything up to this point in, some of the files are slow to load individually

##### Remove ROHverlap segments with homozygous mismatches between samples #####
## plot these results to see extent of mismatch issue
# pdf('../krat_genetics_scripts/figures_output/rohverlap_mismatches.pdf', width = 7, height = 4)
# for(r in 1:nrow(roh.mismatch)){
#   sub <- pos.mismatch[which(pos.mismatch$id.1 == roh.mismatch$id1[r] &
#                             pos.mismatch$id.2 == roh.mismatch$id2[r] &
#                             pos.mismatch$contig == roh.mismatch$contig[r]),]
#   plot(0, 0, xlim = c(0.99*roh.mismatch$start.overlap[r], 1.01 *roh.mismatch$end.overlap[r]), ylim = c(0.5, 1.75),
#        main = paste0(roh.mismatch$id1[r],' & ',roh.mismatch$id2[r]), yaxt = 'n', ylab = '', xlab = 'Chromosome position (bp)')
#     lines(c(0.99*roh.mismatch$start.overlap[r], 1.01 *roh.mismatch$end.overlap[r]), y = c(1,1), lwd = 3)
#     text(mean(roh.mismatch$start.overlap[r], roh.mismatch$end.overlap[r]), 0.75, 
#          labels = paste0('ROHverlap length = ',roh.mismatch$len.overlap[r]), adj = 0)
#     points(sub$pos, rep(1.5, nrow(sub)), pch = 19, col = 'red')
# }
# dev.off()
# 
# ## should probably discard any ROHverlap with at least 1 mismatch - only 17 instances of mismatch in nearly 50k ROHverlap segments
# s <- nrow(plink.overlap)
# for(r in 1:nrow(roh.mismatch)){
#   rm <- rownames(plink.overlap[plink.overlap$id1 == roh.mismatch$id1[r] &
#                         plink.overlap$id2 == roh.mismatch$id2[r] &
#                         plink.overlap$c.index == roh.mismatch$c.index[r] &
#                         plink.overlap$start.overlap == roh.mismatch$start.overlap[r] &
#                         plink.overlap$end.overlap == roh.mismatch$end.overlap[r],])
#   plink.overlap <- plink.overlap[-which(rownames(plink.overlap) == rm),]
# }
# s - nrow(roh.mismatch) == nrow(plink.overlap) ## confirm correct # of ROHverlap segments were removed (TRUE)
# 
# write.csv(plink.overlap, 'corrected_plink_overlap.csv', row.names = FALSE)
plink.overlap <- read.csv('corrected_plink_overlap.csv')

##### >>> plot some overlap patterns to confirm correct calculations #####
## plot for weird P-O pairs with low genetic relatedness but high pedigree probability of P-O relationships
## for 65 contigs with length >= 1e7 (10 Mb)
# pdf('../krat_genetics_scripts/figures_output/parent_offspring_mismatch_4952_4901_4962.pdf', width = 10, height = 6)
# for(c in unique(contigs[contigs$length >= 1e7, 'c.index'])){
#   plink.p <- plink.roh[plink.roh$id == 4952 & plink.roh$c.index == c,]
#   plink.o1 <- plink.roh[plink.roh$id == 4901 & plink.roh$c.index == c,]
#   plink.o2 <- plink.roh[plink.roh$id == 4962 & plink.roh$c.index == c,]
# 
#   overlap.1 <- plink.overlap[which(plink.overlap$id2 == 4952 & plink.overlap$id1 == 4901 &
#                              plink.overlap$c.index == c),]
# 
#   overlap.2 <- plink.overlap[which(plink.overlap$id1 == 4952 & plink.overlap$id2 == 4962 &
#                                      plink.overlap$c.index == c),]
# 
#   plot(0, 0, xlim = c(1, contigs[contigs$c.index == c, 'length']), ylim = c(0.6, 5.4),
#        main = paste0(c), xlab = '', ylab = '', yaxt = 'n')
#     axis(2, at = c(1,2,3,4,5), labels = c('4962\nOverlap','4901\nOverlap','4962','4901','4952'))
#     for(n in 1:nrow(plink.p)){
#       polygon(x = c(plink.p$start[n], plink.p$end[n], plink.p$end[n], plink.p$start[n]),
#               y = c(4.6, 4.6, 5.4, 5.4), border = NA, col = 'darkorchid4')
#     }
#     for(n in 1:nrow(plink.o1)){
#       polygon(x = c(plink.o1$start[n], plink.o1$end[n], plink.o1$end[n], plink.o1$start[n]),
#               y = c(3.6, 3.6, 4.4, 4.4), border = NA, col = 'dodgerblue3')
#     }
#     for(n in 1:nrow(plink.o2)){
#       polygon(x = c(plink.o2$start[n], plink.o2$end[n], plink.o2$end[n], plink.o2$start[n]),
#               y = c(2.6, 2.6, 3.4, 3.4), border = NA, col = 'dodgerblue3')
#     }
#     if(nrow(overlap.1) > 0){
#       for(n in 1:nrow(overlap.1)){
#         polygon(x = c(overlap.1$start[n], overlap.1$end[n], overlap.1$end[n], overlap.1$start[n]),
#                 y = c(1.6, 1.6, 2.4, 2.4), border = NA, col = 'orange3')
#       }
#     }
#     if(nrow(overlap.2) > 0){
#       for(n in 1:nrow(overlap.2)){
#         polygon(x = c(overlap.2$start[n], overlap.2$end[n], overlap.2$end[n], overlap.2$start[n]),
#                 y = c(0.6, 0.6, 1.4, 1.4), border = NA, col = 'orange3')
#       }
#     }
# }
# dev.off()

## plot for full-sib pair whose parents were also full-sibs
## for 65 contigs with length >= 1e7 (10 Mb)
# pdf('../krat_genetics_scripts/figures_output/fullsibs_w_fullsib_parents_4195_4459.pdf', width = 10, height = 4)
# for(c in unique(contigs[contigs$length >= 1e7, 'c.index'])){
#   plink.1 <- plink.roh[plink.roh$id == 4195 & plink.roh$c.index == c,]
#   plink.2 <- plink.roh[plink.roh$id == 4459 & plink.roh$c.index == c,]
#   
#   overlap.1 <- plink.overlap[which(plink.overlap$id2 %in% c(4195, 4459) & plink.overlap$id1 %in% c(4195, 4459) &
#                                      plink.overlap$c.index == c),]
#   
#   plot(0, 0, xlim = c(1, contigs[contigs$c.index == c, 'length']), ylim = c(0.6, 3.4),
#        main = paste0(c), xlab = '', ylab = '', yaxt = 'n')
#     axis(2, at = c(1,2,3), labels = c('Overlap','4459','4195'))
#     if(nrow(plink.1) > 0){
#       for(n in 1:nrow(plink.1)){
#         polygon(x = c(plink.1$start[n], plink.1$end[n], plink.1$end[n], plink.1$start[n]),
#                 y = c(2.6, 2.6, 3.4, 3.4), border = NA, col = 'dodgerblue3')
#       }
#     }
#     if(nrow(plink.2) > 0){
#       for(n in 1:nrow(plink.2)){
#         polygon(x = c(plink.2$start[n], plink.2$end[n], plink.2$end[n], plink.2$start[n]),
#                 y = c(1.6, 1.6, 2.4, 2.4), border = NA, col = 'darkorchid4')
#       }
#     }  
#     if(nrow(overlap.1) > 0){
#       for(n in 1:nrow(overlap.1)){
#         polygon(x = c(overlap.1$start[n], overlap.1$end[n], overlap.1$end[n], overlap.1$start[n]),
#                 y = c(0.6, 0.6, 1.4, 1.4), border = NA, col = 'orange3')
#       }
#     }
# }
# dev.off()

##### Identify regions where ROHs never occur and where ROHs are very common (in these samples) #####
## calculate cumulative contig lengths
contigs$cumsum <- cumsum(as.numeric(contigs$length))

# OUT <- NULL
# # for(c in unique(contigs$c.index)){
#   print(c)
#   sub.rohs <- plink.roh[plink.roh$c.index == c,]
#   locs <- rep(0, contigs[contigs$c.index == c, 'length'])
#   if(nrow(sub.rohs) > 0){
#     for(r in 1:nrow(sub.rohs)){
#       locs[c(sub.rohs$start[r]:sub.rohs$end[r])] <- locs[c(sub.rohs$start[r]:sub.rohs$end[r])] + 1
#     }
#     no.rohs <- which(locs == 0)
#     for(l in 1:length(no.rohs)){
#       if(l == 1){
#         s <- no.rohs[1]
#       } else if(l == length(no.rohs)){
#         e <- no.rohs[l]
#         save <- c(c, s, e)
#         OUT <- rbind(OUT, save)
#       } else if(no.rohs[l+1] == no.rohs[l]+1){
#         next
#       } else if(no.rohs[l+1] > no.rohs[l]+1){
#         e <- no.rohs[l]
#         save <- c(c, s, e)
#         OUT <- rbind(OUT, save)
#         s <- no.rohs[l+1]
#       }
#     }
#   }
# }
# no.roh.regions <- as.data.frame(OUT)
# colnames(no.roh.regions) <- c('c.index','start','end')
# write.csv(no.roh.regions, 'contig_regions_no_rohs.csv', row.names = FALSE)
no.roh.regions <- read.csv('contig_regions_no_rohs.csv')

# OUT <- NULL
# for(c in unique(contigs$c.index)){
#   print(c)
#   sub.rohs <- plink.roh[plink.roh$c.index == c,]
#   locs <- rep(0, contigs[contigs$c.index == c, 'length'])
#   if(nrow(sub.rohs) > 0){
#     for(r in 1:nrow(sub.rohs)){
#       locs[c(sub.rohs$start[r]:sub.rohs$end[r])] <- locs[c(sub.rohs$start[r]:sub.rohs$end[r])] + 1
#     }
#     hi.rohs <- which(locs >= 24) ## where ROHs occur in at least 50% of samples
#     if(length(hi.rohs) > 0){
#       for(l in 1:length(hi.rohs)){
#         if(l == 1){
#           s <- hi.rohs[1]
#         } else if(l == length(hi.rohs)){
#           e <- hi.rohs[l]
#           save <- c(c, s, e)
#           OUT <- rbind(OUT, save)
#         } else if(hi.rohs[l+1] == hi.rohs[l]+1){
#           next
#         } else if(hi.rohs[l+1] > hi.rohs[l]+1){
#           e <- hi.rohs[l]
#           save <- c(c, s, e)
#           OUT <- rbind(OUT, save)
#           s <- hi.rohs[l+1]
#         }
#       }
#     }
#   }
# }
# hi.roh.regions <- as.data.frame(OUT)
# colnames(hi.roh.regions) <- c('c.index','start','end')
# write.csv(hi.roh.regions, 'contig_regions_high_rohs.csv', row.names = FALSE)
hi.roh.regions <- read.csv('contig_regions_high_rohs.csv')

library(rtracklayer)
## read in GFF file
gff <- readGFF('GCF_019054845.1_ASM1905484v1_genomic.gff')
gff <- gff[gff$type == 'gene' & gff$gene_biotype == 'protein_coding',]

## read in contig name key
c.key <- read.table('../assembly_seq_quotes_sample_notes_results/assembly_data/GCF_019054845.1_ASM1905484v1_assembly_report.txt')
c.key <- c.key[c.key$Assigned.Molecule.Location.Type %in% contigs$contig, c(4,6)]
colnames(c.key) <- c('contig','refseq.contig')
contigs <- merge(contigs, c.key, by = 'contig')
contigs <- contigs[order(contigs$length, decreasing = TRUE),]
## add to no.roh.regions and hi.roh.regions
no.roh.regions <- merge(no.roh.regions, contigs[,c(3,5)], by = 'c.index')
hi.roh.regions <- merge(hi.roh.regions, contigs[,c(3,5)], by = 'c.index')


# ### find overlap between no-ROH regions and genes in annotation
# GENE.OVERLAP <- NULL
# for(r in 1:nrow(no.roh.regions)){              ## for each ROH-less region,
#   print(r/nrow(no.roh.regions))
#   s <- no.roh.regions$start[r]
#   e <- no.roh.regions$end[r]
#   sub.gff <- gff[gff$seqid == no.roh.regions$refseq.contig[r],]
#   if(nrow(sub.gff) > 0){
#     ## genes beginning outside of focal ROH region, ending inside
#     if(nrow(sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]) > 0){
#       temp <- sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(no.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#     ## genes beginning inside of a focal ROH region, ending outside
#     if(nrow(sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]) > 0){
#       temp <- sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(no.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#     ## genes completely covering a focal ROH region
#     if(nrow(sub.gff[sub.gff$start < s & sub.gff$end > e,]) > 0){
#       temp <- sub.gff[sub.gff$start < s & sub.gff$end > e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(no.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#     ## genes completely within a focal ROH region
#     if(nrow(sub.gff[sub.gff$start >= s & sub.gff$end <= e,]) > 0){
#       temp <- sub.gff[sub.gff$start >= s & sub.gff$end <= e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(no.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#   }
# }
# no.roh.gene.overlap <- as.data.frame(GENE.OVERLAP)
# nrow(no.roh.gene.overlap)
# 
# ### find overlap between no-ROH regions and genes in annotation
# GENE.OVERLAP <- NULL
# for(r in 1:nrow(hi.roh.regions)){              ## for each ROH-less region,
#   print(r/nrow(hi.roh.regions))
#   s <- hi.roh.regions$start[r]
#   e <- hi.roh.regions$end[r]
#   sub.gff <- gff[gff$seqid == hi.roh.regions$refseq.contig[r],]
#   if(nrow(sub.gff) > 0){
#     ## genes beginning outside of focal ROH region, ending inside
#     if(nrow(sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]) > 0){
#       temp <- sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(hi.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#     ## genes beginning inside of a focal ROH region, ending outside
#     if(nrow(sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]) > 0){
#       temp <- sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(hi.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#     ## genes completely covering a focal ROH region
#     if(nrow(sub.gff[sub.gff$start < s & sub.gff$end > e,]) > 0){
#       temp <- sub.gff[sub.gff$start < s & sub.gff$end > e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(hi.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#     ## genes completely within a focal ROH region
#     if(nrow(sub.gff[sub.gff$start >= s & sub.gff$end <= e,]) > 0){
#       temp <- sub.gff[sub.gff$start >= s & sub.gff$end <= e,]
#       for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#         save <- temp[t, c('Name', 'seqid', 'start', 'end', 'gene_biotype')]   ## save the gene info and the ROH region info
#         save <- cbind(hi.roh.regions[r,], save)
#         colnames(save)[c(2,3,5,7,8)] <- c('region.start','region.end','symbol','gene.start','gene.end')
#         GENE.OVERLAP <- rbind(GENE.OVERLAP, save)
#       }
#     }
#   }
# }
# hi.roh.gene.overlap <- as.data.frame(GENE.OVERLAP)
# nrow(hi.roh.gene.overlap)


#### this bit doesn't work for some reason - not super important rn
# ## add descriptive names to gene info
# gene.desc <- read.table('../assembly_seq_quotes_sample_notes_results/assembly_data/GCF_019054845.1_ASM1905484v1_feature_table.txt', 
#                         header = TRUE, sep = '\t')
# gene.desc <- gene.desc[,c('symbol','name')]
# gene.desc <- gene.desc[which(gene.desc$name != '' & gene.desc$symbol != ''),]
# gene.desc <- gene.desc[!duplicated(gene.desc$symbol),]
# 
# no.roh.gene.overlap <- merge(no.roh.gene.overlap, gene.desc, by = 'symbol', all.x = TRUE)
# hi.roh.gene.overlap <- merge(hi.roh.gene.overlap, gene.desc, by = 'symbol', all.x = TRUE)

# write.csv(no.roh.gene.overlap, 'genes_in_no_roh_regions.csv', row.names = FALSE)
# write.csv(hi.roh.gene.overlap, 'genes_in_high_roh_regions.csv', row.names = FALSE)
no.roh.gene.overlap <- read.csv('genes_in_no_roh_regions.csv')
hi.roh.gene.overlap <- read.csv('genes_in_high_roh_regions.csv')

##### Using random sampling, determine whether protein-coding genes are over-/underrepresented in no/high ROH regions #####
no.roh.lens <- no.roh.regions$end - no.roh.regions$start + 1
no.roh.lens <- no.roh.lens[order(no.roh.lens, decreasing = TRUE)]
no.roh.lens <- no.roh.lens[no.roh.lens >= 50000] ### keep only regions at least 50 kb in length
hi.roh.lens <- hi.roh.regions$end - hi.roh.regions$start + 1
hi.roh.lens <- hi.roh.lens[order(hi.roh.lens, decreasing = TRUE)]
hi.roh.lens <- hi.roh.lens[hi.roh.lens >= 50000] ### keep only regions at least 50 kb in length

### write some data to files so the below can be run in a different script ###
write.table(hi.roh.lens, 'high_roh_region_lengths.txt', sep = '\t', row.names = FALSE, col.names = FALSE)
write.table(no.roh.lens, 'no_roh_region_lengths.txt', sep = '\t', row.names = FALSE, col.names = FALSE)
write.table(contigs, 'contigs_easley.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(gff[,c('seqid','start','end','gene')], 'gff_easley.txt', sep = '\t', row.names = FALSE, quote = FALSE)
### >> next script for this = EASLEY_02_lo_and_hi_roh_region_ID.R
