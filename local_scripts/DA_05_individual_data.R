setwd('/Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/')

library(scales)
library(ghibli)
library(languageR)
library(lme4)
library(piecewiseSEM)

##### Read in data #####
## NGSrelate2 results
ngs.LD <- read.table('ngsrelate_krat_allsamps_LDpruned.res', header = TRUE)
samp.order <- read.table('../vcf_sample_order.txt')
samp.order[,2] <- c(0:47)
## replace indices with sample names in NGS output
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'a', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id1'
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'b', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id2'
ngs.LD <- ngs.LD[order(ngs.LD$a, ngs.LD$b),]

## f(ROH) data
froh <- read.csv('../individual_frohs.csv')
## get rid of bcftools results
froh <- froh[,-(grep('pl.', colnames(froh), fixed = TRUE))]
froh <- froh[,-(grep('gt.', colnames(froh), fixed = TRUE))]

## ROH coordinates
plink.roh <- read.table('../plink_final_settings/krat_plink_roh_phwh_2_phwm_30_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh <- plink.roh[,c(1, 4, 7, 8, 10)]
colnames(plink.roh) <- c('id','contig','start','end','n.snps')
plink.roh$length <- plink.roh$end - plink.roh$start + 1
plink.roh <- plink.roh[plink.roh$length >= 100000,]

## seq stats - below calculated using bcftools stats on krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz
vcf.stats <- read.csv('../../qc_stuff/bcftools_stats/just_sample_stats.txt', check.names = FALSE)
vcf.stats <- vcf.stats[,c('sample','nHets','average depth','nMissing')]
colnames(vcf.stats) <- c('id','n.hets','avg.depth','n.missing')
n.sites <- 11796181 ## # of SNPs in this VCF
vcf.stats$prop.miss <- vcf.stats$n.missing/n.sites
vcf.stats$norm.prop.het <- vcf.stats$n.hets/(n.sites - vcf.stats$n.missing)

## contig information
contigs <- read.table('../contigs_easley.txt', header = TRUE)

## read in pedigree data
ped.dat <- read.csv('../summary_allstats2_withgen_zeroes.csv')
ped.dat$age.at.death <- ped.dat$deathyear - ped.dat$birthyear
all.ped.dat <- ped.dat
ped.dat <- ped.dat[ped.dat$id %in% samp.order$V1,]

## PLINK individual heterozygosity
het <- read.table('LDpruned_plink.het', header = TRUE)
het <- het[,c(1,6)]
colnames(het) <- c('id','F')


##### Look for some correlates of individual fitness #####
## demographic predictor variables
plot(ped.dat$birthyear, ped.dat$off, pch = 19, col = alpha('springgreen4', 0.6))
plot(ped.dat$birthyear, ped.dat$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
plot(ped.dat$birthyear, ped.dat$age.at.death, pch = 19, col = alpha('springgreen4', 0.6))

## somewhat strong relationship between age at death and surviving offspring
plot(ped.dat$age.at.death, ped.dat$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
  abline(lm(ped.dat$off_survive ~ ped.dat$age.at.death))
  summary(lm(ped.dat$off_survive ~ ped.dat$age.at.death)) ## all probably driven by the individuals born in 2005/2006, technical and not biological?
plot(all.ped.dat$age.at.death, all.ped.dat$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
  abline(lm(all.ped.dat$off_survive ~ all.ped.dat$age.at.death))
  summary(lm(all.ped.dat$off_survive ~ all.ped.dat$age.at.death))
temp <- all.ped.dat[all.ped.dat$birthyear <= 2004,]
plot(temp$age.at.death, temp$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
  abline(lm(temp$off_survive ~ temp$age.at.death))
  summary(lm(temp$off_survive ~ temp$age.at.death)) ## relationships look similar with and without those last few years, tho. so maybe useful.

plot(ped.dat$age.at.death, ped.dat$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
  abline(lm(ped.dat$off_survive ~ ped.dat$age.at.death))
  summary(lm(ped.dat$off_survive ~ ped.dat$age.at.death))
plot(ped.dat$age.at.death, ped.dat$off, pch = 19, col = alpha('springgreen4', 0.6))
  abline(lm(ped.dat$off ~ ped.dat$age.at.death))
  summary(lm(ped.dat$off ~ ped.dat$age.at.death))

## consistent with Janna's 2019 paper, looks like offspring produced is driven by lifespan, which is driven by inbreeding (see Fig. 3)
## hold up with genomic inbreeding #s rather than pedigree #s? try doing it right.
froh <- merge(froh, ped.dat[,c('id','birthyear')], by = 'id')
froh <- merge(froh, het, by = 'id')
xlabs <- c('f(ROH)','bin 1 f(ROH)','bin 2 f(ROH)','bin 3 f(ROH)','bin 4 f(ROH)','bin 5 f(ROH)')

## # offspring ~ f(ROH) measures, birthyear as random
x <- 1
pdf('../../krat_genetics_scripts/figures_output/fROH_v_numoffspring.pdf', width = 5, height = 5)
for(c in c(2,8:12)){
  res <- glmer(froh$off ~ froh[,c] + (1|froh$birthyear), family = 'poisson')
  # summary(res)
  rsquared(res)  ## marginal = fixed effects only, conditional = fixed + random effects
  plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(froh[,c]), max(froh[,c])), 
               ylim = c(min(froh$off), max(froh$off)), xlabel = xlabs[x], ylabel = 'Number of offspring')
    points(froh[,c], froh$off, pch = 19, col = alpha('springgreen4', 0.6))
    legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                  paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
           bty = 'n', inset = 0.02)
  x <- x+1
}
dev.off()

## # surviving offspring ~ f(ROH) measures, birthyear as random
x <- 1
pdf('../../krat_genetics_scripts/figures_output/fROH_v_numsurvoffspring.pdf', width = 5, height = 5)
for(c in c(2,8:12)){
  res <- glmer(froh$off_survive ~ froh[,c] + (1|froh$birthyear), family = 'poisson')
  summary(res)
  plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(froh[,c]), max(froh[,c])), 
               ylim = c(min(froh$off_survive), max(froh$off_survive)), xlabel = xlabs[x], ylabel = 'Number of surviving offspring')
  points(froh[,c], froh$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
  x <- x+1
}
dev.off()

## # offspring ~ PLINK F, birthyear as random
pdf('../../krat_genetics_scripts/figures_output/F_v_numoffspring.pdf', width = 5, height = 5)
res <- glmer(froh$off ~ froh$F + (1|froh$birthyear), family = 'poisson')
# summary(res)
rsquared(res)  ## marginal = fixed effects only, conditional = fixed + random effects
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(froh$F), max(froh$F)), 
             ylim = c(min(froh$off), max(froh$off)), xlabel = xlabs[x], ylabel = 'Number of offspring')
points(froh$F, froh$off, pch = 19, col = alpha('springgreen4', 0.6))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
dev.off()

## # surviving offspring ~ PLINK F, birthyear as random
pdf('../../krat_genetics_scripts/figures_output/F_v_numsurvoffspring.pdf', width = 5, height = 5)
res <- glmer(froh$off_survive ~ froh$F + (1|froh$birthyear), family = 'poisson')
# summary(res)
rsquared(res)  ## marginal = fixed effects only, conditional = fixed + random effects
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(froh$F), max(froh$F)), 
             ylim = c(min(froh$off_survive), max(froh$off_survive)), xlabel = xlabs[x], ylabel = 'Number of surviving offspring')
  points(froh$F, froh$off_survive, pch = 19, col = alpha('springgreen4', 0.6))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
dev.off()
