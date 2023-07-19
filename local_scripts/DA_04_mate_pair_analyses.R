setwd('/Users/Avril/Documents/krat_genetics/data/')

library(vcfR)
library(scales)
library(factoextra)
library(adegenet)
library(ghibli)

##### Read in data #####
## sample ROH data (PLINK)
plink.roh <- read.table('plink_final_settings/krat_plink_roh_phwh_2_phwm_30_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom', header = TRUE)
plink.roh <- plink.roh[,c(1, 4, 7, 8, 10)]
colnames(plink.roh) <- c('id','contig','start','end','n.snps')
plink.roh$length <- plink.roh$end - plink.roh$start + 1
plink.roh <- plink.roh[plink.roh$length >= 100000,]

## sample information
ped.dat <- read.csv('../preseq_sample_information/summary_allstats2_withgen_zeroes.csv')
relates <- read.csv('../data/pairwise_samp_relationships.csv')
## purposely sequenced sibling pairs
sib.pairs <- read.csv('/Users/Avril/Documents/krat_genetics/pair_specific_demo_calcs/sibling_pairs.csv')

## below calculated using bcftools stats on krat_final_final_allfiltcontigs_all_samps.recode.vcf.gz (no LD pruning)
vcf.stats <- read.csv('../qc_stuff/bcftools_stats/just_sample_stats.txt', check.names = FALSE)
vcf.stats <- vcf.stats[,c('sample','nHets','average depth','nMissing')]
colnames(vcf.stats) <- c('id','n.hets','avg.depth','n.missing')
n.sites <- 11796181 ## # of SNPs in this VCF
vcf.stats$prop.miss <- vcf.stats$n.missing/n.sites
vcf.stats$norm.prop.het <- vcf.stats$n.hets/(n.sites - vcf.stats$n.missing)

## contig information
contigs <- read.csv('contig_lengths.csv')
contigs <- contigs[order(-contigs$length),]
contigs$c.index <- c(1:nrow(contigs)) ## create contig index, where c.index = 1 == longest contig
tot.len <- sum(contigs$length) ## 2,571,239,112 bp ## length of filtered contigs
## add this contig information to ROH data
plink.roh <- merge(plink.roh, contigs[,c(1,3)], by = 'contig')

## set f(ROH) formatting for plots
froh.lab <- substitute(paste(italic('F')[ROH]))

## read in relatedness info
## ROHverlap (based on non-pruned SNPs)
rohverlap <- read.csv('corrected_plink_overlap.csv')  ## corrected through removal of overlapping segments with homozygous mismatches between individuals
## NGSrelate2 relatedness results (based on LD-pruned SNPs)
ngs.LD <- read.table('./LD_pruned_data_and_results/ngsrelate_krat_allsamps_LDpruned.res', header = TRUE)
samp.order <- read.table('vcf_sample_order.txt')
samp.order[,2] <- c(0:47)
## replace indices with sample names in NGS output
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'a', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id1'
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'b', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id2'
ngs.LD <- ngs.LD[order(ngs.LD$a, ngs.LD$b),]

## read in VCF sample order and genotype information (collected using bcftools query)
samp.order <- read.table('vcf_sample_order.txt')
samp.order <- unlist(samp.order)

##### For all pairs of mates seq'ed, get genetic relatedness and fitness info #####
mates <- relates[relates$mates == 1,]
OUT <- NULL
for(r in 1:nrow(mates)){
  king <- ngs.LD[which((ngs.LD$id1 == mates$id1[r] & ngs.LD$id2 == mates$id2[r]) |
                 (ngs.LD$id1 == mates$id2[r] & ngs.LD$id2 == mates$id1[r])), 'KING']
  rab <- ngs.LD[which((ngs.LD$id1 == mates$id1[r] & ngs.LD$id2 == mates$id2[r]) |
                        (ngs.LD$id1 == mates$id2[r] & ngs.LD$id2 == mates$id1[r])), 'rab']
  rohv <- sum(rohverlap[which((rohverlap$id1 == mates$id1[r] & rohverlap$id2 == mates$id2[r]) |
                        (rohverlap$id1 == mates$id2[r] & rohverlap$id2 == mates$id1[r])), 'len.overlap'])/tot.len
  n.off <- nrow(ped.dat[which((ped.dat$momid == mates$id1[r] & ped.dat$dadid == mates$id2[r]) |
                  (ped.dat$momid == mates$id2[r] & ped.dat$dadid == mates$id1[r])),])
  n.surv <- nrow(ped.dat[which((ped.dat$momid == mates$id1[r] & ped.dat$dadid == mates$id2[r] & ped.dat$deathyear > ped.dat$birthyear) |
                                (ped.dat$momid == mates$id2[r] & ped.dat$dadid == mates$id1[r] & ped.dat$deathyear > ped.dat$birthyear)),])
  n.g.off <- sum(ped.dat[which((ped.dat$momid == mates$id1[r] & ped.dat$dadid == mates$id2[r]) |
                                  (ped.dat$momid == mates$id2[r] & ped.dat$dadid == mates$id1[r])), 'off'])
  n.g.surv <- sum(ped.dat[which((ped.dat$momid == mates$id1[r] & ped.dat$dadid == mates$id2[r]) |
                                 (ped.dat$momid == mates$id2[r] & ped.dat$dadid == mates$id1[r])), 'off_survive'])
  save <- c(mates$id1[r], mates$id2[r], king, rab, rohv, n.off, n.surv, n.g.off, n.g.surv)
  OUT <- rbind(OUT, save)
}
mates <- as.data.frame(OUT)
colnames(mates) <- c('id1','id2','king','rab','rohverlap','n.off','n.surv','n.g.off','n.g.surv')
mates$log.rohverlap <- log(mates$rohverlap)

pdf('../krat_genetics_scripts/figures_output/relatedness_v_fitness_all_pairs.pdf', width = 5, height = 5)
for(r in c(3:5)){
  for(f in c(6:7)){
    plot(mates[,r], mates[,f], xlab = colnames(mates)[r], ylab = colnames(mates)[f],
         pch = 19, col = alpha('darkorchid3', 0.5))
      mod <- lm(mates[,f] ~ mates[,r])
      abline(mod, lty = 1, col = 'grey')
      if(summary(mod)$coefficients[2,4] < 0.10){
        print(paste0('r ',r,' - f ',f))
      }
  }
}
dev.off()


##### For all full-sibling pairs purposefully seq'ed.... doooooo... something #####
# pdf('../krat_genetics_scripts/figures_output/sib_pairs_relatedness_v_fitness.pdf', width = 12, height = 12)
# par(mfrow = c(4,4), mar = c(4.1,4.1,2.1,2.1))
# for(p in 1:nrow(sib.pairs)){
#   s1 <- sib.pairs$sib.1[p]
#   s2 <- sib.pairs$sib.2[p]
#   s1.mates <- mates[which(mates$id1 == s1 | mates$id2 == s1),]
#   s2.mates <- mates[which(mates$id1 == s2 | mates$id2 == s2),]
#   for(r in c(3:5, 10)){
#     for(f in c(6:9)){
#       plot(mates[,r], mates[,f], xlab = colnames(mates)[r], ylab = colnames(mates)[f],
#            pch = 19, col = alpha('darkgrey', 0.5))
#         abline(lm(mates[,f] ~ mates[,r]), lty = 2)
#         points(s1.mates[,r], s1.mates[,f], pch = 17, col = 'darkorchid3')
#         points(s2.mates[,r], s2.mates[,f], pch = 15, col = 'darkturquoise')
#         if(r == 3 & f == 6){
#           legend('topright', legend = paste0(s1,' & ',s2), bty = 'n')
#         }
#     }
#   }
# }
# dev.off()

##### !!! actually, I don't think we can use any grandoffspring measures, because so many of the seq'd individuals were born in 2005 #####
##### !!! .... look into this. might be able to count # offspring produced, but not # surviving for parents born in 2005 #####
# mates <- mates[,-c(8,9)]
# pdf('../krat_genetics_scripts/figures_output/relatedness_v_fitness_all_pairs.pdf', width = 5, height = 5)
# for(r in c(3:5)){
#   for(f in c(6,7)){
#     plot(mates[,r], mates[,f], pch = 19, xlab = colnames(mates)[r], ylab = colnames(mates)[f],
#          col = alpha('darkorchid3', 0.5))
#   }
# }
# dev.off()


##### Identify ROH (no- and hi-regions) shared by pairs with and without surviving offspring #####
mates[mates$n.surv > 0, 'surv.or.no'] <- 1
mates[mates$n.surv == 0, 'surv.or.no'] <- 0
mates$pair.id <- c(1:nrow(mates))

# ## maybe a sliding window approach is best here? 
# ## for each window, assess how many mate pairs have at least 1 bp in a ROHverlap region.
# OUT1 <- NULL        ## for pairs with surviving offspring
# OUT2 <- NULL        ## for pairs without surviving offspring
# w.size <- 100000    ## sliding window size
# s.size <- 50000     ## sliding window step size
# 
# ## split rohverlap data according to the offspring survival categories
# for(r in 1:nrow(mates)){
#   sub.rohverlap <- rohverlap[which((rohverlap$id1 == mates$id1[r] & rohverlap$id2 == mates$id2[r]) |
#                                      rohverlap$id1 == mates$id2[r] & rohverlap$id2 == mates$id1[r]),]
#   sub.rohverlap$pair.id <- r
#   if(mates$surv.or.no[r] == 0){
#     OUT2 <- rbind(OUT2, sub.rohverlap)
#   }
#   if(mates$surv.or.no[r] == 1){
#     OUT1 <- rbind(OUT1, sub.rohverlap)
#   }
# }
# 
# ## calculate depth of ROHverlap for sliding windows for the two categories of mate pairs
# # c <- 300 ### fine value to test
# SURV.OFF.PAIRS <- NULL
# NO.SURV.OFF.PAIRS <- NULL
# WIND.DAT <- NULL
# wind.id <- 1
# 
# for(c in unique(contigs$c.index)){
#   c.len <- contigs[contigs$c.index == c, 'length']
#   sub1 <- OUT1[OUT1$c.index == c,]
#   sub2 <- OUT2[OUT2$c.index == c,]
#   s <- 1
#   e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#   
#   while(e <= c.len){
#     wind.dat <- c(c, wind.id, s, e)
#     print(paste0(c,' - ',e/c.len))
#     WIND.DAT <- rbind(WIND.DAT, wind.dat)
#     ## for pairs with surviving offspring
#     ## comp ROHs beginning outside of window, ending inside
#     if(nrow(sub1[sub1$start < s & sub1$end >= s & sub1$end <= e,]) > 0){
#       temp <- sub1[sub1$start < s & sub1$end >= s & sub1$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         SURV.OFF.PAIRS <- rbind(SURV.OFF.PAIRS, save)
#       }
#     }
#     ## comp ROHs beginning inside of a focal window, ending outside
#     if(nrow(sub1[sub1$start >= s & sub1$start <= e & sub1$end > e,]) > 0){
#       temp <- sub1[sub1$start >= s & sub1$start <= e & sub1$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         SURV.OFF.PAIRS <- rbind(SURV.OFF.PAIRS, save)
#       }
#     }
#     ## comp ROHs completely covering a focal window
#     if(nrow(sub1[sub1$start < s & sub1$end > e,]) > 0){
#       temp <- sub1[sub1$start < s & sub1$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         SURV.OFF.PAIRS <- rbind(SURV.OFF.PAIRS, save)
#       }
#     }
#     ## comp ROHs completely within a focal window
#     if(nrow(sub1[sub1$start >= s & sub1$end <= e,]) > 0){
#       temp <- sub1[sub1$start >= s & sub1$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         SURV.OFF.PAIRS <- rbind(SURV.OFF.PAIRS, save)
#       }
#     }
#     
#     ## for pairs without surviving offspring
#     ## comp ROHs beginning outside of window, ending inside
#     if(nrow(sub2[sub2$start < s & sub2$end >= s & sub2$end <= e,]) > 0){
#       temp <- sub2[sub2$start < s & sub2$end >= s & sub2$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         NO.SURV.OFF.PAIRS <- rbind(NO.SURV.OFF.PAIRS, save)
#       }
#     }
#     ## comp ROHs beginning inside of a focal window, ending outside
#     if(nrow(sub2[sub2$start >= s & sub2$start <= e & sub2$end > e,]) > 0){
#       temp <- sub2[sub2$start >= s & sub2$start <= e & sub2$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         NO.SURV.OFF.PAIRS <- rbind(NO.SURV.OFF.PAIRS, save)
#       }
#     }
#     ## comp ROHs completely covering a focal window
#     if(nrow(sub2[sub2$start < s & sub2$end > e,]) > 0){
#       temp <- sub2[sub2$start < s & sub2$end > e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         NO.SURV.OFF.PAIRS <- rbind(NO.SURV.OFF.PAIRS, save)
#       }
#     }
#     ## comp ROHs completely within a focal window
#     if(nrow(sub2[sub2$start >= s & sub2$end <= e,]) > 0){
#       temp <- sub2[sub2$start >= s & sub2$end <= e,]
#       for(t in 1:nrow(temp)){                                            ## for each overlapping comp ROH,
#         save <- c(wind.id, temp$pair.id[t])
#         NO.SURV.OFF.PAIRS <- rbind(NO.SURV.OFF.PAIRS, save)
#       }
#     }
#     wind.id <- wind.id + 1
#     if(e != c.len){
#       s <- s + s.size
#       e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#     } else{
#       e <- e+1
#     }
#   }
# }
# 
# ## summarize numbers of pairs with ROHs for each window
# wind.dat <- as.data.frame(WIND.DAT)
# colnames(wind.dat) <- c('c.index','wind.id','wind.s','wind.e')
# 
# OUT1 <- NULL
# OUT2 <- NULL
# for(w in 1:max(SURV.OFF.PAIRS[,1])){
#   save <- c(w, length(unique(SURV.OFF.PAIRS[SURV.OFF.PAIRS[,1] == w, 2]))/16)
#   OUT1 <- rbind(OUT1, save)
# }
# for(w in 1:max(NO.SURV.OFF.PAIRS[,1])){
#   save <- c(w, length(unique(NO.SURV.OFF.PAIRS[NO.SURV.OFF.PAIRS[,1] == w, 2]))/24)
#   OUT2 <- rbind(OUT2, save)
# }
# hist(OUT1[,2])
# hist(OUT2[,2])
# 
# ### plot patterns in 1k-window sections, see if anything pops out
# x.lo <- 1
# x.hi <- x.lo + 1000 - 1
# pdf('../krat_genetics_scripts/figures_output/rohs_across_windows.pdf', width = 12, height = 4)
# while(x.hi < max(c(OUT1[,1], OUT2[,1]))){
#   plot(c(OUT1[OUT1[,1] %in% c(x.lo:x.hi), 1], OUT2[OUT2[,1] %in% c(x.lo:x.hi), 1]), 
#        c(OUT1[OUT1[,1] %in% c(x.lo:x.hi), 2], OUT2[OUT2[,1] %in% c(x.lo:x.hi), 2]),
#        main = paste0(x.lo,' - ',x.hi,' windows'), col = 'transparent', xlab = 'Window number', ylim = c(0,1),
#        ylab = 'Proportion of pairs with ROHverlap in window')
#     abline(h = 0.5, lty = 2)
#   
#     # points(OUT1[OUT1[,1] %in% c(x.lo:x.hi), 1], OUT1[OUT1[,1] %in% c(x.lo:x.hi), 2], pch = 19, cex = 0.5, col = 'darkorchid3')
#     # points(OUT2[OUT2[,1] %in% c(x.lo:x.hi), 1], OUT2[OUT2[,1] %in% c(x.lo:x.hi), 2], pch = 19, cex = 0.5, col = 'darkturquoise')
# 
#     lines(OUT1[OUT1[,1] %in% c(x.lo:x.hi), 1], OUT1[OUT1[,1] %in% c(x.lo:x.hi), 2], col = 'darkorchid3')
#     lines(OUT2[OUT2[,1] %in% c(x.lo:x.hi), 1], OUT2[OUT2[,1] %in% c(x.lo:x.hi), 2], col = 'darkturquoise')
#     
#     # polygon(OUT1[OUT1[,1] %in% c(x.lo:x.hi), 1], OUT1[OUT1[,1] %in% c(x.lo:x.hi), 2], col = 'darkorchid3', border = NA)
#     # polygon(OUT2[OUT2[,1] %in% c(x.lo:x.hi), 1], OUT2[OUT2[,1] %in% c(x.lo:x.hi), 2], col = 'darkturquoise', border = NA)    
#     
#     legend('topleft', legend = c('+ surv off (16 pairs)','- surv off (24 pairs)'), lty = 1, col = c('darkorchid3','darkturquoise'))
#     
#   x.lo <- x.lo + 1000
#   x.hi <- x.lo + 1000 - 1
# }
# dev.off()

##### Calculating per-locus, intra-pair genotype similarities --> PCA or something? #####
gts <- read.table('./LD_pruned_data_and_results/krat_LDpruned_extracted_genotypes_numericGTs.txt', sep = '\t')
gts <- gts[,-3]
colnames(gts) <- c('contig','pos',samp.order)
gts <- merge(gts, contigs[,c(1,3)], by = 'contig')
gts <- gts[,c(51,2:50)]

for(c in 1:ncol(gts)){
  gts[,c] <- as.numeric(gts[,c])
}  ## missing genotypes go from . --> NA
gts <- gts[order(gts$c.index, gts$pos),]
rownames(gts) <- c(1:nrow(gts))
gts.with.pos <- gts  ## save a version with contig/pos info retained
gts <- gts[,-c(1,2)]
gts <- as.matrix(gts)

## genotype PCA to see if anything weird falls out
gts.dat <- t(gts)
dim(gts.dat)
gts.dat <- gts.dat[, -(which(colSums(is.na(gts.dat)) > 0))] ## PCA doesn't allow for any NA values
dim(gts.dat)
res.pca <- prcomp(gts.dat, scale = FALSE)

## visualize
fviz_eig(res.pca) ## scree plot
# fviz_pca_ind(res.pca)

## colors points by individual birth year
years <- ped.dat[ped.dat$id %in% rownames(gts.dat), c('id','birthyear')]
years <- years[match(rownames(gts.dat), years$id),]
years$birthyear <- as.factor(years$birthyear)

cols <- sort(unique(years$birthyear))
cols <- as.data.frame(cbind(cols, ghibli_palettes$PonyoMedium))
colnames(cols) <- c('birthyear','colour')

fviz_pca_ind(res.pca, col.ind = years$birthyear, palette = cols$colour, 
             geom = c('point'), addEllipses = FALSE)

### for each pair of mates, calculate per-locus distances and summarize over windows
## run on Easley: EASLEY_01_intrapair_gendist_calcs.R; read in results

# dists <- read.table('LD_pruned_data_and_results/intrapair_gen_distances.txt', header = TRUE)
# colnames(dists) <- c('pair.id','surv.or.no','c.index','wind.id','start','end','n.snps','n.diffs')
# dists <- dists[-which(dists$pair.id == 'V1'),]
# for(c in 1:ncol(dists)){
#   dists[,c] <- as.numeric(dists[,c])
# }
# ## n.snps = number of SNPs with genotypes for both samples; n.diffs = number of allelic differences within that window
# dists$prop.diffs <- dists$n.diffs / (2*dists$n.snps)

## read in PCA-formatted data
dat.pca <- read.table('LD_pruned_data_and_results/pcaformatted_gendists_shared_windowIDs.txt', header = TRUE)
dat.pca <- t(dat.pca)
colnames(dat.pca) <- dat.pca[1,]
dat.pca <- dat.pca[-1,]
rownames(dat.pca) <- c(1:40)

dat.pca <- dat.pca[, -(which(colSums(is.na(dat.pca)) > 0))] ## PCA doesn't allow for any NA values

# ## run PCA
# res.pca <- prcomp(dat.pca, scale = FALSE)
# 
# ## visualize
# fviz_eig(res.pca) ## scree plot
# fviz_pca_ind(res.pca, col.ind = mates$surv.or.no)
# fviz_pca_ind(res.pca, col.ind = mates$n.off)
# fviz_pca_ind(res.pca, col.ind = mates$n.surv)
rm(res.pca) ## takes up a lot of memory

##### DAPC to explore influential windows #####
## DAPC links:
## >> https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

# n.pca <- 32 ## can change this
# n.da <- 1   ## always 1 because k = 2
# 
# ## actually, let's test a few different n.pca values
# pdf('../krat_genetics_scripts/figures_output/LD_pruned/dapc_different_nPCA_postprobs.pdf', width = 6, height = 5)
# for(c in c(30, 32, 34, 36, 38)){
#   n.pca <- c
#   dapc.res <- dapc(dat.pca, grp = mates$surv.or.no, center = TRUE, scale = TRUE, n.pca = n.pca, n.da = n.da)
#   ## n.pca setting really makes a difference here; increasing the value leads to increasing differentiation between the 2 groups
#   scatter(dapc.res) ## only 1 discriminant function when k = 2
#     legend('topleft', legend = paste0('n.pca = ',n.pca), bty = 'n')
#     
#   ## look at posteriors relative to intial fitness assignments
#   temp <- cbind(c(1:40), mates$surv.or.no)
#   
#   par(xpd = FALSE)
#   plot(c(1:40), dapc.res$posterior[,1], col = 'transparent', pch = 19, xlab = 'pair ID', ylab = 'membership probabilities',
#        main = paste0('n.pca = ',n.pca))
#     for(r in 1:nrow(dapc.res$posterior)){
#       lines(c(r, r), c(dapc.res$posterior[r, 1], dapc.res$posterior[r, 2]), col = 'grey')
#     }
#     points(c(1:40), dapc.res$posterior[,2], col = 'red', pch = 19) ## no surviving offspring
#     points(c(1:40), dapc.res$posterior[,1], col = 'blue', pch = 19) ## surviving offspring
#     abline(h = 0.5, lty = 2)
#     par(xpd = TRUE)
#     for(r in 1:nrow(temp)){
#       if(temp[r,2] == 0){
#         points(r, 1.1, pch = 19, col = 'blue')
#       } else{
#         points(r, 1.1, pch = 19, col = 'red')
#       }
#     } ## all match: predicting surv.or.no with discriminant function works
# }
# dev.off()
# 
# n.pca <- 36 ## this is the lowest value that results in no overlap between clusters on discriminant function 1
#             ## 32 also works pretty well -- all of this probably == overfitting though? 
# dapc.res <- dapc(dat.pca, grp = mates$surv.or.no, center = TRUE, scale = TRUE, n.pca = n.pca, n.da = n.da)
# ## check how assignment would work with randomized groups using same # of PCA axes retained;
# ## values closer to 1 = DAPC solution is strongly discriminating and stable;
# ## closer to 0 = weak discrimination or instability of the results
# n.pca <- 39
# dapc.res <- dapc(dat.pca, grp = mates$surv.or.no, center = TRUE, scale = TRUE, n.pca = n.pca, n.da = n.da)
# a.score(dapc.res, n.da = n.da, n.pca = n.pca, plot = TRUE) 
# optim.a.score(dapc.res) ## uh oh, looks like it doesn't work super well. best n.pca = 25, still not good.

n.pca <- 25
n.da <- 1
dapc.res <- dapc(dat.pca, grp = mates$surv.or.no, center = TRUE, scale = TRUE, n.pca = n.pca, n.da = n.da)

scatter(dapc.res) ## only 1 discriminant function when k = 2
legend('topleft', legend = paste0('n.pca = ',n.pca), bty = 'n')

## look at posteriors relative to intial fitness assignments
temp <- cbind(c(1:40), mates$surv.or.no)
par(xpd = FALSE)
plot(c(1:40), dapc.res$posterior[,1], col = 'transparent', pch = 19, xlab = 'pair ID', ylab = 'membership probabilities',
     main = paste0('n.pca = ',n.pca))
  for(r in 1:nrow(dapc.res$posterior)){
    lines(c(r, r), c(dapc.res$posterior[r, 1], dapc.res$posterior[r, 2]), col = 'grey')
  }
  points(c(1:40), dapc.res$posterior[,2], col = 'red', pch = 19) ## no surviving offspring
  points(c(1:40), dapc.res$posterior[,1], col = 'blue', pch = 19) ## surviving offspring
  abline(h = 0.5, lty = 2)
  par(xpd = TRUE)
  for(r in 1:nrow(temp)){
    if(temp[r,2] == 0){
      points(r, 1.1, pch = 19, col = 'blue')
    } else{
      points(r, 1.1, pch = 19, col = 'red')
    }
  }

### check windows with big contributions
head(dapc.res$var.contr) ## contains the contributions of each window along the only discriminant function
contr <- dapc.res$var.contr
contr <- contr[order(contr[,1], decreasing = FALSE),]
hist(contr)
## all contributions in increasing order (takes a long time to plot)
# plot(0, 0, ylim = c(0, max(contr)), xlim = c(0, length(contr)), col = 'transparent', xlab = '')
# for(c in 1:length(contr)){
#   lines(x = c(c, c), y = c(0, contr[c]))
# }

top.10 <- tail(contr, n = 10)[order(tail(contr, n = 10), decreasing = TRUE)]
top.10 <- dat.pca[, colnames(dat.pca) %in% names(top.10)]
mates <- cbind(mates, top.10)

## for each top window, plot distances by group
pdf('../krat_genetics_scripts/figures_output/LD_pruned/dapc_top10_window_distances.pdf', width = 6, height = 5)
for(c in 13:ncol(mates)){
  boxplot(mates[,c] ~ mates$surv.or.no, main = colnames(mates)[c], xlab = 'Surviving offspring or no', 
          ylab = 'Intrapair window distance')
}
dev.off()

## for each window, look at within-pair genotype matching;
## first need to get window coordinates
w.size <- 10000    ## window size
s.size <- 10000    ## window step size
# wind.id <- 1       ## starting window ID #
# WIND.DAT <- NULL
# 
# for(c in unique(contigs$c.index)){
#   print(c/663)
#   c.len <- contigs[contigs$c.index == c, 'length']
#   s <- 1
#   e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#   while(e <= c.len){
#     wind.dat <- c(wind.id, c, s, e)
#     WIND.DAT <- rbind(WIND.DAT, wind.dat)
#     wind.id <- wind.id + 1
#     if(e != c.len){
#       s <- s + s.size
#       e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#     } else{
#       e <- e+1
#     }
#   }
#     if(c == 1){
#       write.table(WIND.DAT, paste0('LD_pruned_data_and_results/window_coords_wsize_',w.size,'.txt'), quote = FALSE, row.names = FALSE)
#       WIND.DAT <- NULL
#     } else{
#       write.table(WIND.DAT, paste0('LD_pruned_data_and_results/window_coords_wsize_',w.size,'.txt'), quote = FALSE, row.names = FALSE,
#                   append = TRUE, col.names=!file.exists(paste0('LD_pruned_data_and_results/window_coords_wsize_',w.size,'.txt')))
#       WIND.DAT <- NULL
#     }
# }
wind.coords <- read.table('LD_pruned_data_and_results/window_coords_wsize_10000.txt', header = TRUE)
colnames(wind.coords) <- c('wind.id','c.index','start','end')

pdf('../krat_genetics_scripts/figures_output/LD_pruned/top_10_dapc_window_genotypes.pdf', width = 7, height = 6)
par(mfrow = c(2,1))
for(c in 13:ncol(mates)){
  wind.id <- colnames(mates)[c]
  c.index <- wind.coords[wind.coords$wind.id == wind.id, 'c.index']
  s <- wind.coords[wind.coords$wind.id == wind.id, 'start']
  e <- wind.coords[wind.coords$wind.id == wind.id, 'end']
  
  ## no surviving offspring
  mates.0 <- mates[mates$surv.or.no == 0,]
  plot(0, 0, xlim = c(s, e),  ylim = c(0, 2), col = 'transparent', ylab = 'intra-pair distance', xlab = 'Chromosomal position', 
       main = paste0('contig ',c.index,' - window ',wind.id,'\nno surviving offspring'), yaxt = 'n')
    axis(2, at = c(0,1,2))
    ## get genotype info and plot by fitness category
    temp <- gts.with.pos[gts.with.pos$c.index == c.index & gts.with.pos$pos >= s & gts.with.pos$pos <= e,]
    temp1 <- temp[,c(1:2)]
    for(r in 1:nrow(mates.0)){
      sub <- temp[,c(2, which(colnames(temp) == mates.0$id1[r]), which(colnames(temp) == mates.0$id2[r]))]
      temp1 <- cbind(temp1, abs(sub[,2] - sub[,3]))
    }
    for(p in 1:nrow(temp1)){
      points(x = temp1$pos[p], y = mean(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE), pch = 19, cex = 1, col = 'red')
      lines(x = c(temp1$pos[p], temp1$pos[p]), 
            y = c(mean(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE) - sd(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE),
                  mean(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE) + sd(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE)),
            col = 'red')
    }
    
  ## surviving offspring
  mates.1 <- mates[mates$surv.or.no == 1,]
  plot(0, 0, xlim = c(s, e),  ylim = c(0, 2), col = 'transparent', ylab = 'intra-pair distance', xlab = 'Chromosomal position', 
       main = paste0('contig ',c.index,' - window ',wind.id,'\nsurviving offspring'), yaxt = 'n')
    axis(2, at = c(0,1,2))
    ## get genotype info and plot by fitness category
    temp <- gts.with.pos[gts.with.pos$c.index == c.index & gts.with.pos$pos >= s & gts.with.pos$pos <= e,]
    for(r in 1:nrow(mates.1)){
      sub <- temp[,c(2, which(colnames(temp) == mates.1$id1[r]), which(colnames(temp) == mates.1$id2[r]))]
      temp1 <- cbind(temp1, abs(sub[,2] - sub[,3]))
    }
    for(p in 1:nrow(temp1)){
      points(x = temp1$pos[p], y = mean(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE), pch = 19, cex = 1, col = 'blue')
      lines(x = c(temp1$pos[p], temp1$pos[p]), 
            y = c(mean(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE) - sd(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE),
                  mean(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE) + sd(unlist(temp1[p, c(3:ncol(temp1))]), na.rm = TRUE)),
            col = 'blue')
    }
}
dev.off()

##### Trying DAPC on shared ROH (hi and no) within 2 fitness groups #####
## for each window, assess whether each mate pair has at least 1 bp in a ROHverlap region.
## add pair ID info to rohverlap data
mp.rohverlap <- rohverlap
for(p in 1:nrow(mates)){
  id1 <- mates$id1[p]
  id2 <- mates$id2[p]
  mp.rohverlap[which((mp.rohverlap$id1 == id1 & mp.rohverlap$id2 == id2) |
                       mp.rohverlap$id1 == id2 & mp.rohverlap$id2 == id1), 'pair.id'] <- mates$pair.id[p]
}
mp.rohverlap <- mp.rohverlap[!is.na(mp.rohverlap$pair.id),]

## output == matrix where each row = 1 mate pair, each column = 1 window
w.size <- 50000    ## sliding window size
s.size <- 25000     ## sliding window step size

## calculate depth of ROHverlap for sliding windows for the two categories of mate pairs
# c <- 300 ### fine value to test
# PAIR.DAT <- matrix(0, nrow = 40)
# WIND.DAT <- NULL
# wind.id <- 1
# 
# for(c in unique(contigs$c.index)){
#   c.len <- contigs[contigs$c.index == c, 'length']
#   sub <- mp.rohverlap[mp.rohverlap$c.index == c,]
#   s <- 1
#   e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
# 
#   while(e <= c.len){
#     wind.dat <- c(wind.id, c, s, e)
#     print(paste0(c/max(contigs$c.index),' - ',e/c.len))
#     WIND.DAT <- rbind(WIND.DAT, wind.dat)
#     ## comp ROHs beginning outside of window, ending inside
#     if(nrow(sub[sub$start < s & sub$end >= s & sub$end <= e,]) > 0){
#       temp <- sub[sub$start < s & sub$end >= s & sub$end <= e,]
#       PAIR.DAT[c(unique(temp$pair.id)), wind.id] <- 1
#     }
#     ## comp ROHs beginning inside of a focal window, ending outside
#     if(nrow(sub[sub$start >= s & sub$start <= e & sub$end > e,]) > 0){
#       temp <- sub[sub$start >= s & sub$start <= e & sub$end > e,]
#       PAIR.DAT[c(unique(temp$pair.id)), wind.id] <- 1
#     }
#     ## comp ROHs completely covering a focal window
#     if(nrow(sub[sub$start < s & sub$end > e,]) > 0){
#       temp <- sub[sub$start < s & sub$end > e,]
#       PAIR.DAT[c(unique(temp$pair.id)), wind.id] <- 1
#     }
#     ## comp ROHs completely within a focal window
#     if(nrow(sub[sub$start >= s & sub$end <= e,]) > 0){
#       temp <- sub[sub$start >= s & sub$end <= e,]
#       PAIR.DAT[c(unique(temp$pair.id)), wind.id] <- 1
#     }
#     wind.id <- wind.id + 1
#     PAIR.DAT <- cbind(PAIR.DAT, rep(0, 40))
#     if(e != c.len){
#       s <- s + s.size
#       e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#     } else{
#       e <- e+1
#     }
#   }
# }
# 
# ## summarize numbers of pairs with ROHs for each window
# wind.dat <- as.data.frame(WIND.DAT)
# rm(WIND.DAT)
# colnames(wind.dat) <- c('wind.id','c.index','wind.s','wind.e')
# write.table(wind.dat, '../data/LD_pruned_data_and_results/window_coords_wsize_50000.txt', quote = FALSE)
wind.dat <- read.table('../data/LD_pruned_data_and_results/window_coords_wsize_50000.txt')

# rownames(PAIR.DAT) <- c(1:40)
# colnames(PAIR.DAT) <- c(1:ncol(PAIR.DAT))
# write.table(PAIR.DAT, '../data/LD_pruned_data_and_results/mate_rohverlap_presence_absence.txt', sep = '\t',
#             quote = FALSE)
pair.dat <- read.table('../data/LD_pruned_data_and_results/mate_rohverlap_presence_absence.txt', check.names = FALSE)

## tryout DAPC
# for(c in c(30, 32, 34, 36, 38)){
#   n.pca <- c
#   dapc.res <- dapc(pair.dat, grp = mates$surv.or.no, center = TRUE, scale = TRUE, n.pca = n.pca, n.da = n.da)
#   ## n.pca setting really makes a difference here; increasing the value leads to increasing differentiation between the 2 groups
#   scatter(dapc.res) ## only 1 discriminant function when k = 2
#     legend('topleft', legend = paste0('n.pca = ',n.pca), bty = 'n')
# 
#   ## look at posteriors relative to intial fitness assignments
#   temp <- cbind(c(1:40), mates$surv.or.no)
# 
#   par(xpd = FALSE)
#   plot(c(1:40), dapc.res$posterior[,1], col = 'transparent', pch = 19, xlab = 'pair ID', ylab = 'membership probabilities',
#        main = paste0('n.pca = ',n.pca))
#     for(r in 1:nrow(dapc.res$posterior)){
#       lines(c(r, r), c(dapc.res$posterior[r, 1], dapc.res$posterior[r, 2]), col = 'grey')
#     }
#     points(c(1:40), dapc.res$posterior[,2], col = 'red', pch = 19) ## no surviving offspring
#     points(c(1:40), dapc.res$posterior[,1], col = 'blue', pch = 19) ## surviving offspring
#     abline(h = 0.5, lty = 2)
#     par(xpd = TRUE)
#     for(r in 1:nrow(temp)){
#       if(temp[r,2] == 0){
#         points(r, 1.1, pch = 19, col = 'blue')
#       } else{
#         points(r, 1.1, pch = 19, col = 'red')
#       }
#     } ## all match: predicting surv.or.no with discriminant function works
# }
# 
# a.score(dapc.res, n.da = n.da, n.pca = n.pca, plot = TRUE)
# optim.a.score(dapc.res) ## uh oh, looks like it doesn't work super well. best n.pca = 25, still not good.

n.pca <- 25
dapc.res <- dapc(pair.dat, grp = mates$surv.or.no, center = TRUE, scale = TRUE, n.pca = n.pca, n.da = n.da)

scatter(dapc.res) ## only 1 discriminant function when k = 2
  legend('topleft', legend = paste0('n.pca = ',n.pca), bty = 'n')

## look at posteriors relative to intial fitness assignments
temp <- cbind(c(1:40), mates$surv.or.no)

par(xpd = FALSE)
plot(c(1:40), dapc.res$posterior[,1], col = 'transparent', pch = 19, xlab = 'pair ID', ylab = 'membership probabilities',
     main = paste0('n.pca = ',n.pca))
  for(r in 1:nrow(dapc.res$posterior)){
    lines(c(r, r), c(dapc.res$posterior[r, 1], dapc.res$posterior[r, 2]), col = 'grey')
  }
  points(c(1:40), dapc.res$posterior[,2], col = 'red', pch = 19) ## no surviving offspring
  points(c(1:40), dapc.res$posterior[,1], col = 'blue', pch = 19) ## surviving offspring
  abline(h = 0.5, lty = 2)
  par(xpd = TRUE)
  for(r in 1:nrow(temp)){
    if(temp[r,2] == 0){
      points(r, 1.1, pch = 19, col = 'blue')
    } else{
      points(r, 1.1, pch = 19, col = 'red')
    }
  }

### check windows with big contributions
head(dapc.res$var.contr) ## contains the contributions of each window along the only discriminant function
contr <- dapc.res$var.contr
contr <- contr[order(contr[,1], decreasing = FALSE),]
hist(contr)
## all contributions in increasing order (takes a long time to plot) -- looks better than the previous approach with genetic distances tho
# plot(0, 0, ylim = c(0, max(contr)), xlim = c(0, length(contr)), col = 'transparent', xlab = '')
# for(c in 1:length(contr)){
#   lines(x = c(c, c), y = c(0, contr[c]))
# }

top.10 <- tail(contr, n = 10)[order(tail(contr, n = 10), decreasing = TRUE)]
top.10 <- dat.pca[, colnames(dat.pca) %in% names(top.10)]
mates <- cbind(mates, top.10)

## for each top window, plot distances by group
pdf('../krat_genetics_scripts/figures_output/LD_pruned/dapc_top10_window_ROHverlap.pdf', width = 6, height = 5) ## this actually doesn't used pruned data but w/e
for(c in 23:ncol(mates)){
  boxplot(mates[,c] ~ mates$surv.or.no, main = colnames(mates)[c], xlab = 'Surviving offspring or no', 
          ylab = 'Intrapair window distance')
}
dev.off()
## yeah, there's nothing here.
