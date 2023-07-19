setwd('/Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/')
`%notin%` <- Negate(`%in%`)
library(scales)
library(nationalparkcolors)
library(lme4)
library(piecewiseSEM)
library(languageR)

# ##### Read in data #####
# ## get contig info table
# contigs <- read.table('../contigs_easley.txt', header = TRUE)
# contigs <- contigs[contigs$length >= 1e6,] ## added Dec. 2, 2022 to only consider contigs >= 1 Mb in length
# 
# ##### >> Filter SnpEff annotations using est-sfs reference == ancestral allele information #####
# ## starting with LD-pruned SNPs with 0 missingness allowed (n = 1,167,410 SNPs)
# annots <- read.table('dspec_n48_snpeff_ann_LDpruned_nomissingness_snpsift.txt', header = TRUE)
# ## clean up columns
# colnames(annots) <- c('chrom','pos','ref','alt','effect','impact','gene','geneid','feature','featureid','biotype','errors')
# annots <- annots[annots$chrom %in% contigs$refseq.contig,]
# # all(annots$impact == annots$gene)     ## true
# # all(annots$impact == annots$geneid)   ## true
# annots$temp <- do.call(rbind, strsplit(annots$impact, split = '|', fixed = TRUE))[,3]
# annots <- annots[,c('chrom','pos','ref','alt','effect','temp','feature','featureid','errors')]
# colnames(annots)[6] <- 'impact'
# annots <- annots[annots$errors == '.',]
# 
# ## read in est-sfs locus input info and output p-values
# input <- read.table('/Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/est_sfs_input_spe_ord_ste_wlocs.txt',
#                     header = TRUE) ## order of alleles goes: A, C, G, T
# pvals <- read.table('/Users/Avril/Documents/programs/est-sfs-release-2.04/Dspec_ord_ste/output-pvalues_spe_ord_ste.txt')
# ## pvals column 3 = probability of major allele being ancestral
# ## get major allele and ancestral probability, combine with locus position information
# nucs <- c('A','C','G','T')
# # OUT <- NULL
# # for(r in 1:nrow(input)){
# #   as <- unlist(strsplit(input[r,3], split = ',', fixed = FALSE)) ## get spectabilis allele counts
# #   nuc <- nucs[which(as == max(as))]                              ## get corresponding allele ID for major allele
# #   if(length(nuc) == 1){                                          ## only keep loci with a major allele (i.e., not 50/50 split)
# #     save <- c(input$c.index[r], input$pos[r], nuc, pvals[r,3])     ## combine info
# #     OUT <- rbind(OUT, save)
# #     if(r %% 1000 == 0){
# #       if(r == 1000){
# #         print(r/nrow(input))
# #         write.table(OUT, 'est-sfs_majoralleles_p-vals.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
# #         OUT <- NULL
# #       } else{
# #         print(r/nrow(input))
# #         write.table(OUT, 'est-sfs_majoralleles_p-vals.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
# #         OUT <- NULL
# #       }
# #     }
# #   }
# # }
# # write.table(OUT, 'est-sfs_majoralleles_p-vals.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
# # rm(OUT)
# est.sfs <- read.table('est-sfs_majoralleles_p-vals.txt', sep = '\t')
# colnames(est.sfs) <- c('c.index','pos','maj.allele','prob')
# ## reformat allele count information in case needed later
# input <- input[,c(1:3)]
# input <- cbind(input[,c(1:2)], do.call(rbind, strsplit(input$spe, split = ',', fixed = TRUE)))
# colnames(input)[3:6] <- nucs
# 
# 
# ##### Read in rest of data #####
# ## pedigree data
# ped.dat <- read.csv('../summary_allstats2_withgen_zeroes.csv')
# ped.dat <- ped.dat[,c('id','inb','off','off_survive','birthyear')]
# 
# ## set f(ROH) formatting for plots
# froh.lab <- substitute(paste(italic('F')[ROH]))
# 
# ## read in ROH coordinates
# rohs <- read.table('../plink_final_settings/krat_plink_roh_phwh_2_phwm_30_phws_50_phzd_50_phzg_1000_phwt_0.05_phzs_50_phzk_100.hom',
#                    header = TRUE)
# rohs <- rohs[,c('FID','CHR','POS1','POS2')]
# colnames(rohs) <- c('id','contig','start','end')
# rohs <- merge(rohs, contigs[,c(1,5)], by = 'contig')
# rohs <- rohs[,c('id','refseq.contig','start','end')]
# rohs <- rohs[rohs$refseq.contig %in% contigs$refseq.contig,]
# 
# ## read in f(ROH) data
# froh <- read.csv('../individual_frohs_1Mbcontigs.csv')
# froh <- froh[,-grep('pl.', colnames(froh), fixed = TRUE)]
# froh <- froh[,-grep('gt.', colnames(froh), fixed = TRUE)]
# 
# ## read in numeric genotypes (0, 1, 2)
# samp.order <- read.table('../vcf_sample_order.txt')
# samp.order <- unlist(samp.order)
# gts <- read.table('krat_LDpruned_extracted_genotypes_numericGTs.txt', sep = '\t')
# gts <- gts[,-3]
# colnames(gts) <- c('contig','pos',samp.order)
# gts <- gts[gts$contig %in% contigs$contig,]
# gts <- merge(gts, contigs[,c(1,5)], by = 'contig') ## add refseq contig IDs to match snpeff output
# gts <- gts[,c(51,2:50)]
# for(c in 2:ncol(gts)){
#   gts[,c] <- as.numeric(gts[,c])
# }
# 
# ## narrow genotypes down to match annotations (i.e., no missingness)
# annots$loc <- paste0(annots$chrom,'-',annots$pos)
# gts$loc <- paste0(gts$refseq.contig,'-',gts$pos)
# length(unique(annots$loc))
# length(unique(gts$loc))
# gts <- gts[gts$loc %in% annots$loc,]
# annots <- annots[annots$loc %in% gts$loc,]
# gts <- gts[,-c(1,2)]
# 
# 
# ##### Check for concordance between reference and inferred ancestral allele #####
# lo <- 0.10
# hi <- 0.90
# nrow(est.sfs[est.sfs$prob >= hi,])  ## 693,679 SNPs
# nrow(est.sfs[est.sfs$prob <= lo,])  ## 137,661 SNPs
# est.sfs <- est.sfs[est.sfs$prob >= hi | est.sfs$prob <= lo,]
# OUT <- NULL
# for(r in 1:nrow(est.sfs)){
#   if(est.sfs$prob[r] <= lo){
#     temp <- input[input$c.index == est.sfs$c.index[r] & input$pos == est.sfs$pos[r], c(3:6)]
#     if(ncol(temp[,temp != 0]) == 2){    ## double-check that locus is biallelic
#       all <- names(temp[order(temp)][3]) ## mostly likely ancestral allele is minor allele
#       save <- c(est.sfs$c.index[r], est.sfs$pos[r], all, est.sfs$prob[r], 0) ## 0 = allele recorded is not major allele
#       OUT <- rbind(OUT, save)
#     } else{
#       print('NOT BIALLELIC')
#     }
#   } else{
#     save <- c(unlist(est.sfs[r,]), 1) ## 1 = allele recorded is major allele
#     OUT <- rbind(OUT, save)
#   }
#   if(r %% 1000 == 0){
#     if(r == 1000){
#       print(r/nrow(est.sfs))
#       write.table(OUT, 'est-sfs_corrected_alleles.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
#       OUT <- NULL
#     } else{
#       print(r/nrow(est.sfs))
#       write.table(OUT, 'est-sfs_corrected_alleles.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
#       OUT <- NULL
#     }
#   }
# }
# write.table(OUT, 'est-sfs_corrected_alleles.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
# rm(OUT)
# est.corr <- read.table('est-sfs_corrected_alleles.txt', sep = '\t')
# colnames(est.corr) <- c('c.index','pos','ancestral','prob','same.as.major')
# 
# ## only keep annotations where reference allele == inferred ancestral allele
# annots <- merge(annots, contigs[,c('c.index','refseq.contig')], by.x = 'chrom', by.y = 'refseq.contig')
# annots$loc <- paste0(annots$c.index,'-',annots$pos)
# est.corr$loc <- paste0(est.corr$c.index,'-',est.corr$pos)
# annots <- annots[annots$loc %in% est.corr$loc,]
# annots <- merge(annots, est.corr[,c('loc','ancestral')], by = 'loc')
# annots[annots$ref == annots$ancestral, 'ref.correct'] <- 1
# annots[annots$ref != annots$ancestral, 'ref.correct'] <- 0
# 
# save.image('PREP_07_data.RData')

##### Load data #####
load('PREP_07_data.RData')

## only keep annotations for loci where the reference allele == inferred ancestral allele
annots <- annots[annots$ref.correct == 1,]
annots$loc <- paste0(annots$chrom,'-',annots$pos)


##### Calculate individual effect totals #####
## 30,000 ft view (just 4 impact categories)
## heterozygous sites
modif.locs <- annots[annots$impact == 'MODIFIER',]
modif.locs <- merge(modif.locs, gts, by = 'loc')
OUT1 <- NULL
for(i in c(14:ncol(modif.locs))){
  OUT1 <- rbind(OUT1, table(modif.locs[,i]))
}

moder.locs <- annots[annots$impact == 'MODERATE',]
moder.locs <- merge(moder.locs, gts, by = 'loc')
OUT2 <- NULL
for(i in c(14:ncol(moder.locs))){
  OUT2 <- rbind(OUT2, table(moder.locs[,i]))
}

low.locs <- annots[annots$impact == 'LOW',]
low.locs <- merge(low.locs, gts, by = 'loc')
OUT3 <- NULL
for(i in c(14:ncol(low.locs))){
  OUT3 <- rbind(OUT3, table(low.locs[,i]))
}

high.locs <- annots[annots$impact == 'HIGH',]
high.locs <- merge(high.locs, gts, by = 'loc')
OUT4 <- NULL
for(i in c(14:ncol(high.locs))){
  OUT4 <- rbind(OUT4, table(high.locs[,i]))
}

all.sums <- as.data.frame(cbind(OUT1, OUT2, OUT3, OUT4))
all.sums$id <- colnames(modif.locs)[c(14:ncol(modif.locs))]
## #s in column names indicate genotype status (0, 1, or 2)
colnames(all.sums) <- c('modif.0','modif.1','modif.2',
                        'moder.0','moder.1','moder.2',
                        'low.0','low.1','low.2',
                        'high.0','high.1','high.2',
                        'id')
table(rowSums(all.sums[,c(1:12)]))
hist(all.sums$high.2)
hist(all.sums$moder.2)
hist(all.sums$low.2)

all.sums <- merge(all.sums, ped.dat, by = 'id')
plot(all.sums$birthyear, all.sums$high.2, pch = 19, col = alpha('dodgerblue3', 0.5))
plot(all.sums$birthyear, all.sums$moder.2, pch = 19, col = alpha('dodgerblue3', 0.5))
plot(all.sums$birthyear, all.sums$low.2, pch = 19, col = alpha('dodgerblue3', 0.5))

#### offspring ~ SnpEff SNPs, birthyear as random
pdf('../../krat_genetics_scripts/figures_output/LD_pruned/snpeff_load_v_fitness.pdf', width = 5, height = 5)
k <- 1
while(k == 1){
## HIGH impact
res <- glmer(all.sums$off ~ all.sums$high.2 + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$high.2), max(all.sums$high.2)),
             ylim = c(min(all.sums$off), max(all.sums$off)), xlabel = 'Alternate homozygous HIGH impact locs',
             ylabel = 'Number of offspring')
  points(all.sums$high.2, all.sums$off, pch = 19, col = alpha('springgreen4', 0.5))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
## MODERATE impact
res <- glmer(all.sums$off ~ all.sums$moder.2 + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$moder.2), max(all.sums$moder.2)),
             ylim = c(min(all.sums$off), max(all.sums$off)), xlabel = 'Alternate homozygous MODERATE impact locs',
             ylabel = 'Number of offspring')
points(all.sums$moder.2, all.sums$off, pch = 19, col = alpha('springgreen4', 0.5))
legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                              paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
       bty = 'n', inset = 0.02)
## LOW impact
res <- glmer(all.sums$off ~ all.sums$low.2 + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$low.2), max(all.sums$low.2)),
             ylim = c(min(all.sums$off), max(all.sums$off)), xlabel = 'Alternate homozygous LOW impact locs',
             ylabel = 'Number of offspring')
  points(all.sums$low.2, all.sums$off, pch = 19, col = alpha('springgreen4', 0.5))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)

#### surviving offspring ~ SnpEff SNPs, birthyear as random
## HIGH impact
res <- glmer(all.sums$off_survive ~ all.sums$high.2 + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$high.2), max(all.sums$high.2)),
             ylim = c(min(all.sums$off_survive), max(all.sums$off_survive)), xlabel = 'Alternate homozygous HIGH impact locs',
             ylabel = 'Number of surviving offspring')
points(all.sums$high.2, all.sums$off_survive, pch = 19, col = alpha('springgreen4', 0.5))
legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                              paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
       bty = 'n', inset = 0.02)
## MODERATE impact
res <- glmer(all.sums$off_survive ~ all.sums$moder.2 + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$moder.2), max(all.sums$moder.2)),
             ylim = c(min(all.sums$off_survive), max(all.sums$off_survive)), xlabel = 'Alternate homozygous MODERATE impact locs',
             ylabel = 'Number of surviving offspring')
points(all.sums$moder.2, all.sums$off_survive, pch = 19, col = alpha('springgreen4', 0.5))
legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                              paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
       bty = 'n', inset = 0.02)
## LOW impact
res <- glmer(all.sums$off_survive ~ all.sums$low.2 + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$low.2), max(all.sums$low.2)),
             ylim = c(min(all.sums$off_survive), max(all.sums$off_survive)), xlabel = 'Alternate homozygous LOW impact locs',
             ylabel = 'Number of surviving offspring')
points(all.sums$low.2, all.sums$off_survive, pch = 19, col = alpha('springgreen4', 0.5))
legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                              paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
       bty = 'n', inset = 0.02)

## HIGH impact
all.sums$alls <- (all.sums$high.2*2)+all.sums$high.1
res <- glmer(all.sums$off ~ all.sums$alls + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$alls), max(all.sums$alls)),
             ylim = c(min(all.sums$off), max(all.sums$off)), xlabel = 'HIGH impact derived alleles',
             ylabel = 'Number of offspring')
  points(all.sums$alls, all.sums$off, pch = 19, col = alpha('springgreen4', 0.5))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
## MODERATE impact
all.sums$alls <- (all.sums$moder.2*2)+all.sums$moder.1
res <- glmer(all.sums$off ~ all.sums$alls + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$alls), max(all.sums$alls)),
             ylim = c(min(all.sums$off), max(all.sums$off)), xlabel = 'MODERATE impact derived alleles',
             ylabel = 'Number of offspring')
points(all.sums$alls, all.sums$off, pch = 19, col = alpha('springgreen4', 0.5))
legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                              paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
       bty = 'n', inset = 0.02)
## LOW impact
all.sums$alls <- (all.sums$low.2*2)+all.sums$low.1
res <- glmer(all.sums$off ~ all.sums$alls + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$alls), max(all.sums$alls)),
             ylim = c(min(all.sums$off), max(all.sums$off)), xlabel = 'LOW impact derived alleles',
             ylabel = 'Number of offspring')
points(all.sums$alls, all.sums$off, pch = 19, col = alpha('springgreen4', 0.5))
legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                              paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
       bty = 'n', inset = 0.02)

#### surviving offspring ~ SnpEff SNPs, birthyear as random
## HIGH impact
all.sums$alls <- (all.sums$high.2*2)+all.sums$high.1
res <- glmer(all.sums$off_survive ~ all.sums$alls + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$alls), max(all.sums$alls)),
             ylim = c(min(all.sums$off_survive), max(all.sums$off_survive)), xlabel = 'HIGH impact derived alleles',
             ylabel = 'Number of surviving offspring')
  points(all.sums$alls, all.sums$off_survive, pch = 19, col = alpha('springgreen4', 0.5))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
## MODERATE impact
all.sums$alls <- (all.sums$moder.2*2)+all.sums$moder.1
res <- glmer(all.sums$off_survive ~ all.sums$alls + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$alls), max(all.sums$alls)),
             ylim = c(min(all.sums$off_survive), max(all.sums$off_survive)), xlabel = 'MODERATE impact derived alleles',
             ylabel = 'Number of surviving offspring')
  points(all.sums$alls, all.sums$off_survive, pch = 19, col = alpha('springgreen4', 0.5))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
## LOW impact
all.sums$alls <- (all.sums$low.2*2)+all.sums$low.1
res <- glmer(all.sums$off_survive ~ all.sums$alls + (1|all.sums$birthyear), family = 'poisson')
plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp, xlim = c(min(all.sums$alls), max(all.sums$alls)),
             ylim = c(min(all.sums$off_survive), max(all.sums$off_survive)), xlabel = 'LOW impact derived alleles',
             ylabel = 'Number of surviving offspring')
  points(all.sums$alls, all.sums$off_survive, pch = 19, col = alpha('springgreen4', 0.5))
  legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
         bty = 'n', inset = 0.02)
k <- k+1
}  
dev.off()


plot(all.sums$high.2, all.sums$off, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off ~ all.sums$high.2))
  summary(lm(all.sums$off ~ all.sums$high.2))
plot(all.sums$moder.2, all.sums$off, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off ~ all.sums$moder.2))
plot(all.sums$low.2, all.sums$off, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off ~ all.sums$low.2))
  summary(lm(all.sums$off ~ all.sums$low.2))
  
plot(all.sums$high.2, all.sums$off_survive, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off_survive ~ all.sums$high.2))
plot(all.sums$moder.2, all.sums$off_survive, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off_survive ~ all.sums$moder.2))
plot(all.sums$low.2, all.sums$off_survive, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off_survive ~ all.sums$low.2))
  summary(lm(all.sums$off_survive ~ all.sums$low.2)) 
  
## response = proportion offspring surviving to 1 year (100% uninteresting!)
plot(all.sums$high.2, all.sums$off_survive/all.sums$off, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off_survive/all.sums$off ~ all.sums$high.2))
plot(all.sums$moder.2, all.sums$off_survive/all.sums$off, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off_survive/all.sums$off ~ all.sums$moder.2))
plot(all.sums$low.2, all.sums$off_survive/all.sums$off, pch = 19, col = alpha('dodgerblue3', 0.5))
  abline(lm(all.sums$off_survive/all.sums$off ~ all.sums$low.2))
  summary(lm(all.sums$off_survive/all.sums$off ~ all.sums$low.2)) 

##### Calculate individual ROH-specific totals #####
annot.gts <- merge(annots, gts, by = 'loc')
annot.gts <- annot.gts[annot.gts$impact != 'MODIFIER',] 
# OUT <- NULL
# for(c in c(14:ncol(annot.gts))){     ## for each individual
#   id <- colnames(annot.gts)[c]       ## get ID
#   print(id)
#   sub.rohs <- rohs[rohs$id == id,]   ## subset to ROHs for that individual
#   sub.rohs$length <- sub.rohs$end - sub.rohs$start + 1
#   for(r in 1:nrow(annot.gts)){
#     if(nrow(sub.rohs[sub.rohs$refseq.contig == annot.gts$chrom[r] &
#                sub.rohs$start <= annot.gts$pos[r] &
#                sub.rohs$end >= annot.gts$pos[r],]) > 0){
#       temp <- sub.rohs[sub.rohs$refseq.contig == annot.gts$chrom[r] &
#                          sub.rohs$start <= annot.gts$pos[r] &
#                          sub.rohs$end >= annot.gts$pos[r],]
#       for(t in 1:nrow(temp)){
#         save <- c(id, annot.gts$chrom[r], annot.gts$pos[r], annot.gts$impact[r], annot.gts[r,c], 1, temp$length[t])
#         OUT <- rbind(OUT, save)
#       }
#     } else{
#       save <- c(id, annot.gts$chrom[r], annot.gts$pos[r], annot.gts$impact[r], annot.gts[r,c], 0, NA)
#       OUT <- rbind(OUT, save)
#     }
#   }
#   if(c == 14){
#     write.table(OUT, 'individual_snpeff_roh_annotations.csv', row.names = FALSE, sep = ',',
#               col.names = !file.exists("individual_snpeff_roh_annotations.csv"), quote = FALSE)
#     OUT <- NULL
#   } else{
#     write.table(OUT, 'individual_snpeff_roh_annotations.csv', row.names = FALSE, sep = ',', quote = FALSE,
#                 col.names = !file.exists("individual_snpeff_roh_annotations.csv"), append = TRUE)
#     OUT <- NULL
#   }
# }

roh.load <- read.csv('individual_snpeff_roh_annotations.csv')
colnames(roh.load) <- c('id','refseq.contig','pos','impact','gt','in.roh','roh.len')

## filtering for loci with multiple impacts (i.e., more than one of each impact type) by removing duplicate lines
dim(roh.load)                           ## 866,592 annotations
dim(roh.load[!duplicated(roh.load),])   ## 511,776 annotations
roh.load <- roh.load[!duplicated(roh.load),]

## for each individual and impact category, total up the number of annotated loci inside and outside of ROH regions
## (genotype counts, not allele counts)
OUT <- NULL
for(i in unique(roh.load$id)){
  print(i)
  for(m in unique(roh.load$impact)){
    sub <- roh.load[roh.load$id == i & roh.load$impact == m,]
    save <- c(i, m, 
              nrow(sub[sub$gt == 0 & sub$in.roh == 0,]), nrow(sub[sub$gt == 0 & sub$in.roh == 1,]),
              nrow(sub[sub$gt == 1 & sub$in.roh == 0,]), nrow(sub[sub$gt == 1 & sub$in.roh == 1,]),
              nrow(sub[sub$gt == 2 & sub$in.roh == 0,]), nrow(sub[sub$gt == 2 & sub$in.roh == 1,]))
    OUT <- rbind(OUT, save)
  }
}
indiv.sum <- as.data.frame(OUT)
colnames(indiv.sum) <- c('id','impact','homref.noroh','homref.roh','het.noroh','het.roh','homalt.noroh','homalt.roh')
## e.g., homref.noroh == # of homozygous reference allele genotypes located outside of ROH regions in an individual

## just some reformatting 
indiv.sum$impact <- ordered(indiv.sum$impact, levels = c('LOW','MODERATE','HIGH'))
for(c in c(3:ncol(indiv.sum))){
  indiv.sum[,c] <- as.numeric(indiv.sum[,c])
}

## within each individual, calculate the proportion of each genotype found inside ROHs (i.e., # loci inside ROH / (outside + inside))
indiv.sum$homref.roh.prop <- indiv.sum$homref.roh/(indiv.sum$homref.noroh + indiv.sum$homref.roh)
indiv.sum$het.roh.prop <- indiv.sum$het.roh/(indiv.sum$het.noroh + indiv.sum$het.roh)
indiv.sum$homalt.roh.prop <- indiv.sum$homalt.roh/(indiv.sum$homalt.noroh + indiv.sum$homalt.roh)

pdf('../../krat_genetics_scripts/figures_output/snpeff_annots_in_out_ROHs.pdf', width = 5, height = 6)
boxplot(indiv.sum$homref.roh.prop ~ indiv.sum$impact) ## homref.roh.prop == proportion of homozygous REF genotypes found inside ROH regions
## plot the raw #s for the homozygous REF genotypes
plot(0,0, xlim = c(0.9, 2.1), ylim = c(0, max(c(indiv.sum$homref.roh, indiv.sum$homref.noroh, indiv.sum$homalt.roh, indiv.sum$homalt.noroh))),
     xlab = '', ylab = '# homozygous REF genotypes', xaxt = 'n')
  axis(1, at = c(1,2), labels = c('Outside ROHs','Inside ROHs'))
  for(i in unique(indiv.sum$id)){
    sub <- indiv.sum[indiv.sum$id == i,]
    h.o <- sub[sub$impact == 'HIGH', 'homref.noroh']
    h.i <- sub[sub$impact == 'HIGH', 'homref.roh']
    m.o <- sub[sub$impact == 'MODERATE', 'homref.noroh']
    m.i <- sub[sub$impact == 'MODERATE', 'homref.roh']
    l.o <- sub[sub$impact == 'LOW', 'homref.noroh']
    l.i <- sub[sub$impact == 'LOW', 'homref.roh']
    points(c(1,2), c(h.o, h.i), pch = 19, cex = 0.75, col = 'darkorchid3')
    lines(c(1,2), c(h.o, h.i), lwd = 0.75, col = 'darkorchid3')
    points(c(1,2), c(m.o, m.i), pch = 19, cex = 0.75, col = 'chocolate2')
    lines(c(1,2), c(m.o, m.i), lwd = 0.75, col = 'chocolate2')
    points(c(1,2), c(l.o, l.i), pch = 19, cex = 0.75, col = 'cyan3')
    lines(c(1,2), c(l.o, l.i), lwd = 0.75, col = 'cyan3')
  }
  legend('topright', legend = c('HIGH','MODERATE','LOW'), col = c('darkorchid3','chocolate2','cyan3'), 
         pch = 19, pt.cex = 0.75, bty = 'n', inset = 0.02, lwd = 0.75)
  
boxplot(indiv.sum$homalt.roh.prop ~ indiv.sum$impact) ## homalt.roh.prop == proportion of homozygous ALT genotypes found inside ROH regions
## plot the raw #s for the homozygous ALT genotypes
plot(0,0, xlim = c(0.9, 2.1), ylim = c(0, max(c(indiv.sum$homalt.roh, indiv.sum$homalt.noroh, indiv.sum$homalt.roh, indiv.sum$homalt.noroh))),
     xlab = '', ylab = '# homozygous ALT genotypes', xaxt = 'n')
  axis(1, at = c(1,2), labels = c('Outside ROHs','Inside ROHs'))
  for(i in unique(indiv.sum$id)){
    sub <- indiv.sum[indiv.sum$id == i,]
    h.o <- sub[sub$impact == 'HIGH', 'homalt.noroh']
    h.i <- sub[sub$impact == 'HIGH', 'homalt.roh']
    m.o <- sub[sub$impact == 'MODERATE', 'homalt.noroh']
    m.i <- sub[sub$impact == 'MODERATE', 'homalt.roh']
    l.o <- sub[sub$impact == 'LOW', 'homalt.noroh']
    l.i <- sub[sub$impact == 'LOW', 'homalt.roh']
    points(c(1,2), c(h.o, h.i), pch = 19, cex = 0.75, col = 'darkorchid3')
    lines(c(1,2), c(h.o, h.i), lwd = 0.75, col = 'darkorchid3')
    points(c(1,2), c(m.o, m.i), pch = 19, cex = 0.75, col = 'chocolate2')
    lines(c(1,2), c(m.o, m.i), lwd = 0.75, col = 'chocolate2')
    points(c(1,2), c(l.o, l.i), pch = 19, cex = 0.75, col = 'cyan3')
    lines(c(1,2), c(l.o, l.i), lwd = 0.75, col = 'cyan3')
  }
  legend('topright', legend = c('HIGH','MODERATE','LOW'), col = c('darkorchid3','chocolate2','cyan3'), 
         pch = 19, pt.cex = 0.75, bty = 'n', inset = 0.02, lwd = 0.75)
dev.off()

indiv.sum <- merge(indiv.sum, froh, by = 'id', all.x = TRUE)
for(m in unique(indiv.sum$impact)){
  sub <- indiv.sum[indiv.sum$impact == m,]
  plot(sub$plink.froh, sub$homref.roh.prop, main = m, pch = 19, col = alpha('dodgerblue3', 0.6),
       xlab = 'f(ROH)', ylab = 'Proportion of REF homozygous sites in ROHs')
    abline(0, 1, lty = 2)
  plot(sub$plink.froh, sub$homalt.roh.prop, main = m, pch = 19, col = alpha('dodgerblue3', 0.6),
       xlab = 'f(ROH)', ylab = 'Proportion of ALT homozygous sites in ROHs')
    abline(0, 1, lty = 2)
}

### plot f(ROH) vs. genotype counts
k <- 1
while(k == 1){
pdf('../../krat_genetics_scripts/figures_output/froh_vs_snpeff_annot_counts.pdf', width = 7, height = 5)
par(mar = c(5.1, 4.1, 4.1, 5.1), xpd = FALSE)
plot(0,0, xlim = c(min(indiv.sum$plink.froh), max(indiv.sum$plink.froh)), ylim = c(0, max(c(indiv.sum$homref.roh, indiv.sum$homref.noroh, indiv.sum$homalt.roh, indiv.sum$homalt.noroh))), xlab = froh.lab, ylab = '# homozygous REF genotypes')
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homref.roh'],
         pch = 19, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homref.roh'],
         pch = 19, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homref.roh'],
         pch = 19, cex = 0.75, col = alpha('cyan3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homref.noroh'],
         pch = 17, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homref.noroh'],
         pch = 17, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homref.noroh'],
         pch = 17, cex = 0.75, col = alpha('cyan3', 0.8))
  par(xpd = TRUE)
  legend('right', inset = -0.2, legend = c('HIGH in','HIGH out','MOD in','MOD out','LOW in','LOW out'), 
         col = c('darkorchid3','darkorchid3','chocolate2','chocolate2','cyan3','cyan3'), pt.cex = 0.75,
         pch = c(19,17,19,17,19,17), bty = 'n')
  
par(mar = c(5.1, 4.1, 4.1, 5.1), xpd = FALSE)
plot(0,0, xlim = c(min(indiv.sum$plink.froh), max(indiv.sum$plink.froh)), ylim = c(0, 40), xlab = froh.lab, ylab = '# homozygous REF genotypes')
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homref.roh'],
         pch = 19, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homref.roh'],
         pch = 19, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homref.roh'],
         pch = 19, cex = 0.75, col = alpha('cyan3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homref.noroh'],
         pch = 17, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homref.noroh'],
         pch = 17, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homref.noroh'],
         pch = 17, cex = 0.75, col = alpha('cyan3', 0.8))
  par(xpd = TRUE)
  legend('right', inset = -0.2, legend = c('HIGH in','HIGH out','MOD in','MOD out','LOW in','LOW out'), 
         col = c('darkorchid3','darkorchid3','chocolate2','chocolate2','cyan3','cyan3'), pt.cex = 0.75,
         pch = c(19,17,19,17,19,17), bty = 'n')
  
par(mar = c(5.1, 4.1, 4.1, 5.1), xpd = FALSE)
plot(0,0, xlim = c(min(indiv.sum$plink.froh), max(indiv.sum$plink.froh)), ylim = c(0, max(c(indiv.sum$homalt.roh, indiv.sum$homalt.noroh, indiv.sum$homalt.roh, indiv.sum$homalt.noroh))), xlab = froh.lab, ylab = '# homozygous ALT genotypes')
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homalt.roh'],
         pch = 19, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homalt.roh'],
         pch = 19, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homalt.roh'],
         pch = 19, cex = 0.75, col = alpha('cyan3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homalt.noroh'],
         pch = 17, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homalt.noroh'],
         pch = 17, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homalt.noroh'],
         pch = 17, cex = 0.75, col = alpha('cyan3', 0.8))
  par(xpd = TRUE)
  legend('right', inset = -0.2, legend = c('HIGH in','HIGH out','MOD in','MOD out','LOW in','LOW out'), 
         col = c('darkorchid3','darkorchid3','chocolate2','chocolate2','cyan3','cyan3'), pt.cex = 0.75,
         pch = c(19,17,19,17,19,17), bty = 'n')
  
par(mar = c(5.1, 4.1, 4.1, 5.1), xpd = FALSE)
plot(0,0, xlim = c(min(indiv.sum$plink.froh), max(indiv.sum$plink.froh)), ylim = c(0, 100), xlab = froh.lab, ylab = '# homozygous ALT genotypes')
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homalt.roh'],
         pch = 19, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homalt.roh'],
         pch = 19, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homalt.roh'],
         pch = 19, cex = 0.75, col = alpha('cyan3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'homalt.noroh'],
         pch = 17, cex = 0.75, col = alpha('darkorchid3', 0.8))
  points(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'homalt.noroh'],
         pch = 17, cex = 0.75, col = alpha('chocolate2', 0.8))
  points(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'homalt.noroh'],
         pch = 17, cex = 0.75, col = alpha('cyan3', 0.8))
  par(xpd = TRUE)
  legend('right', inset = -0.2, legend = c('HIGH in','HIGH out','MOD in','MOD out','LOW in','LOW out'), 
         col = c('darkorchid3','darkorchid3','chocolate2','chocolate2','cyan3','cyan3'), pt.cex = 0.75,
         pch = c(19,17,19,17,19,17), bty = 'n')

dev.off()
k <- k+1
}
  
##### For each individual, calculate the rate of hom ALT genotypes per bp in and out of ROHs #####
indiv.sum$not.froh <- 1-indiv.sum$plink.froh
indiv.sum$alt.hom.per.roh <- indiv.sum$homalt.roh/(indiv.sum$plink.froh * max(contigs$cumsum))
indiv.sum$alt.hom.per.notroh <- indiv.sum$homalt.noroh/(indiv.sum$not.froh * max(contigs$cumsum))
indiv.sum$ref.hom.per.roh <- indiv.sum$homref.roh/(indiv.sum$plink.froh * max(contigs$cumsum))
indiv.sum$ref.hom.per.notroh <- indiv.sum$homref.noroh/(indiv.sum$not.froh * max(contigs$cumsum))

rbPal <- colorRampPalette(c('red','blue'))
indiv.sum$col.froh <- rbPal(25)[as.numeric(cut(indiv.sum$plink.froh,breaks = 25))]

plot(indiv.sum[indiv.sum$impact == 'LOW', 'alt.hom.per.notroh'], indiv.sum[indiv.sum$impact == 'LOW', 'alt.hom.per.roh'],
     pch = 19, col = alpha(indiv.sum$col.froh, 0.5), xlab = 'ALT hom per non-ROH bp', ylab = 'ALT hom per ROH bp',
     main = 'LOW impact')
  abline(0, 1, lty = 2)
plot(indiv.sum[indiv.sum$impact == 'MODERATE', 'alt.hom.per.notroh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'alt.hom.per.roh'],
     pch = 19, col = alpha(indiv.sum$col.froh, 0.5), xlab = 'ALT hom per non-ROH bp', ylab = 'ALT hom per ROH bp',
     main = 'MODERATE impact')
  abline(0, 1, lty = 2)
plot(indiv.sum[indiv.sum$impact == 'HIGH', 'alt.hom.per.notroh'], indiv.sum[indiv.sum$impact == 'HIGH', 'alt.hom.per.roh'],
     pch = 19, col = alpha(indiv.sum$col.froh, 0.5), xlab = 'ALT hom per non-ROH bp', ylab = 'ALT hom per ROH bp',
     main = 'HIGH impact')
  abline(0, 1, lty = 2)
  
plot(indiv.sum[indiv.sum$impact == 'LOW', 'plink.froh'], indiv.sum[indiv.sum$impact == 'LOW', 'alt.hom.per.roh'],
     pch = 19, col = alpha('springgreen4', 0.5), xlab = 'f(ROH)', ylab = 'ALT hom per ROH bp',
     main = 'LOW impact')
plot(indiv.sum[indiv.sum$impact == 'MODERATE', 'plink.froh'], indiv.sum[indiv.sum$impact == 'MODERATE', 'alt.hom.per.roh'],
     pch = 19, col = alpha('springgreen4', 0.5), xlab = 'f(ROH)', ylab = 'ALT hom per ROH bp',
     main = 'MODERATE impact')
plot(indiv.sum[indiv.sum$impact == 'HIGH', 'plink.froh'], indiv.sum[indiv.sum$impact == 'HIGH', 'alt.hom.per.roh'],
     pch = 19, col = alpha('springgreen4', 0.5), xlab = 'f(ROH)', ylab = 'ALT hom per ROH bp',
     main = 'HIGH impact')

##### !! Try these with GLMMs #####
indiv.sum <- merge(indiv.sum, ped.dat[,c('id','birthyear')], by = 'id')
pdf('../../krat_genetics_scripts/figures_output/snpeff_rate_vs_fitness.pdf', width = 5, height = 5)
for(m in unique(indiv.sum$impact)){
  sub <- indiv.sum[indiv.sum$impact == m,]
  
  res <- glmer(off ~ alt.hom.per.roh + (1|birthyear), data = sub, family = 'poisson')
  plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp,
               xlim = c(min(sub$alt.hom.per.roh), max(sub$alt.hom.per.roh)),
               ylim = c(min(sub$off), max(sub$off)), xlabel = 'ALT hom per ROH bp', 
               ylabel = 'Number of offspring', main = m)
    points(sub$alt.hom.per.roh, sub$off, pch = 19, col = alpha('springgreen4', 0.7))
    legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                  paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
           bty = 'n', inset = 0.02)
    
  res <- glmer(off ~ alt.hom.per.notroh + (1|birthyear), data = sub, family = 'poisson')
  plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp,
               xlim = c(min(sub$alt.hom.per.notroh), max(sub$alt.hom.per.notroh)),
               ylim = c(min(sub$off), max(sub$off)), xlabel = 'ALT hom per non-ROH bp', 
               ylabel = 'Number of offspring', main = m)
    points(sub$alt.hom.per.notroh, sub$off, pch = 19, col = alpha('springgreen4', 0.7))
    legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                  paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
           bty = 'n', inset = 0.02)
    
  res <- glmer(off_survive ~ alt.hom.per.roh + (1|birthyear), data = sub, family = 'poisson')
  plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp,
               xlim = c(min(sub$alt.hom.per.roh), max(sub$alt.hom.per.roh)),
               ylim = c(min(sub$off_survive), max(sub$off_survive)), xlabel = 'ALT hom per ROH bp', 
               ylabel = 'Number of surviving offspring', main = m)
    points(sub$alt.hom.per.roh, sub$off_survive, pch = 19, col = alpha('springgreen4', 0.7))
    legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                  paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
           bty = 'n', inset = 0.02)
    
  res <- glmer(off_survive ~ alt.hom.per.notroh + (1|birthyear), data = sub, family = 'poisson')
  plotLMER.fnc(res, linecolor = 'darkgrey', lwd = 2, n = 100, fun = exp,
               xlim = c(min(sub$alt.hom.per.notroh), max(sub$alt.hom.per.notroh)),
               ylim = c(min(sub$off_survive), max(sub$off_survive)), xlabel = 'ALT hom per non-ROH bp', 
               ylabel = 'Number of surviving offspring', main = m)
    points(sub$alt.hom.per.notroh, sub$off_survive, pch = 19, col = alpha('springgreen4', 0.7))
    legend('topright', legend = c(paste0('marg R2 = ',round(rsquared(res)[5], digits = 3)),
                                  paste0('cond R2 = ',round(rsquared(res)[6], digits = 3))),
           bty = 'n', inset = 0.02)
}
dev.off()

hi.froh.ids <- c(4459, 4195) ## indivs with highest f(ROH)
me.froh.ids <- c(5114, 3850) ## indivs with next-highest f(ROH)
pdf('../../krat_genetics_scripts/figures_output/alt_hom_per_bp_inandout_ROHs.pdf', width = 5, height = 7)
for(i in unique(indiv.sum$impact)){
  temp <- indiv.sum[indiv.sum$impact == i,]
  sub <- temp[temp$id %notin% c(hi.froh.ids, me.froh.ids),]
  plot(0, 0, xlim = c(0.5, 2.5), ylim = c(min(c(temp$alt.hom.per.notroh, temp$alt.hom.per.roh)), 
                                            max(c(temp$alt.hom.per.notroh, temp$alt.hom.per.roh))),
       xlab = 'Inside or outside of ROH regions', ylab = 'Rate of ALT hom genotypes per bp',
       xaxt = 'n', main = i)
    axis(1, at = c(1,2), labels = c('Inside','Outside'))
    
    ## inside ROHs
    lo <- mean(temp$alt.hom.per.roh) - sd(temp$alt.hom.per.roh)
    hi <- mean(temp$alt.hom.per.roh) + sd(temp$alt.hom.per.roh)
    polygon(c(0.85, 1.15, 1.15, 0.85), c(lo, lo, hi, hi), border = NA, col = alpha('springgreen4', 0.3))
    points(jitter(rep(1, nrow(sub)), factor = 5), sub$alt.hom.per.roh, col = alpha('springgreen4', 0.7), pch = 19)
    lines(c(0.75, 1.25), c(mean(temp$alt.hom.per.roh), mean(temp$alt.hom.per.roh)), col = 'springgreen4', lwd = 4)
    points(c(1,1), temp[temp$id %in% me.froh.ids, 'alt.hom.per.roh'], bg = 'orange', pch = 24, cex = 1.2)
    points(c(1,1), temp[temp$id %in% hi.froh.ids, 'alt.hom.per.roh'], bg = 'orange', pch = 23, cex = 1.2)
    
    ## outside ROHs
    lo <- mean(temp$alt.hom.per.notroh) - sd(temp$alt.hom.per.notroh)
    hi <- mean(temp$alt.hom.per.notroh) + sd(temp$alt.hom.per.notroh)
    polygon(c(1.85, 2.15, 2.15, 1.85), c(lo, lo, hi, hi), border = NA, col = alpha('springgreen4', 0.3))
    points(jitter(rep(2, nrow(sub)), factor = 5), sub$alt.hom.per.notroh, col = alpha('springgreen4', 0.7), pch = 19)
    lines(c(1.75, 2.25), c(mean(temp$alt.hom.per.notroh), mean(temp$alt.hom.per.notroh)), col = 'springgreen4', lwd = 4)
    points(c(2,2), temp[temp$id %in% me.froh.ids, 'alt.hom.per.notroh'], bg = 'orange', pch = 24, cex = 1.2)
    points(c(2,2), temp[temp$id %in% hi.froh.ids, 'alt.hom.per.notroh'], bg = 'orange', pch = 23, cex = 1.2)
    
    legend('topright', legend = c('f(ROH) ~ 0.28', 'f(ROH) ~ 0.15'), pch = c(23,24), pt.bg = 'orange', inset = 0.05, pt.cex = 1.2)
}
dev.off()


pdf('../../krat_genetics_scripts/figures_output/ref_hom_per_bp_inandout_ROHs.pdf', width = 5, height = 7)
for(i in unique(indiv.sum$impact)){
  temp <- indiv.sum[indiv.sum$impact == i,]
  sub <- temp[temp$id %notin% c(hi.froh.ids, me.froh.ids),]
  plot(0, 0, xlim = c(0.5, 2.5), ylim = c(min(c(temp$ref.hom.per.notroh, temp$ref.hom.per.roh)), 
                                          max(c(temp$ref.hom.per.notroh, temp$ref.hom.per.roh))),
       xlab = 'Inside or outside of ROH regions', ylab = 'Rate of REF hom genotypes per bp',
       xaxt = 'n', main = i)
    axis(1, at = c(1,2), labels = c('Inside','Outside'))
    
    ## inside ROHs
    lo <- mean(temp$ref.hom.per.roh) - sd(temp$ref.hom.per.roh)
    hi <- mean(temp$ref.hom.per.roh) + sd(temp$ref.hom.per.roh)
    polygon(c(0.85, 1.15, 1.15, 0.85), c(lo, lo, hi, hi), border = NA, col = alpha('springgreen4', 0.3))
    points(jitter(rep(1, nrow(sub)), factor = 5), sub$ref.hom.per.roh, col = alpha('springgreen4', 0.7), pch = 19)
    lines(c(0.75, 1.25), c(mean(temp$ref.hom.per.roh), mean(temp$ref.hom.per.roh)), col = 'springgreen4', lwd = 4)
    points(c(1,1), temp[temp$id %in% me.froh.ids, 'ref.hom.per.roh'], bg = 'orange', pch = 24, cex = 1.2)
    points(c(1,1), temp[temp$id %in% hi.froh.ids, 'ref.hom.per.roh'], bg = 'orange', pch = 23, cex = 1.2)
    
    ## outside ROHs
    lo <- mean(temp$ref.hom.per.notroh) - sd(temp$ref.hom.per.notroh)
    hi <- mean(temp$ref.hom.per.notroh) + sd(temp$ref.hom.per.notroh)
    polygon(c(1.85, 2.15, 2.15, 1.85), c(lo, lo, hi, hi), border = NA, col = alpha('springgreen4', 0.3))
    points(jitter(rep(2, nrow(sub)), factor = 5), sub$ref.hom.per.notroh, col = alpha('springgreen4', 0.7), pch = 19)
    lines(c(1.75, 2.25), c(mean(temp$ref.hom.per.notroh), mean(temp$ref.hom.per.notroh)), col = 'springgreen4', lwd = 4)
    points(c(2,2), temp[temp$id %in% me.froh.ids, 'ref.hom.per.notroh'], bg = 'orange', pch = 24, cex = 1.2)
    points(c(2,2), temp[temp$id %in% hi.froh.ids, 'ref.hom.per.notroh'], bg = 'orange', pch = 23, cex = 1.2)
    
    legend('bottomright', legend = c('f(ROH) ~ 0.28', 'f(ROH) ~ 0.15'), pch = c(23,24), pt.bg = 'orange', inset = 0.05, pt.cex = 1.2)
}
dev.off()
##### Not sure what to make of these plots yet #####

##### Take a look at SnpEff annotations for different ROH length bins #####
roh.annots <- roh.load[roh.load$in.roh == 1,] ## 790,045 annots outside of ROHs, 60,611 inside
hist(roh.annots[roh.annots$impact == 'LOW', 'roh.len'])
hist(roh.annots[roh.annots$impact == 'MODERATE', 'roh.len'])
hist(roh.annots[roh.annots$impact == 'HIGH', 'roh.len'])

## for each individual and each impact level, calculate proportion of annotations found in ROHs of different lengths
## set cutoffs for ROH length bins
rohs$length <- rohs$end - rohs$start + 1
c.1 <- 2.5e5
c.2 <- 5e5
c.3 <- 1e6
c.4 <- 2e6

OUT <- NULL
for(i in unique(roh.annots$id)){
  sub <- roh.annots[roh.annots$id == i,]
  sub.rohs <- rohs[rohs$id == i,]
  n.annots <- nrow(sub)
  f.roh <- froh[froh$id == i, 'plink.froh']
  for(m in unique(sub$impact)){
    for(g in unique(sub$gt)){
      temp <- sub[sub$impact == m & sub$gt == g,]
      bin1 <- nrow(temp[temp$roh.len < c.1,])/sum(sub.rohs[sub.rohs$length < c.1, 'length'])
      bin2 <- nrow(temp[temp$roh.len >= c.1 & temp$roh.len < c.2,])/sum(sub.rohs[sub.rohs$length >= c.1 & sub.rohs$length < c.2, 'length'])
      bin3 <- nrow(temp[temp$roh.len >= c.2 & temp$roh.len < c.3,])/sum(sub.rohs[sub.rohs$length >= c.2 & sub.rohs$length < c.3, 'length'])
      bin4 <- nrow(temp[temp$roh.len >= c.3 & temp$roh.len < c.4,])/sum(sub.rohs[sub.rohs$length >= c.3 & sub.rohs$length < c.4, 'length'])
      bin5 <- nrow(temp[temp$roh.len >= c.4,])/sum(sub.rohs[sub.rohs$length >= c.4, 'length'])
      n.bin1 <- nrow(temp[temp$roh.len < c.1,])
      n.bin2 <- nrow(temp[temp$roh.len >= c.1 & temp$roh.len < c.2,])
      n.bin3 <- nrow(temp[temp$roh.len >= c.2 & temp$roh.len < c.3,])
      n.bin4 <- nrow(temp[temp$roh.len >= c.3 & temp$roh.len < c.4,])
      n.bin5 <- nrow(temp[temp$roh.len >= c.4,])
      save <- c(i, f.roh, n.annots, m, g, bin1, bin2, bin3, bin4, bin5, n.bin1, n.bin2, n.bin3, n.bin4, n.bin5,
                sum(sub.rohs[sub.rohs$length < c.1, 'length'])/max(contigs$cumsum),
                sum(sub.rohs[sub.rohs$length >= c.1 & sub.rohs$length < c.2, 'length'])/max(contigs$cumsum),
                sum(sub.rohs[sub.rohs$length >= c.2 & sub.rohs$length < c.3, 'length'])/max(contigs$cumsum),
                sum(sub.rohs[sub.rohs$length >= c.3 & sub.rohs$length < c.4, 'length'])/max(contigs$cumsum),
                sum(sub.rohs[sub.rohs$length >= c.4, 'length'])/max(contigs$cumsum))
      save <- gsub('NaN','NA', save)
      OUT <- rbind(OUT, save)
    }
  }
}  
## values for each individual, impact level, and ROH length size are standardized against the total length of each bin in each individual
roh.bin.annots <- as.data.frame(OUT)
colnames(roh.bin.annots) <- c('id','froh','n.indiv.annots','impact','gt','bin1','bin2','bin3','bin4','bin5','n.bin1','n.bin2','n.bin3','n.bin4','n.bin5',
                              'bin1.froh','bin2.froh','bin3.froh','bin4.froh','bin5.froh')
for(c in grep('bin', colnames(roh.bin.annots))){
  roh.bin.annots[,c] <- as.numeric(roh.bin.annots[,c])
}

cols <- park_palette('Everglades')
# pdf('../../krat_genetics_scripts/figures_output/althom_rate_by_rohbin.pdf', width = 15, height = 9)
par(mfrow = c(3,5))
for(m in unique(roh.bin.annots$impact)){
  sub <- roh.bin.annots[roh.bin.annots$impact == m & roh.bin.annots$gt == 2,]
  plot(sub$bin1.froh, sub$bin1, main = m, pch = 19, col = alpha(cols[1], 0.6))
  plot(sub$bin2.froh, sub$bin2, main = m, pch = 19, col = alpha(cols[2], 0.6))
  plot(sub$bin3.froh, sub$bin3, main = m, pch = 19, col = alpha(cols[3], 0.6))
  plot(sub$bin4.froh, sub$bin4, main = m, pch = 19, col = alpha(cols[4], 0.6))
  plot(sub$bin5.froh, sub$bin5, main = m, pch = 19, col = alpha(cols[5], 0.6))
}
# dev.off()

## uniform y-axes across plots
for(m in unique(roh.bin.annots$impact)){
  temp <- roh.bin.annots[roh.bin.annots$gt == 2,]
  sub <- temp[temp$impact == m,]
  plot(jitter(rep(1, nrow(sub)), factor = 3), sub$bin1, pch = 19, col = alpha(cols[1], 0.6), xlim = c(0.8, 5.2),
       ylim = c(0, max(c(temp$bin1, temp$bin2, temp$bin3, temp$bin4, temp$bin5), na.rm = TRUE)),
       xlab = 'ROH length bin',
       ylab = 'Number of ALT hom genotypes per ROH bin bp', main = m, xaxt = 'n')
    axis(1, at = c(1,2,3,4,5), labels = c('< 250 kb','250-500 kb','500 kb -1 Mb','1-2 Mb','> 2 Mb'))
    points(jitter(rep(2, nrow(sub)), factor = 3), sub$bin2, pch = 19, col = alpha(cols[2], 0.6))
    points(jitter(rep(3, nrow(sub)), factor = 3), sub$bin3, pch = 19, col = alpha(cols[3], 0.6))
    points(jitter(rep(4, nrow(sub)), factor = 3), sub$bin4, pch = 19, col = alpha(cols[4], 0.6))
    points(jitter(rep(5, nrow(sub)), factor = 3), sub$bin5, pch = 19, col = alpha(cols[5], 0.6))
}
## dynamic y-axes across plots
pdf('../../krat_genetics_scripts/figures_output/althom_rate_by_rohbin.pdf', width = 14, height = 4)
par(mfrow = c(1,3))
for(m in unique(roh.bin.annots$impact)){
  temp <- roh.bin.annots[roh.bin.annots$gt == 2,]
  sub <- temp[temp$impact == m,]
  ## Rate of ALT hom genotypes for each length bin
  plot(jitter(rep(1, nrow(sub)), factor = 3), sub$bin1, pch = 19, col = alpha(cols[1], 0.6), xlim = c(0.8, 5.2),
       ylim = c(0, max(c(sub$bin1, sub$bin2, sub$bin3, sub$bin4, sub$bin5), na.rm = TRUE)),
       xlab = 'ROH length bin',
       ylab = 'Number of ALT hom genotypes per ROH bin bp', main = m, xaxt = 'n')
    axis(1, at = c(1,2,3,4,5), labels = c('< 250 kb','250-500 kb','500 kb -1 Mb','1-2 Mb','> 2 Mb'))
    points(jitter(rep(2, nrow(sub)), factor = 3), sub$bin2, pch = 19, col = alpha(cols[2], 0.6))
    points(jitter(rep(3, nrow(sub)), factor = 3), sub$bin3, pch = 19, col = alpha(cols[3], 0.6))
    points(jitter(rep(4, nrow(sub)), factor = 3), sub$bin4, pch = 19, col = alpha(cols[4], 0.6))
    points(jitter(rep(5, nrow(sub)), factor = 3), sub$bin5, pch = 19, col = alpha(cols[5], 0.6))
    if(m == 'HIGH'){
      legend('topleft', legend = 'NOT INFORMATIVE', border = 'red', cex = 3)
    }
  
  # ## Length bin f(ROH) vs rate of ALT hom genotypes
  # plot(0, 0, xlim = c(0, max(c(sub$bin1.froh, sub$bin2.froh, sub$bin3.froh, sub$bin4.froh, sub$bin5.froh), na.rm = TRUE)),
  #      ylim = c(0, max(c(sub$bin1, sub$bin2, sub$bin3, sub$bin4, sub$bin5), na.rm = TRUE)), col = 'transparent',
  #      xlab = 'Bin f(ROH)', ylab = 'Number of ALT hom genotypes per ROH bin bp', main = m)
  #   points(sub$bin1.froh, sub$bin1, pch = 19, col = alpha(cols[1], 0.6))
  #   points(sub$bin2.froh, sub$bin2, pch = 19, col = alpha(cols[2], 0.6))
  #   points(sub$bin3.froh, sub$bin3, pch = 19, col = alpha(cols[3], 0.6))
  #   points(sub$bin4.froh, sub$bin4, pch = 19, col = alpha(cols[4], 0.6))
  #   points(sub$bin5.froh, sub$bin5, pch = 19, col = alpha(cols[5], 0.6))
  #   legend('topright', legend = c('bin 1','bin 2','bin 3','bin 4','bin 5'), col = cols, pch = 19)
  #   if(m == 'HIGH'){
  #     legend('topleft', legend = 'NOT INFORMATIVE', border = 'red', cex = 3)
  #   }
  
  ## Same as previous, but with lines at x = median f(ROH) and y = within-bin range of rates
  plot(0, 0, xlim = c(0, max(c(sub$bin1.froh, sub$bin2.froh, sub$bin3.froh, sub$bin4.froh, sub$bin5.froh), na.rm = TRUE)),
       ylim = c(0, max(c(sub$bin1, sub$bin2, sub$bin3, sub$bin4, sub$bin5), na.rm = TRUE)), col = 'transparent',
       xlab = 'Bin f(ROH)', ylab = 'Number of ALT hom genotypes per ROH bin bp', main = m)
    points(sub$bin1.froh, sub$bin1, pch = 19, col = alpha(cols[1], 0.3))
    points(sub$bin2.froh, sub$bin2, pch = 19, col = alpha(cols[2], 0.3))
    points(sub$bin3.froh, sub$bin3, pch = 19, col = alpha(cols[3], 0.3))
    points(sub$bin4.froh, sub$bin4, pch = 19, col = alpha(cols[4], 0.3))
    points(sub$bin5.froh, sub$bin5, pch = 19, col = alpha(cols[5], 0.3))
    
    lines(c(mean(sub$bin1.froh), mean(sub$bin1.froh)), c(min(sub$bin1), max(sub$bin1)), lwd = 4, col = cols[1])
    lines(c(mean(sub$bin2.froh), mean(sub$bin2.froh)), c(min(sub$bin2), max(sub$bin2)), lwd = 4, col = cols[2])
    lines(c(mean(sub$bin3.froh), mean(sub$bin3.froh)), c(min(sub$bin3), max(sub$bin3)), lwd = 4, col = cols[3])
    lines(c(mean(sub$bin4.froh), mean(sub$bin4.froh)), c(min(sub$bin4), max(sub$bin4)), lwd = 4, col = cols[4])
    lines(c(mean(sub$bin5.froh, na.rm = TRUE), mean(sub$bin5.froh, na.rm = TRUE)), 
          c(min(sub$bin5, na.rm = TRUE), max(sub$bin5, na.rm = TRUE)), lwd = 4, col = cols[5])
    
    legend('topright', legend = c('bin 1','bin 2','bin 3','bin 4','bin 5'), col = cols, pch = 19)
    if(m == 'HIGH'){
      legend('topleft', legend = 'NOT INFORMATIVE', border = 'red', cex = 3)
    }
    
  ## Length bin f(ROH) vs. # of ALT hom genotypes
  plot(0, 0, xlim = c(0, max(c(sub$bin1.froh, sub$bin2.froh, sub$bin3.froh, sub$bin4.froh, sub$bin5.froh), na.rm = TRUE)),
       ylim = c(0, max(c(sub$n.bin1, sub$n.bin2, sub$n.bin3, sub$n.bin4, sub$n.bin5), na.rm = TRUE)), col = 'transparent',
       xlab = 'Bin f(ROH)', ylab = 'Number of ALT hom genotypes', main = m)
    if(m == 'HIGH'){ ## only jitter the results from this category because the range is [0,4]
      points(sub$bin1.froh, jitter(sub$n.bin1, factor = 0.5), pch = 19, col = alpha(cols[1], 0.6))
      points(sub$bin2.froh, jitter(sub$n.bin2, factor = 0.5), pch = 19, col = alpha(cols[2], 0.6))
      points(sub$bin3.froh, jitter(sub$n.bin3, factor = 0.5), pch = 19, col = alpha(cols[3], 0.6))
      points(sub$bin4.froh, jitter(sub$n.bin4, factor = 0.5), pch = 19, col = alpha(cols[4], 0.6))
      points(sub$bin5.froh, jitter(sub$n.bin5, factor = 0.5), pch = 19, col = alpha(cols[5], 0.6))
    } else{
      points(sub$bin1.froh, sub$n.bin1, pch = 19, col = alpha(cols[1], 0.6))
      points(sub$bin2.froh, sub$n.bin2, pch = 19, col = alpha(cols[2], 0.6))
      points(sub$bin3.froh, sub$n.bin3, pch = 19, col = alpha(cols[3], 0.6))
      points(sub$bin4.froh, sub$n.bin4, pch = 19, col = alpha(cols[4], 0.6))
      points(sub$bin5.froh, sub$n.bin5, pch = 19, col = alpha(cols[5], 0.6))
    }
    legend('topleft', legend = c('bin 1','bin 2','bin 3','bin 4','bin 5'), col = cols, pch = 19)
}
dev.off()  
