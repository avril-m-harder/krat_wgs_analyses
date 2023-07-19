library(scales)
`%notin%` <- Negate(`%in%`)
setwd('/Users/Avril/Documents/krat_genetics/qualimap/')

## get contig lengths
lens <- read.table('GCF_019054845.1_ASM1905484v1_assembly_report.txt')
lens <- lens[,c('Assigned.Molecule.Location.Type','Sequence.Length')]
colnames(lens) <- c('contig','len')
ass.len <- sum(lens$len) ## total assembly length

## set minimum contig length
min.l <- 100000
lens <- lens[which(lens$len >= min.l),]

## get sample info
samps <- read.csv('/Users/Avril/Documents/krat_genetics/assembly_seq_quotes_sample_notes_results/reseq_prep/final_reseq_samp_list.csv')
samps <- samps[,c('id','sex')]
samps[samps$sex == 'Male', 'colour'] <- 'dodgerblue'
samps[samps$sex == 'Female', 'colour'] <- 'magenta'

## loop over Qualimap results for all samples (based on BAM files of cleaned reads mapped to entire assembly)
dirs <- list.files()
dirs <- dirs[4:51]
OUT <- NULL
for(i in dirs){
  dat <- read.table(paste0('./',i,'/contig_stats.txt'))
  colnames(dat) <- c('contig','num.mapped.bases','num.base.reads.maybe','mean.cov','sd.cov')
  dat$id <- strsplit(i, split='_', fixed=TRUE)[[1]][2]
  OUT <- rbind(OUT, dat)
}

dat <- merge(OUT, samps, by='id')
dat <- dat[which(dat$contig %in% lens$contig),]
dat <- merge(dat, lens, by='contig')
dat <- dat[order(dat$len, dat$contig),]

j <- 1
for(i in unique(dat$contig)){
  dat[which(dat$contig == i), 'index'] <- j
  j <- j+1
}

plot(dat$index, dat$mean.cov, pch=19, col=alpha(dat$colour, 0.4))

## ID contigs with apparent sex-biased coverage
OUT <- NULL
for(i in unique(dat$contig)){
  sub <- dat[dat$contig == i,]
  m.avg <- mean(sub[which(sub$sex == 'Male'), 'mean.cov'])
  f.avg <- mean(sub[which(sub$sex == 'Female'), 'mean.cov'])
  save <- c(i, m.avg, f.avg)
  OUT <- rbind(OUT, save)
}
sex.cov <- as.data.frame(OUT, row.names=FALSE)
colnames(sex.cov) <- c('contig','m.avg.cov','f.avg.cov')
sex.cov$m.avg.cov <- as.numeric(sex.cov$m.avg.cov)
sex.cov$f.avg.cov <- as.numeric(sex.cov$f.avg.cov)
sex.cov$ratio <- sex.cov$f.avg.cov/sex.cov$m.avg.cov
sex.cov <- sex.cov[order(sex.cov$ratio),]

# pdf('/Users/Avril/Desktop/coverage_sex_ratio.pdf', width=5, height=8)
# par(mfrow=c(3,1))
hist(sex.cov$ratio, breaks=100, main='Coverage sex ratio', xlab='Coverage sex ratio (F/M)')
hist(sex.cov$ratio, breaks=100, xlim=c(0,0.5), ylim=c(0,20), main='Y-chromosome candidate contigs', xlab='Coverage sex ratio (F/M)') ## 0.3 cutoff for Y-chrom? or maybe 0.5 to get rid of that one sus blip?
hist(sex.cov$ratio, breaks=100, xlim=c(1.5,2), ylim=c(0,20), main='X-chromosome candidate contigs', xlab='Coverage sex ratio (F/M)') ## 1.7 cutoff for X-chrom?
# dev.off()

## Y-chrom checks
y.cand <- sex.cov[which(sex.cov$ratio <= 0.3),]
sum(lens[which(lens$contig %in% y.cand$contig), 'len']) ## total length for Y-chrom candidate contigs = 50,870,389 bp
y.dat <- dat[which(dat$contig %in% y.cand$contig),]
j <- 1
for(i in unique(y.dat$contig)){
  y.dat[which(y.dat$contig == i), 'index'] <- j
  j <- j+1
}
# pdf('/Users/Avril/Desktop/Y_chrom_candidate_contigs.pdf', width=8, height=5)
par(xpd=TRUE)
plot(y.dat$index, y.dat$mean.cov, pch=19, col=alpha(y.dat$colour, 0.6), main='Y-chrom candidate contigs',
     xlab='Contig index', ylab='Mean coverage (X)')
  legend('topright', horiz=TRUE, legend=c('Male','Female'), pch=19, col=alpha(c('dodgerblue','magenta'), 0.6), inset=c(0.02,-0.15))
# dev.off()

## X-chrom checks
x.cand <- sex.cov[which(sex.cov$ratio >= 1.7),]
sum(lens[which(lens$contig %in% x.cand$contig), 'len']) ## total length for X-chrom candidate contigs = 116,799,422 bp
x.dat <- dat[which(dat$contig %in% x.cand$contig),]
j <- 1
for(i in unique(x.dat$contig)){
  x.dat[which(x.dat$contig == i), 'index'] <- j
  j <- j+1
}
# pdf('/Users/Avril/Desktop/X_chrom_candidate_contigs.pdf', width=8, height=5)
par(xpd=TRUE)
plot(x.dat$index, x.dat$mean.cov, pch=19, col=alpha(x.dat$colour, 0.6), main='X-chrom candidate contigs',
     xlab='Contig index', ylab='Mean coverage (X)')
  legend('topright', horiz=TRUE, legend=c('Male','Female'), pch=19, col=alpha(c('dodgerblue','magenta'), 0.6), inset=c(0.02,-0.15))
# dev.off()

## Remove X- and Y-chrom candidate contigs
lens <- lens[which(lens$contig %notin% x.cand$contig & lens$contig %notin% y.cand$contig),]
dat <- dat[which(dat$contig %in% lens$contig),]
# pdf('/Users/Avril/Desktop/all_contigs_length_vs_covg.pdf', width=5, height=4)
plot(dat$len, dat$mean.cov, pch=19, cex=0.5, col=alpha('black',0.6), xlab='Contig length', ylab='Contig mean coverage', main='All contigs') ## still a couple of contigs with really high mean coverage
# dev.off()

## Get genome-wide coverage average from Qualimap results (1 genome-wide mean per sample, calculated from autosomal contigs)
OUT <- NULL
for(i in unique(dat$id)){
  sub <- dat[dat$id == i,]
  sub$base.cov <- sub$num.mapped.bases * sub$mean.cov
  save <- c(i, sum(sub$base.cov)/sum(sub$len))
  OUT <- rbind(OUT, save)
}
samp.covs <- as.data.frame(OUT)
samp.covs$V2 <- as.numeric(samp.covs$V2)
colnames(samp.covs) <- c('id','mean.covg')
genome.mean.cov <- mean(samp.covs$mean.covg)

## calculate 1 per-contig mean coverage value by averaging over samples within a contig
contig.cov.means <- NULL
for(i in unique(dat$contig)){
  contig.cov.means <- rbind(contig.cov.means, c(i, mean(dat[which(dat$contig == i), 'mean.cov'])))
}
contig.cov.means <- as.data.frame(contig.cov.means)
colnames(contig.cov.means) <- c('contig','mean.cov')
contig.cov.means$mean.cov <- as.numeric(contig.cov.means$mean.cov)

## Identify high-coverage contigs
## by mean across all samples
high.cov <- unique(contig.cov.means[which(contig.cov.means$mean.cov > 2*genome.mean.cov), 'contig']) ## 37 high-coverage contigs
## and by which contigs have REALLY high depth in >= 1 sample(s)
high.cov <- c(high.cov, unique(dat[which(dat$mean.cov > 3*genome.mean.cov), 'contig']))
high.cov <- unique(high.cov) ## 74 high-coverage contigs

sum(lens[which(lens$contig %in% high.cov), 'len']) ## totaling 65,904,744 bp
sum(lens[which(lens$contig %notin% high.cov), 'len']) ## would leave 2,576,195,715 bp
# pdf('/Users/Avril/Desktop/high_covg_contigs.pdf', width=5, height=4)
plot(dat[which(dat$contig %in% high.cov), 'len'], dat[which(dat$contig %in% high.cov), 'mean.cov'], 
     pch=19, cex=0.5, col=alpha('black',0.6), xlab='Contig length', ylab='Contig mean coverage', main='High-coverage contigs')
# dev.off()

## Identify low-coverage contigs
## by mean across all samples
low.cov <- unique(contig.cov.means[which(contig.cov.means$mean.cov < 0.5*genome.mean.cov), 'contig']) ## 8 low-coverage contigs
## and by which contigs have REALLY high depth in >= 1 sample(s)
low.cov <- c(low.cov, unique(dat[which(dat$mean.cov < (1/5)*genome.mean.cov), 'contig']))
low.cov <- unique(low.cov) ## 12 low-coverage contigs

sum(lens[which(lens$contig %in% low.cov), 'len']) ## totaling 4,956,603 bp
sum(lens[which(lens$contig %notin% low.cov), 'len']) ## would leave 2,637,143,856 bp
# pdf('/Users/Avril/Desktop/low_covg_contigs.pdf', width=5, height=4)
plot(dat[which(dat$contig %in% low.cov), 'len'], dat[which(dat$contig %in% low.cov), 'mean.cov'], 
     pch=19, cex=0.5, col=alpha('black',0.6), xlab='Contig length', ylab='Contig mean coverage', main='Low-coverage contigs')
# dev.off()

## Further filter X-/Y-filtered contigs to get rid of high- and low-coverage contigs
lens <- lens[which(lens$contig %notin% high.cov & lens$contig %notin% low.cov),]
dat <- dat[which(dat$contig %notin% high.cov & dat$contig %notin% low.cov),]
contig.cov.means <- contig.cov.means[which(contig.cov.means$contig %in% dat$contig),]
contig.cov.means <- merge(contig.cov.means, lens, by='contig')

# pdf('/Users/Avril/Desktop/covg_filtered_contigs.pdf', width=5, height=4)
plot(dat$len, dat$mean.cov, pch=19, cex=0.5, col=alpha('black',0.6), xlab='Contig length', ylab='Contig mean coverage', main='Filtered contigs')
plot(contig.cov.means$len, contig.cov.means$mean.cov, pch=19, cex=0.5, col=alpha('black',0.6), xlab='Contig length', ylab='Contig mean coverage', main='Filtered contigs')
# dev.off()

## Post-filtering stats
sum(lens[,'len']) ## 2,571,239,112 bp ()
nrow(lens) ## 663 contigs
hist(dat$sd.cov) ## SDs much lower than unfiltered contig set

## Write list of filtered contigs
write.table(lens$contig, '/Users/Avril/Documents/krat_genetics/contig_filtering/filtered_contigs.txt', 
            quote=FALSE, row.names=FALSE, col.names=FALSE)
## Format contig information to write in BED format for filtering bioinfo files
bed.format <- lens
bed.format$start <- 0
bed.format <- bed.format[,c(1,3,2)]
write.table(bed.format, '/Users/Avril/Documents/krat_genetics/contig_filtering/filtered_contigs.bed',
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

## Using filtered contigs, get depth for samples
OUT <- NULL
for(i in unique(dat$id)){
  temp <- dat[dat$id == i,]
  temp$len.cov <- temp$mean.cov * temp$len
  cov <- sum(temp$len.cov)/sum(temp$len)
  save <- c(i, cov)
  OUT <- rbind(OUT, save)
}
write.table(OUT, '/Users/Avril/Documents/krat_genetics/seq_processing_notes/filtered_contigs_sample_mean_covg.txt',
            sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

## write file of contig-specific coverage information broken down by samples
write.table(dat[,c(1,2,5,6,7)], '/Users/Avril/Documents/krat_genetics/seq_processing_notes/contig_read_depths_by_sample.txt',
            quote = FALSE, row.names = FALSE)
