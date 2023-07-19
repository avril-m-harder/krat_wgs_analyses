### A script to read in MultiQC and samtools flagstat results to track 
### read retention and mapping quality 

setwd('/Users/Avril/Documents/krat_genetics/')
library(scales)

##### MultiQC results #####
### >>> Pre-trim #####
pre.full <- read.table('qc_stuff/fastqc/krat_pretrim/multiqc_data/multiqc_fastqc.txt', sep = '\t', header = TRUE)
pre.full$id <- do.call(rbind, strsplit(pre.full$Sample, split='_'))[,1]
pre.full$read.group <- do.call(rbind, strsplit(pre.full$Sample, split='_'))[,3]
## for each individual and read group, total up original # sequences
OUT <- NULL
for(i in unique(pre.full$id)){
  sub <- pre.full[pre.full$id == i,]
  for(r in unique(sub$read.group)){
    save <- c(i, r, sum(sub[sub$read.group == r, 'Total.Sequences']))
    OUT <- rbind(OUT, save)
  }
}
colnames(OUT) <- c('id','read.group','pretrim.total.seqs')

### >>> Post-trim #####
post.full <- read.table('qc_stuff/fastqc/krat_posttrim/multiqc_data/multiqc_fastqc.txt', sep = '\t', header = TRUE)
post.full$id <- do.call(rbind, strsplit(post.full$Sample, split='_'))[,1]
post.full$read.group <- do.call(rbind, strsplit(post.full$Sample, split='_'))[,3]
OUT1 <- NULL
for(i in unique(post.full$id)){
  sub <- post.full[post.full$id == i,]
  for(r in unique(sub$read.group)){
    save <- c(i, r, sum(sub[sub$read.group == r, 'Total.Sequences']))
    OUT1 <- rbind(OUT, save)
  }
}
colnames(OUT1) <- c('id','read.group','posttrim.total.seqs')

num.reads <- merge(OUT, OUT1, by=c('id','read.group'))

##### !!! # of reads largely unchanged? wtf? #####
### --> actually, this is because TrimGalore did a lot of trimming of shitty bases/
### adapters, and the vast majority of affected reads were retained.

## for each individual, total up the post-trim reads across read groups
OUT <- NULL
for(i in unique(num.reads$id)){
  save <- c(i, sum(as.numeric(num.reads[num.reads$id == i, 'posttrim.total.seqs'])))
  OUT <- rbind(OUT, save)
}
colnames(OUT) <- c('id','posttrim.reads')

##### Samtools flagstat results #####
fns <- list.files('qc_stuff/samtools_flagstat/original_single_sample_files/')
OUT1 <- NULL
for(f in fns){
  dat <- read.table(paste0('qc_stuff/samtools_flagstat/original_single_sample_files/',f), sep = '\t')
  total <- dat[dat$V3 == 'total (QC-passed reads + QC-failed reads)', 1]
  map <- dat[dat$V3 == 'mapped', 1]
  pair <- dat[dat$V3 == 'properly paired', 1]
  id <- strsplit(f, split = '_', fixed = TRUE)[[1]][1]
  save <- c(id, total, map, pair)
  OUT1 <- rbind(OUT1, save)
}
colnames(OUT1) <- c('id','total','mapped','prop.paired')
## combine U of I and Duke totals for double-seq'ed individuals
OUT2 <- NULL
for(i in unique(OUT1[,1])){
  sub <- OUT1[OUT1[,1] == i,]
  if(length(sub) != 4){
    save <- c(i, sum(as.numeric(sub[,2])), sum(as.numeric(sub[,3])), sum(as.numeric(sub[,4])))
    OUT2 <- rbind(OUT2, save)
  } else{
    OUT2 <- rbind(OUT2, sub)
  }
}

samp.num.reads <- merge(OUT, OUT2, by='id')
samp.num.reads$posttrim.reads <- as.numeric(samp.num.reads$posttrim.reads)
samp.num.reads$total <- as.numeric(samp.num.reads$total)
samp.num.reads$mapped <- as.numeric(samp.num.reads$mapped)
samp.num.reads$prop.paired <- as.numeric(samp.num.reads$prop.paired)
samp.num.reads$p.total <- samp.num.reads$total/samp.num.reads$posttrim.reads
samp.num.reads$p.mapped <- samp.num.reads$mapped/samp.num.reads$total
samp.num.reads$p.prop.paired <- samp.num.reads$prop.paired/samp.num.reads$total

##### Combined BAM stats for 10 resequenced samples #####
fns <- list.files('qc_stuff/samtools_flagstat/combined_bams/')
OUT1 <- NULL
for(f in fns){
  dat <- read.table(paste0('qc_stuff/samtools_flagstat/combined_bams/',f), sep = '\t')
  total <- dat[dat$V3 == 'total (QC-passed reads + QC-failed reads)', 1]
  map <- dat[dat$V3 == 'mapped', 1]
  pair <- dat[dat$V3 == 'properly paired', 1]
  id <- strsplit(f, split = '_', fixed = TRUE)[[1]][1]
  save <- c(id, total, map, pair)
  OUT1 <- rbind(OUT1, save)
}
colnames(OUT1) <- c('id','comb.total','comb.mapped','comb.prop.paired')
reseq.check <- merge(samp.num.reads, OUT1, by='id')
all(reseq.check$total == reseq.check$comb.total)
all(reseq.check$mapped == reseq.check$comb.mapped)
all(reseq.check$prop.paired == reseq.check$comb.prop.paired)
### all equivalent, good to go.


##### Read in BCFtools stats results on final VCF for all samples
vcf.stats <- read.csv('qc_stuff/bcftools_stats/just_sample_stats.txt', check.names = FALSE)
hist(vcf.stats$`average depth`)

## total # of sites from bcftools stats:
n.sites <- 11796181
vcf.stats$prop.ref.hom <- vcf.stats$nRefHom/n.sites
vcf.stats$prop.alt.hom <- vcf.stats$nNonRefHom/n.sites
vcf.stats$prop.het <- vcf.stats$nHets/n.sites

hist(vcf.stats$prop.ref.hom)
hist(vcf.stats$prop.alt.hom)
hist(vcf.stats$prop.het)

hist(vcf.stats$nMissing)
vcf.stats$prop.miss <- vcf.stats$nMissing/n.sites
hist(vcf.stats$prop.miss) ## proportion of sites with missing data per sample

## depth : missing data
vcf.stats <- vcf.stats[order(-vcf.stats$prop.miss),]
head(vcf.stats, n = 2) ## 2 samples with concerning amounts of missing data (9% for both)
plot(vcf.stats$`average depth`, vcf.stats$prop.miss, pch = 19, col = alpha('dodgerblue3', 0.6)) ## 4 samples look pretty inconsistent/bad
  points(vcf.stats[c(1:3,5), 'average depth'], 
         vcf.stats[c(1:3,5), 'prop.miss'], pch = 19, col = 'red')
vcf.stats$sample[c(1:3,5)] ## 4 problematic samples

## depth : heterozygosity (positive relationship - not sure what to do about this)
plot(vcf.stats$`average depth`, vcf.stats$prop.het, pch = 19, col = alpha('dodgerblue3', 0.6)) 
  points(vcf.stats[c(1:3,5), 'average depth'], 
         vcf.stats[c(1:3,5), 'prop.het'], pch = 19, col = 'red')
  
plot(vcf.stats$`average depth`, vcf.stats$prop.ref.hom, pch = 19, col = alpha('dodgerblue3', 0.6)) 
  points(vcf.stats[c(1:3,5), 'average depth'], 
         vcf.stats[c(1:3,5), 'prop.ref.hom'], pch = 19, col = 'red')
  
plot(vcf.stats$`average depth`, vcf.stats$prop.alt.hom, pch = 19, col = alpha('dodgerblue3', 0.6)) 
  points(vcf.stats[c(1:3,5), 'average depth'], 
         vcf.stats[c(1:3,5), 'prop.alt.hom'], pch = 19, col = 'red')

## depth : het, accounting for missing data
vcf.stats$norm.prop.het <- vcf.stats$nHets/(n.sites - vcf.stats$nMissing)
plot(vcf.stats$`average depth`, vcf.stats$norm.prop.het, pch = 19, col = alpha('dodgerblue3', 0.6)) 
  points(vcf.stats[c(1:3,5), 'average depth'], 
         vcf.stats[c(1:3,5), 'norm.prop.het'], pch = 19, col = 'red')
  
plot(vcf.stats$prop.het, vcf.stats$norm.prop.het, pch = 19, col = alpha('dodgerblue3', 0.6)) 
  points(vcf.stats[c(1:3,5), 'prop.het'], 
         vcf.stats[c(1:3,5), 'norm.prop.het'], pch = 19, col = 'red')

##### Check effect of missingness on relationship between pedigree and PLINK inbreeding estimates #####
setwd('/Users/Avril/Documents/krat_genetics/data/')
rd <- c('round6')
df.iter <- read.csv(paste0('plink_',rd,'/individual_froh_results_',rd,'.csv'), header = T)

ped.dat <- read.csv('/Users/Avril/Documents/krat_genetics/preseq_sample_information/summary_allstats2_withgen_zeroes.csv')
ped.dat <- ped.dat[ped.dat$id %in% df.iter$id, c('id','inb')]

## for round 5, may want to exclude 4 samples with problematic levels of missingness
excl <- c(5075, 4901, 5018, 5060)

## go with 2 het sites allowed
df.iter <- df.iter[df.iter$phwh == 2,]

for(m in unique(df.iter$phwm)){
  sub <- df.iter[df.iter$phwm == m,]
  sub <- merge(sub, ped.dat, by = 'id')
  plot(sub$inb, sub$froh, pch = 19, col = alpha('dodgerblue3', 0.4), main = paste0('# missing sites = ',m),
       xlab = 'Pedigree inbreeding', ylab = 'PLINK f(ROH)')
    points(sub[sub$id %in% excl, 'inb'], sub[sub$id %in% excl, 'froh'], pch = 19, col = 'red')
}

## doesn't really make a difference in relationship, just shifts f(ROH) values *slightly* for most samples, a lot for 4 samples with high missingness