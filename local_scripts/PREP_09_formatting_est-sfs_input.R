setwd('/Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/')
`%notin%` <- Negate(`%in%`)
library(vcfR)

##### Read in spectabilis VCF and format allele counts #####
## get contig info table
contigs <- read.table('../contigs_easley.txt', header = TRUE)

## for est-sfs, focal species needs to have the same count for each locus
# spe <- read.vcfR('krat_LDpruned_zeromissingness.recode.vcf.gz')
# temp <- as.data.frame(extract.gt(spe, element = 'GT', return.alleles = TRUE))
# rm(spe)
# ## don't care about which individual has which genotype, just allele count totals
# ## split genotypes into alleles
# OUT <- NULL
# for(c in 1:ncol(temp)){
#   print(c)
#   sub <- do.call(rbind, strsplit(temp[,c], "/ | |"))[,c(1,3)]
#   OUT <- cbind(OUT, sub)
# }
# rm(sub)
# ## and count each allele;
# ## order for est-sfs = A, C, G, T
# OUT1 <- NULL
# for(r in 1:nrow(OUT)){
#   if(r %% 10000 == 0){
#     print(r/nrow(OUT))
#     if(r == 10000){
#       colnames(OUT1) <- c('a','c','g','t')
#       write.table(OUT1, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE)
#       OUT1 <- NULL
#     } else{
#       write.table(OUT1, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE, append = TRUE,
#                   col.names=!file.exists('/Users/Avril/Desktop/temp.csv'))
#       OUT1 <- NULL
#     }
#   }
#   a <- length(grep('A', OUT[r,]))
#   c <- length(grep('C', OUT[r,]))
#   g <- length(grep('G', OUT[r,]))
#   t <- length(grep('T', OUT[r,]))
#   save <- c(a,c,g,t)
#   OUT1 <- rbind(OUT1, save)
# }
# write.table(OUT1, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE, append = TRUE,
#             col.names=!file.exists('/Users/Avril/Desktop/temp.csv'))
# OUT1 <- read.csv('/Users/Avril/Desktop/temp.csv', header = TRUE)
# all.cts <- as.data.frame(OUT1)
# all.cts$chrom <- do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,1]
# all.cts$pos <- do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,2]
# rm(OUT, OUT1, temp)
# all.cts <- merge(all.cts, contigs[,c('c.index','contig')], by.x = 'chrom', by.y = 'contig')
# all.cts <- all.cts[,c('c.index','pos','a','c','g','t')]
# write.csv(all.cts, 'dspec_LDpruned_nomissingness_allelecounts.csv', row.names = FALSE)
all.cts <- read.csv('dspec_LDpruned_nomissingness_allelecounts.csv')

##### Read in ordii/stephensi data, check for shared loci, format allele counts #####
## ordii
ord <- read.vcfR('d_ordii_refseqnames_DspecLDprunedSNPs.vcf.gz')
temp <- as.data.frame(extract.gt(ord, element = 'GT', return.alleles = TRUE))
temp$chrom <- paste0(do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,1],'_',do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,2])
temp$pos <- do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,3]
colnames(temp)[1] <- 'gt'
temp$gt.1 <- do.call(rbind, strsplit(temp$gt, split = '/', fixed = TRUE))[,1]
temp$gt.2 <- do.call(rbind, strsplit(temp$gt, split = '/', fixed = TRUE))[,2]
temp <- temp[,-1]
temp[which(temp$gt.1 != temp$gt.2), 'het'] <- 1
temp[which(temp$gt.1 == temp$gt.2), 'het'] <- 0
table(temp$het) ## for ordii, only 22 kb / 2.7 Mb heterozygous
ord <- temp

## stephensi
ste <- read.vcfR('d_stephensi_refseqnames_DspecLDprunedSNPs.vcf.gz')
temp <- as.data.frame(extract.gt(ste, element = 'GT', return.alleles = TRUE))
temp$chrom <- paste0(do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,1],'_',do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,2])
temp$pos <- do.call(rbind, strsplit(rownames(temp), split = '_', fixed = TRUE))[,3]
colnames(temp)[1] <- 'gt'
temp$gt.1 <- do.call(rbind, strsplit(temp$gt, split = '/', fixed = TRUE))[,1]
temp$gt.2 <- do.call(rbind, strsplit(temp$gt, split = '/', fixed = TRUE))[,2]
temp <- temp[,-1]
temp[which(temp$gt.1 != temp$gt.2), 'het'] <- 1
temp[which(temp$gt.1 == temp$gt.2), 'het'] <- 0
table(temp$het) ## for ordii, only 32 kb / 2.7 Mb heterozygous
ste <- temp
rm(temp)

## replace refseq contig names with numeric contig index #'s
ord <- merge(ord, contigs[,c('c.index','refseq.contig')], by.x = 'chrom', by.y = 'refseq.contig')
ord <- ord[,c('c.index','pos','gt.1','gt.2','het')]
ord <- ord[ord$het == 0,]
ste <- merge(ste, contigs[,c('c.index','refseq.contig')], by.x = 'chrom', by.y = 'refseq.contig')
ste <- ste[,c('c.index','pos','gt.1','gt.2','het')]
ste <- ste[ste$het == 0,]

## identify loci covered in all 3 species & subset to those
ord$loc <- paste0(ord$c.index,'-',ord$pos)
ste$loc <- paste0(ste$c.index,'-',ste$pos)
all.cts$loc <- paste0(all.cts$c.index,'-',all.cts$pos)
keep.locs <- c(ord$loc, ste$loc, all.cts$loc)
keep.locs <- names(table(keep.locs)[which(table(keep.locs) == 3)])
ord <- ord[ord$loc %in% keep.locs,]
ste <- ste[ste$loc %in% keep.locs,]
all.cts <- all.cts[all.cts$loc %in% keep.locs,]
all(all.cts$loc == ord$loc)
all(all.cts$loc == ste$loc)
## 901,416 SNPs in final filtered set

## sum alleles for ord and ste
OUT <- NULL
for(r in 1:nrow(ord)){
  if(r %% 10000 == 0){
    if(r == 10000){
      colnames(OUT) <- c('a','c','g','t')
      write.table(OUT, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE)
      OUT <- NULL
    } else{
      write.table(OUT, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE, append = TRUE,
                  col.names=!file.exists('/Users/Avril/Desktop/temp.csv'))
      OUT <- NULL
    }
  }
  a <- length(grep('A', ord[r,c(3)]))
  c <- length(grep('C', ord[r,c(3)]))
  g <- length(grep('G', ord[r,c(3)]))
  t <- length(grep('T', ord[r,c(3)]))
  save <- c(a,c,g,t)
  OUT <- rbind(OUT, save)
}
write.table(OUT, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE, append = TRUE,
            col.names=!file.exists('/Users/Avril/Desktop/temp.csv'))
temp <- read.csv('/Users/Avril/Desktop/temp.csv', header = TRUE)
ord$a <- temp[,1]
ord$c <- temp[,2]
ord$g <- temp[,3]
ord$t <- temp[,4]

OUT <- NULL
for(r in 1:nrow(ste)){
  if(r %% 10000 == 0){
    if(r == 10000){
      colnames(OUT) <- c('a','c','g','t')
      write.table(OUT, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE)
      OUT <- NULL
    } else{
      write.table(OUT, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE, append = TRUE,
                  col.names=!file.exists('/Users/Avril/Desktop/temp.csv'))
      OUT <- NULL
    }
  }
  a <- length(grep('A', ste[r,c(3)]))
  c <- length(grep('C', ste[r,c(3)]))
  g <- length(grep('G', ste[r,c(3)]))
  t <- length(grep('T', ste[r,c(3)]))
  save <- c(a,c,g,t)
  OUT <- rbind(OUT, save)
}
write.table(OUT, file = '/Users/Avril/Desktop/temp.csv', sep = ',', row.names = FALSE, append = TRUE,
            col.names=!file.exists('/Users/Avril/Desktop/temp.csv'))
temp <- read.csv('/Users/Avril/Desktop/temp.csv', header = TRUE)
ste$a <- temp[,1]
ste$c <- temp[,2]
ste$g <- temp[,3]
ste$t <- temp[,4]

##### Prep and check 3 species data and write #####
all.cts <- all.cts[order(all.cts$c.index, all.cts$pos),]
ord <- ord[order(ord$c.index, ord$pos),]
ste <- ste[order(ste$c.index, ste$pos),]
all(all.cts$loc == ord$loc)
all(all.cts$loc == ste$loc)
table(rowSums(all.cts[,c('a','c','g','t')]))
table(rowSums(ord[,c('a','c','g','t')]))
table(rowSums(ste[,c('a','c','g','t')]))

one <- paste0(all.cts$a,',',all.cts$c,',',all.cts$g,',',all.cts$t)
two <- paste0(ord$a,',',ord$c,',',ord$g,',',ord$t)
three <- paste0(ste$a,',',ste$c,',',ste$g,',',ste$t)
final <- cbind(one, two, three)
write.table(final, 'est_sfs_input_spe_ord_ste.txt', sep = '\t', quote = FALSE, row.names = FALSE)
final <- cbind(all.cts$c.index, all.cts$pos, one, two, three)
colnames(final) <- c('c.index','pos','spe','ord','ste')
write.table(final, 'est_sfs_input_spe_ord_ste_wlocs.txt', sep = '\t', quote = FALSE, row.names = FALSE)
