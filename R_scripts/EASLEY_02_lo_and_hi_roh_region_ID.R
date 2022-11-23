setwd('/scratch/avrilh/kratroh_09_R_analyses/')
# setwd('/Users/Avril/Documents/krat_genetics/data/')

##### Read in data #####
hi.roh.lens <- read.table('high_roh_region_lengths.txt', sep = '\t', header = FALSE)
hi.roh.lens <- hi.roh.lens[,1]
no.roh.lens <- read.table('no_roh_region_lengths.txt', sep = '\t', header = FALSE)
no.roh.lens <- no.roh.lens[,1]
contigs <- read.table('contigs_easley.txt', sep = '\t', header = TRUE)
gff <- read.table('gff_easley.txt', sep = '\t', header = TRUE)


## set # of iterations to run
iter <- 1000

# HI.GENE.OVERLAP <- NULL   ## columns will be: iteration # / random region # / random region length / length of gene overlapped
# for(ct in 1:iter){
#   print(ct/iter)
#   region.ct <- 1
#   for(len in 1:length(hi.roh.lens)){
#     # print(paste0(ct,' - ',len/length(hi.roh.lens)))
#     l <- hi.roh.lens[len]
#     loc <- sample(1:max(contigs$cumsum), 1)                     ## randomly select a genome location
#     s <- loc
#     e <- loc + l                                                ## calculate end location
#     contig <- contigs[nrow(contigs[which(contigs$cumsum < loc),])+1, 'refseq.contig']
#     while(contigs[contigs$refseq.contig == contig, 'cumsum'] < e){
#       ## if the end position of the contig is less than the calculated end of the ROH region,
#       loc <- sample(1:max(contigs$cumsum), 1)             ## select a different starting location
#       s <- loc
#       e <- loc + l                                        ## and try again
#       contig <- contigs[nrow(contigs[which(contigs$cumsum < loc),])+1, 'refseq.contig']
#     }
#     sub.gff <- gff[gff$seqid == contig,]
#     if(nrow(sub.gff) > 0){
#       ## determine whether genes overlap that region and save gene length(s), length(s) of overlap, and number of genes overlapped
#       ## genes beginning outside of focal ROH region, ending inside
#       if(nrow(sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]) > 0){
#         temp <- sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]
#         for(t in 1:nrow(temp)){                                                                 ## for each overlapping gene,
#           save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
#           HI.GENE.OVERLAP <- rbind(HI.GENE.OVERLAP, save)
#         }
#       }
#       ## genes beginning inside of a focal ROH region, ending outside
#       if(nrow(sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]) > 0){
#         temp <- sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]
#         for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#           save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
#           HI.GENE.OVERLAP <- rbind(HI.GENE.OVERLAP, save)
#         }
#       }
#       ## genes completely covering a focal ROH region
#       if(nrow(sub.gff[sub.gff$start < s & sub.gff$end > e,]) > 0){
#         temp <- sub.gff[sub.gff$start < s & sub.gff$end > e,]
#         for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#           save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
#           HI.GENE.OVERLAP <- rbind(HI.GENE.OVERLAP, save)
#         }
#       }
#       ## genes completely within a focal ROH region
#       if(nrow(sub.gff[sub.gff$start >= s & sub.gff$end <= e,]) > 0){
#         temp <- sub.gff[sub.gff$start >= s & sub.gff$end <= e,]
#         for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
#           save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
#           HI.GENE.OVERLAP <- rbind(HI.GENE.OVERLAP, save)
#         }
#       }
#     }
#     region.ct <- region.ct + 1
#   }
# }
# 
# HI.RAND.STATS <- NULL
# for(c in 1:iter){
#   sub <- HI.GENE.OVERLAP[HI.GENE.OVERLAP[,1] == c, , drop = FALSE]
#   if(nrow(sub) > 0){
#     save <- c(c, sum(sub[,4]), nrow(sub))
#     HI.RAND.STATS <- rbind(HI.RAND.STATS, save)
#   } else{
#     save <- c(c, 0, 0)
#     HI.RAND.STATS <- rbind(HI.RAND.STATS, save)
#   }
# }
# colnames(HI.RAND.STATS) <- c('iteration','total.gene.len','num.genes')
# write.csv(HI.RAND.STATS, paste0('high_roh_randomization_stats_n',iter,'.csv'), row.names = FALSE)
# 
# HI.RAND.STATS <- read.csv(paste0('high_roh_randomization_stats_n',iter,'.csv'))
# 
# hi.roh.gene.overlap <- read.csv('genes_in_high_roh_regions.csv')
# hi.roh.gene.len <- sum(abs(hi.roh.gene.overlap$gene.start - hi.roh.gene.overlap$gene.end)+1)
# pdf('../krat_genetics_scripts/figures_output/high_roh_region_randomization_histograms.pdf', width = 6, height = 5)
# hist(HI.RAND.STATS[,2], xlab = 'total length of overlapped genes', breaks = 15, main = 'gene overlap with randomly selected regions (hi-ROH)',
#      xlim = c(min(HI.RAND.STATS[,2]), 20000000))
#   abline(v = hi.roh.gene.len, lty = 2, col = 'red')
# dev.off()

##### Running this bit on Easley - more regions to test than for high-ROH regions #####
NO.GENE.OVERLAP <- NULL
for(ct in 1:iter){
  print(ct/iter)
  region.ct <- 1
  for(len in 1:length(no.roh.lens)){
    # print(paste0(ct,' - ',len/length(no.roh.lens)))
    l <- no.roh.lens[len]
    loc <- sample(1:max(contigs$cumsum), 1)                     ## randomly select a genome location
    s <- loc
    e <- loc + l                                                ## calculate end location
    contig <- contigs[nrow(contigs[which(contigs$cumsum < loc),])+1, 'refseq.contig']
    while(contigs[contigs$refseq.contig == contig, 'cumsum'] < e){
      ## if the end position of the contig is less than the calculated end of the ROH region,
      loc <- sample(1:max(contigs$cumsum), 1)             ## select a different starting location
      s <- loc
      e <- loc + l                                        ## and try again
      contig <- contigs[nrow(contigs[which(contigs$cumsum < loc),])+1, 'refseq.contig']
    }
    sub.gff <- gff[gff$seqid == contig,]
    if(nrow(sub.gff) > 0){
      ## determine whether genes overlap that region and save gene length(s), length(s) of overlap, and number of genes overlapped
      ## genes beginning outside of focal ROH region, ending inside
      if(nrow(sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]) > 0){
        temp <- sub.gff[sub.gff$start < s & sub.gff$end >= s & sub.gff$end <= e,]
        for(t in 1:nrow(temp)){                                                                 ## for each overlapping gene,
          save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
          NO.GENE.OVERLAP <- rbind(NO.GENE.OVERLAP, save)
        }
      }
      ## genes beginning inside of a focal ROH region, ending outside
      if(nrow(sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]) > 0){
        temp <- sub.gff[sub.gff$start >= s & sub.gff$start <= e & sub.gff$end > e,]
        for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
          save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
          NO.GENE.OVERLAP <- rbind(NO.GENE.OVERLAP, save)
        }
      }
      ## genes completely covering a focal ROH region
      if(nrow(sub.gff[sub.gff$start < s & sub.gff$end > e,]) > 0){
        temp <- sub.gff[sub.gff$start < s & sub.gff$end > e,]
        for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
          save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
          NO.GENE.OVERLAP <- rbind(NO.GENE.OVERLAP, save)
        }
      }
      ## genes completely within a focal ROH region
      if(nrow(sub.gff[sub.gff$start >= s & sub.gff$end <= e,]) > 0){
        temp <- sub.gff[sub.gff$start >= s & sub.gff$end <= e,]
        for(t in 1:nrow(temp)){                                                 ## for each overlapping gene,
          save <- c(ct, region.ct, l, abs(temp$end[t] - temp$start[t])+1)   ## save the gene info and the ROH region info
          NO.GENE.OVERLAP <- rbind(NO.GENE.OVERLAP, save)
        }
      }
    }
    region.ct <- region.ct + 1
  }
  if(ct == 1){
  	write.table(NO.GENE.OVERLAP, paste0('no_roh_randomregion_geneoverlaps_n',iter,'.csv'), quote = FALSE, row.names = FALSE, sep = ',')
  	NO.GENE.OVERLAP <- NULL
  } else{
  	write.table(NO.GENE.OVERLAP, paste0('no_roh_randomregion_geneoverlaps_n',iter,'.csv'), quote = FALSE, sep = ',', row.names = FALSE, append = TRUE, 
                col.names = !file.exists(paste0('no_roh_randomregion_geneoverlaps_n',iter,'.csv')))
    NO.GENE.OVERLAP <- NULL
  }
}

NO.GENE.OVERLAP <- read.csv(paste0('no_roh_randomregion_geneoverlaps_n',iter,'.csv'))

NO.RAND.STATS <- NULL
for(c in 1:iter){
  sub <- NO.GENE.OVERLAP[NO.GENE.OVERLAP[,1] == c,]
  save <- c(c, sum(sub[,4]), nrow(sub))
  NO.RAND.STATS <- rbind(NO.RAND.STATS, save)
}
colnames(NO.RAND.STATS) <- c('iteration','total.gene.len','num.genes')
write.csv(NO.RAND.STATS, paste0('no_roh_randomization_stats_n',iter,'.csv'), row.names = FALSE)

# NO.RAND.STATS <- read.csv(paste0('no_roh_randomization_stats_n',iter,'.csv'))
# no.roh.gene.overlap <- read.csv('genes_in_no_roh_regions.csv')
# no.roh.gene.len <- sum(abs(no.roh.gene.overlap$gene.start - no.roh.gene.overlap$gene.end)+1)
# hist(NO.RAND.STATS[,2], xlab = 'total gene length', breaks = 20, main = 'gene overlap with no-ROH regions', xlim = c(min(RAND.STATS[,2]), hi.roh.gene.len))
#   abline(v = no.roh.gene.len, lty = 2, col = 'red')

