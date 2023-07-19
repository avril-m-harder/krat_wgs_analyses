setwd('/scratch/avrilh/kratroh_09_R_analyses/')

##### Read in data #####
## sample information
relates <- read.csv('pairwise_samp_relationships.csv')
mates <- relates[relates$mates == 1,]

## contig information
contigs <- read.csv('contig_lengths.csv')
contigs <- contigs[order(-contigs$length),]
contigs$c.index <- c(1:nrow(contigs)) ## create contig index, where c.index = 1 == longest contig
tot.len <- sum(contigs$length) ## 2,571,239,112 bp ## length of filtered contigs

## read in VCF sample order and genotype information (collected using bcftools query)
samp.order <- read.table('vcf_sample_order.txt')
samp.order <- unlist(samp.order)

## ID pairs with and without surviving offspring
mates[mates$n.surv > 0, 'surv.or.no'] <- 1
mates[mates$n.surv == 0, 'surv.or.no'] <- 0
mates$pair.id <- c(1:nrow(mates))


##### Calculating per-locus, intra-pair genotype similarities --> PCA or something? #####
gts <- read.table('krat_LDpruned_extracted_genotypes_numericGTs.txt', sep = '\t')
gts <- gts[,-3]
colnames(gts) <- c('contig','pos',samp.order)
gts <- merge(gts, contigs[,c(1,3)], by = 'contig')
gts <- gts[,c(51,2:50)]
for(c in 1:ncol(gts)){
  gts[,c] <- as.numeric(gts[,c])
}  ## missing genotypes go from . --> NA
gts <- as.matrix(gts)


# ### for each pair of mates, calculate per-locus distances and summarize over windows
# w.size <- 10000    ## window size
# s.size <- 10000    ## window step size
# wind.id <- 1       ## starting window ID #
# WIND.DAT <- NULL
# 
# for(p in 1:nrow(mates)){
#   sub.gts <- gts[,c(1, 2, which(colnames(gts) %in% c(mates$id1[p], mates$id2[p])))]
#   sub.gts <- sub.gts[which(!is.na(sub.gts[,3] & !is.na(sub.gts[,4]))),] ## get rid of loci with missing genotype for either individual
#   for(c in unique(contigs$c.index)){
#     print(paste0(p/40,' - ',c/663))
#     
#     c.len <- contigs[contigs$c.index == c, 'length']
#     temp.gts <- sub.gts[sub.gts[,1] == c, , drop = FALSE]
#     s <- 1
#     e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#     
#     while(e <= c.len){
#       wind.gts <- temp.gts[temp.gts[,2] >= s & temp.gts[,2] <= e,, drop = FALSE]
#       if(nrow(wind.gts) > 0){    ## if there's at least 1 SNP in that window,
#         sum(abs(wind.gts[,4] - wind.gts[,3]))
#         ## pair ID, fitness status, c.index, window ID, window start and end positions, 
#         ## # of SNPs in window, and total # of allelic differences
#         wind.dat <- c(mates$pair.id[p], mates$surv.or.no[p], c, wind.id, s, e, nrow(wind.gts), sum(abs(wind.gts[,4] - wind.gts[,3])))
#       } else{       ## if there are 0 SNPs in that window, save that information
#         wind.dat <- c(mates$pair.id[p], mates$surv.or.no[p], c, wind.id, s, e, 0, NA)
#       }
#       # print(paste0(p/40,' - ',c/663,' - ',e/c.len))
#       WIND.DAT <- rbind(WIND.DAT, wind.dat)
#       wind.id <- wind.id + 1
#       if(e != c.len){
#         s <- s + s.size
#         e <- min(c(s + w.size - 1, contigs[contigs$c.index == c, 'length']))
#       } else{
#         e <- e+1
#       }
#     }
#   }
#   if(mates$pair.id[p] == 1){
#     write.table(WIND.DAT, 'intrapair_gen_distances.txt', quote = FALSE, row.names = FALSE)
#     WIND.DAT <- NULL
#   } else{
#     write.table(WIND.DAT, 'intrapair_gen_distances.txt', quote = FALSE, row.names = FALSE, append = TRUE)
#     WIND.DAT <- NULL
#   }
# }

dists <- read.table('intrapair_gen_distances.txt', header = TRUE)
colnames(dists) <- c('pair.id','surv.or.no','c.index','wind.id','start','end','n.snps','n.diffs')
dists <- dists[-which(dists$pair.id == 'V1'),]
for(c in 1:ncol(dists)){
  dists[,c] <- as.numeric(dists[,c])
}
## n.snps = number of SNPs with genotypes for both samples; n.diffs = number of allelic differences within that window
dists$prop.diffs <- dists$n.diffs / (2*dists$n.snps)

## assign window IDs that are shared across samples
OUT <- NULL
OUT1 <- NULL
univ.wind <- 1
for(c in unique(dists$c.index)){
  print(c/663)
  sub <- dists[dists$c.index == c,]
  for(s in unique(sub$start)){
    ## save updated gen dists table
    temp <- sub[sub$start == s,]
    temp <- temp[order(temp$pair.id),]
    temp$univ.wind <- univ.wind
    OUT <- rbind(OUT, temp)
    ## save shared window value and prop.diffs value
    save <- c(univ.wind, temp$prop.diffs)
    OUT1 <- rbind(OUT1, save)
    
    univ.wind <- univ.wind + 1
  }
  if(c == 1){
    write.table(OUT, 'intrapair_gen_distances_with_shared_windowIDs.txt', quote = FALSE, row.names = FALSE)
    write.table(OUT1, 'pcaformatted_gendists_shared_windowIDs.txt', quote = FALSE, row.names = FALSE)
    
    OUT <- NULL
    OUT1 <- NULL
  } else {
    write.table(OUT, 'intrapair_gen_distances_with_shared_windowIDs.txt', quote = FALSE, row.names = FALSE, append = TRUE, 
                col.names = !file.exists('intrapair_gen_distances_with_shared_windowIDs.txt'))
    write.table(OUT1, 'pcaformatted_gendists_shared_windowIDs.txt', quote = FALSE, row.names = FALSE, append = TRUE,
                col.names = !file.exists('pcaformatted_gendists_shared_windowIDs.txt'))
    
    OUT <- NULL
    OUT1 <- NULL
  }
}
