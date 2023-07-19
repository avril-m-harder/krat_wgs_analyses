setwd('/Users/Avril/Documents/krat_genetics/data/')

library(scales)

## sample ROH data (bcftools/ROH + final PLINK settings)
pl.roh <- read.table('krat_final_allfiltcontigs_all_samps_PL_RG_ONLY.txt', header = FALSE)
colnames(pl.roh) <- c('RG','id','contig','start','end','length','n.snps','quality')
pl.roh <- pl.roh[pl.roh$length >= 100000,]
pl.roh <- pl.roh[,-1]
pl.roh$id <- do.call(rbind, strsplit(pl.roh$id, split = '_'))[,1]

inds <- unique(pl.roh$id) 
samp.dat <- read.csv('../preseq_sample_information/summary_allstats2_withgen_zeroes.csv')

## output format will be:
## id 1 / id 2 / parent-offspring (yes/no) / mates (yes/no) / full sibs (yes/no) / half-sibs (yes/no) / num. off / num. surviving off
## gonna go through and ID ALL relationships first, then only keep the most informative lines (e.g., if a relationship is kept for mating with a half-sib parent and for a half-sib relationship, keep the line describing the first bit)

##### for each individual: #####
OUT <- NULL
save <- rep(NA, 8)
for(i in inds){ ## for each individual
  print(i)
  focal <- samp.dat[samp.dat$id == i,] ## save its info
  
  ## compare against each other individual sequenced and identify all of their relationships
  comps <- inds[which(inds != i)]
  for(c in comps){
    comp <- samp.dat[samp.dat$id == c,]
    save[1] <- i   ## save focal id
    save[2] <- c   ## save comp individual id
    ## check for parent-offspring relationship
    if(focal$momid == c | focal$dadid == c | comp$momid == i | comp$dadid == i){
      save[3] <- 1
    } else{
      save[3] <- 0
    }
    ## check for mate relationship and if found, save offspring information
    if(nrow(samp.dat[(samp.dat$momid == i & samp.dat$dadid == c) | 
                (samp.dat$dadid == i & samp.dat$momid == c),]) > 0){
      sub <- samp.dat[(samp.dat$momid == i & samp.dat$dadid == c) | (samp.dat$dadid == i & samp.dat$momid == c),]
      save[4] <- 1
      save[7] <- length(unique(sub$id))
      save[8] <- length(unique(sub[which(sub$deathyear - sub$birthyear > 0), 'id']))
    } else{
      save[4] <- 0
      save[7] <- 0
      save[8] <- 0
    }
    ## check for full-sib relationship
    if(focal$momid == comp$momid & focal$dadid == comp$dadid & focal$momid != 0 & focal$dadid != 0){
      save[5] <- 1
    } else{
      save[5] <- 0
    }
    ## check for half-sib relationship
    if((focal$momid == comp$momid & focal$dadid != comp$dadid & focal$momid != 0 & comp$momid != 0) |
       (focal$momid != comp$momid & focal$dadid == comp$dadid & focal$dadid != 0 & comp$dadid != 0)){
      save[6] <- 1
    } else{
      save[6] <- 0
    }
    OUT <- rbind(OUT, save)
    save <- rep(NA, 8)
  }
}

colnames(OUT) <- c('id1','id2','parent.off','mates','full.sibs','half.sibs','n.offspring','n.surv.offspring')
relates <- as.data.frame(OUT)
for(c in 1:ncol(relates)){
  relates[,c] <- as.numeric(relates[,c])
}

## get rid of rows with no relationship information
# nrow(relates[rowSums(relates[,c(3:8)]) == 0,])
# nrow(relates[rowSums(relates[,c(3:8)]) > 0,])
## totals to sum
relates <- relates[rowSums(relates[,c(3:8)]) > 0,]

## get rid of redundancy in sample pairs
for(r in 1:nrow(relates)){
  if(relates$id1[r] < relates$id2[r]){
    relates$pair.id[r] <- paste0(relates$id1[r],'-',relates$id2[r])
  } else{
    relates$pair.id[r] <- paste0(relates$id2[r],'-',relates$id1[r])
  }
}
relates <- relates[!duplicated(relates$pair.id),]
relates <- relates[,-ncol(relates)]

write.csv(relates, 'pairwise_samp_relationships.csv', row.names = FALSE)
