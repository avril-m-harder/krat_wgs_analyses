setwd('/Users/Avril/Documents/krat_genetics/preseq_sample_information/')

library(scales)

samp.dat <- read.csv('../preseq_sample_information/summary_allstats2_withgen_zeroes.csv')
ah.samps <- read.table('seqd_sample_ids.txt')
ah.samps <- ah.samps$V1
gl.samps <- read.table('../gina_data/ginas_sample_ids.txt')
gl.samps <- gl.samps$V1

gl.samps[gl.samps %in% ah.samps] ## no overlap
inds <- sort(c(ah.samps, gl.samps))
ind.dat <- samp.dat[samp.dat$id %in% inds,]
ind.dat[ind.dat$id %in% ah.samps, 'seqrd'] <- 1
ind.dat[ind.dat$id %in% gl.samps, 'seqrd'] <- 2

## output format will be:
## id 1 / id 1 sex / id 2 / id 2 sex / 
## parent-offspring (yes/no) / mates (yes/no) / full sibs (yes/no) / half-sibs (yes/no) / 
## num. off / num. surviving off
## gonna go through and ID ALL relationships first, then only keep the most informative lines (e.g., if a relationship is kept for mating with a half-sib parent and for a half-sib relationship, keep the line describing the first bit)

##### for each individual: #####
OUT <- NULL
save <- rep(NA, 12)
for(i in inds){ ## for each individual
  print(i)
  focal <- samp.dat[samp.dat$id == i,] ## save its info
  
  ## compare against each other individual sequenced and identify all of their relationships
  comps <- inds[which(inds != i)]
  for(c in comps){
    comp <- samp.dat[samp.dat$id == c,]
    save[1] <- i   ## save focal id
    if(focal$sex == 'Male'){
      save[2] <- 0
    } else{
      save[2] <- 1
    }
    if(focal$id %in% ah.samps){
      save[3] <- 1
    } else{
      save[3] <- 2
    }
    save[4] <- c   ## save comp individual id
    if(comp$sex == 'Male'){
      save[5] <- 0
    } else{
      save[5] <- 1
    }
    if(c %in% ah.samps){
      save[6] <- 1
    } else{
      save[6] <- 2
    }
    ## check for parent-offspring relationship
    if(focal$momid == c | focal$dadid == c | comp$momid == i | comp$dadid == i){
      save[7] <- 1
    } else{
      save[7] <- 0
    }
    ## check for mate relationship and if found, save offspring information
    if(nrow(samp.dat[(samp.dat$momid == i & samp.dat$dadid == c) | 
                     (samp.dat$dadid == i & samp.dat$momid == c),]) > 0){
      sub <- samp.dat[(samp.dat$momid == i & samp.dat$dadid == c) | (samp.dat$dadid == i & samp.dat$momid == c),]
      save[8] <- 1
      save[11] <- length(unique(sub$id))
      save[12] <- length(unique(sub[which(sub$deathyear - sub$birthyear > 0), 'id']))
    } else{
      save[8] <- 0
      save[11] <- 0
      save[12] <- 0
    }
    ## check for full-sib relationship
    if(focal$momid == comp$momid & focal$dadid == comp$dadid & focal$momid != 0 & focal$dadid != 0){
      save[9] <- 1
    } else{
      save[9] <- 0
    }
    ## check for half-sib relationship
    if((focal$momid == comp$momid & focal$dadid != comp$dadid & focal$momid != 0 & comp$momid != 0) |
       (focal$momid != comp$momid & focal$dadid == comp$dadid & focal$dadid != 0 & comp$dadid != 0)){
      save[10] <- 1
    } else{
      save[10] <- 0
    }
    OUT <- rbind(OUT, save)
    save <- rep(NA, 12)
  }
}

colnames(OUT) <- c('id1','id1.sex','id1.seqrd','id2','id2.sex','id2.seqrd','parent.off','mates','full.sibs','half.sibs','n.offspring','n.surv.offspring')
relates <- as.data.frame(OUT)
for(c in 1:ncol(relates)){
  relates[,c] <- as.numeric(relates[,c])
}

## get rid of rows with no relationship information
## totals to sum
relates <- relates[rowSums(relates[,c(7:12)]) > 0,]

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

write.csv(relates, 'n97_pairwise_samp_relationships.csv', row.names = FALSE)
write.csv(ind.dat, 'n97_sample_demo_information.csv', row.names = FALSE)

## some plots
plot(table(ind.dat$birthyear), lwd = 3, xlab = 'Birth year', ylab = 'Count')
table(ind.dat$sex) ## 52 females, 45 males
hist(ind.dat[ind.dat$sex == 'Female', 'off'], xlab = '# offspring')
hist(ind.dat[ind.dat$sex == 'Female', 'off_survive'], xlab = '# offspring surviving')
hist(ind.dat[ind.dat$sex == 'Female', 'inb'], xlab = 'Inbreeding', breaks = 20)
