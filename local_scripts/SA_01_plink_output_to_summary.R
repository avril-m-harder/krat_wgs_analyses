library(scales)

## set round and data set
rd <- 'round6'

### Use Globus to download *.hom files for round, add to appropriate directory, then run below
setwd(paste0('/Users/Avril/Documents/krat_genetics/data/plink_',rd,'/',rd,'_hom_files/'))
fns <- list.files()
  
##### Summarize all coordinates in 1 file #####
z <- 1
for(fn in fns){
  print(z/length(fns))
  ## read in data
  dat <- read.table(fn, header=TRUE)
  
  if(nrow(dat) > 0){
    ## split file name to get parameter settings
    phwh.ns <- do.call(rbind, strsplit(fn, split='phwh_', fixed=TRUE))[,2]
    a <- unique(do.call(rbind, strsplit(phwh.ns, split='_', fixed=TRUE))[,1])
  
    phwm.ns <- do.call(rbind, strsplit(fn, split='phwm_', fixed=TRUE))[,2]
    b <- unique(do.call(rbind, strsplit(phwm.ns, split='_', fixed=TRUE))[,1])
  
    phws.ns <- do.call(rbind, strsplit(fn, split='phws_', fixed=TRUE))[,2]
    c <- unique(do.call(rbind, strsplit(phws.ns, split='_', fixed=TRUE))[,1])
  
    phzd.ns <- do.call(rbind, strsplit(fn, split='phzd_', fixed=TRUE))[,2]
    d <- unique(do.call(rbind, strsplit(phzd.ns, split='_', fixed=TRUE))[,1])
  
    phzg.ns <- do.call(rbind, strsplit(fn, split='phzg_', fixed=TRUE))[,2]
    e <- unique(do.call(rbind, strsplit(phzg.ns, split='_', fixed=TRUE))[,1])
  
    phwt.ns <- do.call(rbind, strsplit(fn, split='phwt_', fixed=TRUE))[,2]
    f <- unique(do.call(rbind, strsplit(phwt.ns, split='_', fixed=TRUE))[,1])
  
    phzs.ns <- do.call(rbind, strsplit(fn, split='phzs_', fixed=TRUE))[,2]
    g <- unique(do.call(rbind, strsplit(phzs.ns, split='_', fixed=TRUE))[,1])
  
    phzk.ns <- do.call(rbind, strsplit(fn, split='phzk_', fixed=TRUE))[,2]
    h <- unique(do.call(rbind, strsplit(phzk.ns, split='_', fixed=TRUE))[,1])
    h <- gsub('.hom', '' ,h)
    
    ## prep data to keep and write
    dat <- dat[,c(1,4,7,8,10)]
    dat$phwh <- a
    dat$phwm <- b
    dat$phws <- c
    dat$phzd <- d
    dat$phzg <- e
    dat$phwt <- f
    dat$phzs <- g
    dat$phzk <- h
    
    if(fn == fns[1]){
      write.table(dat, paste0('/Users/Avril/Documents/krat_genetics/data/plink_',rd,'/PLINK_all_coordinates_',rd,'.txt'),
                  quote = FALSE, row.names = FALSE, sep='\t', col.names = FALSE)
    } else{
      write.table(dat, paste0('/Users/Avril/Documents/krat_genetics/data/plink_',rd,'/PLINK_all_coordinates_',rd,'.txt'),
                  append = TRUE, quote = FALSE, row.names = FALSE, sep='\t', col.names = FALSE)
    }
    z <- z+1
  }
}


##### Summarize individual f(ROH) data #####
## contig info
contigs <- read.csv('../../contig_lengths.csv')
tot.len <- sum(contigs$length) ## 2,571,239,112 bp

plink.out <- read.table(paste0('../PLINK_all_coordinates_',rd,'.txt'))
colnames(plink.out) <- c('id','contig','start','end','n.snps','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk')

OUT <- NULL
write.table(OUT, paste0('../individual_froh_results_',rd,'.txt'),
            sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

z <- 1
for(a in unique(plink.out$phwh)){
  for(b in unique(plink.out$phwm)){
    for(c in unique(plink.out$phws)){
      for(d in unique(plink.out$phzd)){
        for(e in unique(plink.out$phzg)){
          for(f in unique(plink.out$phwt)){
            for(g in unique(plink.out$phzs)){
              for(h in unique(plink.out$phzk)){
                sub <- plink.out[plink.out$phwh == a &
                                   plink.out$phwm == b & plink.out$phws == c & plink.out$phzd == d &
                                   plink.out$phzg == e & plink.out$phwt == f & plink.out$phzs == g &
                                   plink.out$phzk == h,]
                if(nrow(sub) > 0){
                  for(i in unique(sub$id)){
                    temp <- sub[sub$id == i,]
                    temp$length <- temp$end - temp$start + 1
                    froh <- sum(temp$length)/tot.len
                    save <- c(i, a, b, c, d, e, f, g, h, froh)
                    OUT <- rbind(OUT, save)
                  }
                }
                write.table(OUT, paste0('../individual_froh_results_',rd,'.txt'),
                            sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
                OUT <- NULL
                z <- z+1
                print(z/length(fns))
              }
            }
          }
        }
      }
    }
  }
}

OUT <- read.table(paste0('../individual_froh_results_',rd,'.txt'))
colnames(OUT) <- c('id','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk','froh')
froh.stats <- as.data.frame(OUT)
froh.stats <- froh.stats[,c('id','froh','phwh','phwm','phws','phzd','phzg','phwt','phzs','phzk')]
write.csv(froh.stats, paste0('../individual_froh_results_',rd,'.csv'), row.names = FALSE)
froh.stats <- read.csv(paste0('../individual_froh_results_',rd,'.csv'))

