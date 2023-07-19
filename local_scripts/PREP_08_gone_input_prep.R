setwd('/Users/Avril/Documents/krat_genetics/data/')
`%notin%` <- Negate(`%in%`)

## get contig info table
contigs <- read.table('contigs_easley.txt', header = TRUE)

## sample names
samps <- read.table('vcf_sample_order.txt')
samps <- samps$V1

## pedigree data
ped.dat <- read.csv('summary_allstats2_withgen_zeroes.csv')
ped.dat <- ped.dat[,c('id','inb','off','off_survive','birthyear')]
ped.dat <- ped.dat[ped.dat$id %in% samps,]

## relationship data
relates <- read.csv('pairwise_samp_relationships.csv')

##### Manual changes to pedigree-defined relationships based on genetic information (taken from DA_03.R) #####
## ** see genetic_pedigree_mismatches.pptx and genetic_pedigree_mismatches_NOTES.txt for details **
## low probability of both parents correct (0.5525 for 4915)
relates[which((relates$id1 == 4915 & relates$id2 == 4795) |
                (relates$id2 == 4915 & relates$id1 == 4795)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4915 & relates$id2 == 4796) |
                (relates$id2 == 4915 & relates$id1 == 4796)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4915 & relates$id2 == 4952) |
                (relates$id2 == 4915 & relates$id1 == 4952)), 'half.sibs'] <- 0

## all 3 genetic estimates are low, even though pedigree probz aren't awful
relates[which((relates$id1 == 4952 & relates$id2 == 4795) |
                (relates$id2 == 4952 & relates$id1 == 4795)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4952 & relates$id2 == 4796) |
                (relates$id2 == 4952 & relates$id1 == 4796)), 'half.sibs'] <- 0

## low probability for shared parent (father), low genetic estimates, f(ROH)verlap inconclusive
relates[which((relates$id1 == 4910 & relates$id2 == 4897) |
                (relates$id2 == 4910 & relates$id1 == 4897)), 'half.sibs'] <- 0

## all 3 genetic estimates are low
relates[which((relates$id1 == 4901 & relates$id2 == 4943) |
                (relates$id2 == 4901 & relates$id1 == 4943)), 'half.sibs'] <- 0
relates[which((relates$id1 == 4901 & relates$id2 == 4962) |
                (relates$id2 == 4901 & relates$id1 == 4962)), 'half.sibs'] <- 0

relates[which((relates$id1 == 4976 & relates$id2 == 5054) |
                (relates$id2 == 4976 & relates$id1 == 5054)), 'full.sibs'] <- 0


## Questionable reassignments because marginal posterior probz are pretty high for these two sets in the pedigree;
## fROHverlap low here, too
relates[which((relates$id1 == 4962 & relates$id2 == 4943) |
                (relates$id2 == 4962 & relates$id1 == 4943)), 'full.sibs'] <- 0

##### Need to generate multiple lists of randomly selected individuals that are not full-sibs nor parent-offspring #####
## Ideally, also not half-sibs, but idk that we're gonna have enough samples for that.
## May want to create a bunch of sample sets from birthyear == 2005 to see how contemporaneous samples need to be /
## effect of sampling individuals across multiple birthyears
imp.relates <- relates[which(relates$parent.off == 1 | relates$full.sibs == 1 | relates$half.sibs == 1),]

## start by generating sets of individuals born in 2005
b.2005 <- ped.dat[ped.dat$birthyear == 2005,]

## what is the maximum sample size that can be taken from 2005 where individuals are not related? looks like ~14
temp <- imp.relates[which(imp.relates$id1 %in% b.2005$id & imp.relates$id2 %in% b.2005$id),]
samps <- b.2005$id[b.2005$id %notin% c(temp$id1, temp$id2)] ## 8 individuals do not have documented 
                                                            ## relationships to others born in 2005
sort(table(c(temp$id1, temp$id2))) ## most only have 1 or 2 relationships, 4 have 4, and 1 has 6

OUT <- NULL
OUT1 <- NULL
for(c in 1:250){
  print(c)
  tbd <- as.numeric(names(sort(table(c(temp$id1, temp$id2)))))
  samps <- b.2005$id[b.2005$id %notin% c(temp$id1, temp$id2)]
  while(length(tbd) >= 1){
    id <- sample(tbd, 1)
    samps <- c(samps, id)
    ids <- unique(c(temp[which(temp$id1 == id | temp$id2 == id), 'id1'], temp[which(temp$id1 == id | temp$id2 == id), 'id2']))
    ids <- ids[which(ids %notin% id)]
    tbd <- tbd[tbd %notin% c(id,ids)]
  }
  if(length(samps) < 23){
    save <- c(c, length(samps))
    OUT <- rbind(OUT, save)
    OUT1 <- rbind(OUT1, sample(samps, 12, replace = FALSE))
  }
}

## OUT1 == lists of 12 ~unrelated individuals born in 2005 (each row = 1 set)

full.names <- read.table('full_sample_names.txt')
full.names <- cbind(full.names, do.call(rbind, strsplit(full.names$V1, split = '_', fixed = TRUE))[,1])
colnames(full.names) <- c('full','id')

# OUT1[1,]
# write.table(full.names[full.names$id %in% OUT1[1,], 'full'], '/Users/Avril/Desktop/go_test_sample_list.txt', quote = FALSE,
#             row.names = FALSE)

for(c in 1:100){
  write.table(full.names[full.names$id %in% OUT1[c,], 'full'], paste0('/Users/Avril/Desktop/temp/',c,'.txt'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
}
