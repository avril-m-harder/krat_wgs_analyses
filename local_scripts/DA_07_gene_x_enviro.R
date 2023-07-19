#### Attempting to incorporate both genomic and environmental predictor variables to explain individual (female) and pair fitness (lol)
library(scales)
library(ghibli)
library(lme4)
library(lmtest)
setwd('/Users/Avril/Documents/krat_genetics/data/')

##### 1. Read in data #####
##### 1A. Remote sensing/PRISM data #####
#### cell assignments
cells <- read.csv('seq_indivs_mound_locations.csv')

#### TC data 
## landscape-level (i.e., 1 value per metric per year for entire set of occupied/occ + adjacent (buffered) cells)
##### !!! look at summarize_cells_with_mounds.R and summarize here what the different lags mean if anything turns up interesting #####
l.6mo.buffered <- read.csv('/Users/Avril/Documents/krat_remote_sensing/intermediate_data/occ_and_adj_cells_6monthinterval_means.csv')
l.12mo.buffered <- read.csv('/Users/Avril/Documents/krat_remote_sensing/intermediate_data/occ_and_adj_cells_12monthinterval_means_seasoneq.csv')

## cell-level (i.e., 1 value per cell/focal cell + adjacent cells per year/lag combo)
## explanation of seasons in "explanation_of_month_and_weather_based_seasons.png"
c.3mo.buffered <- read.csv('/Users/Avril/Documents/krat_remote_sensing/intermediate_data/archive/cells_w_mounds_sumstats_incl_adjacents.csv')
c.3mo.focalonly <- read.csv('/Users/Avril/Documents/krat_remote_sensing/intermediate_data/archive/cells_w_mounds_sumstats.csv')

#### PRSIM data
## read in monthly data
d.dat <- read.csv('/Users/Avril/Documents/krats/krat_data_and_paper2/prism_data_analysis/22_11_10_daily_data.csv', header=TRUE)

## read in annual data (no lags)
y.dat <- read.csv('/Users/Avril/Documents/krats/krat_data_and_paper2/prism_data_analysis/22_11_10_annual_data.csv', header = TRUE)


##### 1B. Pedigree/genomic data #####
#### all pedigree data
ped.dat <- read.csv('../preseq_sample_information/summary_allstats2_withgen_zeroes.csv')
ped.dat$age.at.death <- ped.dat$deathyear - ped.dat$birthyear
ped.dat[ped.dat$birthyear > 2005, 'age.at.death'] <- NA ## we can't really know when these individuals died

#### relatedness data
## ROHverlap (based on non-pruned SNPs)
rohverlap <- read.csv('corrected_plink_overlap.csv')  ## corrected through removal of overlapping segments with homozygous mismatches between individuals
## NGSrelate2 relatedness results (based on LD-pruned SNPs)
ngs.LD <- read.table('./LD_pruned_data_and_results/ngsrelate_krat_allsamps_LDpruned.res', header = TRUE)
samp.order <- read.table('vcf_sample_order.txt')
samp.order[,2] <- c(0:47)
## replace indices with sample names in NGS output
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'a', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id1'
ngs.LD <- merge(ngs.LD, samp.order, by.x = 'b', by.y = 'V2')
colnames(ngs.LD)[ncol(ngs.LD)] <- 'id2'
ngs.LD <- ngs.LD[order(ngs.LD$a, ngs.LD$b),]

#### mate pair data
relates <- read.csv('pairwise_samp_relationships.csv')
mates <- relates[relates$mates == 1,]

#### sibling data
sib.pairs <- read.csv('/Users/Avril/Documents/krat_genetics/pair_specific_demo_calcs/sibling_pairs.csv')

#### individual inbreeding data
## f(ROH) data
froh <- read.csv('individual_frohs.csv')
## get rid of bcftools results
froh <- froh[,-(grep('pl.', colnames(froh), fixed = TRUE))]
froh <- froh[,-(grep('gt.', colnames(froh), fixed = TRUE))]

##### 1C. Divide data for females out by years in which offspring were produced #####
## include 0's for years when females were alive but no offspring were detected
fem.froh <- froh[froh$sex == 'Female',]  ## n = 27
OUT <- NULL
for(f in unique(fem.froh$id)){
  save <- unlist(fem.froh[fem.froh$id == f, c('id','plink.froh','norm.prop.het'),])
  temp <- ped.dat[ped.dat$momid == f,]
  by <- ped.dat[ped.dat$id == f, 'birthyear']
  dy <- ped.dat[ped.dat$id == f, 'deathyear']
  for(y in by:dy){
    if(y == by){
      age.0 <- 1
    } else{
      age.0 <- 0
    }
    save1 <- c(save, y, age.0,
               length(unique(temp[temp$birthyear == y, 'id'])))
    OUT <- rbind(OUT, save1)
  }
}
## not calculating surviving offspring because too many females were born in 2005 to reasonably evalute that
fem.dat <- as.data.frame(OUT)
colnames(fem.dat) <- c('id','froh','het','repro.year','age.0','n.off') ## repro year = years in which n.off were produced;
                                                                       ## age.0 tells whether mom is/isn't age 0 in that year
table(fem.dat[fem.dat$age.0 == 1, 'n.off']) ## double-checking
fem.dat <- fem.dat[-which(fem.dat$age.0 == 1),]
fem.dat <- fem.dat[,-5]

##### 2. Individual fitness #####
##### 2A. Female perspective: comparing explanatory power across 2 levels of models #####
##### >> I. # offspring ~ f(ROH) #####
## using random effect of repro year
mod <- lmer(fem.dat$n.off ~ fem.dat$froh + (1|fem.dat$repro.year), REML = TRUE)
summary(mod)
ranef(mod)
fixef(mod)
hist(resid(mod))

## using annual PRISM data
temp <- merge(fem.dat, y.dat, by.x = 'repro.year', by.y = 'year')
cols <- c(6:ncol(temp))
for(c in cols){
  red.mod <- lm(temp$n.off ~ temp$froh)
  summary(red.mod)
  ful.mod <- lm(temp$n.off ~ temp$froh + temp[,c])
  summary(ful.mod)
  res <- lrtest(ful.mod, red.mod)
  if(res[2,5] < 0.05){
    print(colnames(temp)[c])
  }
} 

## using remote-sensing data
##### !!! would need additional rs years here if preliminary data look interesting #####
colnames(cells)[4] <- 'repro.year'
temp1 <- merge(fem.dat, cells[,c('id','cell.num','repro.year')], by = c('id', 'repro.year'))
colnames(temp1)[2] <- 'year'

## landscape-level, 12-month
temp2 <- merge(temp1, l.12mo.buffered, by = 'year')
for(i in unique(temp$interval)){
  cols <- grep('mean', colnames(temp), fixed = TRUE)
  temp <- temp2[temp2$interval == i,]
  for(c in cols){
    red.mod <- lm(temp$n.off ~ temp$froh)
    summary(red.mod)
    ful.mod <- lm(temp$n.off ~ temp$froh + temp[,c])
    summary(ful.mod)
    res <- lrtest(ful.mod, red.mod)
    if(res[2,5] < 0.05){
      print(paste0(i,' - ',colnames(temp)[c]))
    }
  }   ## wetness, brightness, temp.k show up (szneq for all 3, without equalization for brightness and temp.k)
}
cols <- which(colnames(temp) %in% c('wetness.mean.szneq','brightness.mean.szneq','temp.k.mean.szneq'))
for(c in cols){
  temp <- temp2[temp2$interval == 3,]
  red.mod <- lm(temp$n.off ~ temp$froh)
  print(summary(red.mod))
  ful.mod <- lm(temp$n.off ~ temp$froh + temp[,c])
  plot(temp[,c], temp$n.off, main = colnames(temp)[c])
  print(summary(ful.mod))
  res <- lrtest(ful.mod, red.mod)
  print(res)
  if(res[2,5] < 0.05){
    print(colnames(temp)[c])
  }
}

##### >> II. # offspring ~ heterozygosity #####
## using annual PRISM data
temp <- merge(fem.dat, y.dat, by.x = 'repro.year', by.y = 'year')
cols <- grep('mean', colnames(temp), fixed = TRUE)
for(c in cols){
  red.mod <- lm(temp$n.off ~ temp$het)
  summary(red.mod)
  ful.mod <- lm(temp$n.off ~ temp$het + temp[,c])
  summary(ful.mod)
  res <- lrtest(ful.mod, red.mod)
  if(res[2,5] < 0.05){
    print(colnames(temp)[c])
  }
} 


##### 3. Pair fitness #####
