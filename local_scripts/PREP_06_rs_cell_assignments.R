library(sp)
library(raster)
library(RStoolbox)
library(viridis)
library(ggplot2)
library(scales)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(TeachingDemos)
source('/Users/Avril/Documents/krat_remote_sensing/krat_remote_sensing_scripts/c2l2readMeta.R')
`%notin%` <- Negate(`%in%`)


##### 1. Read in data #####
##### 1A. Read in population (Peter's) data #####
pop.dat <- read.csv('/Users/Avril/Documents/krat_genetics/preseq_sample_information/KRATP.csv')
## lat/long headings reversed in original file
colnames(pop.dat)[7:8] <- c('long','lat')

##### 1B. Read in more population (Janna's) data #####
j.dat <- read.csv('/Users/Avril/Documents/krats/krat_data_and_paper2/summary_allstats2_withgen_zeroes.csv')
ages <- j.dat[,c('id','birthyear','deathyear')]
ages$age <- ages$deathyear - ages$birthyear
ages <- ages[,c(1,4)]

## read in 1 cropped scene to get background image and cell #s (rows, columns)
ls5.stack <- brick('/Users/Avril/Documents/krat_remote_sensing/C2L2_cropped_landsat45tm_scenes/LT05_L2SP_035038_20020622_20200905_02_T1_CROPPED.grd')
raster::plotRGB(ls5.stack, r=3, g=2, b=1, scale=ls5.stack@data@max[c(3,2,1)], margins=FALSE)

##### 1C. Read in manual mound:cell assignments (n=26) #####
### >> also needed to re-name unique mounds that were all labeled 'R2' in the original database 
man.ass <- read.csv('/Users/Avril/Documents/krat_remote_sensing/sorting_out_mound_names/2021_first_round/manual_mound_cell_assignments.csv')
## rename 'R2' mounds to match their actual locations (unique names R2.1-R2.6)
r2 <- man.ass[man.ass$terr=='R2',]
for(i in 1:nrow(r2)){
  pop.dat[which(pop.dat$lat == r2$lat[i] & pop.dat$long == r2$long[i]), 'terr'] <- r2$new.db.name[i]
}

##### 1D. Read in mound coordinates #####
mnd.locs <- read.csv('/Users/Avril/Documents/krat_remote_sensing/intermediate_data/mound_GPS_coords_n188.csv')
mnd.locs$ID <- 1:nrow(mnd.locs) ## add column for later matching up with TC extracted pixel values
mnd.cells <- mnd.locs[,c(4,5)] ## save ID and db.name for saving cell names later
mnd.locs <- mnd.locs[c('long','lat','database.name')] ## rearrange/rename to match man.locs below
colnames(mnd.locs)[3] <- 'terr'
coordinates(mnd.locs) <- c('long','lat') ## converts to SpatialPointsDataFrame object for plotting
proj4string(mnd.locs) <- CRS("+proj=longlat +datum=WGS84") ## still in degrees
mnd.locs <- sp::spTransform(mnd.locs, CRS("+proj=utm +zone=12 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")) ## now in meters
### manual mound:cell assignments
man.ass <- read.csv('/Users/Avril/Documents/krat_remote_sensing/sorting_out_mound_names/2021_first_round/manual_mound_cell_assignments.csv')
man.cells <- man.ass[,c('terr','new.db.name','cell')]
man.cells[!is.na(man.cells$new.db.name), 'terr'] <- man.cells[!is.na(man.cells$new.db.name), 'new.db.name']
man.cells$ID <- c((max(mnd.cells$ID)+1):(max(mnd.cells$ID)+nrow(man.cells)))
colnames(man.cells)[1] <- 'database.name'
man.cells <- man.cells[,c('database.name','ID')]
man.locs <- as.data.frame(xyFromCell(ls5.stack, man.ass$cell))
man.locs$terr <- man.ass$terr
coordinates(man.locs) <- c('x','y') ## already in meters
proj4string(man.locs)<- CRS("+proj=utm +zone=12 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
all.locs <- rbind(mnd.locs, man.locs)

##### 1E. Read in mound:cell key #####
mc.key <- read.csv(paste0('/Users/Avril/Documents/krat_remote_sensing/C2L2_tc_output_tables/C2L2_mnd_manual_cloudcheck.csv'))
mc.key <- mc.key[,c('database.name','cell.num')]
mc.key <- mc.key[!duplicated(mc.key),]
colnames(mc.key)[1] <- 'terr'

## set extent for all analyses
lo.x <- 663500
hi.x <- 665600
lo.y <- 3497000
hi.y <- 3499750
ext <- extent(lo.x, hi.x, lo.y, hi.y)

##### 1F. Read in list of sequenced individuals #####
seq.samps <- read.table('/Users/Avril/Documents/krat_genetics/data/vcf_sample_order.txt')
seq.samps <- seq.samps$V1

##### 2. Adult location:reproduction #####
## For each adult individual, record mound location and fitness measures for
## each year of that individual's life. Also record juvenile mound location.
## Assumptions:
#### (1) Juvenile (offspring==1) location is natal mound.
#### (2) Adult (offspring==0) location(s) are reproductive mounds for their
####      respective years.
### Requires that juvenile mound != adult mounds.
## ** Not sure that this is any better than my original approach, just based
## ** on assigning individual to mounds as juveniles and nothing else.
OUT <- NULL
for(i in unique(seq.samps)){
  sub <- pop.dat[pop.dat$id == i,]                     ## subset to individual,
  sub <- sub[,c(2,3,6,11)]                             ## keep ID, lifestage, mound, year,
  sub$sex <- j.dat[j.dat$id == i, 'sex']               ## add sex information,
  OUT <- rbind(OUT, sub)
}

locs <- OUT

## add cell numbers to data
locs <- merge(locs, mc.key, by='terr', all.x = TRUE)

write.csv(locs, '/Users/Avril/Documents/krat_genetics/data/seq_indivs_mound_locations.csv', row.names = FALSE)

### for all samples
OUT <- NULL
for(i in unique(pop.dat$id)){
  sub <- pop.dat[pop.dat$id == i,]                     ## subset to individual,
  sub <- sub[,c(2,3,6,11)]                             ## keep ID, lifestage, mound, year,
  sub$sex <- j.dat[j.dat$id == i, 'sex']               ## add sex information,
  OUT <- rbind(OUT, sub)
}

locs <- OUT

## add cell numbers to data
locs <- merge(locs, mc.key, by='terr', all.x = TRUE)
write.csv(locs, '/Users/Avril/Documents/krat_genetics/data/all_indivs_mound_locations.csv', row.names = FALSE)
