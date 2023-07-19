###########################################################################
###                          Samarth Mathur, PhD                        ###
###                       The Ohio State University                     ###
###                                                                     ###
###     Date Created: 07/11/22                  Last Modified: 07/11/22 ###
###########################################################################
###########################################################################
###                   ROH_sensitivity.R  		                            ###
###########################################################################

#### PREREQUISITES #####
# load packages
library(ggplot2)
library(plyr)
library(dplyr)
library(DataCombine)
library(sensitivity)
library(reshape2)
`%notin%` <- Negate(`%in%`)

## color settings (1 per parameter)
library(NatParksPalettes)
pal <- natparks.pals(name="DeathValley", n=7, type="discrete")
cols <- cbind(pal[c(1:3,5:7)], c('phwh','phwm','phws','phzg','phwt','phzs'))
colnames(cols) <- c('col','Variables')

# Good intro about sensitivity analysis: 
#https://complementarytraining.net/simple-sensitivity-analysis-with-r/


##### Empirical Data #####
setwd('/Users/Avril/Documents/krat_genetics/data/')

### loop over rounds/iterations (i.e., different parameter setting combinations) and analyze results 

## define rounds to be analyzed
rds <- c('round5')

## for round 5, want to exclude 4 samples with problematic levels of missingness
# excl <- c(5075, 4901, 5018, 5060)

for(rd in rds){
  df.iter <- read.csv(paste0('plink_',rd,'/individual_froh_results_',rd,'.csv'), header = T)
  # df.iter <- df.iter[df.iter$id %notin% excl,]
  
  ## write a table of all the parameter settings
  param.file <- file(paste0('plink_',rd,'/param_settings_',rd,'.txt'))
  for(c in 3:ncol(df.iter)){
    writeLines(c(paste0(colnames(df.iter[c]),': '), unlist(unique(df.iter[,c]))), param.file)
  }
  close(param.file)
  
  sink(paste0('plink_',rd,'/param_settings_',rd,'.txt'))
  for(c in 3:ncol(df.iter)){
    cat(paste0(colnames(df.iter[c]),':\n'))
    cat(unlist(unique(df.iter[,c])))
    cat('\n\n')
  }
  sink()
  
  ## phzk and phzd not varied, remove from analysis
  df.iter <- df.iter[, -c(which(colnames(df.iter) %in% c('phzk','phzd')))]
  
  # Mean and SD for each ind (calculated over all settings combos with ROH calls for each ind)
  mean.iter1 <- ddply(df.iter, "id", summarise, meanCall=mean(froh), sdCall=sd(froh))
    
  # Standardized Rank Regression Coefficients
  inds <- unique(df.iter$id)
  
  srcOriginal1 <- NULL
  srcError1 <- NULL
  srcBias1 <- NULL
  for (ind in inds)
  {
    df <- df.iter[which(df.iter$id ==ind),]
    targets <- df[,2]
    ## detect predictors for each round based on which are variable
    CS <- NULL
    for(p in c(3:ncol(df))){
      if(length(unique(df[,p])) > 1){
        CS <- c(CS, p)
      }
    }
    predictors <- df[,c(CS)]
    srcInd <- sensitivity::src(predictors, targets, rank = T, logistic = TRUE, nboot = 100, conf = 0.95)
    df3 <- c(ind,as.array(t(srcInd$SRC))[1,])
    df4 <- c(ind,as.array(t(srcInd$SRC))[3,])
    df5 <- c(ind,as.array(t(srcInd$SRC))[2,])
    srcOriginal1 <- rbind(srcOriginal1,df3)
    srcError1 <- rbind(srcError1,df4)
    srcBias1 <- rbind(srcBias1,df5)
  }
  
  srcOriginal1 <- as.data.frame(srcOriginal1)
  srcError1 <- as.data.frame(srcError1)
  srcBias1 <- as.data.frame(srcBias1)
  
  colnames(srcOriginal1)[1] <- "id"
  colnames(srcBias1)[1] <- "id"
  
  # Plot SRCs for each variable
  origMelt1 <- melt(data = srcOriginal1, id.vars="id",
                    measure.vars = colnames(srcOriginal1)[colnames(srcOriginal1) %notin% c('id')],
                    variable.name = "Variables",
                    value.name = "SRC")
  
  ### base R figure
  origMelt1 <- merge(origMelt1, cols, by = 'Variables')
  origMelt1$SRC <- as.numeric(origMelt1$SRC)
  max <- length(unique(origMelt1$Variables))
  text.size <- 1.25
  shrink <- 400 ## higher # here ==> narrower x-direction spread for points
  alph <- 0.4
  pt.size <- 0.8
  xmin <- 0.85
  xmax <- max + 0.15
  
  if(length(unique(origMelt1$Variables)) > 2){
    pdf(paste0('../krat_genetics_scripts/figures_output/sensitivity_analysis/',rd,'_paramvsSRC.pdf'), width = 6.25, height = 5)
  } else{
    pdf(paste0('../krat_genetics_scripts/figures_output/sensitivity_analysis/',rd,'_paramvsSRC.pdf'), width = 3, height = 5)
    xmin <- 0.65
    xmax <- max + 0.35
  }
  x <- 1
  plot(0,0, xlim = c(xmin, xmax), ylim = c(-1, 1), xaxt = 'n', xlab = 'Parameter', ylab = 'Standardized regression coefficient', col = 'transparent', cex.axis = text.size, cex.lab = text.size, main = rd)
    axis(1, at = c(1:max), labels = unique(origMelt1$Variables), cex.axis = text.size)
    abline(h = 0, lty = 2)
    for(v in unique(origMelt1$Variables)){
      sub <- origMelt1[origMelt1$Variables == v,]
      f <- sample(c(-100:100), nrow(sub))
      f <- f/shrink+x
      points(f, sub$SRC, col = alpha(sub$col[1], alph), pch = 19, cex = pt.size)
      arrows(x0 = x, x1 = x, y0 = (mean(sub$SRC, na.rm = TRUE) - sd(sub$SRC, na.rm = TRUE)),
             y1 = (mean(sub$SRC, na.rm = TRUE) + sd(sub$SRC, na.rm = TRUE)),
             lwd = 2, col = 'black', code=3, angle=90, length=0)
      points(x, mean(sub$SRC), pch = 23, col = 'black', bg = sub$col[1], cex = pt.size + 0.5, lwd = 2)
      
      x <- x+1
    }
  dev.off()
  
  # SRCs vs true value
  rohSRC1 <- inner_join(srcOriginal1, mean.iter1 ,by="id")
  
  # Plot each variable vs true FROH
  rohMelt1 <- melt(data = rohSRC1, id.vars=c("id","meanCall"),
                   measure.vars = colnames(rohSRC1)[colnames(rohSRC1) %notin% c('id','meanCall','sdCall')],
                   variable.name = "Variables",
                   value.name = "SRC")
  
  rohMelt1 <- merge(rohMelt1, cols, by = 'Variables')
  alph <- 0.8
  pt.size <- 0.8
  txt.size <- 1.25
  
  pdf(paste0('../krat_genetics_scripts/figures_output/sensitivity_analysis/',rd,'_callfROH_vs_SRC.pdf'), width = 7, height = 5)
  par(mai = c(1.02,0.82,0.82,1.42))
  plot(rohMelt1$meanCall, rohMelt1$SRC, ylim = c(-1, 1), col = 'transparent', pch = 19, xlab = substitute(paste('Mean called ',italic('F')[ROH])), ylab = 'Standardized regression coefficient', cex.axis = txt.size, cex.lab = txt.size,
       main = paste0(rd))
    abline(h = 0, lty = 2)
    points(rohMelt1$meanCall, rohMelt1$SRC, col = alpha(rohMelt1$col, alph), pch = 19, cex = pt.size)
    for(v in unique(rohMelt1$Variables)){
      sub <- rohMelt1[rohMelt1$Variables == v,]
      vals <- loess.smooth(sub$meanCall, sub$SRC, span = 0.75,
                           family = c("gaussian"), col = sub$col[1])
      lines(vals$x, vals$y, col = sub$col[1], lwd = 3)
    }
    par(xpd = TRUE)
    legend('right', lwd = 3, col = unique(rohMelt1$col), legend = unique(rohMelt1$Variables), bty = 'n', inset = -0.30, cex = txt.size)
  dev.off()
  
  ## make a zoomed-in figure to focus on those far-left points
  pdf(paste0('../krat_genetics_scripts/figures_output/sensitivity_analysis/',rd,'_callfROH_vs_SRC_zoomed_in.pdf'), width = 7, height = 5)
  par(mai = c(1.02,0.82,0.82,1.42))
  sub <- rohMelt1[order(rohMelt1$meanCall),]
  sub <- sub[c(1:(nrow(sub)-(4*length(unique(rohMelt1$Variables))))),] ## get rid of 4 samples * # variables
  plot(sub$meanCall, sub$SRC, ylim = c(-1, 1), col = 'transparent', pch = 19, xlab = substitute(paste('Mean called ',italic('F')[ROH])), ylab = 'Standardized regression coefficient', cex.axis = txt.size, cex.lab = txt.size,
       main = paste0(rd))
    abline(h = 0, lty = 2)
    points(sub$meanCall, sub$SRC, col = alpha(sub$col, alph), pch = 19, cex = pt.size)
    for(v in unique(sub$Variables)){
      temp <- sub[sub$Variables == v,]
      vals <- loess.smooth(temp$meanCall, temp$SRC, span = 0.75,
                           family = c("gaussian"), col = temp$col[1])
      lines(vals$x, vals$y, col = temp$col[1], lwd = 3)
    }
    par(xpd = TRUE)
    legend('right', lwd = 3, col = unique(sub$col), legend = unique(sub$Variables), bty = 'n', inset = -0.30, cex = txt.size)
  dev.off()
}

