### Converting VCF to genotype format suitable for reading in to 'related' package
setwd('/Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/')

library(vcfR)

## read in LD-pruned VCF file
vcf <- read.vcfR('krat_final_final_allfiltcontigs_all_samps_LDpruned.recode.vcf.gz')
gts <- extract.gt(vcf, element = 'GT', return.alleles = TRUE)
gts <- t(gts)


## format for 'related': column 1 = individual IDs, 1 column per allele (i.e., 2 columns per SNP); 
## missing data = 0 (so alleles have to be letters I guess?); no header row
for(r in 1:nrow(gts)){
  # gts[r,] <- gsub('.', '0/0', gts[r,], fixed = TRUE) ## convert missing data to zeroes
  # gts[r,] <- gsub('|', '/', gts[r,], fixed = TRUE)
  gts[r,] <- gsub('/', ' ', gts[r,], fixed = TRUE)
  rownames(gts)[r] <- unlist(strsplit(rownames(gts)[r], split = '_'))[1]
}
table(gts[1,])

write.table(gts, 'related_formatted_genotypes.txt', col.names = FALSE, row.names = TRUE, quote = FALSE, sep = ' ')

#

















# OUT <- NULL
# write.table(OUT, 'related_formatted_genotypes.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
# for(c in 1:ncol(gts)){
#   sub <- gts[,c]                                            ## get all genotypes for sample
#   id <- unlist(strsplit(colnames(gts)[c], split = '_'))[1]  ## get sample ID
#   print(id)
#   sub <- gsub('|', '/', sub, fixed = TRUE)                  ## make genotype formatting consistent
#   save <- NULL
#   save <- c(save, id)
#   for(s in 1:length(sub)){                                  ## for each locus,
#     # print(s/length(sub))
#     if(sub[s] == '.'){
#       save <- c(save, 0, 0)
#     } else{
#       ## save each allele separately
#       save <- c(save, unlist(strsplit(sub[s], split = '/'))[1], unlist(strsplit(sub[s], split = '/'))[2])
#     }
#   }
#   OUT <- rbind(OUT, save)
#   write.table(OUT, 'related_formatted_genotypes.txt', append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
#   OUT <- NULL
# }
