# cd /Users/Avril/Documents/krat_genetics/data/
# cp krat_plink_rohverlap_GTs_ALL_SAMPLES.txt krat_plink_rohverlap_ALL_SAMPLES_numericGTs.txt

cd /Users/Avril/Documents/krat_genetics/data/LD_pruned_data_and_results/
cp krat_LDpruned_extracted_genotypes.txt krat_LDpruned_extracted_genotypes_numericGTs.txt


sed -i '' 's/\.\/\./0/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
sed -i '' 's/0\/1/0/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
sed -i '' 's/0\|1/0/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
echo "hets gone"
sed -i '' 's/0\/0/1/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
sed -i '' 's/0\|0/1/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
echo "ref homs gone"
sed -i '' 's/1\/1/2/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
sed -i '' 's/1\|1/2/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
echo "ref hets gone"
sed -i '' 's/\t\.\t/\t0\t/g' krat_LDpruned_extracted_genotypes_numericGTs.txt
