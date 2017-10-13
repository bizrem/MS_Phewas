#REM 13/10/17
#Extract all SNPs, merge them(as extracted from different chr separately), qc appropriately, recode according to ref allele and generate grs using the meta-analysis snps
#Non-HLA

dataDir="/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/"
incluDir="/mnt/storage/home/bizrem/ms_phewas/data/grs_generation/input"
outputDir="/mnt/storage/home/bizrem/ms_phewas/data/grs_generation/output"

### Extraction of the SNP - note that using bgenix changes it to bestguess

#bgenix installed locally by me
for chr in {01..22}; do
bgenix \
-g ${dataDir}/dosage_bgen/data.chr$chr.bgen -i ${dataDir}/dosage_bgen/data.chr$chr.bgen.bgi -incl-rsids ${incluDir}/ms_meta_snps.txt >  ${outputDir}/ms_chr$chr.bgen
done

#Merge all extracted SNPs from chr together into one file
#get freq files as need freq for being able to force the reference allele later on - deciding if the correct strand ie the one with the risk allele has been genotyped in ukbiboank
#double check that the chr_files.txt is in appropriate directory

#Remove consent withdrawn, related individuals, Non-europeans, sex-mismatch, recommended exclusions
#generating txt files of both and then removing them using plink command

qcDir="/mnt/storage/private/mrcieu/research/UKBIOBANK_Array_Genotypes_500k_HRC_Imputation/data/derived/"

cat ${qcDir}/standard_exclusions/data.combined_recommended.txt ${qcDir}/ancestry/data.non_white_british.txt ${qcDir}/relateds_exclusions/data.minimal_relateds.txt ${qcDir}/relateds_exclusions/data.highly_relateds.txt  > ${incluDir}/ukbb_exclusions.txt

#Need to move the bgen files into bed bim fam format so the old plink can work with them
#Removing people at this stage too
#Plink is installed on my personnal space

for chr in {01..22}; do
plink2 \
--bgen ${outputDir}/ms_chr$chr.bgen \
--sample ${dataDir}/sample-stats/data.chr01.sample \
--remove  ${incluDir}/ukbb_exclusions.txt \
--make-bed \
--out ${outputDir}/ms_chr$chr
done


plink \
 --bfile ${outputDir}/ms_chr01 \
 --merge-list ${outputDir}/chr_files.txt \
 --freq \
 --make-bed \
 --out ${outputDir}/ms_snps_all_chr_exclusions

 #forcing specific reference allele to be the risk allele
#give one file (.raw) with all snps coded as dosages the right way round
#make sure your SNP list with risk alleles is in your home directory

plink \
 --bfile ${outputDir}/ms_snps_all_chr_exclusions \
 --reference-allele ${incluDir}/ms_risk_alleles.txt \
 --recode A \
 --freq \
 --make-bed \
 --out ${outputDir}/ms_snps_all_chr_exclusions_RA

 #Make score
 #need to convert txt file with weights(odds ratios) into a .raw file

 plink \
  --bfile ${outputDir}/ms_all_chr_exclusions_RA \
  --score ${incluDir}/ms_effect.raw \
  --out ${outputDir}/ms_grs
