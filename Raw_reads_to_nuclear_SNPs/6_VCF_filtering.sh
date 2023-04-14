# 7_VCF_filtering.sh

## Since different parts of the pipeline require different resources, the master script here is split in four scripts

######################################################## Beginning script 1

#!/bin/bash -l
#SBATCH -J proc
#SBATCH --account=project_2001443
#SBATCH -t 48:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=12G

module load bioconda/3
source activate my_seqdata

module load biokit
cd /scratch/project_2001443/vcf/filt

##
## 0. Pre-filtering: normalisation, indels, non-SNPs, read imbalance, decomposition ---------------
##

vt normalize -n -r ../../../ref/Formica_hybrid_v1_wFhyb_Sapis.fa ../all_samples.vcf.gz | bgzip -c > all_samples.normalized.vcf.gz

bcftools filter --threads 8 -Oz -s+ --SnpGap 2 all_samples.normalized.vcf.gz > all_samples.normalized.SnpGap_2.vcf.gz && \

bcftools filter --threads 8 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ all_samples.normalized.SnpGap_2.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 8 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 & SAF>=1 && SAR>=1' all_samples.normalized.SnpGap_2.NonSNP.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 8 -O z -f PASS all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 8 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir /scratch/project_2001443 --max-mem 4G -O z > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \

bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz




##
## 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
##

bcftools filter --threads 8 --include 'QUAL >= 30 && TYPE="snp"' -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 8 --min-alleles 2 --max-alleles 2 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz"
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz




##
## 2. Prep for max DP filtering -------------------------------------------------------------------
##

### a. Correct field IDs (MNV issue)

# Extract header from VCF
bcftools view -h all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz

# Check fraction missing data
vcftools --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz --missing-indv --out all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.MISSING_DATA

echo "Average missing data over all individuals:"
cut -f5 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.MISSING_DATA.imiss | grep -v 'F_' | awk '{sum+=$1;} END{print sum/NR;}'


### b. Compute mean DP per individual
vcftools --depth --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz --out all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader


### c. Compute individual DP thresholds
module load r-env
scriptR minMaxDP.R

# minMaxDP.R script copied below, for the record

      # dat = read.table("all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth", h = T)
      # dat$MIN_DEPTH = round((dat$MEAN_DEPTH-.5)/2) # min depth is half of mean depth
      # dat$MAX_DEPTH = round((dat$MEAN_DEPTH)*2) # max depth is twice mean depth
      # write.table(dat, "all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth.MinMax", col.names = T, row.names = F, quote = F, sep = "\t")


### d. Create results folder
mkdir logs
mkdir individualDP


######################################################## End script 1




##
## 3. Filtering on sequencing depth ---------------------------------------------------------------
##

######################################################## Beginning script 2

#!/bin/bash -l
#SBATCH -J splitDP
#SBATCH --account=project_2001443
#SBATCH -o logs/splitDP_%a.out
#SBATCH -e logs/splitDP_%a.err
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-90
#SBATCH --ntasks 1
#SBATCH --mem=2G

module load biokit
cd /scratch/project_2001443/vcf/filt

# Get file name and sample ID
ind=$(sed -n "$SLURM_ARRAY_TASK_ID"p ind.list)
mkdir individualDP/${ind}
cd individualDP/${ind}

# Extract single individual from VCF and apply DP filters

dpmin=2 # hard-coded, we will apply the same min filter accross all individuals afterwards
dpmax=$(grep $ind ../../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.idepth.MinMax | cut -f5)

echo "Processing individual $ind with $dpmin < DP <= $dpmax"
bcftools view -s ${ind} -Ou ../../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz | bcftools filter -e "FORMAT/DP < ${dpmin} | FORMAT/DP >= ${dpmax}" --set-GTs . -Oz > ${ind}.vcf.gz &&
bcftools index -t ${ind}.vcf.gz

######################################################## End script 2




##
## 4. HWE testing ---------------------------------------------------------------------------------
##

######################################################## Beginning script 3

#!/bin/bash -l
#SBATCH -J filt3
#SBATCH --account=project_2001443
#SBATCH -t 48:00:00
#SBATCH -p small
#SBATCH --ntasks 4
#SBATCH --mem=8G


cd /scratch/project_2001443/vcf/filt/individualDP
module load biokit
module load r-env


### Combine all samples
ls */*gz > female.list
bcftools merge -l female.list -Oz --threads 4 > ../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.vcf.gz


### Test HWE in females
cd ..
vcftools --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.vcf.gz --hardy --out FEMALES_hwe

### Extract sites with heterozygote excess, in R
Rscript hwe_filtering.R

# hwe_filtering.R copied below, for the record

      # # get the filename to read the hwe test
      # filename = "FEMALES_hwe" # file name of reference table
      # # print the command line options
      # print(paste("file name with hwe output:", filename))
      #
      # # Read the HW results
      # hwvalues <- read.table(paste(filename,".hwe",sep=""), header=T, stringsAsFactors = F)
      # # Get distribution of the p-values
      # print(paste("he excess", sum(hwvalues$P_HET_EXCESS<0.01)))
      # print(paste("he deficit", sum(hwvalues$P_HET_DEFICIT<0.01)))
      # print(paste("he overall", sum(hwvalues$P_HWE<0.01)))
      #
      # # remove only the sites with excess of heterozygotes
      # indextoremove = which(hwvalues$P_HET_EXCESS<0.01)
      # # Create BED file with the sites that fail the HW test for excess of heterozygotes
      # # BED file has three entries:
      # # chromosome, position-1, position (it is position-1 because it is assumed that the first base is 0)
      # position = hwvalues[indextoremove, 2]
      # bedmatrix = matrix(c(hwvalues[indextoremove, 1], format(position-1, scientific = F), format(position, scientific = F)), ncol=3, byrow=F)
      # write("#Sites with heterozygosity excess", file=paste("nohwe_excess_",filename,".bed",sep=""))
      # write(t(bedmatrix), file=paste("nohwe_excess_",filename,".bed",sep=""), ncolumns = 3, append=T)


### Remove these sites in both sex-specific VCFs

vcftools --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.vcf.gz --exclude-bed nohwe_excess_FEMALES_hwe.bed --recode --recode-INFO-all --stdout | bgzip >  all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz




##
## 5. Min DP filtering  ---------------------------------------------------------------------------
##

# Individual genotypes not reaching the minimal cutoff dpmin are tagged as missing (--set-GTs .)

# DP 8
dpmin=8
bcftools filter --threads 4 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.indDP.hwe.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP${dpmin}.hwe.vcf.gz


##
## 6. Filter based on missing data ----------------------------------------------------------------
##

### FILTERING BASED ON ALLELIC NUMBER (AN)

# 90 diploid samples = 180 alleles expected (the name "allele" is misleading here, as all sites are biallelic)
# 10% missing data = 162 alleles min.


# DP8
inputfile=all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.vcf.gz
bcftools view --threads 4 -Oz -i 'AN >= 162' ${inputfile} > ${inputfile%.vcf.gz}.AN10percMiss.vcf.gz


# Index VCFs, count number of sites left overall & percentage of missing data
for i in *percMiss.vcf.gz ; do
  echo $i ;
  bcftools index -ft $i
  bcftools index -n $i
  vcftools --gzvcf ${i} --missing-indv --out ${i%.vcf.gz}
done

######################################################## End script 3



######################################################## Beginning script 4

# Thin with 20kb distances, exclude unmapped regions (Scaffold00) and the social chromosome (Scaffold 03),
# exclude a collaborative sample prepared alongside our data but not belonging to this study (110-FaquH),
# and run minor allele frequency filter

#start interactive mode
sinteractive --account project_2001443 --mem 4000

module load biokit
FULLVCF=/scratch/project_2001443/vcf/filt/all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.minDP8.hwe.AN10percMiss.vcf.gz
OUTVCF=$SCRATCH/vcf/filtered_modified_vcf/all_samples.DP8.hwe.AN10.mac2.noScaff0003.thin20kb.vcf.gz

vcftools --gzvcf $FULLVCF --not-chr Scaffold00 --not-chr Scaffold03 --remove-indv 110-FaquH --thin 20000 --mac 2 --recode --stdout | bgzip > $OUTVCF
bcftools index -t $OUTVCF

######################################################## End script 4
### THIS IS THE END.




