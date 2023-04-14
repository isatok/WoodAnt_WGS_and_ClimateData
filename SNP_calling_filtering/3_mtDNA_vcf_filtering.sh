#!/bin/bash -l
#SBATCH -J proc
#SBATCH --account=project_2001443
#SBATCH -t 12:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=6G


#####
##### 0. Pre-filtering -----------
#####

module load bioconda/3
module load biokit

source activate /projappl/project_2001443/bioconda3_env/my_seqdata/ 

cd /scratch/project_2001443/mt_WFhyb_Sapis/mtDNA_0422_ms/vcf


vt normalize -n -r /scratch/project_2001443/mt_WFhyb_Sapis/mtDNA_0422_ms/ref/mtDNA_wFhyb_Sapis.fa mt_F005_C1000_all_samples.vcf.gz | bgzip -c > all_samples.normalized.vcf.gz

bcftools filter --threads 1 -Oz -s+ --SnpGap 2 all_samples.normalized.vcf.gz > all_samples.normalized.SnpGap_2.vcf.gz && \

bcftools filter --threads 1 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ all_samples.normalized.SnpGap_2.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 1 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 & SAF>=1 && SAR>=1' all_samples.normalized.SnpGap_2.NonSNP.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 1 -O z -f PASS all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 1 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir /scratch/project_2001443 --max-mem 4G -O z > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \

bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz

#####
##### 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
#####

bcftools filter --threads 1 --include 'QUAL >= 30 && TYPE="snp"' -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 1 --min-alleles 2 --max-alleles 2 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'"
gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz


#####
##### 2. Filtering on sequencing depth--------------------------------------------------
#####


##
## Prep for max DP filtering
##

### a. Correct field IDs

# Extract header from VCF
bcftools view -h all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz

bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.fixedHeader.vcf.gz #248 SNPs

### This is the end.