#!/bin/bash -l
#SBATCH -J comp_sort
#SBATCH -o /scratch/project_2001443/mt_WFhyb_Sapis/mtDNA_0422_ms/vcf/logs/compress_sort_%a.out
#SBATCH -e /scratch/project_2001443/mt_WFhyb_Sapis/mtDNA_0422_ms/vcf/logs/compress_sort_%a.err
#SBATCH --account=project_2001443
#SBATCH -t 02:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=4G

module load biokit
cd /scratch/project_2001443/mt_WFhyb_Sapis/mtDNA_0422_ms/vcf

# Compress
/appl/soft/bio/samtools/htslib-1.9/bgzip mt_F005_C1000_all_samples_raw.vcf

# Sort
bcftools sort -m 1G -O z -o mt_F005_C1000_all_samples.vcf.gz -T ./tmp_sort mt_F005_C1000_all_samples_raw.vcf.gz

# Index
/appl/soft/bio/samtools/htslib-1.9/tabix -p vcf mt_F005_C1000_all_samples.vcf.gz

# Count the number of sites #928
gunzip -c mt_F005_C1000_all_samples.vcf.gz | grep -v "^#" | cut -f 1 | uniq -c > mt_F005_C1000_site_count_raw.tab
