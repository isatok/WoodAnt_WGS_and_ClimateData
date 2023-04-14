# 6_SNP_calling.sh


cd /scratch/project_2001443/vcf




###
### Create additional input files
###

# 1. BAM list w/out sample 121 (excluded due to low quality)
ls /scratch/project_2001443/bam/nodupl_RG_clip/*.bam > bam.tmp
grep -v 121 bam.tmp > bam.list ; rm bam.tmp

# 2. Split ref in 50 kb regions to speed up the analysis
source activate freebayes
fasta_generate_regions.py Formica_hybrid_v1.fa.fai 50000 > Formica_hybrid_v1_50kb_regions.tmp




###
### Full SNP calling
###

module load bioconda
source activate freebayes # version: v1.3.1

cd /scratch/project_2001443/vcf
REF=/scratch/project_2001443/reference_genome
RES=/scratch/project_2001443/vcf

/appl/soft/bio/bioconda/miniconda3/envs/freebayes/bin/freebayes-puhti \
  -time 72 \
  -regions $REF/Formica_hybrid_v1_50kb_regions.txt \
  -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -L $RES/bam.list \
  -k --genotype-qualities --skip-coverage 36000 \
  -out $RES/raw/all_samples_raw.vcf




###
### Compress & sort
###

#!/bin/bash -l
#SBATCH -J comp_sort
#SBATCH -o /scratch/project_2001443/vcf/raw/compress_sort.out
#SBATCH -e /scratch/project_2001443/vcf/raw/compress_sort.err
#SBATCH --account=project_2001443
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=8G

module load biokit
cd /scratch/project_2001443/vcf/raw
# Compress
/appl/soft/bio/samtools/htslib-1.9/bgzip all_samples_raw.vcf

# Sort
bcftools sort -m 1G -O z -o all_samples.vcf.gz -T ./tmp_sort all_samples_raw.vcf.gz

# Index
/appl/soft/bio/samtools/htslib-1.9/tabix -p vcf all_samples.vcf.gz
bcftools index -n all_samples.vcf.gz

gunzip -c all_samples.vcf.gz | grep -v "^#" | cut -f 1 | uniq -c > site_count_per_chromosome.tab




### END
