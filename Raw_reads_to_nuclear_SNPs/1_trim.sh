# 1_trim.sh

#!/bin/bash -l
#SBATCH -J trim
#SBATCH -o trim/logs/trim_%j.out
#SBATCH -e trim/logs/trim_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 12:00:00
#SBATCH -p small
#SBATCH --array=1-71
#SBATCH --ntasks 4
#SBATCH --mem=4G

cd /scratch/project_2001443/fastq

# Get file name
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list)
echo $file

# Load required modules
module load biokit
module load trimmomatic/0.39

# Trim reads
trimmomatic PE -threads 4 -phred33 \
raw/R1/$file"_R1_001.fastq.gz" raw/R2/$file"_R2_001.fastq.gz" \
trim/$file"_1.trim.pair.fq.gz" trim/unpair/$file"_1.trim.unpair.fq.gz" \
trim/$file"_2.trim.pair.fq.gz" trim/unpair/$file"_2.trim.unpair.fq.gz" \
ILLUMINACLIP:/scratch/project_2001443/fastq/adapters/general_adapters.fa:2:30:10:2:keepBothReads LEADING:10 TRAILING:10 S
LIDINGWINDOW:4:15 MINLEN:50

# Quality control
fastqc trim/$file"_1.trim.pair.fq.gz" -o trim/fastqc && \
fastqc trim/$file"_2.trim.pair.fq.gz" -o trim/fastqc

### END
