# 4_clip_overlaps_RG.sh


# Create the input file list
cd /scratch/project_2001443/bam/nodupl_RG
ls *bam | sort -nk1,1 > input.list ; cd ..

# Create the output folder
mkdir /scratch/project_2001443/bam/nodupl_RG_clip

# Prep RG metadata: sort as input list, extract only relevant column and swap underscores with dashes
sort -nk1,1 ID_to_ID_species.txt | cut -f2 | sed 's/_/-/' > sample_RG.txt


###
### sbatch script below ------------------------------------------------------------------------
###



#!/bin/bash -l
#SBATCH -J clipOverlapRG
#SBATCH -o /scratch/project_2001443/bam/logs/clipOverlap_RG_%j.out
#SBATCH -e /scratch/project_2001443/bam/logs/clipOverlap_RG_%j.err
#SBATCH --account=project_2001443
#SBATCH -t 04:00:00
#SBATCH -p small
#SBATCH --array=1-71
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4GB

# Load modules
module load biokit
module load bioconda/3
source activate my_seqdata

# Define directories
cd /scratch/project_2001443/bam/nodupl_RG
FINALDIR=/scratch/project_2001443/bam/nodupl_RG_clip

# Get file & sample ID
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p input.list)
sample=${file%_wRG*}

echo "### Processing $file"

###
### Correct read groups
###

myRG=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2001443/bam/sample_RG.txt)

java -Xmx4G -jar /appl/soft/bio/picard/picard-tools-2.21.4/picard.jar AddOrReplaceReadGroups \
    I=$file \
    O=$FINALDIR/${sample}"_wRG.bam" \
    RGID=$myRG \
    RGPL=illumina \
    RGLB=1 \
    RGSM=$myRG \
    RGPU=1 \
    TMP_DIR=/scratch/project_2001443/tmp


###
### Clip overlaps
###

# Reads are sorted by coordinates (see "samtools sort" command after mapping)
bam clipOverlap --in $FINALDIR/${sample}"_wRG.bam" --out $FINALDIR/${sample}"_wRG_clip.bam" --stats --params


###
### Index
###

samtools index $FINALDIR/$sample"_wRG_clip.bam"



### END
