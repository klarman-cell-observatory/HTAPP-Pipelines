#! /bin/bash

#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l os=RedHat7
#$ -l h_vmem=3G
#$ -pe smp 12
#$ -binding linear:12
#$ -l h_rt=3:00:00
#$ -e mkfastq-6.0.2.err
#$ -o mkfastq-6.0.2.out
#$ -R y
source /broad/software/scripts/useuse
reuse -q .bcl2fastq2-v2.20.0
reuse Python-2.7
reuse UGER
export PATH=/broad/xavierlab_datadeposit/cporter/bin/cellranger-6.0.2:$PATH
cellranger mkfastq --run=/ahg/regev_nextseq/Data02/171117_NB501337_0302_AHMCCTBGX3 --samplesheet=/ahg/regevdata/projects/MantonSkin/fastq/10x_liver004_005/10x_liver004_005.csv --output-dir=/ahg/regevdata/projects/MantonSkin/fastq/10x_liver004_005
