#!/bin/bash

#$ -o $HOME/logs/lgssm_dac_comparison.$TASK_ID.stdout
#$ -e $HOME/logs/lgssm_dac_comparison.$TASK_ID.stderr
#$ -l h_vmem=16G,h_rt=24:00:00
#$ -S /bin/bash

. /etc/profile
module add R
time R CMD BATCH lgssm_dac_comparison_buster.R $SGE_TASK_ID
