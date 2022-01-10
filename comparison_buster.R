#!/bin/bash

#$ -o $HOME/logs/lgssm_crossover.R$TASK_ID.stdout
#$ -e $HOME/logs/lgssm_crossover.R$TASK_ID.stderr
#$ -l h_vmem=16G,h_rt=72:00:00
#$ -S /bin/bash

. /etc/profile
module add R
time R CMD BATCH lgssm_crossover.R $SGE_TASK_ID
