#!/bin/bash

#$ -o $HOME/d4f_nsmc_stpf/logs/lgssm_crossover.R$TASK_ID.stdout
#$ -e $HOME/d4f_nsmc_stpf/logs/lgssm_crossover.R$TASK_ID.stderr
#$ -l h_vmem=16G,h_rt=24:00:00
#$ -S /bin/bash

. /etc/profile
module add R
time R CMD BATCH d4f_nsmc_stpf/lgssm_crossover.R $SGE_TASK_ID
