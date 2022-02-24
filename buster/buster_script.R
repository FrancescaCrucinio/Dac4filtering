#!/bin/bash

#$ -o $HOME/run_dac/logs/algo_comparison.R$TASK_ID.stdout
#$ -e $HOME/run_dac/logs/algo_comparison.R$TASK_ID.stderr
#$ -l h_vmem=4G,h_rt=1:00:00
#$ -S /bin/bash

. /etc/profile
module add R
time R CMD BATCH --no-save --no-restore run_dac/algo_comparison.R $SGE_TASK_ID
