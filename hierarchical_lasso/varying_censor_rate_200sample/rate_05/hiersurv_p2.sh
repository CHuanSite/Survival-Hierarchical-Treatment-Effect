#! /bin/sh
#$ -cwd
#$ -t 1107-1306
#$ -M hchen130@jhu.edu
module load R
Rscript hiersurv_p2.R
