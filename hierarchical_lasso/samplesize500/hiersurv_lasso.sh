#! /bin/sh
#$ -cwd
#$ -t 1107-1606
#$ -M hchen130@jhu.edu
module load R
Rscript hiersurv_lasso.R

