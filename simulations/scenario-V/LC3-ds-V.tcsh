#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/lc-sims

#$ -cwd

#set environment variables

#output files
#$ -o LC3-ds5.out
#$ -e LC3-ds5.err

#request memory

#request nodes
#$ -pe local 40

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/dissert/lc-sims/LC3-ds-V.R

#clean up temp files?

