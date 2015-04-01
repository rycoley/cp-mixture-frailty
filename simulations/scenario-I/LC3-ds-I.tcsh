#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/lc-sims

#$ -cwd

#set environment variables

#output files
#$ -o LC3-ds1.out
#$ -e LC3-ds1.err

#request memory

#request nodes
#$ -pe local 35

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/dissert/lc-sims/LC3-ds-I.R

#clean up temp files?

