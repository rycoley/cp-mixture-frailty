#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/lc-sims

#$ -cwd

#set environment variables

#output files
#$ -o cpe-ds5.out
#$ -e cpe-ds5.err

#request memory

#request nodes
#$ -pe local 20

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH  /home/bst/student/rcoley/dissert/lc-sims/cpe-ds-V.R cpe-ds-V.Rout

#clean up temp files?

