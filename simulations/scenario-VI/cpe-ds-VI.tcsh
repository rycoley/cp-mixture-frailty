#!/bin/tcsh
#

#set working directory. any temporary files go here.
cd /home/bst/student/rcoley/dissert/cp-sims

#$ -cwd

#set environment variables

#output files
#$ -o cpe-ds6.out
#$ -e cpe-ds6.err

#request memory


#request nodes
#$ -pe local 20

# notifications
#$ -m abe
#$ -M ryc@jhu.edu

#execute command
R CMD BATCH --vanilla /home/bst/student/rcoley/dissert/cp-sims/cpe-ds-VI.R cpe-ds-VI.Rout

#clean up temp files?

