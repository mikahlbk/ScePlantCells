#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -q  long 				# Specify queue
#$ -N  run_JAN18_Animate16   		 	# Specify job name

mkdir Animate16
./program Animate16
