#!/bin/csh

#$ -M mbkuhn@math.ucr.edu
#$ -m abe
#$ -q  long 				# Specify queue
#$ -N  run_NOV13   		 	# Specify job name

mkdir Animate2
./program Animate2
