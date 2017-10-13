#!/bin/csh

#$ -M mbkuhn@math.ucr.edy 	# Email address for job notifications
#$ -m  abe					# Send mail when job begins, ends and aborts
#$ -q  long 				# Specify queue
#$ -N  run_October13	 	# Specify job name

./program
