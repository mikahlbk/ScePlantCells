#!/bin/csh

#$ -q  debug 				# Specify queue
#$ -N  run_Oct25   		 	# Specify job name

#setenv PATH /afs/crc.nd.edu/user/a/awhitake/PlantCells/ScePlantCells:$PATH
mkdir Animate_debug
./program Animate_debug
