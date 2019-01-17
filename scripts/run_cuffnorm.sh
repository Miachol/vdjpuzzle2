#!/bin/bash

#PBS -P va1

#PBS -q hugemem

#PBS -l walltime=48:00:00
#PBS -l ncpus=7
#PBS -l mem=500G

#PBS -l wd

set -x -e -o pipefail # echo on, command fails causes script to exit, pipes fail

# if PATH_PARAM has been passed to the script, then set PATH
# this is because sometimes PATH gets overwritten on slave nodes, even when using -V
if [ ! -z ${PATH_PARAM+x} ]; then
	export PATH=$PATH_PARAM
fi

cuffnorm --use-sample-sheet -o CuffNorm $ANNOTATION $CUFFOUTPUT/sample_sheet.txt
