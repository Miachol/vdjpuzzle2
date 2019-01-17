#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)


#PBS -P va1

#PBS -q normal

#PBS -l walltime=24:00:00
#PBS -l ncpus=1
#PBS -l mem=20G

#PBS -l wd

set -x # echo on, command fails causes script to exit, pipes fail

#parameters
#parameter one needs to be the ABSOLUTE path where cell sequences are located WITHOUT /
if [[ -n "$P3" ]]; then
	param1=$P1
	param2=$P2
	param3=$P3
	param4=$P4
	param5=$P5
	param6=$P6
else
	param1=$1
	param2=$2
	param3=$3
	param4=$4
	param5=$5
	param6=$6
fi

# if PATH_PARAM has been passed to the script, then set PATH
# this is because sometimes PATH gets overwritten on slave nodes, even when using -V
if [ ! -z ${PATH_PARAM+x} ]; then
	export PATH=$PATH_PARAM
fi


CELL_PATH=$param1
CHAIN_ARRAY=($param5)
CHAIN_PREFIX_ARRAY=($param6)

#unzip all files in a dir

echo "P1: $param1
P2: $param2
P3: $param3
P4: $param4
P5: $param5
P6: $param6"

index=0
for chain in "${CHAIN_ARRAY[@]}"
do
	migmap -S $param2 -R ${chain//C} --by-read --data-dir=$CONDA_PREFIX/share/igblast $CELL_PATH/${chain}.fa $CELL_PATH/${chain}.out
	tail -n+2 $CELL_PATH/${chain}.out > $CELL_PATH/${chain}.tmp
	cut -f1 $CELL_PATH/${chain}.tmp -d " " > $CELL_PATH/reads_${CHAIN_PREFIX_ARRAY[$index]}.txt
	cut -c 2- $CELL_PATH/reads_${CHAIN_PREFIX_ARRAY[$index]}.txt | xargs -n 1 $SAMTOOLS faidx $CELL_PATH/${chain}.fa > $param3/matching_reads_${CHAIN_PREFIX_ARRAY[$index]}_$param4.fa

	index=$((index+1))
done