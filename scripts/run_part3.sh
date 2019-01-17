#!/bin/bash
#This script execut the RNA-Seq alingnamet pipeline sequentially, multithreading is implemented within the scripts (Trimmomatic, tophat and cufflinks)
#This script should be executed in parallel on different cells indicating the right directory ($1)

#PBS -P va1

#PBS -q normal

#PBS -l walltime=48:00:00
#PBS -l ncpus=8
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
	param7=$P7
	param8=$P8
	param9=$P9
	param10=$P10
	param11=$p11
else
	param1=$1
	param2=$2
	param3=$3
	param4=$4
	param5=$5
	param6=$6
	param7=$7
	param8=$8
	param9=$9
	param10=${10}
	param11=${11}
fi

# if PATH_PARAM has been passed to the script, then set PATH
# this is because sometimes PATH gets overwritten on slave nodes, even when using -V
if [ ! -z ${PATH_PARAM+x} ]; then
	export PATH=$PATH_PARAM
fi

CELL_PATH=$param1

FNAME1=`find -L ${CELL_PATH} -name "*fastq.gz" | egrep ".+_(R1_001|R1|1)\.fastq\.gz" | grep -v "PAIRED" | xargs basename` # ${CELL_PATH}/${param2}1.fastq.gz"
FNAME2=`find -L ${CELL_PATH} -name "*fastq.gz" | egrep ".+_(R2_001|R2|2)\.fastq\.gz" | grep -v "PAIRED" | xargs basename` #"${CELL_PATH}/${param2}2.fastq.gz"
Q1=${CELL_PATH}/$FNAME1
Q2=${CELL_PATH}/$FNAME2
Q3=$param4/VDJ_p3_$param2
Q4=$param4
CHAIN_ARRAY=($param6)
CHAIN_PREFIX_ARRAY=($param7)

echo "P1: $param1
P2: $param2
P3: $param3
P4: $param4
P5: $param5
P6: $param6
P7: $param7
P8: $param8
P9: $param9
P10: $param10
P11: $param11"

rm $Q3/overlapping_reads*
for prefix in "${CHAIN_PREFIX_ARRAY[@]}"
do
	rm $Q3/out1${prefix}.fastq
	rm $Q3/out2${prefix}.fastq
done
mkdir -p $Q3
mkdir -p $Q3/out

if [ "$param5" -ge 1 ]; then
	Q1="${CELL_PATH}/PAIRED_${FNAME1}"
	Q2="${CELL_PATH}/PAIRED_${FNAME2}"
fi

if [ "$param11" -ge 1 ]; then
	echo "Trim galore"
	filename1="${FNAME1%.*}"
	filename1="${filename1%.*}"
	filename2="${FNAME2%.*}"
	filename2="${filename2%.*}"
	Q1="${CELL_PATH}/{filename1}_val_1.fq.gz"
	Q2="${CELL_PATH}/{filename2}_val_2.fq.gz"
fi

for chain in "${CHAIN_ARRAY[@]}"
do
	bowtie2 --no-unal -p $param8 -k 1 --np 0 --rdg 1,1 --rfg 1,1 -x $Q4/assembled${param9}_genome/${chain} -1 $Q1 -2 $Q2 --al-conc $Q3/reads_${chain}_%.fastq -S $Q4/${chain}.sam
done

for chain in "${CHAIN_ARRAY[@]}"
do
	Trinity --left $Q3/reads_${chain}_1.fastq --right $Q3/reads_${chain}_2.fastq --seqType fq --max_memory 10G --output $Q3/trinity_out_dir
	mv $Q3/trinity_out_dir/Trinity.fasta $Q3/${chain}.fa
	rm -rf $Q3/trinity_out_dir
done

mkdir -p $param4/summary
echo "Running MIGMAP"
for chain in "${CHAIN_ARRAY[@]}"
do
	if [[ $param10 -ge 1 ]] ; then
		migmap -S $param3 -R ${chain//C} --data-dir=$CONDA_PREFIX/share/igblast --details fr1nt,cdr1nt,fr2nt,cdr2nt,fr3nt,fr4nt $Q3/${chain}.fa $param4/summary/${chain//C}_$param2
	else
		migmap -S $param3 -R ${chain//C} --data-dir=$CONDA_PREFIX/share/igblast $Q3/${chain}.fa $param4/summary/${chain//C}_$param2
	fi
done
