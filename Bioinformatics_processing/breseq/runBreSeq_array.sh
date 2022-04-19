#!/bin/bash

#use the current working directory and current modules
#$ -cwd -V

#$ -b y -l mem_free=24G -P jravel-lab -q threaded.q -pe thread 8 -N MU01B_breseq -j y -o /local/scratch/mfrance/logs/ -e /local/scratch/mfrance/logs/

#setting the number of jobs to be executed
#$ -t 1-29

cd /local/scratch/mfrance/MU01B/08_breSeq

infile=`sed -n "$SGE_TASK_ID p" Mapping_instructions.csv`

refdir=/local/scratch/mfrance/MU01B/08_breSeq/reference
readdir=/local/groupshare/ravel/mfrance/small_projects/MU01B/08_breSeq

ref="$(cut -d ',' -f 1 <<< "$infile" )"
mg1="$(cut -d ',' -f 2 <<< "$infile" )"
mg2="$(cut -d ',' -f 3 <<< "$infile" )"

mg2=$(echo $mg2|tr -d '\r')

mkdir ${ref}_match
cd ${ref}_match
/home/mfrance/software/breseq-0.33.2/bin/breseq -j 8 -n ${ref}_match -r $refdir/${ref}.gbf $readdir/${mg1}_R1.fq.gz $readdir/${mg1}_R2.fq.gz $readdir/${mg1}_unpaired.fq.gz

rm -r 0*
rm -r data

cd ..

mkdir ${ref}_unmatch
cd  ${ref}_unmatch
/home/mfrance/software/breseq-0.33.2/bin/breseq -j 8 -n ${ref} -r $refdir/${ref}.gbf $readdir/${mg2}_R1.fq.gz $readdir/${mg2}_R2.fq.gz $readdir/${mg2}_unpaired.fq.gz

rm -r 0*
rm -r data


