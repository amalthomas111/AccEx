#!/bin/bash
#echo "$#"
if [ "$#" -ne 1 ]; then
    printf "\n"
    echo "bash script.sh indexbamfile"
    echo "---------------------------------------------------------------------"
    exit
fi

bam=$1
samtools idxstats ${bam} |awk -v name=${bam} \
        '{OFS="\t"}BEGIN {total=0;reads=0} {if(length($1)<=5){total += $3;if($1=="chrM")reads=$3;}} END {print name,total/1000000,reads/1000000,reads*100/total}'
