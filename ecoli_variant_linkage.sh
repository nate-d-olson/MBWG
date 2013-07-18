#!/usr/bin/sh
for bam in *F.bam;
do
	sam=$(echo $bam | sed 's/bam/sam/')
	csv=$(echo $bam | sed 's/bam/csv/')
	samtools view -o $sam $bam
	python SAM_sequencingError_parser.py $sam Ecoli_16S_consensus.fasta $csv
done

