#!/bin/sh

#code for testing proportions with down sampling 
for bam in *bam;
	do
	if [[ $bam == *coli* ]]; then 
	f="Ecoli_16S_consensus.fasta"
	else 
	f="Lmono_16S_consensus.fasta"
	fi
	for i in $(seq 1 9); do
		for j in 0.0001 0.001 0.1;do
			frac=$(($i * $j)) 
			pup=$(echo $bam | sed "s/.bam/-$frac.pileup/")             
			java -jar ~/bin/GenomeAnalysisTK.jar -T Pileup -I $bam -R $f -dfrac $frac -o $pup;
		done            	
	done
done
	
