#!/bin/sh

# code for testing required proportions for biologically variant positions 
# by: Nate Olson
# organization: NIST
# written 03/11/2013
#
# requirements:
# GenomeAnalysisTK.jar
# code written assuming the requirements are in ~/bin/
# execute in the same folder as the bam files being downsampled
# along with the reference fasta files
 
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
	
