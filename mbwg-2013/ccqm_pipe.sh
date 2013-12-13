#!/bin/sh

#
# pipeline for processing data from CCQM study
# by: Nate Olson
# organization: NIST
# written 03/11/2013
#
# requirements:
# TMAP for mapping
# samtools, picard, and gatk for processing sam files
# samtools for mpileup
# scripts for processing mpileup files

# code written assuming the requirements are in ~/bin/

# Functions
index_reference() {
	~/bin/tmap index -f $1
	~/bin/samtools faidx $1
}

process_fastq_bam() {
	#takes fastq as input and processes whole through to mpileup and gatk unified genotyper
	
	fq=$1
	org=$2
	method=$3
	lab=$4
	rep=$5
	se=$6
	st=$7

	# determining which reference to use
	if [[ $org == "Ecoli" ]]; then
		f=$Ecoli_ref
	else
		f=$Lmono_ref
	fi
	
	#mapping reference datasets
	tmap mapall -f $f -r $fq -n 2 -v -Y -u -o 0 stage1 map4 >$f.$fq.TMAP.sam

	#modifying sequence alignment
	samtools view -bSh -o $f.$fq.TMAP.bam $f.$fq.TMAP.sam
	java -Xmx4g -jar ~/bin/AddOrReplaceReadGroups.jar I= $f.$fq.TMAP.bam O= $f.$fq.TMAP.h.bam RGID=$lab RGLB=$method RGPU=$rep RGPL=$se RGSM=$org 
	java -Xmx4g -jar ~/bin/SortSam.jar I= $f.$fq.TMAP.h.bam O= $f.$fq.TMAP.h.sorted.bam SO=coordinate

	#do not remove replicates for amplicon datasets
	if [[ $method == 'pcr' ]]; then 
		mv $f.$fq.TMAP.h.sorted.bam $f.$fq.TMAP.h.sorted.dedup.bam
	else
		java -Xmx4g -jar ~/bin/MarkDuplicates.jar REMOVE_DUPLICATES=TRUE I=$f.$fq.TMAP.h.sorted.bam O=$f.$fq.TMAP.h.sorted.dedup.bam M=$f.$fq.TMAP.metric.txt
		rm $f.$fq.TMAP.h.sorted.bam $f.$fq.TMAP.metric.txt
	fi

	samtools index $f.$fq.TMAP.h.sorted.dedup.bam
	java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $f.$fq.TMAP.h.sorted.dedup.bam -T RealignerTargetCreator -o $f.$fq.TMAP.h.sorted.bam.intervals
	bam=$(echo $fq | sed 's/fastq/bam/')
	java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $f.$fq.TMAP.h.sorted.dedup.bam -T IndelRealigner -targetIntervals $f.$fq.TMAP.h.sorted.bam.intervals -o $bam
	
	samtools index $bam
	#GATK variant calling
	vcf=$(echo $fq | sed 's/fastq/vcf/')
	java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $bam -T UnifiedGenotyper -glm SNP -out_mode EMIT_ALL_SITES -o $vcf

	#generating mpileup file
	pup=$(echo $fq | sed 's/fastq/pileup/')
	samtools mpileup -f $f $bam > $pup

	#removing excess files
	rm  $f.$fq.TMAP.sam $f.$fq.TMAP.bam $f.$fq.TMAP.h.sorted.dedup.bam $f.$fq.TMAP.h.sorted.dedup.bam.bai\
	 $f.$fq.TMAP.h.sorted.bam.intervals $f.$fq.TMAP.h.bam

}

get_metadata_and_process_dataset (){
	# splits lines of code in metadata file and assigns appropriate variables
	# required variables
	# 1. data_set.file
	# 2. organism
	# 3. PCR or shotgun
	# 4. sample ID - lab that performed the sequencing
	# 5. sequencing unit - replicate number (1 for most)
	# 6. sequencing platform
	
	local meta_line=$1
	data_set= $(echo $meta_line | cut -d ',' -f1)
	org=      $(echo $meta_line | cut -d ',' -f2)
	method=   $(echo $meta_line | cut -d ',' -f3)
	lab=      $(echo $meta_line | cut -d ',' -f4)
	rep=      $(echo $meta_line | cut -d ',' -f5)
	se=       $(echo $meta_line | cut -d ',' -f6)
	
	echo "Working on:"	
	echo "$data_set $org $method $lab $rep $se $st $trim"
	
	# modifying datafile names for analysis 
	cp $data_set "$org-$se-$lab-$rep-$st-$trim.fastq"
	fastq="$org-$se-$lab-$rep-$st-$trim.fastq"
	
	# map fastq dataset to reference
	echo "processing fastq ..."
	process_fastq_bam $fastq $org $method $lab $rep $se

	# process mpileup
	echo "performing pileup counts ..."
	pup=$(echo $fastq | sed 's/fastq/pileup/')
	perl PileUpParse.pl $pup $pup.csv
}


# indexing and assigning reference variables
# assumes references are E. coli and L. monocytogenes
if [[ $2 == *coli* ]]; then 
	Ecoli_ref=$2 
	Lmono_ref=$3
else 
	Ecoli_ref=$3 
	Lmono_ref=$2
fi
index_reference $Ecoli_ref
index_reference $Lmono_ref

# Read metadata input file command line argument
while read line; 
	do
	get_metadata_and_process_dataset $line	
done < $1


#cleaning up intermediate files generated during analysis
rm *.tmap.* *.dict *.bai *.fai *.idx


