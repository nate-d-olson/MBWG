#!/bin/sh

#
# pipeline for processing data from CCQM study for variant call analysis
# by: Nate olson
# organisation: NIST
# written 05/15/2013
# this is a revised version to include bwa as an additional mapper and samtools for snp calling
# Requirements:
# TMAP and bwa for mapping
# samtools, picard, and gatk for processing sam files
# samtools for mpileup
# samtools and gatk for snp calling
# scripts for processing mpileup files

# code written for requirements in ~/bin folder

# Functions
index_reference() {
	tmap index -f $1
	samtools faidx $1
	bwa index $1
	java -Xmx4g -jar ~/bin/CreateSequenceDictionary.jar R=$1 O=$1.dict	
}

mapping(){
	#takes fastq and reference fasta as input and maps using tmap and bwasw
	f=$1 # the reference fasta is the first argument passed
	fq=$2 # the read set in fastq format is the second argument passed
	
	prefix2=$(echo $fq | sed 's/.fastq//')
	#mapping datasets using tmap
	tmap mapall -f $f -r $fq -n 4 -v -Y -u -o 0 stage1 map4 > $prefix2-TMAP.bam #$prefix2-TMAP.sam
	#samtools view -bSh -o $prefix2-TMAP.bam $prefix2-TMAP.sam
	#mapping datasets using bwa bwasw
	bwa bwasw -t 4 $f $fq > $prefix2-bwa.sam
	samtools view -bSh -o $prefix2-bwa.bam $prefix2-bwa.sam
}

bam_processing(){
	# take bam files as input adds full read ground, then sorts and indexes

	f=$1 # the reference fasta is the first argument passed
	bam=$2 # the sam file is the second argument passed
	lab=$3
	method=$4
	rep=$5
	se=$6
	org=$7

	prefix3=$(echo $bam | sed 's/.bam/-basic/')	
	#modifying sequence alignment
	java -Xmx4g -jar ~/bin/AddOrReplaceReadGroups.jar I= $bam O= $prefix3.bam RGID=$lab RGLB=$method RGPU=$rep RGPL=$se RGSM=$org 
	echo "SortingSam ..."	
	java -Xmx4g -jar ~/bin/SortSam.jar I= $prefix3.bam O= $prefix3.sort.bam SO=coordinate
	samtools index $prefix3.sort.bam
	mv $prefix3.sort.bam $prefix3.bam
}

bam_refiner(){
	f=$1
	bam=$2
	method=$3

	prefix4=$(echo $bam | sed 's/-basic.bam/-refine/')
	#do not remove replicates for amplicon datasets
	if [[ $method != 'pcr' ]]; then 
		java -Xmx4g -jar ~/bin/MarkDuplicates.jar REMOVE_DUPLICATES=TRUE I=$bam O=$prefix4.bam M=$prefix4.metric.txt

	else
		cp $bam $prefix4.bam	
	fi
	
	samtools index $prefix4.bam
	java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $prefix4.bam -T RealignerTargetCreator -o $prefix4.intervals
	java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $prefix4.bam -T IndelRealigner -targetIntervals $prefix4.intervals -o $prefix4.realigned.bam
	mv $prefix4.realigned.bam $prefix4.bam
	rm $prefix4.intervals $prefix4.metric.txt
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
	
	prefix1=$(echo $fq | sed 's/.fastq//')
	# determining which reference to use
	if [[ $org == "Ecoli" ]]; then
		f=$Ecoli_ref
	else
		f=$Lmono_ref
	fi

	echo $fq

	mapping $f $fq

	for bam in $prefix1-TMAP.bam $prefix1-bwa.bam;
	do
		bam_processing $f $bam $lab $method $rep $se $org
		bam_ref=$(echo $bam | sed 's/.bam/-basic.bam/')
		bam_refiner $f $bam_ref $method
	done

	for bam in $prefix1-TMAP-basic.bam $prefix1-TMAP-refine.bam \
			$prefix1-bwa-basic.bam $prefix1-bwa-refine.bam;
	do
		samtools index $bam		
		prefix5=$(echo $bam | sed 's/.bam//')
		#GATK variant calling
		java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -R $f -I $bam -T UnifiedGenotyper -glm SNP -o $prefix5-gatk.vcf

		#samtools variant calling
		samtools mpileup -uf $f $bam > $prefix5-sam.pileup
		bcftools view -vcg $prefix5-sam.pileup > $prefix5-sam.vcf
	done
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
	# 7. strain source
	# 8. trim to specific genus (T or F)
	
	local meta_line=$1
	data_set=$(echo $meta_line | cut -f1 -d ',')
	org=$(echo $meta_line | cut -d ',' -f2)
	method=$(echo $meta_line | cut -d ',' -f3)
	lab=$(echo $meta_line | cut -d ',' -f4)
	rep=$(echo $meta_line | cut -d ',' -f5)
	se=$(echo $meta_line | cut -d ',' -f6)
	st=$(echo $meta_line | cut -d ',' -f7)
	trim=$(echo $meta_line | cut -d ',' -f8)
	
	echo "Working on:"	
	echo "$data_set $org $method $lab $rep $se $st $trim"
	
	# modifying datafile names for analysis 
	# will I want to move the analysis into a specific folder?
	cp $data_set "$org-$se-$lab-$rep-$st-$trim.fastq"
	fastq="$org-$se-$lab-$rep-$st-$trim.fastq"
	
	# map fastq dataset to reference
	echo "processing fastq ..."
	process_fastq_bam $fastq $org $method $lab $rep $se
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


