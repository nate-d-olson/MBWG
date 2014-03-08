<<<<<<< HEAD
MBWG
====

Code used in analysing data generated during CCQM MBWG Microbial Identity Studies
=======
#_CCQM MBWG Microbial Identity 2013_

_Description:  The code in this repository is still in development.  International inter-laboratory study sequencing 16S rRNA from two genomic reference materials._

## List of scripts used to analyze the results.

_Scripts used in analyzing data from the 2013 CCQM MBWG Microbial Identity Working Group.
Publication in preparation, the reference will be added when published. Note that the scripts in this repository were developed specifically to analyze the results of the MBWG study and therefore will require appropriate modifications for application to other data sets._

_List of scripts with brief description, list of requirements, and usage:_
-_ccqm_pipe.sh_
	-_Description:_
		_Main script used in the study. Maps data sets to a reference sequence and generates variant calls._
	-_Requirements:_ 
		1. _TMAP (mapping reads to reference, https://git hub.com/iontorrent/TMAP)_
		2. _samtools (processing sam files, http://samtools.sourceforge.net/)_
		3. _Picard (processing sam files, http://picard.sourceforge.net/)_
		4. _GenomeAnalysisToolKit (processing sam files, http://www.broadinstitute.org/gatk/)_
	-_Usage:_
		_bash ccqm_pipe.sh metadatafile.txt org1_16S_reference.fasta org2_16S_reference.fasta_
		_Notes:_
		-_metadatafile.txt: Includes a list of the sequence data files in fastq format along with the associated meta data. See ccqm_metadata.txt for example_
		-_org(1/2)_16S_reference.fasta: The consensus 16S rRNA reference sequences the reads are mapped to, use Ns for biologically variant positions._ 

		_Description of ccqm_pipe metadata file
		comma separated file in the following format
		data_set,org,method,lab,rep,se,st,trim
		data_set = name of the fastq file with the sequencing reads
		org = template organism - used to determine which reference sequence to map to
		method = whether the sequencing data were amplicons (e.g. 454) or shotgun (e.g. Ion Torrent)
		rep = lab platform replicate number
		se = sequencing platform
		st = strain source
		trim = if reads not identified to the correct genus level were removed based on RDP *the use of this data was discontinued early in the study_

-_CCQM_coverage.R_
	-_Description: code for generating coverage plot (Figure 2)_
	-_Requirements: R packages ggplot2 and reshape2_
	-_Usage: can execute from R IDE or commandline, will need to change the variable indicating the location of the parsed pileup files generated with ccqm_pipe.sh_

-_ccqm_biovar_proportions.R_
	-_Description: code for Bayesian proportions analysis and figure_
	-_Requirements: R packages ggplot2, reshape2, and plyr_
	-_Usage: can run from commandline for IDE, make sure to change pileup_csv variable to the location of parsed pileup files._

-_ccqm_downsampling_proportion.sh_
	-_Description: code for testing required proportions for biologically variant positions by generating pileup files from downsampled bam files._
	-_Requirements: GenomeAnalysisTK.jar_
	-_Usage: execute in the same folder as the bam files being downsampled along with the reference fasta files._

-_subsample_proportions_analysis.R_
	-_Description: subsampled proportion analysis and plots for subsampled biologically variant positions._
	-_Requirements: R packages ggplot2, reshape2, and plyr_
	-_Usage: can run from commandline for IDE, make sure to change pileup_csv variable to the location of parsed pileup files from ccqm_downsampling_proportions._

-_SAM_parser.py_
	-_Description: Converts pileup file to a tab separated file with individual bases as rows in the document._
	-_Requirements: python package Biopython_
	-_Usage: from commandline "pileup_parse.py input.pileup output.tsv"_

-_lmono_variant_linkage_figure.r_
	-_Description: Generates full length sequence figure and data table for maximum likelihood analysis (lmono_variant_linkage_analysis.r). _
	-_Requirements: R packages reshape2 and ggplot2_
	-_Usage: Can run from commandline or within an R IDE._

-_lmono_variant_linkage_analysis.r_
	-_Description: Determination of the individual full length 16S rRNA gene copies_
	-_Requirements: R packages reshape2 and ggplot2 and files parsed by SAM_parser_
	-_Usage: Can run from commandline or within an R IDE._

-_ecoli_variant_linkage_analysis.R_
	-_Description: Code for generating data and figure for full length sequence analysis_
	-_Requirements: R packages reshape2 and ggplot2 and files parsed by SAM_parser_
	-_Usage: Can run from commandline or within an R IDE._

-_contam_analysis.py_
	-_Description:_
		_Script used to parse results of the taxonomic classification of the reads using the RDP Bayesian classifier (http://rdp.cme.msu.edu/)._
	-_Requirements:_
		_Pandas python package_
	-_Usage:_
		_python contam_analysis.py_
		_Notes:_
		-_Requires a metadata file: ccqm_metadata_rdp_analysis.txt (see repository for example), requires as input results from the RDP classifier as implemented in QIIME (http://qiime.org/)_
>>>>>>> 4c8c845d416efeafb29505c3b9e38a46db443a9a
