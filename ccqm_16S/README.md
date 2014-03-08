## To Do:
Write script generate links to binaries
change ccqm-16S-bioinformatics/data to results
need to rethink document structure with stat/ pub/ and bioinformatics directories
pick up documentation at - fastq_stats.sh

# CCQM MBWG Microbial Identity 2013  
Description: International interlaboratory study sequencing 16s rRNA from two genomic reference materials.  Manuscript is in prep.  

Scripts and pipeline presented have been run on **Mac session info** and **Linux session info**

### Requirements  
Unix based operating system is required for a number of the pipeline dependencies.  

#### Third Party Software and Packages
**Software** see links below for installation procedures  

- sratoolkit [http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=software)  
- BWA [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)  
- TMAP [https://github.com/nh13/TMAP](https://github.com/nh13/TMAP)  
- samtools [http://sourceforge.net/projects/samtools](http://sourceforge.net/projects/samtools)  
- Picard [http://picard.sourceforge.net/](http://picard.sourceforge.net/)
- GATK [http://www.broadinstitute.org/gatk](http://www.broadinstitute.org/gatk) 
- mothur [http://www.mothur.org/wiki/Download_mothur](http://www.mothur.org/wiki/Download_mothur)

**Programming languages and Packages**  

* Python  
  *  biopython  
* R  
  * ggplot2  
  * reshape2  
  * plyr 
* Bash 
  
## Procedure for reproducing sequence analysis  
1. [Data gathering](#dg)
2. [conservative bases](#cb)
3. [biologically variable positions](#bv)
4. [indiviudual gene copy sequences](#vs)

### General Notes
Directory structure  
  
    ccqm_16S/    
        ccqm_16S_stats/ 
           data/                   -- data files generated as part of statistical analysis
           
        ccqm_16S_pub/  
           tables/                 -- Rnw scripts for publication tables
           supplemental/           -- documents included as supplemental material for publication
           figures/                -- pdf of publication figures
           
        ccqm_16S_bioinformatics/
            bin/                   -- executables for third party software
            fastq-data/            -- raw sequence data
            resources/             -- additional input files required to run pipelines
            results/               -- output from bioinformatic pipelines
            scripts/               -- scripts written for bioinformatic analysis
            src/                   -- source code for third party software
            trace-data/            -- trace data obtained from GenBank
Users need to obtain third party software from link listed above and move the required compiled binaries to the *ccqm_16S_bioinformatics/bin* directory or create links.

1. AddOrReplaceReadGroups.jar
2. bcftools
3. bwa
4. CreateSequenceDictionary.jar
5. GenomeAnalysisTK.jar
6. mothur
7. prinseq-lite.pl
8. samtools
9. SortSam.jar
10. tmap



### [Data Gathering](id:dg)
* Sanger seqeuncing data - archived in the _GenBank Trace Archive_
	- Go to [http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?view=search](http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?view=search)  
	- Input search query string: `CENTER_NAME = "NIST" or CENTER_NAME = "LGC" or CENTER_NAME = "ISP" or CENTER_NAME = "ATCC"` (should return 495 total sequences)  
	- Select **fasta** and **qual**
	- Save query results  
	- Move download tar file to *trace-data* directory 
	- Uncompress tar file and generate fastq files for individual datasets using *trace_to_fastq.sh*.  Note *trace_to_fastq.sh* script requires the *fastaqualt_to_fastq.py* script  
* Next generation sequencing data - archived in the _GenBank SRA_   
	* *sra_to_fastq.sh* - script will download and dump datasets as fastq for next generation sequencing datasets generated in the study.
	* Additional steps included to split the NMIA data by barcodes using *mothur*.  The dataset is "dumped" as an sff file which is split by barcode converted into fasta and qual.  The *fastaqual_to_fastq.py* script is used to convert the files into fastq files for *E. coli* and *L. monocytogenes*.
	
### [Conservative Bases](id:cb)
* pipeline comparison - *ccqm_pipeline_comparison.sh*
	* requirements: in *resources directory*: *metadata.txt*, *Ecoli_16S_consensus.fasta*, and *Lmono_16S_consensus*; software requirements noted in script header. 	

### [Biologically Variable Positions](id:bv)

### [Indiviudual Gene Copy Sequences](id:vs)

