################################################################################
## Script for parsing SAM file 
##
## Created by: Nate Olson
## October 18, 2012
## Modified July, 2013
##
## parsing based on SAM format specifications v1.4-r985, September 7, 2011
################################################################################


#!/usr/bin/env python

__doc__ = """
SYNOPSIS

    pileup_parse.py input.pileup output.tsv

DESCRIPTION

    Converts pileup file to a tab separated
    file with individual bases as rows in the document.
    Only positions where the reference is amgibious are written to file.

USAGE
    Input file must be labeled .pileup for script to function
    
    output filter - ***need to incorporate this function***
        -ft Filter type - use to state whether to retain or exclude
                          specified filters
        -F filter - specific filters to keep or remove can be any
                    character in a pileup base or qual file, range of
                    base positions in relation to the reference as well
                    as any of the following:
                        location: start middle end
                        read direction: forward reverse

AUTHOR

        Nate Olson <nolson@nist.gov>

"""

__author__ = "Nate Olson"
__version__ = '0.1'

#importing required modules
import sys
import re
from Bio import SeqIO

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] == "help" or sys.argv[1] == "-h":
            print help()
        elif len(sys.argv) < 2:
            print("Input and Output files not specified" % help())
        else:
		sam_file = sys.argv[1]
		ref_file = sys.argv[2]
		output_file = sys.argv[3]
		sam_expand(sam_file, ref_file, output_file)
        
    else:
        print("No arguments were provided see documentation below:\n\n %s"  %help())

def help():
    return globals()['__doc__']


def sam_expand(sam_file, ref_file, output_file):
	# reading sam input file
	sam = open(sam_file, 'r')

	#reading reference fasta file
	ref = {}
	for seq_record in SeqIO.parse(ref_file,"fasta"):
		ref[seq_record.id] = seq_record.seq

	# creating output file
	samTOcsv = open(output_file,'w')

	# write to file column headers
	samTOcsv.write('read_name,refseq_name,ref_value, map_qual,cigar_value,base_call,quality_score,read_position, ref_position\n')


	# iterating through each line of the sam file
	for line in sam:
	    if not line.startswith("@"): # skip elements that start with @
		lineSplit = line.split('\t')
		# skip elements with a flag score (second element in the array) of 4 (reads not mapped)
		if lineSplit[2] != 4 and int(lineSplit[4]) > 10:
		    # expand cigar (element 5 of sub arrays)
		    # array of letters in the cigar
		    letters = re.findall('[A-Z]',lineSplit[5])
		    # array of numbers in the cigar as intigers 
		    numbers = [int(x) for x in re.findall(r'\d+', lineSplit[5])]
		    expCigar = ''
		    for i in range(0,len(letters)):
		            expCigar = expCigar + numbers[i]*letters[i]

		    # writing to output file individual line for each base
		    # split cigar, sequence, and qual into arrays
		    cigar = list(expCigar) # cigar
		    sequence = list(lineSplit[9]) # sequence
		    qual = list(lineSplit[10]) # qual
			    
		    #reference sequence
		    ref_seq = ref[lineSplit[2]]
		    # interate through each position of each of the reads creating a line for each
		    # setting initial ref position value
		    # for offsetting sequence and qual looping position for deletion in the read
		    seq_qual_pos = 0
		    ref_pos = int(lineSplit[3]) - 1
		    ref_position = ref_pos

		    #iterating through the cigar
		    for i in range(0,len(cigar)):
			# setting values and increasing counts based on cigar value
			if cigar[i] == 'D':
			    # setting values
			    seq_value = '*'
		  	    qual_value = '*'
			    ref_value = ref_seq[ref_pos]
			    ref_position = ref_pos
			    # increasing counts
			    ref_pos = ref_pos + 1 
					               	
			elif cigar[i] == "I":
			    # setting values
			    seq_value = sequence[seq_qual_pos]
			    qual_value = ord(qual[seq_qual_pos]) - 33
			    ref_value = '*'
			    if ref_position % 1 != 0: 
					''' incrementing by 0.05 to incorporate 
					    insertions, approach only works 
					    for 20 insertions '''
				ref_position = ref_position + 0.05
		            else:
				ref_position = ref_pos + 0.05
			    # increasing counts
			    seq_qual_pos = seq_qual_pos + 1
					
			elif cigar[i] == "M":
			    # setting values 
			    seq_value = sequence[seq_qual_pos]
			    qual_value = ord(qual[seq_qual_pos]) - 33
			    ref_value = ref_seq[ref_pos]
			    ref_position = ref_pos
			    # increasing counts
			    ref_pos = ref_pos + 1
			    seq_qual_pos = seq_qual_pos + 1
					
			elif cigar[i] == "S":
	    		    seq_qual_pos = seq_qual_pos + 1
			    continue
			else:
			    continue

			# writing values to output file - only outputs positions where the reference is an ambiguous base
			if ref_value not in ["A","C","G","T","*"] and cigar[i] == "M":
				samTOcsv.write(lineSplit[0] + ',' + 
						lineSplit[2]  + ',' + 
						ref_value + ',' + 
						lineSplit[4] + ',' +  
						cigar[i] + ',' + 
						seq_value + ',' + 
						str(qual_value) + ',' + 
						str(seq_qual_pos) + ',' + 
						str(ref_position) + '\n')
				int(seq_qual_pos)
	samTOcsv.close()
	sam.close()


if __name__ == '__main__':
    main()
    sys.exit()
