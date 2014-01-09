import HTSeq
import numpy
from matplotlib import pyplot
import itertools
import csv
import sys,argparse
import os.path

#########################################################################################
#                                                                                       #
#       HTSeq analysis.                                                                 #
#                                                                                       #
#       This script will take BAM and GTF input and analyse the coverage of genes       #
#       in the transcriptome data. Using functions provided by HTSeq python module      #
#                                                                                       #
#       HTSeq                                                                           #
#       http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html                   #
#                                                                                       #
#                                                                                       #
#                                                                                       #
#               (c) Ernst Thuer 2014                                                    #
#########################################################################################


# arguments for commandline input and help
####################################################
parser = argparse.ArgumentParser(description='Read BAM and GTF files, to produce readcount')
parser.add_argument('-bam',
                    dest='filename',
                    required = True,
                    help='Input a BAM file containing RNASeq data',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-gtf',
                    dest='gtf',
                    required = True,
                    help='Input a GTF file containing annotations',
                    metavar = 'FILE',
                    #type=lambda x: is_valid_file(parser,x)
                    )

parser.add_argument('-o',
                    dest='output',
                    required = False,
                    default='HTSeq_counts.csv',
                    help='Output a .csv file containing Gene names and counts, defaults to HTSeq_counts.csv',
                    metavar = 'FILE',
                    #type=argparse.FileType('w')
                    )

parser.add_argument('-feat',
                    dest='feature',
                    required = False,
                    default='CDS',
                    help='Feature that is targeted, from the GTF File. Usually CDS or exons',
                    metavar = 'FILE',
                    #type=argparse.FileType('w')
                    )

args = parser.parse_args()

#####################################################

#create empty count dict
counts = {}


# open input files with HTSeq reader 

try:
	bam_reader = HTSeq.BAM_Reader( "%s" %(args.filename) )
	gtf_file = HTSeq.GFF_Reader( "%s" %(args.gtf), end_included=True )
except:
	print "File not found exception"

#create an empty Genomic Array

try:
        exons = HTSeq.GenomicArrayOfSets( "auto", stranded=False ) #, typecode='O' 
except:
	print "HTseq array not loaded properly"
	

#populate the exons array with distinct feature type,  CDS would be a good default
	
for feature in gtf_file:
        if feature.type == args.feature:
                exons[feature.iv] += feature.name
                counts[ feature.name ] = 0
                
# analyse the bam input for intersect

             
for alnmt in bam_reader:
    if alnmt.aligned:
       intersection_set = None
       for iv2, step_set in exons[ alnmt.iv ].steps():
           if intersection_set is None:
                try:
                      intersection_set = step_set.copy()
                except:
                        print "intersection set failed"
           else:
              intersection_set.intersection_update( step_set )
       if len( intersection_set ) == 1:
               counts[ list(intersection_set)[0] ] += 1

# write to file 

with open ("%s" %(args.output),'w') as output:
        output_csv=csv.writer(output,delimiter='\t')
        for name in sorted(counts.keys() ) :
                output_csv.writerow([name,counts[name]])


# clean up
#args.filename.close()
#args.gtf.close()
#args.output.flush()
#args.output.close()
