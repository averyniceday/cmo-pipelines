#!/usr/bin/env python

import argparse
import merge_mafs
import os
import shutil
import subprocess
import sys

# ------------------------------------------------------------------------------------------------------------
def generate_annotate_maf_call(annotator_jar, input_maf, output_maf):
    annotate_maf_call = "java -jar %s -f %s -o %s -i uniprot" % (annotator_jar, input_maf, output_maf)
    return annotate_maf_call

def annotate_maf_call(annotator, input_maf, output_maf):
    annotate_maf_call = generate_annotate_maf_call(annotator, input_maf, output_maf)
    annotate_maf_status = subprocess.call(annotate_maf_call, shell = True)
    if annotate_maf_status != 0:
        print "Failed to annotate %s." % input_maf
        exit(1)
    print "Successfully annotated %s, annotated MAF written to %s" % (input_maf, output_maf) 
# ------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--annotator-jar", help = "path to annotator jar", required = True)
    parser.add_argument("-i", "--input-mafs", help = "comma-delimited directories (distinguished by center) to converted", required = True)
    parser.add_argument("-o", "--output-maf", help = "temp directory for processing", required = True)
        
    args = parser.parse_args()
    annotator_jar = args.annotator_jar
    input_mafs = args.input_mafs
    output_maf = args.output_maf

    if not os.path.exists(annotator):
            print "%s does not exist. Exiting..." % (annotator) 
            exit(1)

    merged_maf = output_maf + ".tmp"
    merge_mafs.run(input_mafs, merged_maf) 
    annotate_maf(annotator_jar, merged_maf, output_maf)

if __name__ == '__main__':
    main()
