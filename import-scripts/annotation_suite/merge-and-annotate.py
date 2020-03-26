#!/usr/bin/env python

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import merge_mafs

# ------------------------------------------------------------------------------------------------------------
def annotate_maf(input_maf, output_maf):
    annotate_maf_call = generate_annotator_call(lib, input_maf, output_maf)
    subprocess.call(annotate_maf_call)

def generate_annotator_call(lib, root_directory, temp_directory):
    annotator_call = 'java -jar ' + annotator_jar + ' -f ' + destination_directory + '/data_mutations_extended.txt ' + '-o ' + destination_directory + '/' + ANNOTATED_MAF_FILE_PATTERN + ' -i mskcc'
    return annotator_call
# ------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lib", help = "directory containing required scripts (vcf2vcf, vcf2maf, annotator)", required = True)
    parser.add_argument("-i", "--input-mafs", help = "comma-delimited directories (distinguished by center) to converted", required = True)
    parser.add_argument("-o", "--output-maf", help = "temp directory for processing", required = True)
        
    args = parser.parse_args()
    lib = args.lib
    input_mafs = args.input_mafs
    output_maf = args.output_maf

    merge_mafs.run(input_mafs, output_maf + ".tmp") 
    annotate_maf(output_maf + ".tmp", output_maf)

if __name__ == '__main__':
    main()
