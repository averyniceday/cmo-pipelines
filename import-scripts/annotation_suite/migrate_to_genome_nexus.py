#!/usr/bin/env python

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys

# ------------------------------------------------------------------------------------------------------------
# Functions for generating executable commands

def generate_vcf2vcf_call(lib, root_directory, temp_directory):
    vcf2vcf_call = 'python ' + lib + '/generate-clinical-subset.py --study-id=' + cancer_study_id + ' --clinical-file=' + source_directory + '/data_clinical.txt --filter-criteria="PATIENT_ID=' + patient_list + '" --subset-filename=' + destination_directory + "/subset_file.txt"
    return vcf2vcf_call


        print 'python standardize_mutation_data.py --input-directory [path/to/input/directory] --output-directory [path/to/output/directory] --center [default name for center] --sequence-source [WGS | WXS] --extensions [comma-delimited list of extensions]'

def generate_vcf2maf_call(lib, root_directory, temp_directory):
    vcf2maf_call = 'python ' + lib + '/vcf2maf.py --study-id=' + cancer_study_id + ' --clinical-file=' + source_directory + '/data_clinical.txt --filter-criteria="PATIENT_ID=' + patient_list + '" --subset-filename=' + destination_directory + "/subset_file.txt"
    return vcf2vcf_call

def generate_annotator_call(lib, root_directory, temp_directory):
    annotator_call = 'java -jar ' + annotator_jar + ' -f ' + destination_directory + '/data_mutations_extended.txt ' + '-o ' + destination_directory + '/' + ANNOTATED_MAF_FILE_PATTERN + ' -i mskcc'
    return annotator_call

# ------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lib", help = "directory containing required scripts (vcf2vcf, vcf2maf, annotator)", required = True)
    parser.add_argument("-r", "--root-directory", help = "root directory with vcfs to be converted/reannotated", required = True)
    parser.add_argument("-t", "--temp-directory", help = "temp directory for processing", required = True)
    parser.add_argument('-c', '--center', help = "Center name", required = False)
    parser.add_argument('-s', '--sequence-source', help = "Sequencing source (standard MAF field = 'Sequencing_Source'), e.g., WXS or WGS", required = False)

    args = parser.parse_args()
    lib = args.lib
    root_directory = args.root_directory
    temp_directory = args.temp_directory

    call_vcf2vcf(root_directory, lib)
    call_vcf2maf(root_directory, lib)
    call_annotator(root_directory, lib)
    merge_all_mafs(root_directory, lib)

if __name__ == '__main__':
    main()
