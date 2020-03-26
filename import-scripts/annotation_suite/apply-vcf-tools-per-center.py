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
# being called on a per file basis
# DOES REQUIRE samtools for a step -- is that step needed?
# does Tom ever use vcf2vcf and if so, does he provide bam files/fasta 

    'help!' => \$help,
    'man!' => \$man,
    'input-vcf=s' => \$input_vcf,
    'output-vcf=s' => \$output_vcf,
    'vcf-tumor-id=s' => \$vcf_tumor_id,
    'vcf-normal-id=s' => \$vcf_normal_id,
    'new-tumor-id=s' => \$new_tumor_id,
    'new-normal-id=s' => \$new_normal_id,
    'tumor-bam=s' => \$tumor_bam,
    'normal-bam=s' => \$normal_bam,
    'ref-fasta=s' => \$ref_fasta,
    'add-header=s' => \$add_header,
    'add-info=s' => \$add_info,
    'retain-info=s' => \$retain_info,
    'retain-format=s' => \$retain_format,
    'remap-chain=s' => \$remap_chain,
    'add-filters!' => \$add_filters
def generate_vcf2vcf_call(lib, root_directory, temp_directory):
    vcf2vcf_call = 'perl ' + lib + '/generate-clinical-subset.py --study-id=' + cancer_study_id + ' --clinical-file=' + source_directory + '/data_clinical.txt --filter-criteria="PATIENT_ID=' + patient_list + '" --subset-filename=' + destination_directory + "/subset_file.txt"
    return vcf2vcf_call


# ------------------------------------------------------------------------------------------------------------
# general call
'python standardize_mutation_data.py --input-directory [path/to/input/directory] --output-directory [path/to/output/directory] --center [default name for center] --sequence-source [WGS | WXS] --extensions [comma-delimited list of extensions]'

# assuming it does exactly what cyriac's script does minus the annotation
# one center per call because it's organized with all center files in one directory
# confirm whether sequence source always the same for all vcfs/mafs per center
# if not, pain in the butt - have to map sequence source per vcf/maf
#
# example call
'python vcf2maf.py --input-directory centerA --output-directory tmp/centerA --center centerA --sequence-source <?> --extensions .vcf, .txt/.maf'

def generate_vcf2maf_call(lib, root_directory, temp_directory):
    vcf2maf_call = 'python ' + lib + '/vcf2maf.py --study-id=' + cancer_study_id + ' --clinical-file=' + source_directory + '/data_clinical.txt --filter-criteria="PATIENT_ID=' + patient_list + '" --subset-filename=' + destination_directory + "/subset_file.txt"
    return vcf2vcf_call
# ------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------
def generate_merge_mafs_call(lib, input_mafs, output_maf):
    merge_mafs_call = 'python merge_mafs.py -i %s -o %s' % (input_mafs, output_maf)
    return merge_mafs_call     

# ------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lib", help = "directory containing required scripts (vcf2vcf, vcf2maf, annotator)", required = True)
    parser.add_argument("-i", "--input-directory", help = "comma-delimited directories (distinguished by center) to converted", required = True)
    parser.add_argument("-t", "--temp-directory", help = "temp directory for processing", required = True)
    parser.add_argument('-c', '--center', help = "Center name", required = False)
    parser.add_argument('-s', '--sequence-source', help = "Sequencing source (standard MAF field = 'Sequencing_Source'), e.g., WXS or WGS", required = False)

    sequencing_source - differs my panels
    args = parser.parse_args()
    lib = args.lib
    input_directories = args.input_directories
    temp_directory = args.temp_directory

    // sample id name either file name or the column name
    
    for directory in input_directories.split(","):
        # in subdir: temp_dir/center/normalized_vcfs1,2,3...
        call_vcf2vcf(directory, temp_directory, lib)
            for vcf in directory:
                # parse vcf and find out the sample id
                # unpaired samples - tumor sample id
                if sample name is TUMOR:
                    filename is samplename
                if SAMPLE_ID:
                    sample id is sample


        # in subdir: temp_dir/center/mafs1,2,3 (normalized vcfs gone)
        call_vcf2maf(directory, temp_directory, lib)
    
    for directory in temp_directory:
        maf_list = []
        for maf in directory:
            maf_list.append(maf)
    
    merge_all_mafs(lib, maf_list, output_maf)    
    call_annotator(output_maf, lib)
   
vcf2maf -> calls vcf2vcf
maf2maf -> maf2vcf 

if __name__ == '__main__':
    main()
