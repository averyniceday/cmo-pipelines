#!/usr/bin/env python

import argparse
import merge_mafs
import os
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

def generate_vcf2maf_call(lib, root_director
# ------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lib", help = "directory containing required scripts (vcf2vcf, vcf2maf, annotator)", required = True)
    parser.add_argument("-i", "--input-directory", help = "comma-delimited directories (distinguished by center) to converted", required = True)
    parser.add_argument("-t", "--temp-directory", help = "temp directory for processing", required = True)
    parser.add_argument('-c', '--center', help = "Center name", required = False)
    
    # ignore for now
    # parser.add_argument('-s', '--sequence-source', help = "Sequencing source (standard MAF field = 'Sequencing_Source'), e.g., WXS or WGS", required = False)

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
