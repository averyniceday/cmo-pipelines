#!/usr/bin/env python

import argparse
import merge_mafs
import os
import subprocess
import sys

def setup_working_directory(temp_directory, center):
    working_directory = os.path.join(temp_directory, center)
    if os.path.isdir(working_directory):
        shutil.rmtree(working_directory)
    shutil.copytree(input_directory, working_directory)
    return working_directory
        
def get_mafs(input_directory):
    return [os.path.join(input_directory, filename] for filename in os.listdir(input_directory) if filename.endswith(".maf")]

def get_vcfs(input_directory):
    return [os.path.join(input_directory, filename] for filename in os.listdir(input_directory) if filename.endswith(".vcf")]
# ------------------------------------------------------------------------------------------------------------
# Functions relating to vcf2vcf call
FAILED_VCF2VCF = []

def generate_vcf2vcf_call(lib, input_maf, output_maf, vcf_tumor_id):
    vcf2vcf_call = "perl %s/vcf2vcf.pl --input-vcf %s --output-vcf %s --vcf-tumor-id %s --ref-fasta ~/.vep/homo_sapiens/95_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa" % (lib, input_maf, output_maf, vcf_tumor_id)
    return vcf2vcf_call
# given a directory, run vcf2vcf on each *.vcf file

def call_vcf2vcf(lib, input_directory):
    vcfs = get_vcfs(input_directory)
    for input_vcf in vcfs:
        # write output vcf to temp directory for later steps
        # prevents current directory from getting too messy
        output_vcf = "standardized." + input_vcf 
        vcf_tumor_id = get_vcf_tumor_id(input_vcf)
        vcf2vcf_call = generate_vcf2vcf_call(lib, input_vcf, output_vcf, vcf_tumor_id)
        vcf2vcf_status = subprocess.call(vcf2vcf_call, shell = True)
        # store failed vcf2vcf files for later logging
        if vcf2vcf_status != 0:
            FAILED_VCF2VCF.append(input_vcf)

def get_vcf_tumor_id(input_vcf):
    filename = input_vcf.rstrip(".vcf")
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



def generate_vcf2maf_call(lib, root_directory, temp_directory):
    vcf2maf_call =
    return vcf2maf_call
# ------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--lib", help = "directory containing required scripts (vcf2vcf, vcf2maf, annotator)", required = True)
    parser.add_argument("-i", "--input-directory", help = "directory with vcfs/mafs to process", required = True)
    parser.add_argument("-t", "--temp-directory", help = "temp directory for processing", required = True)
    parser.add_argument('-c', '--center', help = "Center name to be added", required = False)
    

    args = parser.parse_args()
    center = args.center
    lib = args.lib
    input_directory = args.input_directory
    temp_directory = args.temp_directory
    
    if not os.path.isdir(input_directory):
        print "Specified input directory %s does not exist, exiting..." % input_directory
        exit(1)
    if not os.path.isdir(temp_directory):
        print "Specified temp directory %s does not exist, exiting..." % temp_directory
        exit(1)

    # copy everything into a working directory and pre-clean
    working_directory = setup_working_directory(temp_directory, center)
    preclean_vcf.run(working_directory)

    # run vcf2vcf on all vcfs in directory
    # outputs are "standardized." + <vcf_name> + ".vcf" e.g standardized.10678314.vcf
    # delete all non-startswith "standardized" vcfs since we don't need them anymore
    call_vcf2vcf(lib, working_directory)

    # if vcf2maf does not work on mafs then...
    # run some maf standardizer on all mafs in the directory
    standardize_maf(lib, working_directory)
    
    # run vcf2maf on all vcfs in directory
    # remove all vcfs
    call_vcf2maf(lib, working_directory)
    
    mafs = get_mafs(working_directory)
    merge_mafs.run(','.join(mafs), center + "_maf.txt")
    
for vcf in directory
        # parse vcf and find out the sample id
        # unpaired samples - tumor sample id
        if sample name is TUMOR:
            filename is samplename
        if SAMPLE_ID:
            sample id is sample

if __name__ == '__main__':
    main()
