import os
import sys
import argparse
# Script for mergining MAFs
# Creates one aggregate header taken from column names acorss MAFs
# Writes out rows while mapping to header
# NAs and equivalents replaced with ''
# - sequenced samples header???

def process_datum(val):
    """ Cleans up datum. """
    try:
        vfixed = val.strip()
    except AttributeError:
        vfixed = ''
    if vfixed in ['', 'NA', 'N/A', None]:
        return ''
    return vfixed

def get_header(data_file):
    """
        Returns header from MAF
        Assumes file is tab-delimited
        Assumes header is first non-commented line
    """
    header = []
    with open(data_file, "r") as header_source:
        for line in header_source:
            if not line.startswith("#"):
                header = map(str.strip, line.rstrip().split('\t'))
                break
    return header

def merge_headers(data_filenames):
    header = []
    for fname in data_filenames:
        new_header = [hdr for hdr in get_header(fname) if hdr not in header]
        header.extend(new_header)
    return header

# does order of MAF columns matter (e.g HUGO symbol at front?)
def merge_input_mafs(input_mafs, output_mafs):
    merged_header = merge_headers(input_mafs) # .sort(enforced_order?)
    sequenced_samples = set()
    

    for input_maf in input_mafs:
        with open(input_maf, "r") as maf:
            file_header = get_header(input_maf)
            for line in maf:
                # add comment processing later (sequenced samples header?)
                if line.startswith("#"):
                    continue
                # map row values to current header columns
                mapped_row = dict(zip(file_header, map(lambda x: process_datum(x), line.rstrip("\n").split("\t"))))
                if mapped_row.get(TUMOR_SAMPLE_BARCODE):
                    sequenced_samples.add(mapped_row.get(TUMOR_SAMPLE_BARCODE))
                # full record for the merged MAF (if mapped_row does not contain a column in merged_header, it is blank)
                normalized_row = map(lambda x: mapped_row.get(x, ''), merged_header)
                rows_to_write.append("\t".join(normalized_row))
                     
    with open(output_maf, "w") as output_maf:
        output_maf.write("#sequenced_samples: %s" % (" ".join(sequenced_samples)))
        output_maf.write("\t".join(merged_header) + "\n")
        output_maf.write("\n".join(rows_to_write) + "\n")

def run(input_mafs, output_maf):
    input_maf_list = map(str.strip, input_mafs.split(','))
    for input_maf in input_maf_list:
        if not os.path.isfile(input_maf):
            print input_maf + " cannot be found. Exiting..."
            exit(1)
        
    merge_input_mafs(input_maf_list, output_maf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-mafs", help = "comma-deliited list of MAFs to merge", required = True)
    parser.add_argument("-o", "--output-maf", help = "output filename for merged MAF", required = True)

    args = parser.parse_args()
    input_mafs = args.input_mafs
    output_maf = args.output_maf
    
    run(input_mafs, output_maf)
