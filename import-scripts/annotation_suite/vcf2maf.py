# imports
import sys
import os
import optparse
import re

# ---------------------------- GLOBALS ----------------------------
MAF_HEADER = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", 
        "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", 
        "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", 
        "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1", 
        "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2", "Verification_Status", 
        "Validation_Status", "Mutation_Status", "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score", 
        "BAM_File", "Sequencer", "t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count"]
VCF_HEADER = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]
DEFAULT_NCBI_BUILD = 'GRCh37'
DEFAULT_STRAND = '+'
DEFAULT_VALIDATION_STATUS = 'Unknown'
DEFAULT_VERIFICATION_STATUS = 'Unknown'
DEFAULT_MUTATION_STATUS = 'Somatic'

MERGED_MAF_DATA_KEYS = set()

VALID_VARIANT_CLASSIFICATIONS = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", 
        "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", 
        "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region"]
VALID_VARIANT_TYPES = ["SNP", "DNP", "TNP", "ONP", "DEL", "INS"]
MUTATION_FILTER = ["3'Flank", "3'UTR", "5'Flank", "5'UTR", "IGR", "Intron", "RNA", "Silent"]
MUTATION_STATUS_FILTER = ["LOH", "Wildtype"]
VALID_CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']

# --------------------- CGI DATA GLOBALS -------------------------
CGI_VARIANT_CLASS_FILTER = ["INTRON", "TSS-UPSTREAM", "UTR5", "UTR3", "UTR", "SPAN5", "SPAN3", "SPAN", "SYNONYMOUS", "IGR", "NO-CHANGE", "UPSTREAM"]
CGI_VARIANT_CLASS_EXTRA_HANDLING = ["UTR", "SPAN", "FRAMESHIFT"]
CGI_INDEL_VARIANT_CLASSES = ["INSERT", "DELETE", "INSERT+", "DELETE+", "FRAMESHIFT"]
# NO DIRECT MAPPING FOR: UTR, SPAN, FRAMESHIFT - These require additional data either from other variant classes 
# if datum contains comma-delimited variant classes (i.e., NO-CHANGE,DELETE,FRAMESHIFT for one record)
# or value of variant type 
# > If only UTR is present then use IGR as default
# > If only SPAN is present then use IGR as default
# > If only FRAMESHIFT is present and variant type is INS/DEL then use Frame_Shift_Ins, Frame_Shift_Del.
#   Otherwise, if variant type is not INS/DEL then use Missense_Mutation
CGI_VARIANT_CLASS_MAP = {"INTRON":"Intron", 
        "TSS-UPSTREAM":"5'Flank",
        "UPSTREAM":"5'Flank",
        "UTR5":"5'UTR",
        "UTR3":"3'UTR",
        "SPAN5":"5'UTR",
        "SPAN3":"3'UTR",
        "SYNONYMOUS":"Silent",
        "MISSTART":"Translation_Start_Site",
        "DONOR":"Splice_Site",
        "ACCEPTOR":"Splice_Site",
        "DISRUPT":"Splice_Site",
        "NO-CHANGE":"Silent",
        "MISSENSE":"Missense_Mutation",
        "NONSENSE":"Nonsense_Mutation",
        "NONSTOP":"Nonstop_Mutation",
        "DELETE":"In_Frame_Del",
        "INSERT":"In_Frame_Ins",
        "DELETE+":"Frame_Shift_Del",
        "INSERT+":"Frame_Shift_Ins"}
# ----------------------------------------------------------------


def get_standardized_sample_id(sample_id):
        try:
                standardized_sample_id = get_target_standardized_sample_id(sample_id)
        except AttributeError:
                standardized_sample_id = get_other_standardized_sample_id(sample_id)            
        return standardized_sample_id

def get_target_standardized_sample_id(sample_id):
        try:
                m = re.search('(TARGET-[0-9]*-[0-9A-Z]*-[0-9]*)[A-Z.]*-[0-9]*', sample_id)
                standardized_sample_id = m.group(1)
        except AttributeError:
                m = re.search('(TARGET-[0-9]*-[0-9A-Z]*-[0-9A-Z.]*)', sample_id)
                standardized_sample_id = m.group(1) + "-01"
        return standardized_sample_id

def get_other_standardized_sample_id(sample_id):
        try:
                m = re.search('([HTMCP|BLGSP]+-\d*-\d*-\d*-[0-9A-Z]*-\d*)', sample_id)
                standardized_sample_id = m.group(1)
        except AttributeError:
                if len(sample_id.split('-')) == 1:
                        standardized_sample_id = sample_id + '-01'
                else:
                        print 'Could not resolve unknown sample id pattern: ' + sample_id
                        return sample_id
        return standardized_sample_id


def get_header(filename):
        """ Returns the file header. """
        data_file = open(filename, 'rU')
        filedata = [x for x in data_file.readlines() if not x.startswith('#')]
        header = map(str.strip, filedata[0].split('\t'))
        data_file.close()
        return header

def process_datum(value):
        try:
                vfixed = value.strip()
        except AttributeError:
                vfixed = ''
        if vfixed in ['', 'NA', 'N/A', None, '.', '?']:
                return ''
        return vfixed

def get_maf_key(maf_data):
        sample_id = maf_data['Tumor_Sample_Barcode']
        chromosome = maf_data['Chromosome']
        start_pos = maf_data['Start_Position']
        end_pos = maf_data['End_Position']
        ref_allele = maf_data['Reference_Allele']
        alt_allele = maf_data['Tumor_Seq_Allele2']
        return (sample_id, chromosome, start_pos, end_pos, ref_allele, alt_allele)


def resolve_tumor_seq_allele(data, ref_allele):
        tum_seq_allele1 = ''
        tum_seq_allele2 = ''
        if 'mutated_to_allele' in data.keys():
                return data['mutated_to_allele']

        tum_seq_allele1_columns = ['Tumor_Seq_Allele1', 'TumorSeq_Allele1']
        tum_seq_allele2_columns = ['Tumor_Seq_Allele2', 'TumorSeq_Allele2']
        for i,column in enumerate(tum_seq_allele1_columns):             
                if column in data.keys():
                        tum_seq_allele1 = process_datum(data[column])
                        tum_seq_allele2 = process_datum(data[tum_seq_allele2_columns[i]])

                        # if at least one is not empty then exit for-loop
                        if tum_seq_allele1 != '' or tum_seq_allele2 != '':
                                break

        # if temp tum seq alleles are not equal then set tum_seq_allele that doesn't equal ref_allele
        # otherwise there is probably something wrong with the data
        if tum_seq_allele1 == '' and tum_seq_allele2 == '':
                return ''
        if tum_seq_allele1 == ref_allele:
                return tum_seq_allele2
        else:
                return tum_seq_allele1


def resolve_variant_type(data, ref_allele, tum_seq_allele):
        for column in ['Variant_Type', 'VariantType', 'mut_type', 'mutation_type']:
                if column in data.keys():
                        variant_type = process_datum(data[column].upper())
                        break
        if variant_type == '1':
                return 'SNP'
        elif variant_type == '2':
                return 'INS'
        elif variant_type == '3':
                return 'DEL'

        # if variant type is empty or not valid then try resolving it based on ref allele and tum seq allele values
        if variant_type not in VALID_VARIANT_TYPES:             
                if len(ref_allele) > len(tum_seq_allele):
                        variant_type = 'DEL'
                elif len(ref_allele) < len(tum_seq_allele):
                        variant_type = 'INS'
                elif len(ref_allele) == len(tum_seq_allele):
                        if len(ref_allele) == 1:
                                variant_type = 'SNP'
                        elif len(ref_allele) == 2:
                                variant_type = 'DNP'
                        elif len(ref_allele) == 3:
                                variant_type = 'TNP'
                        else:
                                variant_type = 'ONP'
        if variant_type == '':
                print 'Could not salvage variant type from alleles'
                print 'Ref allele=' + ref_allele, 'Tumor Allele=' + tum_seq_allele
        # if variant type still not in VALID_VARIANT_TYPES then must be SUB/SNV
        if variant_type == 'SNV':
                variant_type = 'SNP'
        elif variant_type.upper() == 'SUB':
                variant_type = 'ONP'
        return variant_type

def resolve_complex_variant_classification(data, variant_class_list, variant_type):     
        # MISSTART variant classes take precedence over other CGI variant classes
        if 'MISSTART' in variant_class_list:
                return CGI_VARIANT_CLASS_MAP['MISSTART']

        # if length of variant class list is 1 and variant class is in list of filtered CGI variants then return
        # direct mapping of variant class or IGR by default, as in the cases of SPAN and UTR variant classes
        if len(variant_class_list) == 1 and variant_class_list[0] in CGI_VARIANT_CLASS_FILTER:
                return CGI_VARIANT_CLASS_MAP.get(variant_class_list[0], 'IGR')

        # check if any CGI variant classes that we normally filter are present in given variant_class_list
        filtered_cgi_var_classes = [var_class for var_class in variant_class_list if var_class in CGI_VARIANT_CLASS_FILTER]
        if len(filtered_cgi_var_classes) > 0:
                # if filtered CGI variant classes are present then return the first one that can be directly mapped to 
                # a standard variant classification
                for var_class in filtered_cgi_var_classes:
                        if var_class in CGI_VARIANT_CLASS_MAP.keys():
                                return CGI_VARIANT_CLASS_MAP[var_class]

        # Map var classes to standard variant classifications or use orig var class name if not found.
        # Variant classes that do not map directly include FRAMESHIFT, SPAN, UTR
        variant_class_candidates = map(lambda x: CGI_VARIANT_CLASS_MAP.get(x, x), variant_class_list)

        # splice sites take precedence over other variant classifications
        if 'Splice_Site' in variant_class_candidates:
                return 'Splice_Site'

        # if variant type is INS/DEL or INSERT/DELETE/FRAMESHIFT/INSERT+/DELETE+ in input variant_class_list then 
        # return either In_Frame_Ins/Del or Frame_Shift_Ins/Del
        indel_variant_classes = [var_class for var_class in variant_class_list if var_class in CGI_INDEL_VARIANT_CLASSES]
        if variant_type in ['INS', 'DEL'] and len(indel_variant_classes) > 0:
                for var_class in indel_variant_classes:
                        if var_class == 'FRAMESHIFT':
                                return 'Frame_Shift_' + variant_type.title()
                        else:
                                return CGI_VARIANT_CLASS_MAP[var_class]

        # if indel variant classes present but variant_type is not INS or DEL - variant type will need to be fixed later
        for var_class in indel_variant_classes:
                if var_class != 'FRAMESHIFT':
                        return CGI_VARIANT_CLASS_MAP[var_class]

        # if variant class is not an indel then variant class is either Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation, or NO-CHANGE/Silent
        # if variant type is INS/DEL then return Frame_Shift_Ins/Del if Nonsense_Mutation or Nonstop_Mutation in variant class candidates
        # or In_Frame_Ins/Del if Missense_Mutation in candidates - Nonsense/Nonstop mutations take precedence over Missense
        if variant_type in ['INS', 'DEL']:
                if 'Nonsense_Mutation' in variant_class_candidates or 'Nonstop_Mutation' in variant_class_candidates:
                        return 'Frame_Shift_' + variant_type.title()
                elif 'Missense_Mutation' in variant_class_candidates:
                        return 'In_Frame_' + variant_type.title()

        # if variant type not INS/DEL then must be SNP/SNV/SUB
        # Return Nonsense_Mutation, Nonstop_Mutation, or Missense_Mutation with Nonsense/Nonstop taking precedence
        for var_class in ['Nonstop_Mutation', 'Nonsense_Mutation','Missense_Mutation']:
                if var_class in variant_class_candidates:
                        return var_class

        # if variant class is empty but variant type is SNP then return missense mutation
        if variant_type in ['SNP', 'DNP', 'TNP', 'ONP']:
                return 'Missense_Mutation'

        # if variant class can't be resolve by this point then arbitrarily return first in list that can be mapped directly
        return variant_class_candidates[0]


def resolve_variant_classification(data, variant_type, ref_allele, alt_allele):
        variant_class = ''
        for column in ['Variant_Classification', 'class', 'Transcript architecture around variant']:
                if column in data.keys():
                        variant_class = process_datum(data[column])
                        break
        # if variant classification is valid then return, else try to resolve value
        if variant_class in VALID_VARIANT_CLASSIFICATIONS:
                return variant_class
        
        # if empty string then assume missense or indel - let annotator resolve correct variant class
        if variant_class == '':
                # if indel then determine whether in frame or out of frame
                if variant_type in ['INS', 'DEL']:
                        if variant_type == 'INS':
                                in_frame_variant = (len(alt_allele) % 3 == 0)
                        else:
                                in_frame_variant = (len(ref_allele) % 3 == 0)
                        if in_frame_variant:
                                variant_class = 'In_Frame_' + variant_type.title()
                        else:
                                variant_class = 'Frame_Shift_' + variant_type.title()
                else:
                        # let annotator figure it out from the VEP response
                        variant_class = 'Missense_Mutation'
        else:
                variant_class_list = re.split('[\\|, ]', variant_class)
                variant_class = resolve_complex_variant_classification(data, variant_class_list, variant_type)
        return variant_class

def resolve_start_position(data):
        start_pos = ''
        for column in ['Start_Position', 'Start_position', 'POS', 'chromosome_start']:
                if column in data.keys():
                        start_pos = process_datum(data[column])
                        break
        return start_pos

def resolve_end_position(data, start_pos, variant_type, ref_allele):
        end_pos = ''
        for column in ['End_Position', 'End_position', 'chromosome_end']:
                if column in data.keys():
                        end_pos = process_datum(data[column])
                        break
        # if insertion then end pos is start pos + 1
        if variant_type == 'INS' or ref_allele == '-':
                try:
                        end_pos = str(int(start_pos)+1)
                except ValueError:
                        print data
                        sys.exit(2)

        # resolve end pos from ref allele length if empty string
        if end_pos == '':
                end_pos = str(int(start_pos) + len(ref_allele) - 1)

        return end_pos

def resolve_ref_allele(data):
        ref = ''
        for col in ['Reference_Allele', 'reference_genome_allele']:
                if col in data.keys():
                        ref = data[col]
        return ref

def resolve_variant_allele_data(data, maf_data):        
        ref_allele = resolve_ref_allele(data)
        tum_seq_allele = resolve_tumor_seq_allele(data, ref_allele)
        variant_type = resolve_variant_type(data, ref_allele, tum_seq_allele)

        # fix ref allele and tum seq allele for INS or DEL variant types
        if variant_type == 'INS':
                ref_allele = '-'
        elif variant_type == 'DEL':
                tum_seq_allele = '-'

        variant_class = resolve_variant_classification(data, variant_type, ref_allele, tum_seq_allele)
        # fix variant type just in case it was missed before
        if variant_class.endswith('Ins') and variant_type != 'INS':
                variant_type = 'INS'
        elif variant_class.endswith('Del') and variant_type != 'DEL':
                variant_type = 'DEL'

        # resolve start and end positions
        start_pos = resolve_start_position(data)
        end_pos = resolve_end_position(data, start_pos, variant_type, ref_allele)

        maf_data['Variant_Classification'] = variant_class
        maf_data['Variant_Type'] = variant_type
        maf_data['Reference_Allele'] = ref_allele
        maf_data['Tumor_Seq_Allele2'] = tum_seq_allele
        maf_data['Start_Position'] = start_pos
        maf_data['End_Position'] = end_pos

        # resolve tumor seq allele 1 - default is set as ref allele
        maf_data['Tumor_Seq_Allele1'] = ref_allele
        if 'tumour_genotype' in data.keys():
                tum_seq_allele1 = data['tumour_genotype'].split('/')[0]
                if tum_seq_allele1 != '':
                        maf_data['Tumor_Seq_Allele1'] = tum_seq_allele1

        return maf_data

def resolve_counts_data(data, maf_data):
        # get defaults from expected maf columns first
        t_depth = process_datum(data.get('t_depth',''))
        t_ref_count = process_datum(data.get('t_ref_count',''))
        t_alt_count = process_datum(data.get('t_alt_count',''))
        n_depth = process_datum(data.get('n_depth', ''))
        n_ref_count = process_datum(data.get('n_ref_count',''))
        n_alt_count = process_datum(data.get('n_alt_count',''))
        if 'Tumor_ReadCount_Total' in data.keys():
                # tumor counts
                t_depth = process_datum(data.get('Tumor_ReadCount_Total',''))
                t_ref_count = process_datum(data.get('Tumor_ReadCount_Ref',''))
                t_alt_count = process_datum(data.get('Tumor_ReadCount_Alt',''))
                # normal counts
                n_depth = process_datum(data.get('Normal_ReadCount_Total',''))
                n_ref_count = process_datum(data.get('Normal_ReadCount_Ref',''))
                n_alt_count = process_datum(data.get('Normal_ReadCount_Alt',''))
        if 'TTotCovVal' in data.keys():
                # tumor counts
                t_depth = process_datum(data.get('TTotCovVal',''))
                t_alt_count = process_datum(data.get('TVarCovVal',''))
                if not t_depth in ['', 'None'] and not t_alt_count in ['', 'None']:
                        t_ref_count = str(int(t_depth) - int(t_alt_count))
                # normal counts
                n_depth = process_datum(data.get('NTotCovVal',''))
                n_alt_count = process_datum(data.get('NVarCovVal',''))
                if not n_depth in ['', 'None'] and not n_alt_count in ['', 'None']:
                        n_ref_count = str(int(n_depth)-int(n_alt_count))
        elif 'Ref_Allele_Coverage' in data.keys():
                # tumor counts
                t_ref_count = process_datum(data.get('Ref_Allele_Coverage',''))
                t_alt_count = process_datum(data.get('Tumor_Seq_Allele1_Coverage',''))
                if t_ref_count != '' and t_alt_count != '':
                        t_depth = str(int(t_ref_count) + int(t_alt_count))
                # normal counts
                n_ref_count = process_datum(data.get('Normal_Ref_Allele_Coverage', ''))
                n_alt_count = process_datum(data.get('Normal_Seq_Allele1_Coverage', ''))
                if n_ref_count != '' and n_alt_count != '':
                        n_depth = str(int(n_ref_count) + int(n_alt_count))
        if 'total_read_count' in data.keys():
                # tumor counts
                t_depth = process_datum(data.get('total_read_count',''))
                t_alt_count = process_datum(data.get('mutant_allele_read_count',''))
                if not t_depth in ['', 'None'] and not t_alt_count in ['', 'None']:
                        t_ref_count = str(int(t_depth) - int(t_alt_count))


        maf_data['t_depth'] = t_depth
        maf_data['t_ref_count'] = t_ref_count
        maf_data['t_alt_count'] = t_alt_count
        maf_data['n_depth'] = n_depth
        maf_data['n_ref_count'] = n_ref_count
        maf_data['n_alt_count'] = n_alt_count

        return maf_data

def resolve_hugo_symbol(data):
        hugo_symbol = ''
        for column in ['Hugo_Symbol', 'HugoSymbol', 'Gene Symbol']:
                if column in data.keys():
                        hugo_symbol = process_datum(data[column].split('|')[0])
                        break
        if hugo_symbol == '':
                hugo_symbol = 'Unknown'
        return hugo_symbol

def resolve_match_normal_sample_barcode(data):
        barcode = ''
        for column in ['Match_Normal_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'matched_sample_id']:
                if column in data.keys():
                        barcode = process_datum(data[column])
                        break
        return barcode

def resolve_chromosome(data):
        chromosome = ''
        for column in ['Chromosome', 'chromosome', 'CHROM']:
                if column in data.keys():
                        chromosome = process_datum(data[column]).replace('chr','')
                        break   
        return chromosome.split('_')[0]

def resolve_mutation_status(data):
        mutation_status = process_datum(data.get('Mutation_Status', DEFAULT_MUTATION_STATUS))
        if mutation_status == '':
                mutation_status = DEFAULT_MUTATION_STATUS
        return mutation_status

def resolve_center_name(data, center_name):
        center = process_datum(data.get('Center', ''))
        if center == '':
                center = center_name
        return center

def resolve_sequence_source(data, sequence_source):
        # convert sequence strategy to sequence source values if column present in data
        if 'sequencing_strategy' in data.keys():
                if data['sequencing_strategy'] == '1':
                        return 'WGS'
                elif data['sequencing_strategy'] == '3':
                        return 'WXS'
        # fall back on parsing sequence source column or fall on default
        seq_source = process_datum(data.get('Sequence_Source', ''))
        if seq_source == '':            
                seq_source = sequence_source
        return seq_source

def init_maf_record():
        # init maf record
        maf_data = dict(zip(MAF_HEADER, ['' for column in MAF_HEADER]))

        # set defaults
        maf_data['Strand'] = DEFAULT_STRAND
        maf_data['Validation_Status'] = DEFAULT_VALIDATION_STATUS
        return maf_data 

def resolve_ncbi_build(data):
        ncbi_build = ''
        for col in ['NCBI_Build', 'assembly_version']:
                if col in data.keys():
                        ncbi_build = process_datum(data.get(col, ''))
                        break
        if ncbi_build == '':
                ncbi_build = DEFAULT_NCBI_BUILD
        return ncbi_build

def resolve_dbsnp_rs(data):
        dbsnp_rs = ''
        for column in ['dbSNP_RS', 'dbSNP rsID']:
                if column in data.keys():
                        dbsnp_rs = process_datum(data[column])
                        break
        return dbsnp_rs

def resolve_verification_status(data):
        ver_status = ''
        for col in ['Verification_Status', 'verification_status']:
                if col in data.keys():
                        ver_status = data[col]
        if ver_status == '1':
                ver_status = 'Verified'
        elif ver_status == '':
                ver_status = DEFAULT_VERIFICATION_STATUS
        return ver_status

def resolve_sequencer(data):
        sequencer = ''
        for col in ['Sequencer', 'platform']:
                if col in data.keys():
                        sequencer = process_datum(data.get(col,''))
                        break
        if sequencer == '60':
                sequencer = 'Illumina HiSeq 2000'
        return sequencer

def resolve_validation_method(data):
        validation_method = ''
        for col in ['Validation_Method', 'verification_platform']:
                if col in data.keys():
                        validation_method = process_datum(data.get(col,''))
                        break
        if validation_method == '60':
                validation_method = 'Illumina HiSeq (RNAseq)'
        elif validation_method == '6':
                validation_method = '454 sequencing'
        elif validation_method == '67':
                validation_method = 'Ion Torrent PGM'
        elif validation_method == '-888':
                validation_method = 'NA'
        return validation_method

def resolve_match_norm_seq_alleles(data, maf_data):
        norm_allele1 = ''
        for col in ['Match_Norm_Seq_Allele1', 'mutated_to_allele']:
                if col in data.keys():
                        norm_allele1 = process_datum(data.get(col,''))
                        break

        norm_allele2 = ''
        for col in ['Match_Norm_Seq_Allele2', 'control_genotype']:
                if col in data.keys():
                        norm_allele2 = process_datum(data.get(col,''))
                        if col == 'control_genotype':
                                norm_allele2 = norm_allele2.split('/')[1]
                        break
        maf_data['Match_Norm_Seq_Allele1'] = norm_allele1
        maf_data['Match_Norm_Seq_Allele2'] = norm_allele2
        return maf_data

def create_maf_record(data, center_name, sequence_source):
        # init maf record
        maf_data = init_maf_record()

        # set easy to resolve values    
        try:
                for col in ['Tumor_Sample_Barcode', 'analyzed_sample_id']:
                        if col in data.keys():
                                sid_column = col
                                break
                maf_data['Tumor_Sample_Barcode'] = get_standardized_sample_id(data[sid_column])
        except AttributeError:
                print data
                sys.exit(2)
        maf_data['Matched_Norm_Sample_Barcode'] = resolve_match_normal_sample_barcode(data)     
        maf_data['Hugo_Symbol'] = resolve_hugo_symbol(data)
        maf_data['Entrez_Gene_Id'] = process_datum(data.get('Entrez_Gene_Id','').split('|')[0])
        maf_data['NCBI_Build'] = resolve_ncbi_build(data)       
        maf_data['Chromosome'] = resolve_chromosome(data)
        maf_data['dbSNP_RS'] = resolve_dbsnp_rs(data)
        maf_data['Sequencing_Phase']  = process_datum(data.get('Sequencing_Phase', ''))
        maf_data['Sequence_Source'] = resolve_sequence_source(data, sequence_source)
        maf_data['Validation_Method'] = resolve_validation_method(data)
        maf_data['Center'] = resolve_center_name(data, center_name)
        maf_data['Verification_Status'] = resolve_verification_status(data)
        maf_data['Mutation_Status'] = resolve_mutation_status(data)
        maf_data['Center'] = resolve_center_name(data, center_name)
        maf_data['Sequencer'] = resolve_sequencer(data)

        if maf_data['Verification_Status'] == 'Verified':
                maf_data['Validation_Status'] = 'Valid'


        # resolve counts, variant alleles, and norm alleles data
        resolve_counts_data(data, maf_data)
        resolve_variant_allele_data(data, maf_data)
        resolve_match_norm_seq_alleles(data, maf_data)
        return maf_data

def get_vcf_sample_id(filename):
        sample_id = get_standardized_sample_id(re.split('[_ .]', os.path.basename(filename))[0])
        return sample_id

def resolve_vcf_allele(vcf_data):
        ref_allele = ''
        alt_allele = ''

        if vcf_data['ALT'] in ['<DEL>', '<DUP>', '<INV>', '<TRA>']:
                if vcf_data['REF'] == 'N' or vcf_data['REF'] == '':
                        ref_allele = get_consensus_sequence_from_vcf_info(vcf_data['INFO'])
                else:
                        ref_allele = vcf_data['REF']

                if ref_allele != 'N' and ref_allele != '':
                        if vcf_data['ALT'] == '<DEL>':
                                alt_allele = '-'                
                        if vcf_data['ALT'] == '<INV>':
                                alt_allele = ref_allele[::-1]
                        elif vcf_data['ALT'] == '<DUP>':
                                alt_allele = ref_allele*2
        else:
                ref_allele = process_datum(vcf_data['REF'].split(',')[0])
                alt_allele = process_datum(vcf_data['ALT'].split(',')[0])
                if ref_allele == '':
                        ref_allele = '-'
                if alt_allele == '':
                        alt_allele = '-'
        return ref_allele,alt_allele

def get_consensus_sequence_from_vcf_info(vcf_info):
        sequence = ''
        if 'CONSENSUS' in vcf_info:
                for data in vcf_info.split(';'):
                        if data.startswith('CONSENSUS'):
                                sequence = data.split('=')[1]
        return sequence

def resolve_vcf_variant_type(ref_allele, alt_allele):
        variant_type = ''
        # first check if indel
        if ref_allele == '-' or len(ref_allele) < len(alt_allele):
                variant_type = 'INS'
        elif alt_allele == '-' or len(alt_allele) < len(ref_allele):
                variant_type = 'DEL'
        else:
                # check whether variant type is type of polymorphism
                if len(ref_allele) == len(alt_allele):
                        if len(ref_allele) == 1:
                                variant_type = 'SNP'
                        elif len(ref_allele) == 2:
                                variant_type = 'DNP'
                        elif len(ref_allele) == 3:
                                variant_type = 'TNP'
                        else:
                                variant_type = 'ONP'
        # if variant type is still empty then report
        if variant_type == '':
                print 'Could not salvage variant type from alleles'
                print 'Ref allele=' + ref_allele, 'Alt allele=' + alt_allele
        return variant_type


def resolve_vcf_counts_data(vcf_data, maf_data):
        info_data = vcf_data['INFO'].split(';')

        n_alt_count = ''
        n_ref_count = ''
        n_depth = ''
        t_alt_count = ''
        t_ref_count = ''
        t_depth = ''    
        for data in info_data:
                if data.startswith('NORMALT'):
                        n_alt_count = process_datum(data.split('=')[1])
                elif data.startswith('NORMREF'):
                        n_ref_count = process_datum(data.split('=')[1])
                elif data.startswith('TUMALT'):
                        t_alt_count = process_datum(data.split('=')[1])
                elif data.startswith('TUMREF'):
                        t_ref_count = process_datum(data.split('=')[1])
        if n_ref_count != '' and n_alt_count != '':
                n_depth = str(int(n_ref_count) + int(n_alt_count))
        if t_ref_count != '' and t_alt_count != '':
                t_depth = str(int(t_ref_count) + int(t_alt_count))
        maf_data['t_depth'] = t_depth
        maf_data['t_ref_count'] = t_ref_count
        maf_data['t_alt_count'] = t_alt_count
        maf_data['n_depth'] = n_depth
        maf_data['n_ref_count'] = n_ref_count
        maf_data['n_alt_count'] = n_alt_count
        return maf_data

def resolve_vcf_variant_allele_data(vcf_data, maf_data):
        # get ref allele and alt allele values
        ref_allele,alt_allele = resolve_vcf_allele(vcf_data)
        start_pos = resolve_start_position(vcf_data)
        variant_type = ''
        end_pos = ''
        variant_class = ''

        if (ref_allele != 'N' or ref_allele != '') and alt_allele != '':
                # indels from vcf need to be shifted by one nucleotide and start position needs to be incremented by one
                if ref_allele[0] == alt_allele[0] and ref_allele != alt_allele:
                        # shift ref allele and alt allele by one nucleotide, set as '-' if len == 1
                        if len(ref_allele) == 1:
                                ref_allele = '-'
                        else:
                                ref_allele = ref_allele[1:]
                        if len(alt_allele) == 1:
                                alt_allele = '-'
                        else:
                                alt_allele = alt_allele[1:]
                        
                        # fix start position value
                        if start_pos != '':
                                start_pos = str(int(start_pos) + 1)

                # resolve variant type, end position, and variant class
                variant_type = resolve_vcf_variant_type(ref_allele, alt_allele)
                end_pos = resolve_end_position(vcf_data, start_pos, variant_type, ref_allele)
                variant_class = resolve_variant_classification(vcf_data, variant_type, ref_allele, alt_allele)

        maf_data['Variant_Classification'] = variant_class
        maf_data['Variant_Type'] = variant_type
        maf_data['Reference_Allele'] = ref_allele
        maf_data['Tumor_Seq_Allele1'] = ref_allele
        maf_data['Tumor_Seq_Allele2'] = alt_allele
        maf_data['Start_Position'] = start_pos
        maf_data['End_Position'] = end_pos
        return maf_data


def create_maf_record_from_vcf(sample_id, center_name, sequence_source, vcf_data, is_germline_data):
        """
                Creates MAF record from VCF data.
        """
        # init maf record
        maf_data = init_maf_record()

        # set easy to resolve values
        maf_data['Tumor_Sample_Barcode'] = sample_id
        maf_data['Center'] = center_name
        maf_data['Hugo_Symbol'] = 'Unknown'
        maf_data['Chromosome'] = resolve_chromosome(vcf_data)
        maf_data['Sequence_Source'] = resolve_sequence_source(vcf_data, sequence_source)
        maf_data['Verification_Status'] = DEFAULT_VERIFICATION_STATUS
        if is_germline_data:
                maf_data['Mutation_Status'] = "GERMLINE"
        else:
                maf_data['Mutation_Status'] = DEFAULT_MUTATION_STATUS
                        
        

        # resolve counts and variant allele data
        resolve_vcf_counts_data(vcf_data, maf_data)
        resolve_vcf_variant_allele_data(vcf_data, maf_data)
        return maf_data

def is_valid_vcf_maf_record(maf_record):
        return not ((maf_record['Reference_Allele'] == 'N' or maf_record['Reference_Allele'] == '') and maf_record['Tumor_Seq_Allele2'] == '')


def extract_vcf_data_from_file(filename, center_name, sequence_source):
        """
                Load data from a VCF file.
                VCF header:     CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR           
        """
        rejected_vars_by_varclass = 0
        rejected_vars_by_mutstatus = 0
        rejected_vars_by_inc_data = 0

        sample_id = get_vcf_sample_id(filename)
        is_germline_data = ('germline' in filename)
        data_file = open(filename, 'rU')
        data_reader = [line for line in data_file.readlines() if not line.startswith('#')] # ignore commented lines
        maf_data = []
        vcf_data_errors = []
        print 'Loading data from file:', filename
        for line in data_reader:
                vcf_data = dict(zip(VCF_HEADER, map(str.strip, line.split('\t'))))
                maf_record = create_maf_record_from_vcf(sample_id, center_name, sequence_source, vcf_data, is_germline_data)
                if not is_valid_vcf_maf_record(maf_record):
                        vcf_data_errors.append(vcf_data)
                if maf_record_already_processed(maf_record):
                        continue
                maf_data.append(maf_record)
        data_file.close()
        print_data_loading_summary(filename, len(maf_data), rejected_vars_by_varclass, rejected_vars_by_mutstatus, rejected_vars_by_inc_data, vcf_data_errors)
        return maf_data 


def extract_maf_data_from_file(filename, center_name, sequence_source):
        rejected_vars_by_varclass = 0
        rejected_vars_by_mutstatus = 0
        rejected_vars_by_inc_data = 0
        records_loaded = 0

        header = get_header(filename)
        data_file = open(filename, 'rU')
        data_reader = [line for line in data_file.readlines() if not line.startswith('#')][1:]
        maf_data = []
        print 'Loading data from file:', filename
        for line in data_reader:
                add_record = True
                if map(str.strip, line.split('\t')) == header:
                        print 'Double header in file! Skipping this line'
                        continue
                data = dict(zip(header, map(str.strip, line.split('\t'))))
                maf_record = create_maf_record(data, center_name, sequence_source)
                # if maf_record_already_processed(maf_record):
                #       continue

                # determine whether record should be skipped or not
                # if maf_record['Variant_Classification'] in MUTATION_FILTER:
                #       add_record = False
                #       rejected_vars_by_varclass += 1
                # elif maf_record['Reference_Allele'] == '' and maf_record['Tumor_Seq_Allele1'] == '' and maf_record['Tumor_Seq_Allele2'] == '':
                #       add_record = False
                #       rejected_vars_by_inc_data += 1
                # elif maf_record['Mutation_Status'] in MUTATION_STATUS_FILTER:
                #       add_record = False
                #       rejected_vars_by_mutstatus += 1

                # add record if not rejected
                if add_record:
                        maf_data.append(maf_record)
                records_loaded += 1
        data_file.close()
        print_data_loading_summary(filename, records_loaded, rejected_vars_by_varclass, rejected_vars_by_mutstatus, rejected_vars_by_inc_data, [])
        return maf_data


def print_data_loading_summary(filename, records_loaded, rejected_vars_by_varclass, rejected_vars_by_mutstatus, rejected_vars_by_inc_data, vcf_data_errors):
        print '\tTotal records loaded:', records_loaded
        if rejected_vars_by_varclass > 0:
                print '\tTotal records filtered by variant classification:', rejected_vars_by_varclass
        if rejected_vars_by_mutstatus > 0:
                print '\tTotal records rejected due to LOH or Wildtype Mutation Status:', rejected_vars_by_mutstatus
        if rejected_vars_by_inc_data > 0:
                print '\tTotal records rejected due to incomplete allele data:', rejected_vars_by_inc_data
        if len(vcf_data_errors) > 0:
                print '\tNumber of VCR data errors: ', len(vcf_data_errors)
        print


def get_mod_val(total_recs):
        if total_recs > 100000:
                mod_val = 20000
        elif total_recs > 50000 and total_recs <= 100000:
                mod_val = 10000
        elif total_recs > 10000 and total_recs <= 50000:
                mod_val = 5000
        else:
                mod_val = 1000
        return mod_val


def write_standardized_mutation_file(maf_data, output_filename):
        output_file = open(output_filename, 'w')
        output_file.write('\t'.join(MAF_HEADER))
        for data in maf_data:
                formatted_data = map(lambda x: process_datum(data.get(x, '')), MAF_HEADER)
                output_file.write('\n' + '\t'.join(formatted_data))
        output_file.close()
        print 'Standardized MAF written to: ' + output_filename


def generate_maf_from_input_data(input_directory, output_directory, extensions, center_name, sequence_source, has_vcf_data):
        print '\nLoading data from input directory:', input_directory
        print '\n\tSearching for files with extensions:', ', '.join(extensions) 

        for filename in os.listdir(input_directory):
                MERGED_MAF_DATA_KEYS = set()
                extract_data = False
                for ext in extensions:
                        if filename.endswith(ext) and 'normal' not in filename:
                                extract_data = True 
                                break
                if extract_data:
                        output_filename = os.path.join(output_directory, filename + '.temp')
                        if has_vcf_data and '.vcf' in filename:
                                maf_data = extract_vcf_data_from_file(os.path.join(input_directory, filename), center_name, sequence_source)                            
                        else:
                                maf_data = extract_maf_data_from_file(os.path.join(input_directory, filename), center_name, sequence_source)
                        write_standardized_mutation_file(maf_data, output_filename)


def maf_record_already_processed(maf_data):
        maf_key = get_maf_key(maf_data)
        if maf_key in MERGED_MAF_DATA_KEYS:
                return True
        else:
                MERGED_MAF_DATA_KEYS.add(maf_key)
                return False

def is_valid_record(maf_data):
        if not maf_data['Chromosome'] in VALID_CHROMOSOMES:
                return False
        return True

def usage():
        print 'python standardize_mutation_data.py --input-directory [path/to/input/directory] --output-directory [path/to/output/directory] --center [default name for center] --sequence-source [WGS | WXS] --extensions [comma-delimited list of extensions]'
        sys.exit(1)

def main():
        parser = optparse.OptionParser()
        parser.add_option('-i', '--input-directory', action = 'store', dest = 'inputdir', help = "input data directory")
        parser.add_option('-o', '--output-directory', action = 'store', dest = 'outputdir', help = "output data directory")
        parser.add_option('-c', '--center', action = 'store', dest = 'centername', help = "name of center (standard MAF field = 'Center')")
        parser.add_option('-s', '--sequence-source', action = 'store', dest = 'seqsource', help = "Sequencing source (standard MAF field = 'Sequencing_Source'), e.g., WXS or WGS")
        parser.add_option('-x', '--extensions', action = 'store', dest = 'extensions', help = "File extensions to look for (e.g., .vcf, .maf), comma-delimited list")

        (options, args) = parser.parse_args()
        input_directory = options.inputdir
        center_name = options.centername
        output_directory = options.outputdir
        sequence_source = options.seqsource
        extensions = options.extensions

        if not input_directory or not output_directory or not center_name or not extensions:
                print 'Input directory, output filename, extension(s) to search for, and center name must be provided.'
                usage()
        if not sequence_source:
                sequence_source = ''
        else:
                sequence_source = sequence_source.upper()

        has_vcf_data = ('vcf' in extensions)
        generate_maf_from_input_data(input_directory, output_directory, extensions.split(','), center_name, sequence_source, has_vcf_data)

if __name__ == '__main__':
        main()
