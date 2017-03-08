#!/usr/bin/python
# -*- coding: UTF-8 -*-
# get_CAI_from_fasta.py
# written on 2/3/2017 by Xiangchen Li

# import time
# import argparse
from sys import exit
from scipy import stats
from Bio import SeqIO


def read_wi_table(wi_file, wi_column):
    """
    Reads in the relative adaptiveness file
    """
    tmp_wi_dict = {}
    line_number = 0
    with open(wi_file, 'r') as infile:
        for line in infile:
            line_number += 1
            if line_number == 1:
                continue  # skip header
            line = line.strip("\n").split("\t")
            codon = line[0]
            wi = line[wi_column]
            tmp_wi_dict[codon] = wi
    return tmp_wi_dict
    

def get_cai(input_file, wi_dict, nucleotide_list, exclude, cpg_codon_list, debug):
    """
    Function that reads through a set of coding sequences
    in fasta format and gets the geometric mean of the wi
    values for each codon in the sequence (The CAI).
    """
    num_seqs = 0
    tmp_seq_id_list = []
    tmp_cai_list = []
    fasta_seqs = SeqIO.parse(open(input_file), 'fasta')
    skipped_codon_count = 0
    for seq in fasta_seqs:
        num_seqs += 1
        seq_string = str(seq.seq)
        codon_count = len(seq_string) / 3
        index_series = range(0, codon_count * 3, 3)
        codon_list = []
        wi_list = []
        for i in index_series:
            # some codons will have ambiguous nucleotides, so skip those
            good = 1
            codon = seq_string[i:i+3].upper()
            for nucleotide in codon:
                if nucleotide in nucleotide_list:
                    continue
                else:
                    good = 0
            if good == 1:
                # now optionally exclude CpG codons based on the -exclude_cpg_codons argument
                if exclude: 
                    # either skip the CpG codons
                    if exclude == 'yes':
                        if codon in cpg_codon_list:
                            # skip the CpG codons
                            skipped_codon_count += 1
                            continue
                        else:
                            # include non CpG codons
                            codon_list.append(codon)
                    # or skip all other codons
                    elif exclude == 'only':
                        if codon not in cpg_codon_list:
                            # skip the CpG codons
                            skipped_codon_count += 1
                            continue
                        else:
                            # include non CpG codons
                            codon_list.append(codon)
                    else:
                        exit("\nError, {} is an unexpected value given for -exclude_cpg_codons".format(exclude))
                else:
                    # Default: include all codons if running normally
                    codon_list.append(codon)
        for codon in codon_list:
            wi_list.append(float(wi_dict[codon]))
        if len(wi_list) > 5:
            cai = stats.gmean(wi_list)
        else:
            cai = 'NA'
        # record the results in two lists for output
        tmp_cai_list.append(cai)
        tmp_seq_id_list.append(seq.id)
        if debug:
            print("\n------------")
            print("gene = {0}".format(seq.id))
            print("codon_list:")
            print(codon_list)
            print("wi_list:")
            print(wi_list)
            print("CAI = {0}".format(cai))
    if exclude:
        if exclude == 'yes':
            print("\n{0} CpG codons were not considered in the calculation of CAI".format(skipped_codon_count))
        else:
            print("\n{0} non-CpG codons were not considered in the calculation of CAI".format(skipped_codon_count))
    return tmp_seq_id_list, tmp_cai_list


def output(output_file, seq_id_list, cai_list):
    """
    Output the results as a table
    """
    with open(output_file, 'w') as out:
        out.write('contig\tcai')
        for i in range(len(seq_id_list)):
            out_string = "{}\t{}".format(seq_id_list[i], cai_list[i])
            out.write('\n' + out_string)


def calculate_cai(input_file, output_file, wi_file, wi_column, exclude):
    # program_name = 'get_CAI_from_fasta.py'
    # last_update = '2/3/2017'
    # author = 'Xiangchen Li'
    # version_number = '1.0'
    # print("\nRunning Program {0}...".format(program_name))
    # version_string = '{0} version {1} Last Updated {2} by {3}'.format(program_name, version_number,
    #                                                                   last_update, author)
    # description = '''
    # Description:
    # This script calculates the codon adaptation index (CAI) for a set of sequences
    # based on an input set of relative adaptiveness values (wi) for each codon.
    # Relative adaptiveness values can be generated using calculate_relative_adaptiveness.py.
    #
    # Outputs a table of sequence IDs from the fasta file as rows linked with their
    # respective CAI values.
    # '''
    # additional_program_info = '''
    # Additional Program Information:
    #
    # See example of relative adaptiveness table at bottom of script.
    # '''
    # start_time = time.time()  # Keep tracking of how long the script takes to run
    # # Set Up Argument Parsing
    # # Create argument parser that will automatically return help texts from global variables above
    # parser = argparse.ArgumentParser(description=description, epilog=additional_program_info)
    # parser.add_argument('-i', required=False, dest='input', help='The the input fasta file you want CAI values for')
    # parser.add_argument('-o', required=True, dest='out', help='The desired name for the output file.')
    # parser.add_argument('-w', required=True, dest='wi',
    #                     help='A table of relative adaptiveness values (wi) for each codon')
    # parser.add_argument('-c', required=False, default=3, dest='col',
    #                     help='The column number where the wi values are in the relative adaptiveness table. '
    #                          'Codons are assumed to be in the first column. Default = 3rd column.')
    # parser.add_argument('-exclude_cpg_codons', required=False, default=None, dest='exclude_cpg',
    #                     help='Argument for whether to exclude CpG codons in the calculation of CAI. '
    #                          'Default is to include all codons. '
    #                          'Change to "yes" to exclude CpG codons from the calculation. '
    #                          'Change to "only" to base the CAI measures only on CpG codons.')
    # args = parser.parse_args()
    # Assign Arguments
    debug = False
    # infile_name = args.input
    # outfile_name = args.out
    # wi_file = args.wi
    # wi_column = int(args.col) - 1
    nucleotide_list = ['A', 'T', 'C', 'G']
    # exclude = args.exclude_cpg
    cpg_codon_list = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',  # serine
                      'CCT', 'CCC', 'CCA', 'CCG',                # proline
                      'ACT', 'ACC', 'ACA', 'ACG',                # threonine
                      'GCT', 'GCC', 'GCA', 'GCG',                # alanine
                      'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']  # arginine
    if exclude == 'yes':
        print("\nNote, calculating CAI with amino acids coded for by CpG codons excluded.")
    elif exclude == 'only':
        print("\nNote, calculating CAI based only on amino acids coded for by CpG codons.")
    elif not exclude:
        print("\nCalculating CAI normally with all amino acids included.")
    else:
        print("\nCalculating CAI normally with all amino acids included.")
    wi_dict = read_wi_table(wi_file, wi_column)
    (seq_id_list, cai_list) = get_cai(input_file, wi_dict, nucleotide_list,
                                      exclude, cpg_codon_list, debug)
    output(output_file, seq_id_list, cai_list)
    # return time to run
    # time = time.time() - start_time
    # print('\nTime took to run: {}'.format(time))
    # example of input wi table
    # these were calculated from a set of the top 5%
    # highest expressed genes in an A.millepora RNA-seq data
    # example_relative_adaptiveness_table = """
    # Example Relative Adaptiveness Table:
    # codon	aa	wi	rscu
    # GCA	Ala	0.956834532374	1.33
    # GCC	Ala	0.575539568345	0.8
    # GCG	Ala	0.338129496403	0.47
    # GCT	Ala	1.0	1.39
    # AGA	Arg	1.0	1.98
    # AGG	Arg	0.570707070707	1.13
    # CGA	Arg	0.484848484848	0.96
    # CGC	Arg	0.338383838384	0.67
    # CGG	Arg	0.222222222222	0.44
    # CGT	Arg	0.414141414141	0.82
    # AAC	Asn	0.834862385321	0.91
    # AAT	Asn	1.0	1.09
    # GAC	Asp	0.652892561983	0.79
    # GAT	Asp	1.0	1.21
    # TGC	Cys	0.754385964912	0.86
    # TGT	Cys	1.0	1.14
    # CAA	Gln	1.0	1.06
    # CAG	Gln	0.88679245283	0.94
    # GAA	Glu	1.0	1.21
    # GAG	Glu	0.652892561983	0.79
    # GGA	Gly	1.0	1.56
    # GGC	Gly	0.519230769231	0.81
    # GGG	Gly	0.339743589744	0.53
    # GGT	Gly	0.711538461538	1.11
    # CAC	His	0.754385964912	0.86
    # CAT	His	1.0	1.14
    # ATA	Ile	0.478571428571	0.67
    # ATC	Ile	0.671428571429	0.94
    # ATT	Ile	1.0	1.4
    # CTA	Leu	0.344827586207	0.5
    # CTC	Leu	0.475862068966	0.69
    # CTG	Leu	0.793103448276	1.15
    # CTT	Leu	0.965517241379	1.4
    # TTA	Leu	0.565517241379	0.82
    # TTG	Leu	1.0	1.45
    # AAA	Lys	1.0	1.09
    # AAG	Lys	0.834862385321	0.91
    # ATG	Met	1.0	1.0
    # TTC	Phe	0.639344262295	0.78
    # TTT	Phe	1.0	1.22
    # CCA	Pro	1.0	1.63
    # CCC	Pro	0.39263803681	0.64
    # CCG	Pro	0.282208588957	0.46
    # CCT	Pro	0.779141104294	1.27
    # AGC	Ser	0.698529411765	0.95
    # AGT	Ser	0.882352941176	1.2
    # TCA	Ser	1.0	1.36
    # TCC	Ser	0.558823529412	0.76
    # TCG	Ser	0.389705882353	0.53
    # TCT	Ser	0.882352941176	1.2
    # ACA	Thr	1.0	1.53
    # ACC	Thr	0.490196078431	0.75
    # ACG	Thr	0.346405228758	0.53
    # ACT	Thr	0.777777777778	1.19
    # TGG	Trp	1.0	1.0
    # TAC	Tyr	0.960784313725	0.98
    # TAT	Tyr	1.0	1.02
    # GTA	Val	0.467153284672	0.64
    # GTC	Val	0.569343065693	0.78
    # GTG	Val	0.890510948905	1.22
    # GTT	Val	1.0	1.37
    # """
