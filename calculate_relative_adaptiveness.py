#!/usr/bin/python
# -*- coding: UTF-8 -*-
# calculate_relative_adaptiveness.py
# written on 2/3/2017 by Xiangchen Li


import time
import argparse
from sys import exit


def read_rscu():    # set up a list of the codons
    """
    Function to read in the RSCU values.
    Output is a dictionary linking amino acids with 
    two parallel lists, their codons, and the codon's 
    RSCU values extracted from the input file.
    """
    line_number = 0
    codon_list = []
    aa_list = []
    with open(infile_name, 'r') as infile:
        for line in infile:
            line_number += 1
            line = line.strip("\n").split("\t")
            if line_number == 1:
                # gather the codons from the input file from the first line
                codon_list = line[1:]
                # set up parallel list of amino acids
                aa_list = []
                for i in codon_list:
                    aa = genetic_code[i]
                    aa_list.append(aa)
                tmp_rscu_dict = {}
                # set up dictionary to store the set of codons and rscu values for each amino acid
                for i in aa_list:
                    tmp_rscu_dict[i] = [[], []]
            if line_number == 2:
                rscu_list = line[1:]
                for i in range(len(rscu_list)):
                    codon = codon_list[i]
                    aa = aa_list[i]
                    rscu = rscu_list[i]
                    if rscu == 0:
                        rscu = 0.5
                    # this is to prevent w values of zero, which would reduce CAI to zero
                    # (see Behura and Severson 2013; Biological Reviews 88:49-61)
                    tmp_rscu_dict[aa][0].append(codon)
                    tmp_rscu_dict[aa][1].append(rscu)
    if debug:
        print("\nResults from gathering RSCU data:")
        aas = list(tmp_rscu_dict.keys())
        aas.sort()
        for i in aas:
            print(i)
            print(tmp_rscu_dict[i])
    return tmp_rscu_dict


def calculate_w(tmp_rscu_dict):
    """Function to calculate relative
    adaptiveness from the assembled rscu data"""
    # organize the amino acids and alphabetize them
    with open(outfile_name, 'w') as out:
        header = "codon\taa\twi\trscu"
        out.write(header)
        aa_list = list(tmp_rscu_dict.keys())
        aa_list.remove('STOP')
        aa_list.sort()
        for i in aa_list:
            rscus = tmp_rscu_dict[i][1]
            codons = tmp_rscu_dict[i][0]
            xmax = float(max(rscus))
            if len(codons) != len(rscus):
                exit("Error, don't have an rscu for each codon for {0}".format(i))
            for c in range(len(codons)):
                codon = codons[c]
                ru = float(rscus[c])
                w = ru / xmax
                out_string = "\n{0}\t{1}\t{2}\t{3}".format(codon, i, w, ru)
                out.write(out_string)


if __name__ == '__main__':
    program_name = 'calculate_relative_adaptiveness.py'
    last_update = '2/3/2017'
    author = 'Xiangchen Li'
    version_number = '1.0'
    print("\nRunning Program {0}...".format(program_name))

    version_string = '{0} version {1} Last Updated {2} by {3}'.format(program_name, version_number,
                                                                      last_update, author)
    description = '''
    Description:
    This program calculates the relative adaptiveness (w) of each codon
    from the RSCU values from a concatenated set of highly expressed genes.

    For each codon i that codes for an amino acid j, the relative
    adaptiveness wij is equal to the ratio of the RSCU for that
    codon to that of the most abundant synonymous codon:

    wij = RSCUi / RSCUimax

    (This is equal to the ratio of their counts)
    '''
    additional_program_info = '''
    Additional Program Information:
    To use this first get a set of highly expressed genes.
    Then extract and concatenate their coding sequences in to a single fasta.
    Calculate RSCU values from the concatenated sequence using get_RSCU_for_fasta.py.
    eg:
    get_RSCU_for_fasta.py -i top5_hiExpressed_contigs_concatenated.fasta -o top5_RSCU.txt

    Use the RSCU values as input for this script.
    '''
    start_time = time.time()
    # keeps track of how long the script takes to run
    # Set Up Argument Parsing
    # create argument parser that will automatically return help texts from global variables above
    parser = argparse.ArgumentParser(description=description, epilog=additional_program_info)
    parser.add_argument('-i', required=True, dest='input',
                        help='The the input file name with RSCUs for each codon from a set of highly expressed genes')
    parser.add_argument('-o', required=True, dest='output',
                        help='The the desired output file name. '
                             'Each codon will be grouped in a table with its '
                             'amino acid and its relative adaptiveness value')
    args = parser.parse_args()
    # Assign Arguments
    infile_name = args.input
    outfile_name = args.output
    debug = True
    # set up genetic code
    genetic_code = {"TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
                    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
                    "TAT": "Tyr", "TAC": "Tyr", "TAA": "STOP", "TAG": "STOP",
                    "TGT": "Cys", "TGC": "Cys", "TGA": "STOP", "TGG": "Trp",
                    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
                    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
                    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
                    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
                    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
                    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
                    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
                    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
                    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
                    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
                    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
                    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}
    rscu_dict = read_rscu()
    calculate_w(rscu_dict)
    # return time to run
    time = time.time() - start_time
    print('\nTime took to run: {0}'.format(str(time)))
