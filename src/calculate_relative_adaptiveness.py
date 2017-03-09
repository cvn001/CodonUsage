#!/usr/bin/python
# -*- coding: UTF-8 -*-
# calculate_relative_adaptiveness.py
# written on 2/3/2017 by Xiangchen Li


from sys import exit

from src.global_items import genetic_code


def read_rscu(input_file, debug):    # set up a list of the codons
    """
    Function to read in the RSCU values.
    Output is a dictionary linking amino acids with 
    two parallel lists, their codons, and the codon's 
    RSCU values extracted from the input file.
    """
    line_number = 0
    codon_list = []
    aa_list = []
    with open(input_file, 'r') as infile:
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
                rscu_dict = {}
                # set up dictionary to store the set of codons and rscu values for each amino acid
                for i in aa_list:
                    rscu_dict[i] = [[], []]
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
                    rscu_dict[aa][0].append(codon)
                    rscu_dict[aa][1].append(rscu)
    if debug:
        print("\nResults from gathering RSCU data:")
        aas = list(rscu_dict.keys())
        aas.sort()
        for i in aas:
            print(i)
            print(rscu_dict[i])
    return rscu_dict


def wi(rscu_dict, output_file):
    """Function to calculate relative
    adaptiveness from the assembled rscu data"""
    # organize the amino acids and alphabetize them
    with open(output_file, 'w') as out:
        header = "codon\taa\twi\trscu"
        out.write(header)
        aa_list = list(rscu_dict.keys())
        # aa_list.remove('STOP')
        aa_list.sort()
        for i in aa_list:
            rscus = rscu_dict[i][1]
            codons = rscu_dict[i][0]
            xmax = float(max(rscus))
            if len(codons) != len(rscus):
                exit("Error, don't have an rscu for each codon for {0}".format(i))
            for c in range(len(codons)):
                codon = codons[c]
                ru = float(rscus[c])
                w = ru / xmax
                out_string = "\n{0}\t{1}\t{2}\t{3}".format(codon, i, w, ru)
                out.write(out_string)


def calculate_wi(input_file, output_file):
    # program_name = 'calculate_relative_adaptiveness.py'
    # last_update = '2/3/2017'
    # author = 'Xiangchen Li'
    # version_number = '1.0'
    # print("\nRunning Program {0}...".format(program_name))
    #
    # version_string = '{0} version {1} Last Updated {2} by {3}'.format(program_name, version_number,
    #                                                                   last_update, author)
    # description = '''
    # Description:
    # This program calculates the relative adaptiveness (w) of each codon
    # from the RSCU values from a concatenated set of highly expressed genes.
    #
    # For each codon i that codes for an amino acid j, the relative
    # adaptiveness wij is equal to the ratio of the RSCU for that
    # codon to that of the most abundant synonymous codon:
    #
    # wij = RSCUi / RSCUimax
    #
    # (This is equal to the ratio of their counts)
    # '''
    # additional_program_info = '''
    # Additional Program Information:
    # To use this first get a set of highly expressed genes.
    # Then extract and concatenate their coding sequences in to a single fasta.
    # Calculate RSCU values from the concatenated sequence using get_RSCU_for_fasta.py.
    # eg:
    # get_RSCU_for_fasta.py -i top5_hiExpressed_contigs_concatenated.fasta -o top5_RSCU.txt
    #
    # Use the RSCU values as input for this script.
    # '''
    debug = False
    rscu_dict = read_rscu(input_file, debug)
    wi(rscu_dict, output_file)
