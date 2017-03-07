#!/usr/bin/python
# -*- coding: UTF-8 -*-
# get_RSCU_from_fasta.py
# written on 2/3/2017 by Xiangchen Li

import time
import argparse
import numpy as np
from Bio import SeqIO
from collections import defaultdict


# use the genetic code dictionary above to build
# a new one to store counts in
def create_count_dict():
    """Function to make a clean dictionary
    to keep counts of codons in"""
    tmp_codon_list = list(genetic_code.keys())
    tmp_codon_list.sort()
    counts_dict = defaultdict()
    for i in genetic_code.keys():
        counts_dict[i] = 0
    return counts_dict


# we also need the reverse dictionary, where each list of synonymous
# codons is paired with its amino acid
def build_inverse_code():
    """
    Function to reverse the genetic code entered above
    so that amino acids key to lists of synonymous codons
    that code for them
    """
    aa_list = []
    rev_code = defaultdict()
    for i in codon_list:
        aa = genetic_code[i]
        if aa not in aa_list:
            aa_list.append(aa)
            rev_code[aa] = [i]
        else:
            rev_code[aa].append(i)
    return aa_list, rev_code


def get_counts():
    """
    Function to read the fasta file and
    count the codon use
    """
    gene_list = []
    data_dict = defaultdict()
    total_ambiguous = 0  # keep count of ambiguous codons that have 'N's in them
    total_codons = 0
    fas_seqs = SeqIO.parse(open(infile_name), 'fasta')
    for seqRecord in fas_seqs:
        seq = str(seqRecord.seq)
        seq_id = seqRecord.id
        gene_list.append(seq_id)
        codon_indices = range(0, len(seq), 3)
        seq_codons = []  # list of the codons in this sequence
        data_dict[seq_id] = create_count_dict()  # set up empty counts dictionary for this gene
        for i in codon_indices:
            total_codons += 1
            codon = seq[i:(i + 3)].upper()
            # skip ambiguous codons that have an N in them
            if "N" in codon:
                total_ambiguous += 1
                continue
            # skip ambiguous codons that have other characters
            try:
                genetic_code[codon]
            except KeyError:
                print("Warning, codon {0} is not in standard genetic code and will be ignored".format(codon))
                total_ambiguous += 1
                continue
            seq_codons.append(codon)
            data_dict[seq_id][codon] += 1
    print("\n{0} out of {1} Codons were ambiguous and not "
          "scored in the RSCU calculations".format(total_ambiguous, total_codons))
    return gene_list, data_dict


def get_aa_totals(aa_list, codon_count_dict):
    """Function to get the total number of times an amino
    acid appears in a sequence by totaling the counts of its
    synonymous codons"""
    aa_totals = defaultdict()
    for aa in aa_list:
        aa_totals[aa] = 0
    for codon in codon_list:
        codon_counts = codon_count_dict[codon]
        aa = genetic_code[codon]
        aa_totals[aa] += codon_counts
    return aa_totals


def calculate_rscu(aa_list, rev_code, gene_list, data_dict):
    """Go through each gene and convert the counts data to RSCU"""
    if interesting in aa_list:
        print("\n-----------------------------------")
        print("Results for amino acid of interest: {0}".format(interesting))
    rscu_dict = defaultdict()
    for gene in gene_list:
        rscu_dict[gene] = defaultdict()
        codon_count_dict = data_dict[gene]
        aa_totals = get_aa_totals(aa_list, codon_count_dict)
        for codon in codon_list:
            codon_count = codon_count_dict[codon]
            aa = genetic_code[codon]
            number_synonymous = float(len(rev_code[aa]))
            aa_count = float(aa_totals[aa])
            if aa_count == 0:
                rscu_dict[gene][codon] = "0.00"
                continue
            expected = aa_count / number_synonymous
            # rather than count missing codons as zero, all codons get a minimum count of 1
            # if codonCount == 0:
            #     codonCount += 1
            rscu = float(codon_count) / expected
            if codon == interesting:
                print("\n-----------------------------------")
                print("Results for codon of interest: {0}".format(interesting))
                print("codon Count = {0}".format(codon_count))
                print("amino acid = {0}".format(aa))
                print("amino acid count = {0}".format(aa_count))
                print("expected = {0}".format(expected))
                print("------------------------------------\n")
            if aa == interesting:
                print("count for codon {0} = {0}".format(codon, codon_count))
                print("amino acid count = {0}".format(aa_count))
                print("expected = {0}".format(expected))
            rscu_dict[gene][codon] = str(np.round(rscu, decimals=2))
    if interesting in aa_list:
        print("------------------------------------\n")
    return rscu_dict


def organize_codons(rev_code, aa_list):
    aa_list.sort()
    sorted_codon_list = []
    sorted_aa = []
    for aa in aa_list:
        s_codons = rev_code[aa]
        for sc in s_codons:
            sorted_codon_list.append(sc)
            sorted_aa.append(aa)
    return sorted_codon_list, sorted_aa


def output(rscu_dict, gene_list, sorted_codon_list, sorted_aa):
    """Function to output the data as a table"""
    header_list = ['contig']
    header_2_list = ['contig']
    for i in sorted_codon_list:
        header_list.append(i)
    for i in sorted_aa:
        header_2_list.append(i)
    if debug:
        print(len(header_list))
        print(header_list)
        print(len(header_2_list))
        print(header_2_list)
    with open(outfile_name, 'w') as out:
        header_1 = "\t".join(header_list)
        # header_2 = "\t".join(header_2_list)
        out.write(header_1)
        # out.write("\n" + header2)
        for gene in gene_list:
            rscu_values = []
            for codon in sorted_codon_list:
                rscu_values.append(rscu_dict[gene][codon])
            data_list = [gene] + rscu_values
            out.write("\n" + "\t".join(data_list))


if __name__ == '__main__':
    program_name = 'get_RSCU_from_fasta.py'
    last_updated = '2/3/2017'
    by = 'Xiangchen Li'
    version_number = '1.0'
    print("\nRunning Program {0}...".format(program_name))
    version_string = '{0} version {1} Last Updated {2} by {3}'.format(program_name, version_number,
                                                                      last_updated, by)
    description = '''
    Description:
    This program reads through a fasta file of coding sequences.
    It divides the sequence into codons assuming the reading frame begins with the
    first nucleotide.
    It outputs a table with each sequence ID linked to the RSCU (relative synonymous codon usage)
    for each codon and the species name (66 column table).
    '''
    additional_program_info = '''
    Additional Program Information:
    Note the fasta must be RNA.

    This script has been tested against output from dnaSP.

    Ambiguous codons that include any uncertain nucleotide
    such as N, W or Y are ignored in the analysis.
    '''
    start_time = time.time()  # keeps track of how long the script takes to run
    # Set Up Argument Parsing
    # create argument parser that will automatically return help texts from global variables above
    parser = argparse.ArgumentParser(description=description,
                                     epilog=additional_program_info)
    parser.add_argument('-i', required=True, dest='input', help='The the input file name')
    parser.add_argument('-o', required=True, dest='output', help='The the desired output file name')
    parser.add_argument('-interest', required=False, default='non', dest='interest',
                        help='Enter a codon or amino acid (three letter notation) you are interested '
                             'in to view its count data individually in standard output.')
    args = parser.parse_args()
    # Usage
    # get_RSCU_from_fasta.py -i sequences.fasta -o codon_RSCUs.txt
    # Assign Arguments
    infile_name = args.input
    outfile_name = args.output
    interesting = args.interest
    debug = False
    # Species = args.spp
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
    # set up a list of the codons
    codon_list = list(genetic_code.keys())  # global variable of all codons
    codon_list.sort()
    (all_aa_list, all_rev_code) = build_inverse_code()
    (all_gene_list, all_data_dict) = get_counts()
    all_rscu_dict = calculate_rscu(all_aa_list, all_rev_code, all_gene_list, all_data_dict)
    (the_sorted_codon_list, all_sorted_aa) = organize_codons(all_rev_code, all_aa_list)
    output(all_rscu_dict, all_gene_list, the_sorted_codon_list, all_sorted_aa)
    # return time to run
    time = time.time() - start_time
    print('\nTime took to run: {0}'.format(time))
