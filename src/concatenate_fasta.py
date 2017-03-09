#!/usr/bin/env python
# concatenate_fasta.py
# written on 3/7/2017 by Xiangchen Li

# import time
# import argparse
from Bio import SeqIO


def read_sub(subset_file):
    """
    Reads in the subset file
    """
    sub_list = []
    with open(subset_file, 'r') as infile:
        for line in infile:
            sub_list.append(line.strip("\n"))
    print("\nLooking for {} sequences to concatenate".format(len(sub_list)))
    return sub_list


def read_file(input_file, output_file, nucleotide_list, subset_file):
    """
    Function to read in a file as a list of lists
    """
    # Iterate through the seqs
    seqs_num = 0
    with open(output_file, 'w') as out:
        out.write('>{0}\n'.format(output_file))
        fas_seqs = SeqIO.parse(open(input_file), 'fasta')
        if subset_file == "None":
            for seq in fas_seqs:
                seqs_num += 1
                new_seq = str(seq.seq)
                # replace any non-canonical nucleotides with N
                edit_seq = ''
                for i in new_seq:
                    if i in nucleotide_list:
                        edit_seq += i
                    else:
                        edit_seq += 'N'
                out.write(edit_seq)
        else:
            sub_list = read_sub(subset_file)
            for seq in fas_seqs:
                if seq.id in sub_list:
                    seqs_num += 1
                    new_seq = str(seq.seq)
                    # Replace any non-canonical nucleotides with N
                    edit_seq = ''
                    for i in new_seq:
                        if i in nucleotide_list:
                            edit_seq += i
                        else:
                            edit_seq += 'N'
                    out.write(edit_seq)
            if seqs_num == len(sub_list):
                print("\nNice, found all {0} of the sequences from the subset file".format(seqs_num))
            else:
                print("\nFound {0} of the {1} sequences in the subset file".format(str(seqs_num), str(len(sub_list))))


def concatenate(input_file, output_file, subset_file=None):
    # program_name = 'concatenate_fasta.py'
    # last_updated = '3/7/2017'
    # author = 'Xiangchen Li'
    # version_number = '1.0'
    # print("\nRunning Program {0}...".format(program_name))
    # version_string = '{0} version {1} Last Updated {2} by {3}'.format(program_name, version_number,
    #                                                                   last_updated, author)
    # description = '''
    # Description:
    # This script concatenates all sequences in a fasta file into one single long sequence.
    #
    # '''
    # additional_program_info = '''
    # Additional Program Information:
    # The purpose of this is to run an entire set of coding sequences as one long sequence in dnaSP.
    # If the argument -sub is supplied then only a subset of sequences given as a table will be concatenated.
    # Otherwise the entire set of sequences are concatenated.
    # '''
    # start_time = time.time()  # keeps track of how long the script takes to run
    # Set Up Argument Parsing
    # Create argument parser that will automatically return help texts from global variables above
    # parser = argparse.ArgumentParser(description=description, epilog=additional_program_info)
    # parser.add_argument('-i', '--input', required=False, dest='input', help='The the input file')
    # parser.add_argument('-o', '--output', required=True, dest='out', help='The desired name for the output file')
    # parser.add_argument('-sub', required=False, default="NONE", dest='sub',
    #                     help='A subset of sequence names to pull and concatenate')
    # args = parser.parse_args()

    # Assign Arguments
    # input_file_name = args.input
    # output_file_name = args.out
    # use_subset_file = args.sub
    nucleotide_list = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'N']
    read_file(input_file, output_file, nucleotide_list, subset_file)
    # Return time to run
    # time = time.time() - start_time
    # print('\nTime took to run: {0}'.format(time))
