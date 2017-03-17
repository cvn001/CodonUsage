#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to get effective number of codons (ENC) from fasta
#               Input is a Seq format object (Bio.SeqIO)
# Created by galaxy on 2017/3/10 0010 17:37

# import os
# from Bio import SeqIO
from collections import defaultdict
from src.global_items import *


def degenerated():
    """
    This function is used to load degenercy AAs.
    :return:
    """
    degeneracy_dict = defaultdict()
    codon_dict = defaultdict()
    for amino_acid, codons in degenerated_codons.items():
        degeneracy_class = len(codons)
        degeneracy_dict[amino_acid] = degeneracy_class
        for codon in codons:
            codon_dict[codon] = amino_acid
    return degeneracy_dict, codon_dict


def read_seq(seq, codon_aa_dict):
    aa_dict = defaultdict(int)
    seq_codon_dict = defaultdict(int)
    max_seq = len(seq)
    query_codons = [seq[i:i + 3] for i in range(0, max_seq, 3)]
    for each_codon in query_codons:
        aa = codon_aa_dict[each_codon]
        aa_dict[aa] += 1
        seq_codon_dict[each_codon] += 1
    return aa_dict, seq_codon_dict


def enc_calculation(aa_dict, codon_dict, seq_codon_dict, degeneracy_dict):
    totb = defaultdict(int)
    numaa = defaultdict(int)
    error_t = False
    fold = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    for j in degeneracy_dict.keys():
        deg = degeneracy_dict[j]
        if j == 'STOP':  # skip STOP codon
            continue
        if aa_dict[j] <= 1:  # aa_dict[j]: nnaa + i
            bb = 0
        else:
            s2 = 0
            for x in codon_dict.keys():
                if codon_dict[x] != j:
                    continue
                if seq_codon_dict[x] == 0:
                    k2 = 0.0
                else:
                    k2 = pow((seq_codon_dict[x] / aa_dict[j]), 2)
                s2 += k2
            bb = (aa_dict[j] * s2 - 1.0) / (aa_dict[j] - 1.0)  # homozygosity
        if bb > 0.0000001:
            totb[deg] += bb
            numaa[deg] += 1
        fold[deg] += 1
    enc_tot = fold[1]
    for z in range(2, 9):
        if fold[z]:
            if numaa[z] and totb[z] > 0:
                averb = totb[z] / numaa[z]
            elif z == 3 and numaa[2] > 0 and numaa[4] > 0 and fold[z] == 1:
                averb = (totb[2] / numaa[2] + totb[4] / numaa[4]) * 0.5
            else:
                error_t = True
                break
            enc_tot += fold[z] / averb
    if error_t:
        result = 0
    elif enc_tot <= 61:
        result = enc_tot
    else:
        result = 61.00
    return result


def get_enc(query_seq):
    (degeneracy_dict, codon_dict) = degenerated()
    seq_string = str(query_seq).upper().replace('U', 'T')
    (aa_dict, seq_codon_dict) = read_seq(seq_string, genetic_code)
    result = enc_calculation(aa_dict, codon_dict, seq_codon_dict, degeneracy_dict)
    return result


# if __name__ == '__main__':
#     my_path = os.getcwd()
#     seq_file = os.path.join(my_path, 'test.fna')
#     for each_record in SeqIO.parse(seq_file, 'fasta'):
#         seq_enc = get_enc(each_record.seq)
#         print(seq_enc)
