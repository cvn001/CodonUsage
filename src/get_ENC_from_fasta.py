#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to get effective number of codons (ENC) from fasta
#               Input is a Seq format object (Bio.SeqIO)
# Created by galaxy on 2017/3/10 0010 17:37

import os
from Bio import SeqIO
from collections import defaultdict
from src.global_items import genetic_code, degenerated_codons


def degenerated():
    degeneracy_dict = defaultdict(list)
    codon_dict = defaultdict()
    for amino_acid, codons in degenerated_codons.items():
        degeneracy_class = len(codons)
        degeneracy_dict[amino_acid] = [codons, degeneracy_class]
        for codon in codons:
            codon_dict[codon] = codons
    return degeneracy_dict, codon_dict


def read_seq(seq, codon_dict, degeneracy_dict):
    max_seq = len(seq)
    query_codons = [seq[i:i + 3] for i in range(0, max_seq, 3)]
    codon_num = len(query_codons)
    counts = defaultdict(int)
    seq_degeneracy_dict = defaultdict(list)
    aa_dict = defaultdict(int)
    for codon in query_codons:
        amino_acid = genetic_code[codon]
        if amino_acid in degeneracy_dict:
            aa_degeneracy = degeneracy_dict[amino_acid][1]
            seq_degeneracy_dict[aa_degeneracy].append(amino_acid)
            counts[codon] += 1
            aa_dict[amino_acid] = 1
    data = defaultdict(float)
    for codon in query_codons:
        if codon in codon_dict:
            totals = sum(counts[deg] for deg in codon_dict[codon])
            frequency = counts[codon] / totals
        else:
            frequency = 'NA'
        data[codon] = frequency
    return data, seq_degeneracy_dict, codon_num, aa_dict


def calculator(data, degeneracy_dict, seq_degeneracy_dict, codon_num, aa_dict):
    aa_freq_dict = defaultdict(float)
    for amino_acid in aa_dict.keys():
        aa_codons = degeneracy_dict[amino_acid][0]
        aa_s = 0
        for each_codon in aa_codons:
            codon_frequency = data[each_codon]
            aa_s += codon_frequency ** 2
        aa_freq = (aa_s * codon_num - 1) / (codon_num - 1)
        aa_freq_dict[amino_acid] = aa_freq
    deg_dict = {2: 9, 3: 1, 4: 5, 6: 3}
    f_dict = defaultdict()
    for each_deg in deg_dict.keys():
        if each_deg != 3:
            f_sum = 0
            # deg_aa_dict = defaultdict()
            # for aa in seq_degeneracy_dict[each_deg]:
            #     deg_aa_dict[aa] = 1
            # # deg_aa_num = len(deg_aa_dict.keys())
            for each_aa in seq_degeneracy_dict[each_deg]:
                f_sum += aa_freq_dict[each_aa] / len(seq_degeneracy_dict[each_deg])
            f_dict[each_deg] = f_sum
    f_dict[3] = (f_dict[2] + f_dict[4]) / 2
    enc = 2 + 9 / f_dict[2] + 1 / f_dict[3] + 5 / f_dict[4] + 3 / f_dict[6]
    return enc


def get_enc(query_seq, precision=2):
    (degeneracy_dict, codon_dict) = degenerated()
    seq_string = str(query_seq).upper().replace('U', 'T')
    (data, seq_degeneracy_dict, codon_num, aa_dict) = read_seq(seq_string, codon_dict, degeneracy_dict)
    enc = calculator(data, degeneracy_dict, seq_degeneracy_dict, codon_num, aa_dict)
    return round(enc, precision)


if __name__ == '__main__':
    my_path = os.getcwd()
    seq_file = os.path.join(my_path, 'test.fna')
    for each_record in SeqIO.parse(seq_file, 'fasta'):
        seq_enc = get_enc(each_record.seq)
        print(seq_enc)
