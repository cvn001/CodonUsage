#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to get effective number of codons (ENC) from fasta
#               Input is a Seq format object (Bio.SeqIO)
# Created by galaxy on 2017/3/10 0010 17:37

from Bio import SeqIO
from collections import defaultdict
from src.global_items import genetic_code, degenerated_codons


def degenerated():
    # prepare the dict of degenerated codons
    degeneracy_dict = defaultdict(list)
    codon_dict = defaultdict()
    all_dict = defaultdict()
    for amino_acid, codons in degenerated_codons.items():
        degeneracy_class = len(codons)
        all_dict.setdefault(degeneracy_class, []).append(amino_acid)
        degeneracy_dict[amino_acid] = codons
        for codon in codons:
            codon_dict[codon] = codons
    return degeneracy_dict, codon_dict, all_dict


def read_seq(seq, codon_dict, degeneracy_dict):
    max_seq = len(seq)
    query_codons = [seq[i:i + 3] for i in range(0, max_seq, 3)]
    counts = defaultdict(int)
    amino_acid_dict = defaultdict(int)
    seq_degeneracy_dict = defaultdict(int)
    for codon in query_codons:
        amino_acid = genetic_code[codon]
        aa_degeneracy = degeneracy_dict[amino_acid][0]
        seq_degeneracy_dict[aa_degeneracy] += 1
        counts[codon] += 1
    # actual calculation of frequencies
    data = defaultdict(float)
    for codon in query_codons:
        if codon in codon_dict:
            amino_acid = genetic_code[codon]
            totals = sum(counts[deg] for deg in codon_dict[codon])
            amino_acid_dict[amino_acid] = totals
            frequency = float(counts[codon]) / totals
        else:
            frequency = 'NA'
        data[codon] = frequency
    return data, amino_acid_dict, seq_degeneracy_dict


def calculator(data, degeneracy_dict, amino_acid_dict, seq_degeneracy_dict, all_dict):
    aa_dict = defaultdict(float)
    for amino_acid in degeneracy_dict.keys():
        aa_codons = degeneracy_dict[amino_acid][0]
        aa_s = 0
        for each_codon in aa_codons:
            codon_frequency = data[each_codon]
            aa_s += codon_frequency ** 2
        total_codons = amino_acid_dict[amino_acid]
        aa_f = (total_codons * aa_s - 1) / (total_codons - 1)
        aa_dict[amino_acid] = aa_f
    degeneracy_list = [2, 3, 4, 6]
    f_dict = defaultdict()
    for each_deg in degeneracy_list:
        aa_num = seq_degeneracy_dict[each_deg]
        each_deg_aa_list = all_dict[each_deg]
        f_mean = 0
        for aa in each_deg_aa_list:
            f_mean += aa_dict[aa] / aa_num
        f_dict[each_deg] = f_mean
    enc = 2 + 9 / f_dict[2] + 1 / f_dict[3] + 5 / f_dict[4] + 3 / f_dict[6]
    return enc


def get_enc(query_seq, precision=2):
    (degeneracy_dict, codon_dict, all_dict) = degenerated()
    seq_string = str(query_seq).upper().replace('T', 'U')
    (data, amino_acid_dict, seq_degeneracy_dict) = read_seq(seq_string, codon_dict, degeneracy_dict)
    seq_enc = calculator(data, degeneracy_dict, amino_acid_dict, seq_degeneracy_dict, all_dict)
    return round(seq_enc, precision)


