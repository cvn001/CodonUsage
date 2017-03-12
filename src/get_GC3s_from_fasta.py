#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script used to get G+C content
#               3rd position of synonymous codons (GC3s)
#               from each fasta format sequence
#               Input is a Seq format object (Bio.SeqIO)
# Created by galaxy on 17-3-12 4:12pm


def read_seq(seq):
    seq_string = str(seq).upper()
    max_seq = len(seq_string)
    seq_codons = [seq_string[i:i + 3] for i in range(0, max_seq, 3)]
    return seq_codons


def get_gc3s(seq, precision=2):
    exclude_codons = ['ATG', 'TGG', 'TGG', 'TAA', 'TAG', 'TGA']
    seq_codons = read_seq(seq)
    gc_content = 0
    ex = 0
    for each_codon in seq_codons:
        s = each_codon[2]
        if each_codon not in exclude_codons:
            if 'G' == s or s == 'C':
                gc_content += 1
        else:
            ex += 1
    seq_gc3s = gc_content / (len(seq_codons) - ex)
    return round(seq_gc3s, precision)
