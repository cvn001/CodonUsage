#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to calculate similarity index (SiD)
# Created by Xiangchen Li on 2017/3/19 21:15

from collections import defaultdict
from src.global_items import genetic_code


def get_sid(virus_rscu_file, host_rscu_file):
    for pass_codon in ["TAG", "TAA", "TGA", "ATG", "TGG"]:
        del genetic_code[pass_codon]
    virus_rscu_dict = defaultdict()
    with open(virus_rscu_file, 'r') as f1:
        for each_line in f1.readlines()[1:]:
            v_list = each_line.strip().split('\t')
            v_codon = v_list[0]
            v_rscu = v_list[1]
            virus_rscu_dict[v_codon] = float(v_rscu)
    host_rscu_dict = defaultdict()
    with open(host_rscu_file, 'r') as f2:
        for each_line in f2.readlines()[1:]:
            h_list = each_line.strip().split('\t')
            h_codon = h_list[0]
            h_rscu = h_list[1]
            host_rscu_dict[h_codon] = float(h_rscu)
    aa = 0
    bb = 0
    cc = 0
    for codon in genetic_code.keys():
        aa += virus_rscu_dict[codon] * host_rscu_dict[codon]
        bb += pow(virus_rscu_dict[codon], 2)
        cc += pow(host_rscu_dict[codon], 2)
    """
    R(A,B) is defined as the cosine value of the angle included
    between the A and B spatial vectors, and represents the degree of
    similarity between the virus and host overall codon usage patterns.

    D(A,B) represents the potential effect of the overall codon usage
    of the host on that of virus, and its value ranges from 0 to 1.0.
    """
    rr = aa / pow(bb * cc, 0.5)     # rr -> R(A,B)
    dd = (1 - rr) / 2       # dd -> D(A,B)
    return dd
