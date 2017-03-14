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
    # prepare the dict of degenerated codons
    degeneracy_dict = defaultdict(list)
    codon_dict = defaultdict()
    for amino_acid, codons in degenerated_codons.items():
        degeneracy_class = len(codons)
        degeneracy_dict[amino_acid] = [codons, degeneracy_class]
        for codon in codons:
            codon_dict[codon] = codons
    return degeneracy_dict, codon_dict


def read_seq(seq, codon_dict, degeneracy_dict):
    # 读取序列的长度
    max_seq = len(seq)
    # 分离出所有的密码子形成一个数组
    query_codons = [seq[i:i + 3] for i in range(0, max_seq, 3)]
    codon_num = len(query_codons)
    print(codon_num)
    # 初始化一个记录
    counts = defaultdict(int)
    # amino_acid_dict = defaultdict(int)
    seq_degeneracy_dict = defaultdict(list)
    aa_dict = defaultdict(int)
    for codon in query_codons:
        # 获取每个codon所对应的amino acid
        amino_acid = genetic_code[codon]
        # 获取此氨基酸对应的兼并性
        aa_degeneracy = degeneracy_dict[amino_acid][1]
        # 统计该兼并性在u序列中出现的次数
        seq_degeneracy_dict[aa_degeneracy].append(amino_acid)
        # 统计该codon在序列出现的次数
        counts[codon] += 1
        # 统计出现的氨基酸
        aa_dict[amino_acid] = 1
    # actual calculation of frequencies
    data = defaultdict(float)
    for codon in query_codons:
        if codon in codon_dict:
            # amino_acid = genetic_code[codon]
            # 该氨基酸所对应的所有同义密码子的数量的总和
            totals = sum(counts[deg] for deg in codon_dict[codon])
            # amino_acid_dict[amino_acid] = totals
            # 该codon出现的次数除以其对应的氨基酸的所有同义密码子出现的次数
            frequency = counts[codon] / totals
        else:
            frequency = 1.00
        # print('{0}\t{1}\n'.format(genetic_code[codon], str(totals)))
        data[codon] = frequency
    return data, seq_degeneracy_dict, codon_num, aa_dict


def calculator(data, degeneracy_dict, seq_degeneracy_dict, codon_num, aa_dict):
    aa_freq_dict = defaultdict(float)
    # deg_dict = defaultdict(int)
    for amino_acid in aa_dict.keys():
        aa_codons = degeneracy_dict[amino_acid][0]  # 该氨基酸对应的密码子
        # aa_deg = degeneracy_dict[amino_acid][1]     # 该氨基酸对应的简并性
        aa_s = 0
        # deg_dict[aa_deg] += 1
        # 计算该氨基酸所对应的所有同义密码子的频率平方之和
        for each_codon in aa_codons:
            codon_frequency = data[each_codon]
            aa_s += codon_frequency ** 2
        # 该氨基酸在序列中出现的所有同义密码子的数量，也就是该氨基酸出现的次数
        # total_codons = amino_acid_dict[amino_acid]
        aa_freq = (aa_s * codon_num - 1) / (codon_num - 1)
        aa_freq_dict[amino_acid] = aa_freq
    # 兼并性列表
    degeneracy_list = [2, 3, 4, 6]
    f_dict = defaultdict()
    for each_deg in degeneracy_list:
        # aa_num = len(seq_degeneracy_dict[each_deg])
        f_mean = 0
        aa_dict = defaultdict()
        for aa in seq_degeneracy_dict[each_deg]:
            aa_dict[aa] = 1
        deg_aa_num = len(aa_dict.keys())
        for each_aa in aa_dict.keys():
            f_mean += aa_freq_dict[each_aa] / deg_aa_num
        f_dict[each_deg] = f_mean
    # 最后计算得到该序列的enc值
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
        print(each_record.id)
        print(seq_enc)
