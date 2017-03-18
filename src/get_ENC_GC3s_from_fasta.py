#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to calculate:
#               1. GC3s (s);
#               2. ENC;
#               3. ENC expected: (2 + s + 29/(s^2 + (1 − s)^2))
# Created by Xiangchen Li on 2017/3/18 17:10

from Bio import SeqIO
from collections import defaultdict
from src.get_GC3s_from_fasta import get_gc3s
from src.get_ENC_from_fasta import get_enc


def get_enc_gc3s(input_file, output_file, precision=2):
    seq_result_dict = defaultdict()
    for seq_record in SeqIO.parse(input_file, 'fasta'):
        seq_string = str(seq_record.seq)
        seq_id = str(seq_record.id)
        seq_s = get_gc3s(seq_string, precision)
        seq_enc = get_enc(seq_string, precision)
        expected_enc = 2 + seq_s + 29 / (pow(seq_s, 2) + pow((1-seq_s), 2))
        seq_result_dict[seq_id] = [seq_s, seq_enc, round(expected_enc, precision)]
    with open(output_file, 'w') as f1:
        header = 'ID\tGC3s\tENC\tENC_expected\n'
        f1.write(header)
        for seq, result in seq_result_dict.items():
            result_line = '{0}\t{1}\t{2}\t{3}\n'.format(seq, str(result[0]),
                                                        str(result[1]), str(result[2]))
            f1.write(result_line)
