#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to calculate ENC (Nc) and GC3s from each genome.
# Created by galaxy on 2017/3/2 0002 21:38

import os
import sys
import re
import subprocess


def enc_expected_fun(s):
    enc_expected = 2 + s + 29 / (s ** 2 + (1 - s) ** 2)
    return enc_expected


def run_codonw(each_type, fasta_file, out_file):
    devnull = open(os.devnull, 'w')
    cmd = '{0} {1} {2} -noblk -nomenu -nowarn -silent -human -fasta -enc'.format(each_type,
                                                                                 fasta_file,
                                                                                 out_file)
    try:
        subprocess.call(cmd, shell=True, stdout=devnull, stderr=devnull)
    except OSError:
        sys.exit(1)


def main_process(input_file, output_file):
    run_list = [input_file, output_file]
    gc3s_bin = os.path.join('src', 'gc3s')
    run_codonw(gc3s_bin, run_list[0], run_list[1])
    tmp_result_file = run_list[1]
    final_result_lines = ''
    with open(tmp_result_file, 'r') as f1:
        result_lines = f1.readlines()
        tmp_header = result_lines[0].strip()
        result_header = re.sub(r'\s+', '\t', tmp_header) + '\tENC_expected\n'
        a_line = result_lines[1].strip()
        b_line = re.sub(r'\s+', '\t', a_line)
        b_list = b_line.strip().split('\t')
        gc3s = float(b_list[2])
        enc_expected = enc_expected_fun(gc3s)
        result_line = '{0}\t{1}\t{2}\n'.format(b_list[1], b_list[2], str(enc_expected))
        final_result_lines += result_header + result_line
    final_result_file = run_list[1]
    os.remove(tmp_result_file)
    with open(final_result_file, 'w') as f2:
        f2.write(final_result_lines)


def get_enc_gc3s(input_file, output_file):
    main_process(input_file, output_file)


