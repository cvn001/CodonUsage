#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: Integrating all results
# Created by galaxy on 2017/3/9 0009 16:48

import os
from collections import defaultdict


def integrate_enc(output_dir):
    enc_dir = os.path.join(output_dir, 'ENC_GC3s')
    integrate_dir = os.path.join(output_dir, 'integrate_result')
    if not os.path.exists(integrate_dir):
        os.makedirs(integrate_dir)
    result_dict = defaultdict()
    for root, dirs, files in os.walk(enc_dir):
        for each_file in files:
            fname = os.path.splitext(each_file)[0]
            f_path = os.path.join(enc_dir, each_file)
            with open(f_path, 'r') as f1:
                result_line = f1.readlines()[1].strip()
                result_dict[fname] = result_line
    result_file = os.path.join(integrate_dir, 'all_ENC_GC3s.txt')
    with open(result_file, 'w') as f2:
        header = 'Strain\tNC\tGC3s\tENC\n'
        f2.write(header)
        for strain, result in result_dict.items():
            f2.write('{0}\t{1}\n'.format(strain, result))
    message = 'ENC results have been integrated.'
    return message
