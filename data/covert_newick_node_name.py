#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: 
# Created by galaxy on 2017/3/9 0009 13:53

import os
from ete3 import Tree
from collections import defaultdict


my_path = os.getcwd()
organisms_file = os.path.join(my_path, 'organisms.txt')
strain_dict = defaultdict()
with open(organisms_file, 'r') as f1:
    for each_line in f1.readlines():
        a_list = each_line.strip().split('\t')
        strain_id = a_list[1]
        strain_name = a_list[0]
        strain_dict[strain_id] = strain_name
tmp_tree = os.path.join(my_path, '19_species_tree.nwk')
t = Tree(tmp_tree, format=1)
for node in t:
    a = node.name
    b = strain_dict[a]
    node.name = b
converted_tree = os.path.join(my_path, '19_species_tree_converted.nwk')
t.write(outfile=converted_tree, format=1)
