#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: There are some global items.
# Created by galaxy on 2017/3/8 0008 22:10

genetic_code = {"TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
                "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
                "TAT": "Tyr", "TAC": "Tyr", "TAA": "STOP", "TAG": "STOP",
                "TGT": "Cys", "TGC": "Cys", "TGA": "STOP", "TGG": "Trp",
                "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
                "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
                "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
                "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
                "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
                "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
                "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
                "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
                "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
                "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
                "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
                "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}

nucleotide_list = ['A', 'T', 'C', 'G']

cpg_codon_list = ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',  # serine
                  'CCT', 'CCC', 'CCA', 'CCG',  # proline
                  'ACT', 'ACC', 'ACA', 'ACG',  # threonine
                  'GCT', 'GCC', 'GCA', 'GCG',  # alanine
                  'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']  # arginine