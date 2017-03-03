#!/usr/bin/python
# -*- coding: UTF-8 -*-
# Introduction: This script is used to calculate ENC (Nc) and GC3s from each genome.
# Created by galaxy on 2017/3/2 0002 21:38

import os
import sys
import time
import shutil
import subprocess
import logging.handlers
import traceback
from multiprocessing import Pool, cpu_count
from argparse import ArgumentParser


def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))


def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog="main_process.py")
    parser.add_argument("-i", "--indir", dest="indirname", action="store", default=None,
                        help="Input directory name")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store", default=None,
                        help="Output directory")
    parser.add_argument("-t", "--threads", type=int, dest="threads", default=cpu_count(),
                        help="How many threads will be used? [default all]")
    return parser.parse_args()


def run_codonw(each_type, fasta_file, out_file, blk_file):
    devnull = open(os.devnull, 'w')
    cmd = '{0} {1} {2} {3} -noblk -nomenu -nowarn -silent -human -fasta -enc'.format(each_type, fasta_file,
                                                                                     out_file, blk_file)
    try:
        subprocess.call(cmd, shell=True, stdout=devnull, stderr=devnull)
    except OSError:
        logger.info('Try to run codonw but failed, please check the dependency.'.format(args.part))
        logger.error(last_exception())
        sys.exit(1)


def multi_process(each_type):
    run_list = []
    type_dir = os.path.join(args.outdirname, each_type)
    if not os.path.exists(type_dir):
        os.makedirs(type_dir)
    else:
        shutil.rmtree(type_dir)
        os.makedirs(type_dir)
    out_dir = os.path.join(type_dir, 'out')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    blk_dir = os.path.join(type_dir, 'blk')
    if not os.path.exists(blk_dir):
        os.makedirs(blk_dir)
    for root, dirs, files in os.walk(args.indirname):
        for each_file in files:
            strain_id = os.path.splitext(each_file)[0]
            fasta_file = os.path.join(args.indirname, each_file)
            out_file = os.path.join(out_dir, strain_id + '.out')
            blk_file = os.path.join(blk_dir, strain_id + '.blk')
            run_list.append([fasta_file, out_file, blk_file])
    p = Pool(args.threads)
    for each_item in run_list:
        p.apply_async(run_codonw, args=(each_type, each_item[0], each_item[1], each_item[2]))
    p.close()
    p.join()


def main():
    type_list = ['enc', 'gc3s']
    for each_type in type_list:
        logger.info('Processing {0} Calculating'.join(each_type.upper()))
        multi_process(each_type)


if __name__ == '__main__':
    # Run as script
    # Parse command-line
    args = parse_cmdline()
    # Set up logging
    logger = logging.getLogger('main_process.py: %s' % time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    if not os.path.exists(args.outdirname):
        os.makedirs(args.outdirname)
    else:
        shutil.rmtree(args.outdirname)
        os.makedirs(args.outdirname)
    logger.info("Output directory: %s", args.outdirname)
    main()
    # Report that we've finished
    logger.info("All jobs have been done: %s.", time.asctime())
    logger.info("Total time taken: %.2fs", (time.time() - t0))
